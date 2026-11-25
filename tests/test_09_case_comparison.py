"""Verification of cycle metrics and multi-case comparison against closed-form results.

References: ideal Otto-cycle efficiency (Heywood, 2018, Eq. 5.31), analytical
Wiebe burn fractions, and the Celik et al. (2008) GCI bound on a power-law triplet.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from combusta.synthetic import _ScalarNasa, _ScalarSliderCrank

from combusta import (
    DRY_AIR,
    ISO_OCTANE,
    CycleTrace,
    EngineCase,
    NASAGasMixture,
    OperatingPoint,
    PhysicalBoundsError,
    SliderCrankGeometry,
    compare_cases,
    compute_cycle_metrics,
    compute_gci,
    uncertainty_from_gci,
)

GEOM = SliderCrankGeometry(bore_m=0.086, stroke_m=0.086, conrod_m=0.145, compression_ratio=12.0)
THETA = np.linspace(-180.0, 180.0, 7201)
P_IVC = 1.0e5
T_IVC = 330.0
A_WIEBE = 5.0
M_WIEBE = 2.0


def _constant_gamma_gas(gamma: float) -> NASAGasMixture:
    a1 = gamma / (gamma - 1.0)
    coeffs = (a1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    return NASAGasMixture("constant-gamma", DRY_AIR.molar_mass_kg_per_mol, coeffs, coeffs)


def _wiebe_rate(theta: float, start: float, duration: float, q_total: float) -> float:
    u = (theta - start) / duration
    if u <= 0.0 or u >= 1.0:
        return 0.0
    norm = 1.0 - math.exp(-A_WIEBE)
    return q_total * A_WIEBE * (M_WIEBE + 1.0) / duration * u**M_WIEBE * math.exp(-A_WIEBE * u ** (M_WIEBE + 1.0)) / norm


def _wiebe_crossing(fraction: float, start: float, duration: float) -> float:
    norm = 1.0 - math.exp(-A_WIEBE)
    u = (-math.log(1.0 - fraction * norm) / A_WIEBE) ** (1.0 / (M_WIEBE + 1.0))
    return start + u * duration


def _synthetic_case(
    name: str,
    *,
    gas: NASAGasMixture,
    start: float,
    duration: float,
    q_total: float,
) -> EngineCase:
    """Adiabatic single-zone cycle dP/dθ = [(γ-1) dQ/dθ - γ P dV/dθ] / V with γ = γ(T)."""
    kin = GEOM.evaluate(THETA)
    v_ivc = float(kin.volume_m3[0])
    mass = P_IVC * v_ivc / (gas.specific_gas_constant * T_IVC)

    crank, scalar_gas = _ScalarSliderCrank(GEOM), _ScalarNasa(gas)  # plain-float forward model

    def rhs(th: float, y: np.ndarray) -> list[float]:
        v, dv = crank(th)
        t = float(y[0]) * v / (mass * gas.specific_gas_constant)
        g = scalar_gas.gamma(t)
        return [((g - 1.0) * _wiebe_rate(th, start, duration, q_total) - g * float(y[0]) * dv) / v]

    sol = solve_ivp(rhs, (THETA[0], THETA[-1]), [P_IVC], t_eval=THETA, rtol=1e-11, atol=1e-6, max_step=0.05)
    assert sol.success
    p = sol.y[0]
    trace = CycleTrace(
        crank_angle_deg=THETA,
        pressure_pa=p,
        temperature_k=p * kin.volume_m3 / (mass * gas.specific_gas_constant),
    )
    op = OperatingPoint(
        rpm=2000.0,
        fuel=ISO_OCTANE,
        fuel_mass_per_cycle_kg=q_total / ISO_OCTANE.lhv_j_per_kg if q_total > 0.0 else None,
        trapped_mass_kg=mass,
    )
    return EngineCase(name, GEOM, op, trace, source="synthetic")


@pytest.mark.strict_audit
class TestCycleMetrics:
    def test_near_instantaneous_combustion_reaches_otto_efficiency(self) -> None:
        gamma = 1.35
        case = _synthetic_case("otto", gas=_constant_gamma_gas(gamma), start=-1.0, duration=2.0, q_total=800.0)
        metrics = compute_cycle_metrics(case)
        eta_otto = 1.0 - GEOM.compression_ratio ** (1.0 - gamma)
        assert metrics["eta_indicated_gross"] == pytest.approx(eta_otto, rel=2.0e-3)
        assert metrics["imep_gross"] == pytest.approx(metrics["imep_net"], rel=1e-12)
        w = metrics["imep_gross"] * GEOM.displaced_volume_m3
        assert metrics["indicated_power_gross"] == pytest.approx(w * 2000.0 / 120.0, rel=1e-12)
        m_fuel = case.operating_point.fuel_mass_per_cycle_kg
        assert m_fuel is not None
        assert metrics["isfc_gross"] == pytest.approx(m_fuel / w, rel=1e-12)

    def test_wiebe_burn_fractions_recovered(self) -> None:
        start, duration = -10.0, 40.0
        case = _synthetic_case("wiebe", gas=DRY_AIR, start=start, duration=duration, q_total=900.0)
        metrics = compute_cycle_metrics(case, combustion_window_deg=(-20.0, 40.0))
        for key, frac in (("ca10", 0.1), ("ca50", 0.5), ("ca90", 0.9)):
            assert metrics[key] == pytest.approx(_wiebe_crossing(frac, start, duration), abs=0.02)
        assert metrics["heat_release_fraction"] == pytest.approx(1.0, abs=2.0e-3)
        assert any("frozen composition" in n for n in metrics.notes)

    def test_motored_trace_has_zero_work_and_no_combustion(self) -> None:
        case = _synthetic_case("motored", gas=DRY_AIR, start=0.0, duration=10.0, q_total=0.0)
        metrics = compute_cycle_metrics(case)
        # Reversible adiabatic compression-expansion: zero net work, peak pressure at TDC
        assert abs(metrics["imep_gross"]) < 50.0  # Pa, against IMEP scale ~1e6 Pa when fired
        assert metrics["theta_p_max"] == pytest.approx(0.0, abs=1e-3)
        assert metrics.burn is None and math.isnan(metrics["ca50"])
        assert any("No combustion detected" in n for n in metrics.notes)

    def test_inconsistent_polytrope_is_not_mistaken_for_combustion(self) -> None:
        kin = GEOM.evaluate(THETA)
        p = P_IVC * (kin.volume_m3[0] / kin.volume_m3) ** 1.30
        op = OperatingPoint(rpm=1500.0, fuel=ISO_OCTANE, trapped_mass_kg=4.0e-4)
        trace = CycleTrace(THETA, p, temperature_k=np.full_like(THETA, 600.0))
        metrics = compute_cycle_metrics(EngineCase("polytrope", GEOM, op, trace, source="synthetic"))
        v_ratio = float(kin.volume_m3[0] / np.min(kin.volume_m3))
        assert metrics["p_max"] == pytest.approx(P_IVC * v_ratio**1.30, rel=1e-12)
        assert metrics.burn is None


@pytest.fixture(scope="module")
def early_late() -> tuple[EngineCase, EngineCase]:
    early = _synthetic_case("early", gas=DRY_AIR, start=-12.0, duration=35.0, q_total=900.0)
    late = _synthetic_case("late", gas=DRY_AIR, start=3.0, duration=35.0, q_total=900.0)
    return early, late


@pytest.mark.unit
class TestCaseComparison:

    def test_deltas_follow_combustion_phasing(self, early_late) -> None:
        cmp = compare_cases(list(early_late), combustion_window_deg=(-25.0, 50.0))
        assert cmp.case_names == ("early", "late")
        assert cmp.row("ca50").deltas[1] == pytest.approx(15.0, abs=0.03)
        assert cmp.row("p_max").deltas[1] < 0.0
        assert cmp.row("eta_indicated_gross").deltas[1] < 0.0  # retarded combustion loses expansion work
        assert cmp.row("ca50").relative_deltas()[0] == 0.0
        assert all(s is None for s in cmp.row("p_max").significant)

    def test_significance_against_numerical_uncertainty(self, early_late) -> None:
        dp = compare_cases(list(early_late), combustion_window_deg=(-25.0, 50.0)).row("p_max").deltas[1]
        small = {"early": {"p_max": abs(dp) / 10}, "late": {"p_max": abs(dp) / 10}}
        large = {"early": {"p_max": abs(dp)}, "late": {"p_max": abs(dp)}}
        sig = compare_cases(list(early_late), uncertainties=small, combustion_window_deg=(-25.0, 50.0))
        non = compare_cases(list(early_late), uncertainties=large, combustion_window_deg=(-25.0, 50.0))
        assert sig.row("p_max").significant == (None, True)
        assert non.row("p_max").significant == (None, False)  # |Δ| < sqrt(2)|Δ|

    def test_reports_flag_synthetic_data(self, early_late) -> None:
        cmp = compare_cases(list(early_late), baseline="late", combustion_window_deg=(-25.0, 50.0))
        md = cmp.to_markdown()
        assert cmp.contains_synthetic and "Synthetic verification data" in md
        assert "late (baseline)" in md and "CA50" in md
        df = cmp.to_dataframe()
        assert len(df) == 2 * len(cmp.rows)
        assert set(df["case"]) == {"early", "late"}

    def test_input_validation(self, early_late) -> None:
        early, _ = early_late
        with pytest.raises(ValueError):
            compare_cases([early])
        with pytest.raises(ValueError):
            compare_cases([early, early])
        with pytest.raises(KeyError):
            compare_cases(list(early_late), uncertainties={"nope": {"p_max": 1.0}})


@pytest.mark.unit
class TestUncertaintyFromGCI:
    def test_gci_bounds_true_error_on_power_law(self) -> None:
        # phi(h) = 10 + 2 h^2, h = 1, 1.5, 2.25 (p = 2, r = 1.5): exact continuum value 10
        h = (1.0, 1.5, 2.25)
        phi = [10.0 + 2.0 * x**2 for x in h]
        gci = compute_gci(phi[0], phi[1], phi[2], h1=h[0], h2=h[1], h3=h[2])
        u = uncertainty_from_gci(gci)
        assert u == pytest.approx(1.25 * abs(phi[0] - phi[1]) / (1.5**2 - 1.0), rel=1e-10)
        assert abs(phi[0] - 10.0) < u

    def test_divergent_triplet_is_rejected(self) -> None:
        with pytest.warns(UserWarning, match="Divergent"):
            gci = compute_gci(10.0, 12.0, 12.5, r21=1.5, r32=1.5)  # eps21 > eps32: R = 4
        assert gci.convergence_regime == "DIVERGENT"
        with pytest.raises(PhysicalBoundsError):
            uncertainty_from_gci(gci)
