"""Verification of fuel-agnostic thermodynamics, x_b-dependent heat release, and EngineCase.

Reference values come from independent literature (NIST-JANAF / ATcT tables,
Heywood 2018 App. D), not from the code under test.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from combusta import (
    HYDROGEN,
    ISO_OCTANE,
    METHANE,
    N_DODECANE,
    JET_A1_SURROGATE,
    CycleTrace,
    EngineCase,
    ForteThermoData,
    OperatingPoint,
    PhysicalBoundsError,
    SliderCrankGeometry,
    blend_burned_unburned,
    compute_heat_release,
    element_moles,
    make_product_mixture,
    make_reactant_mixture,
    product_moles,
    reactant_moles,
)
from combusta.thermo.nasa_polynomials import SPECIES_DB

R_U = 8.314462618


@pytest.mark.unit
class TestSpeciesThermo:
    @pytest.mark.parametrize(
        ("species", "cp_ref", "hf_ref_kj", "hf_tol_kj"),
        [
            ("H2", 28.836, 0.0, 0.01),  # NIST-JANAF
            ("H2O", 33.590, -241.826, 0.01),  # NIST-JANAF
            ("CH4", 35.69, -74.60, 0.05),  # NIST cp; ATcT / GRI-Mech 3.0 enthalpy of formation
            ("CO2", 37.12, -393.51, 0.01),  # NIST-JANAF
        ],
    )
    def test_cp_and_formation_enthalpy_at_298k(
        self, species: str, cp_ref: float, hf_ref_kj: float, hf_tol_kj: float
    ) -> None:
        sp = SPECIES_DB[species]
        w = sp.molar_mass_kg_per_mol
        assert float(sp.cp(298.15)) * w == pytest.approx(cp_ref, rel=2.0e-3)
        assert float(sp.enthalpy(298.15)) * w / 1e3 == pytest.approx(hf_ref_kj, abs=hf_tol_kj)

    @pytest.mark.parametrize("species", ["H2", "CH4"])
    def test_polynomial_continuity_at_1000k(self, species: str) -> None:
        sp = SPECIES_DB[species]
        t_lo, t_hi = np.nextafter(1000.0, 0.0), 1000.0
        assert float(sp.cp(t_lo)) == pytest.approx(float(sp.cp(t_hi)), rel=1e-5)
        assert float(sp.enthalpy(t_lo)) == pytest.approx(float(sp.enthalpy(t_hi)), rel=1e-5, abs=50.0)


@pytest.mark.unit
class TestStoichiometry:
    @pytest.mark.parametrize(
        ("fuel", "afr_ref"),
        [(HYDROGEN, 34.3), (METHANE, 17.2), (ISO_OCTANE, 15.1), (N_DODECANE, 15.0)],
    )
    def test_stoichiometric_afr_matches_literature(self, fuel, afr_ref: float) -> None:
        # Heywood (2018) App. D tabulates (A/F)s to three significant figures.
        assert fuel.stoichiometric_afr == pytest.approx(afr_ref, abs=0.06)

    @pytest.mark.parametrize(("fuel", "x_ref"), [(HYDROGEN, 0.295), (METHANE, 0.095)])
    def test_stoichiometric_fuel_volume_fraction(self, fuel, x_ref: float) -> None:
        moles = reactant_moles(fuel, 1.0)
        x_fuel = moles[fuel.vapor_species] / sum(moles.values())
        assert x_fuel == pytest.approx(x_ref, abs=1.5e-3)

    @pytest.mark.parametrize("fuel", [HYDROGEN, METHANE, ISO_OCTANE, N_DODECANE, JET_A1_SURROGATE])
    @pytest.mark.parametrize("lam", [1.0, 1.4, 2.5])
    def test_element_conservation_lean(self, fuel, lam: float) -> None:
        reactants = reactant_moles(fuel, lam, include_fuel_vapor=False)
        reactants["FUEL"] = 1.0
        before = element_moles(reactants, fuel)
        after = element_moles(product_moles(fuel, lam))
        for element in before:
            assert after[element] == pytest.approx(before[element], rel=1e-12, abs=1e-12)

    def test_rich_hydrogen_leaves_excess_h2_and_conserves_elements(self) -> None:
        lam = 0.8
        products = product_moles(HYDROGEN, lam)
        assert products["O2"] == 0.0
        assert products["H2"] == pytest.approx(1.0 - lam, rel=1e-9)
        before = element_moles(reactant_moles(HYDROGEN, lam))
        after = element_moles(products)
        for element in before:
            assert after[element] == pytest.approx(before[element], rel=1e-12, abs=1e-12)

    def test_rich_hydrocarbon_is_rejected(self) -> None:
        with pytest.raises(PhysicalBoundsError, match="water-gas shift"):
            product_moles(METHANE, 0.9)

    def test_stoichiometric_hydrogen_products(self) -> None:
        moles = product_moles(HYDROGEN, 1.0)
        assert moles["O2"] == 0.0
        x_h2o = moles["H2O"] / sum(moles.values())
        # H2 + 0.5 (air) -> H2O + inerts: 1 / (1 + 0.5 * (1 - x_O2) / x_O2)
        x_o2 = 0.20946
        assert x_h2o == pytest.approx(1.0 / (1.0 + 0.5 * (1.0 - x_o2) / x_o2), rel=1e-12)


@pytest.mark.unit
class TestBlendedProperties:
    def test_limits_and_monotonic_blend(self) -> None:
        unburned = make_reactant_mixture(HYDROGEN, 2.0)
        burned = make_product_mixture(HYDROGEN, 2.0)
        t = np.full(5, 1500.0)
        x_b = np.linspace(0.0, 1.0, 5)
        props = blend_burned_unburned(t, x_b, unburned, burned)
        assert props.gamma[0] == pytest.approx(float(unburned.gamma(1500.0)), rel=1e-12)
        assert props.gamma[-1] == pytest.approx(float(burned.gamma(1500.0)), rel=1e-12)
        assert np.all(np.diff(props.gamma) < 0.0)  # water-laden products have lower gamma
        assert np.allclose(props.cp_j_per_kg_k - props.cv_j_per_kg_k, props.r_specific_j_per_kg_k)

    def test_rejects_out_of_range_burn_fraction(self) -> None:
        mix = make_reactant_mixture(HYDROGEN, 2.0)
        with pytest.raises(PhysicalBoundsError):
            blend_burned_unburned(np.array([800.0]), np.array([1.2]), mix, mix)


def _wiebe(theta: np.ndarray, start: float, duration: float, a: float = 5.0, m: float = 2.0) -> np.ndarray:
    """Wiebe burn fraction normalized to reach exactly 1 at the end of combustion."""
    u = np.clip((theta - start) / duration, 0.0, 1.0)
    return (1.0 - np.exp(-a * u ** (m + 1.0))) / (1.0 - math.exp(-a))


@pytest.mark.strict_audit
class TestVariableCompositionHeatRelease:
    """Code verification by the method of manufactured solutions.

    The forward model integrates the single-zone pressure equation
    dP/dθ = [(γ-1) dQ/dθ - γ P dV/dθ] / V with γ = γ(T, x_b(θ)), T = P V / (m R(x_b)),
    and a prescribed Wiebe x_b(θ). The inverse must recover dQ/dθ and x_b(θ).
    This verifies the algorithm only; it is not a validation against engine data.
    """

    def test_recovers_manufactured_heat_release_and_burn_fraction(self) -> None:
        geom = SliderCrankGeometry(bore_m=0.086, stroke_m=0.086, conrod_m=0.145, compression_ratio=12.0)
        lam = 2.0
        unburned = make_reactant_mixture(HYDROGEN, lam)
        burned = make_product_mixture(HYDROGEN, lam)
        theta = np.linspace(-150.0, 120.0, 5401)
        kin = geom.evaluate(theta)
        v_ivc = float(kin.volume_m3[0])
        m = 1.0e5 * v_ivc / (unburned.specific_gas_constant * 330.0)
        m_fuel = m / (1.0 + lam * HYDROGEN.stoichiometric_afr)
        q_total = m_fuel * HYDROGEN.lhv_j_per_kg
        start, duration = -15.0, 45.0

        def x_b_of(th: float) -> float:
            return float(_wiebe(np.array([th]), start, duration)[0])

        def dq_of(th: float) -> float:
            u = min(max((th - start) / duration, 0.0), 1.0)
            if u <= 0.0 or u >= 1.0:
                return 0.0
            return q_total * 5.0 * 3.0 / duration * u**2 * math.exp(-5.0 * u**3) / (1.0 - math.exp(-5.0))

        def rhs(th: float, y: np.ndarray) -> list[float]:
            k = geom.evaluate(np.array([th]))
            v, dv = float(k.volume_m3[0]), float(k.dv_dtheta_m3_per_deg[0])
            xb = x_b_of(th)
            r_mix = (1 - xb) * unburned.specific_gas_constant + xb * burned.specific_gas_constant
            t = float(y[0]) * v / (m * r_mix)
            g = float(blend_burned_unburned(np.array([t]), np.array([xb]), unburned, burned).gamma[0])
            return [((g - 1.0) * dq_of(th) - g * float(y[0]) * dv) / v]

        sol = solve_ivp(rhs, (theta[0], theta[-1]), [1.0e5], t_eval=theta, rtol=1e-11, atol=1e-6, max_step=0.05)
        assert sol.success
        p = sol.y[0]
        x_b_true = _wiebe(theta, start, duration)
        r_true = (1 - x_b_true) * unburned.specific_gas_constant + x_b_true * burned.specific_gas_constant
        t_trace = p * kin.volume_m3 / (m * r_true)

        variable = compute_heat_release(
            theta, p, kin.volume_m3,
            temperature_k=t_trace,
            dv_dtheta_m3_per_deg=kin.dv_dtheta_m3_per_deg,
            mixture=unburned,
            burned_mixture=burned,
            combustion_window_deg=(start - 5.0, start + duration + 5.0),
        )
        frozen = compute_heat_release(
            theta, p, kin.volume_m3,
            temperature_k=t_trace,
            dv_dtheta_m3_per_deg=kin.dv_dtheta_m3_per_deg,
            mixture=unburned,
        )
        assert variable.mass_fraction_burned is not None
        assert 2 <= variable.property_iterations <= 20
        assert np.max(np.abs(variable.mass_fraction_burned - x_b_true)) < 2.0e-3

        q_var = float(variable.cumulative_gross_j[-1])
        q_frozen = float(frozen.cumulative_gross_j[-1])
        err_var = abs(q_var - q_total) / q_total
        err_frozen = abs(q_frozen - q_total) / q_total
        assert err_var < 2.0e-3
        assert err_frozen > 10.0 * err_var  # frozen unburned gamma biases the energy

    def test_single_mixture_path_is_unchanged(self) -> None:
        theta = np.linspace(-180.0, 180.0, 721)
        geom = SliderCrankGeometry(bore_m=0.08, stroke_m=0.09, conrod_m=0.15, compression_ratio=16.0)
        kin = geom.evaluate(theta)
        p = 1.0e5 * (kin.volume_m3[0] / kin.volume_m3) ** 1.35
        res = compute_heat_release(theta, p, kin.volume_m3, gamma=1.35, dv_dtheta_m3_per_deg=kin.dv_dtheta_m3_per_deg)
        assert res.mass_fraction_burned is None
        assert res.property_iterations == 0
        assert np.max(np.abs(res.dq_net_j_per_deg)) < 1e-4


@pytest.mark.unit
class TestEngineCase:
    @staticmethod
    def _geometry() -> SliderCrankGeometry:
        return SliderCrankGeometry(bore_m=0.086, stroke_m=0.086, conrod_m=0.145, compression_ratio=12.0)

    def _forte_like(self, theta_offset_deg: float = 0.0) -> ForteThermoData:
        theta = np.linspace(-180.0, 180.0, 361)
        kin = self._geometry().evaluate(theta + theta_offset_deg)
        p = 1.0e5 * (kin.volume_m3[0] / kin.volume_m3) ** 1.33
        return ForteThermoData(
            crank_angle_deg=theta,
            crank_angle_rad=np.deg2rad(theta),
            pressure_pa=p,
            volume_m3=kin.volume_m3,
            source_units={"pressure": "MPa", "volume": "cm^3"},
        )

    def test_forte_adapter_and_volume_consistency(self) -> None:
        op = OperatingPoint(rpm=2000.0, fuel=HYDROGEN, lambda_air=2.0, trapped_mass_kg=5.0e-4)
        case = EngineCase.from_forte_thermo("open chamber", self._forte_like(), geometry=self._geometry(), operating_point=op)
        assert case.source == "cfd" and not case.is_synthetic
        assert case.metadata["unit:pressure"] == "MPa"
        assert case.volume_consistency_error() < 1e-12

        shifted = EngineCase.from_forte_thermo(
            "shifted", self._forte_like(theta_offset_deg=1.0), geometry=self._geometry(), operating_point=op
        )
        assert shifted.volume_consistency_error() > 1e-2  # a 1 deg crank offset is detected

    def test_temperature_from_ideal_gas(self) -> None:
        op = OperatingPoint(rpm=2000.0, fuel=HYDROGEN, lambda_air=2.0, trapped_mass_kg=5.0e-4)
        thermo = self._forte_like()
        case = EngineCase.from_forte_thermo("c", thermo, geometry=self._geometry(), operating_point=op)
        r = make_reactant_mixture(HYDROGEN, 2.0).specific_gas_constant
        assert np.allclose(case.temperature_k(), thermo.pressure_pa * thermo.volume_m3 / (5.0e-4 * r))

    def test_fuel_energy_and_validation(self) -> None:
        op = OperatingPoint(rpm=1500.0, fuel=HYDROGEN, lambda_air=2.5, fuel_mass_per_cycle_kg=1.0e-5)
        assert op.fuel_energy_per_cycle_j == pytest.approx(1.0e-5 * 119.96e6)
        assert op.equivalence_ratio == pytest.approx(0.4)
        with pytest.raises(PhysicalBoundsError):
            OperatingPoint(rpm=-1.0, fuel=HYDROGEN)
        with pytest.raises(PhysicalBoundsError):
            OperatingPoint(rpm=1500.0, fuel=HYDROGEN, lambda_air=0.0)
        theta = np.linspace(0.0, 10.0, 11)
        with pytest.raises(ValueError):
            CycleTrace(crank_angle_deg=theta[::-1], pressure_pa=np.ones(11))
        with pytest.raises(PhysicalBoundsError):
            CycleTrace(crank_angle_deg=theta, pressure_pa=-np.ones(11))
        with pytest.raises(ValueError):
            EngineCase("x", self._geometry(), op, CycleTrace(theta, np.ones(11)), source="unknown")  # type: ignore[arg-type]
