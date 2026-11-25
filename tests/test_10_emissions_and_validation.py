"""Verification of emission post-processing and ASME V&V 20 validation utilities.

Checks use conservation laws, literature constants (NIST molar masses), and
closed-form statistics rather than values produced by the code under test.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from combusta import (
    DRY_AIR,
    HYDROGEN,
    ISO_OCTANE,
    SPECIES_MOLAR_MASS,
    CycleTrace,
    EngineCase,
    OperatingPoint,
    PhysicalBoundsError,
    SliderCrankGeometry,
    combustion_efficiency,
    compare_cases,
    compute_cycle_metrics,
    compute_gci,
    emission_index,
    exhaust_composition,
    experimental_uncertainty,
    input_uncertainty,
    numerical_uncertainty_from_gci,
    sensitivity_coefficient,
    species_mass_from_mole_fraction,
    specific_emission,
    uncertainty_from_gci,
    validate,
)

GEOM = SliderCrankGeometry(bore_m=0.086, stroke_m=0.086, conrod_m=0.145, compression_ratio=12.0)


@pytest.mark.unit
class TestExhaust:
    @pytest.mark.parametrize(
        ("species", "molar_mass_g"),
        [("NO", 30.006), ("NO2", 46.006), ("NOx", 46.006), ("CO", 28.010), ("CO2", 44.010), ("H2O", 18.015)],
    )
    def test_molar_masses_match_nist(self, species: str, molar_mass_g: float) -> None:
        assert SPECIES_MOLAR_MASS[species] * 1e3 == pytest.approx(molar_mass_g, abs=1.5e-3)

    @pytest.mark.parametrize("fuel", [HYDROGEN, ISO_OCTANE])
    @pytest.mark.parametrize("lam", [1.0, 1.6, 2.5])
    def test_exhaust_mass_equals_fuel_plus_air(self, fuel, lam: float) -> None:
        m_f = 2.0e-5
        ex = exhaust_composition(fuel, lam, m_f)
        m_exhaust = ex.total_moles * ex.mean_molar_mass_kg_per_mol
        assert m_exhaust == pytest.approx(m_f * (1.0 + lam * fuel.stoichiometric_afr), rel=1e-6)

    def test_stoichiometric_hydrogen_water_fraction(self) -> None:
        x_o2 = 0.20946
        expected = 1.0 / (1.0 + 0.5 * (1.0 - x_o2) / x_o2)
        assert exhaust_composition(HYDROGEN, 1.0, 1e-5).water_mole_fraction == pytest.approx(expected, rel=1e-12)

    def test_ppm_to_mass_wet_and_dry(self) -> None:
        kwargs = dict(fuel=HYDROGEN, lambda_air=2.0, fuel_mass_kg=1.0e-5)
        ex = exhaust_composition(HYDROGEN, 2.0, 1.0e-5)
        wet = species_mass_from_mole_fraction(100e-6, "NOx", **kwargs)
        assert wet == pytest.approx(100e-6 * ex.total_moles * 46.0055e-3, rel=1e-9)
        dry = species_mass_from_mole_fraction(100e-6, "NOx", basis="dry", **kwargs)
        assert dry / wet == pytest.approx(1.0 - ex.water_mole_fraction, rel=1e-12)
        with pytest.raises(KeyError):
            species_mass_from_mole_fraction(1e-6, "XYZ", **kwargs)
        with pytest.raises(PhysicalBoundsError):
            species_mass_from_mole_fraction(1.5, "NOx", **kwargs)

    def test_specific_emission_and_index_definitions(self) -> None:
        assert specific_emission(2.0e-6, 500.0) * 3.6e9 == pytest.approx(2.0e-6 / 500.0 * 3.6e9)
        assert emission_index(2.0e-6, 1.0e-5) * 1e3 == pytest.approx(200.0)
        with pytest.raises(PhysicalBoundsError):
            specific_emission(1e-6, 0.0)

    def test_combustion_efficiency_heywood(self) -> None:
        m_f = 1.0e-5
        assert combustion_efficiency(HYDROGEN, m_f, {"H2": 0.02 * m_f, "NOx": 1e-7}) == pytest.approx(0.98, rel=1e-12)
        m_co, m_hc = 2.0e-7, 1.0e-7
        expected = 1.0 - (m_co * 10.10e6 + m_hc * ISO_OCTANE.lhv_j_per_kg) / (m_f * ISO_OCTANE.lhv_j_per_kg)
        assert combustion_efficiency(ISO_OCTANE, m_f, {"CO": m_co, "HC": m_hc}) == pytest.approx(expected, rel=1e-12)


def _motored_like_case(name: str, emissions: dict[str, float]) -> EngineCase:
    theta = np.linspace(-180.0, 180.0, 1441)
    kin = GEOM.evaluate(theta)
    mass = 1.0e5 * kin.volume_m3[0] / (DRY_AIR.specific_gas_constant * 330.0)
    x = np.clip((theta + 10.0) / 40.0, 0.0, 1.0)
    p = 1.0e5 * (kin.volume_m3[0] / kin.volume_m3) ** 1.35 * (1.0 + 1.5 * (1.0 - np.exp(-5.0 * x**3)))
    trace = CycleTrace(theta, p, temperature_k=p * kin.volume_m3 / (mass * DRY_AIR.specific_gas_constant))
    op = OperatingPoint(rpm=2000.0, fuel=HYDROGEN, fuel_mass_per_cycle_kg=1.0e-5, trapped_mass_kg=mass)
    return EngineCase(name, GEOM, op, trace, source="synthetic", emissions_kg_per_cycle=emissions)


@pytest.mark.unit
class TestEmissionMetrics:
    def test_metrics_include_specific_emissions_and_efficiency(self) -> None:
        case = _motored_like_case("a", {"NOx": 3.0e-8, "H2": 1.0e-7})
        m = compute_cycle_metrics(case)
        w = m["imep_gross"] * GEOM.displaced_volume_m3
        assert m["is_NOx"] == pytest.approx(3.0e-8 / w, rel=1e-12)
        assert m["ei_NOx"] == pytest.approx(3.0e-3, rel=1e-12)
        assert m["combustion_efficiency"] == pytest.approx(0.99, rel=1e-12)

    def test_comparison_rows_and_units(self) -> None:
        a = _motored_like_case("open", {"NOx": 3.0e-8})
        b = _motored_like_case("prechamber", {"NOx": 1.5e-8})
        cmp = compare_cases([a, b])
        row = cmp.row("is_NOx")
        assert row.display_unit == "g/kWh"
        assert row.relative_deltas()[1] == pytest.approx(-0.5, rel=1e-12)
        assert "Emission index NOx" in cmp.to_markdown()

    def test_negative_emission_mass_rejected(self) -> None:
        with pytest.raises(PhysicalBoundsError):
            _motored_like_case("bad", {"NOx": -1.0})


@pytest.mark.unit
class TestValidationVV20:
    def test_scalar_validation_inside_and_outside_uncertainty(self) -> None:
        inside = validate(10.5, 10.0, u_num=0.1, u_input=0.2, u_d=0.2)
        assert float(inside.u_val) == pytest.approx(0.3, rel=1e-12)
        assert float(inside.comparison_error) == pytest.approx(0.5)
        assert not bool(inside.error_exceeds_uncertainty)  # |E| = 0.5 <= 2 * 0.3
        assert "not distinguishable" in inside.summary()

        outside = validate(10.5, 9.0, u_num=0.1, u_input=0.2, u_d=0.2)
        assert bool(outside.error_exceeds_uncertainty)
        assert float(outside.model_error_lower) == pytest.approx(1.5 - 0.6)
        assert float(outside.model_error_upper) == pytest.approx(1.5 + 0.6)

    def test_pointwise_trace_validation(self) -> None:
        d = np.linspace(1.0, 2.0, 11)
        s = d + np.where(d > 1.5, 0.3, 0.0)
        res = validate(s, d, u_d=np.full(11, 0.05))
        assert np.array_equal(res.error_exceeds_uncertainty, d > 1.5)
        assert "45.5 %" in res.summary()

    def test_experimental_uncertainty_type_a_and_b(self) -> None:
        samples = np.array([1.0, 2.0, 3.0, 4.0, 5.0])  # s^2 = 2.5, s^2 / n = 0.5
        assert float(experimental_uncertainty(samples=samples)) == pytest.approx(math.sqrt(0.5), rel=1e-12)
        u = experimental_uncertainty(samples=samples, systematic_standard_uncertainty=math.sqrt(0.5))
        assert float(u) == pytest.approx(1.0, rel=1e-12)
        u2 = experimental_uncertainty(standard_deviation=math.sqrt(2.5), n_samples=5)
        assert float(u2) == pytest.approx(math.sqrt(0.5), rel=1e-12)
        with pytest.raises(ValueError):
            experimental_uncertainty(samples=np.array([1.0]))

    def test_central_difference_is_exact_for_quadratic(self) -> None:
        def model(x: float) -> float:
            return 3.0 * x**2 + 2.0 * x + 1.0

        c = sensitivity_coefficient(model(2.1), model(1.9), 0.1)
        assert float(c) == pytest.approx(14.0, rel=1e-12)

    def test_input_uncertainty_propagation(self) -> None:
        u = input_uncertainty([2.0, 3.0], [0.1, 0.2])
        assert float(u) == pytest.approx(math.sqrt(0.4), rel=1e-12)

    def test_numerical_uncertainty_from_gci(self) -> None:
        h = (1.0, 1.5, 2.25)
        phi = [10.0 + 2.0 * x**2 for x in h]
        gci = compute_gci(phi[0], phi[1], phi[2], h1=h[0], h2=h[1], h3=h[2])
        assert numerical_uncertainty_from_gci(gci) == pytest.approx(uncertainty_from_gci(gci) / 1.15, rel=1e-12)
        with pytest.warns(UserWarning, match="Divergent"):
            divergent = compute_gci(10.0, 12.0, 12.5, r21=1.5, r32=1.5)
        with pytest.raises(PhysicalBoundsError):
            numerical_uncertainty_from_gci(divergent)
