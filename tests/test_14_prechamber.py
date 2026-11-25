"""Pre-chamber two-volume analysis: orifice flow, open-system heat release, jet metrics.

Verification:

- orifice flow against closed-form isentropic values (critical pressure ratio 0.5283
  and choked flow function 0.6847 for gamma = 1.4), continuity at the critical ratio,
  zero flow at equal pressures;
- the open-system energy balance reproduces isentropic blowdown of a rigid zone;
- a manufactured two-volume cycle (forward integration with prescribed Wiebe heat
  release in both volumes) is inverted from traced masses and from the orifice model,
  recovering both heat releases and the discharge coefficient.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from combusta import HYDROGEN, PhysicalBoundsError, SliderCrankGeometry, make_reactant_mixture
from combusta.prechamber import (
    PreChamberGeometry,
    analyze_pre_chamber,
    critical_pressure_ratio,
    estimate_discharge_coefficient,
    orifice_mass_flow,
    pre_chamber_metrics,
    zone_apparent_heat_release,
)
from combusta.synthetic import manufactured_pre_chamber_cycle
from combusta.thermo.nasa_polynomials import DRY_AIR

GEOM = SliderCrankGeometry(bore_m=0.086, stroke_m=0.086, conrod_m=0.145, compression_ratio=12.0)
PC = PreChamberGeometry(volume_m3=1.0e-6, orifice_count=6, orifice_diameter_m=1.2e-3, discharge_coefficient=0.7)
MIX = make_reactant_mixture(HYDROGEN, 2.0)


@pytest.mark.unit
class TestOrificeFlow:
    def test_critical_ratio_air(self) -> None:
        assert float(critical_pressure_ratio(1.4)) == pytest.approx(0.528282, rel=1e-6)

    def test_choked_flow_function(self) -> None:
        """Choked: m_dot sqrt(R T0) / (C_D A p0) = sqrt(gamma) (2/(gamma+1))^((gamma+1)/(2(gamma-1)))."""
        flow = orifice_mass_flow(1.0e6, 1000.0, 1.0e5, gamma=1.4, r_specific=287.0,
                                 flow_area_m2=1.0e-6, discharge_coefficient=1.0)
        expected = 1.0e-6 * 1.0e6 / math.sqrt(287.0 * 1000.0) * math.sqrt(1.4) * (2.0 / 2.4) ** 3.0
        assert float(flow) == pytest.approx(expected, rel=1e-12)
        assert math.sqrt(1.4) * (2.0 / 2.4) ** 3.0 == pytest.approx(0.684731, rel=1e-5)

    def test_choked_flow_independent_of_downstream(self) -> None:
        kw = dict(gamma=1.35, r_specific=300.0, flow_area_m2=2.0e-6, discharge_coefficient=0.8)
        a = orifice_mass_flow(5.0e6, 1500.0, 1.0e6, **kw)
        b = orifice_mass_flow(5.0e6, 1500.0, 2.0e6, **kw)
        assert float(a) == pytest.approx(float(b), rel=1e-14)

    def test_continuous_at_critical_ratio_and_zero_at_equal_pressure(self) -> None:
        kw = dict(gamma=1.3, r_specific=300.0, flow_area_m2=1.0e-6, discharge_coefficient=1.0)
        rc = float(critical_pressure_ratio(1.3))
        below = orifice_mass_flow(1.0e6, 900.0, 1.0e6 * rc * (1 - 1e-9), **kw)
        above = orifice_mass_flow(1.0e6, 900.0, 1.0e6 * rc * (1 + 1e-9), **kw)
        assert float(below) == pytest.approx(float(above), rel=1e-7)
        assert float(orifice_mass_flow(1.0e6, 900.0, 1.0e6, **kw)) == 0.0

    def test_subsonic_value(self) -> None:
        """Unchoked flow for p_T/p_0 = 0.9, gamma = 1.4 from the isentropic formula."""
        r, g = 0.9, 1.4
        expected = r ** (1 / g) * math.sqrt(2 * g / (g - 1) * (1 - r ** ((g - 1) / g)))
        flow = orifice_mass_flow(1.0e5, 300.0, 0.9e5, gamma=g, r_specific=287.0,
                                 flow_area_m2=1.0, discharge_coefficient=1.0)
        assert float(flow) == pytest.approx(1.0e5 / math.sqrt(287.0 * 300.0) * expected, rel=1e-12)

    def test_rejects_reversed_states(self) -> None:
        with pytest.raises(ValueError, match="swap"):
            orifice_mass_flow(1.0e5, 300.0, 2.0e5, gamma=1.4, r_specific=287.0,
                              flow_area_m2=1.0e-6, discharge_coefficient=0.7)

    def test_geometry_validation(self) -> None:
        assert PC.flow_area_m2 == pytest.approx(6 * math.pi * 1.2e-3**2 / 4)
        with pytest.raises(PhysicalBoundsError):
            PreChamberGeometry(volume_m3=1e-6, orifice_count=4, orifice_diameter_m=1e-3, discharge_coefficient=1.2)


@pytest.mark.unit
def test_rigid_zone_isentropic_blowdown() -> None:
    """dQ = 0 for a rigid zone whose state follows p ~ m^gamma while losing mass (constant gamma)."""
    from combusta.thermo.nasa_polynomials import NASAGasMixture

    # A gas with constant cp: NASA coefficients (a1, 0, 0, 0, 0, a6, a7).
    const = NASAGasMixture(name="const-cp", molar_mass_kg_per_mol=DRY_AIR.molar_mass_kg_per_mol,
                           low_coeffs=(3.5, 0, 0, 0, 0, 0, 0), high_coeffs=(3.5, 0, 0, 0, 0, 0, 0))
    gamma = 1.4
    v = 1.0e-6
    m = np.linspace(1.0e-5, 0.6e-5, 201)
    theta = np.linspace(0.0, 20.0, 201)
    p = 5.0e6 * (m / m[0]) ** gamma
    t = p * v / (m * const.specific_gas_constant)
    dm = np.gradient(m, theta)
    dp = np.gradient(p, theta)
    dq = zone_apparent_heat_release(p, dp, np.full_like(p, v), np.zeros_like(p), t, dm,
                                    const.enthalpy(t), const)
    scale = float(np.max(np.abs(v * dp / (gamma - 1.0))))
    assert float(np.max(np.abs(dq[2:-2]))) < 1e-4 * scale


@pytest.fixture(scope="module")
def cycle():
    return manufactured_pre_chamber_cycle(GEOM, PC, HYDROGEN, lambda_air=2.0, rpm=2000.0)


@pytest.mark.strict_audit
class TestManufacturedPreChamber:
    def analyze(self, cycle, **kw):
        return analyze_pre_chamber(
            cycle.crank_angle_deg, cycle.pre_chamber_pressure_pa, cycle.main_chamber_pressure_pa,
            cycle.main_chamber_volume_m3, PC, rpm=cycle.rpm, mixture_pc=MIX, mixture_mc=MIX,
            main_chamber_dv_dtheta_m3_per_deg=cycle.main_chamber_dv_dtheta_m3_per_deg, **kw,
        )

    def test_forward_model_conserves_mass_and_releases_pc_charge_share(self, cycle) -> None:
        total = cycle.pre_chamber_mass_kg + cycle.main_chamber_mass_kg
        np.testing.assert_allclose(total, cycle.trapped_mass_kg, rtol=1e-12)
        v_c = GEOM.displaced_volume_m3 / (GEOM.compression_ratio - 1.0)
        share = cycle.q_pc_total_j / (cycle.q_pc_total_j + cycle.q_mc_total_j)
        assert share == pytest.approx(PC.volume_m3 / (PC.volume_m3 + v_c), rel=1e-12)

    @pytest.mark.parametrize("mode", ["masses", "temperatures", "orifice"])
    def test_recovers_both_heat_releases(self, cycle, mode) -> None:
        kw = {
            "masses": dict(pre_chamber_mass_kg=cycle.pre_chamber_mass_kg,
                           main_chamber_mass_kg=cycle.main_chamber_mass_kg),
            "temperatures": dict(pre_chamber_temperature_k=cycle.pre_chamber_temperature_k,
                                 main_chamber_temperature_k=cycle.main_chamber_temperature_k),
            "orifice": dict(trapped_mass_kg=cycle.trapped_mass_kg),
        }[mode]
        result = self.analyze(cycle, **kw)
        assert result.mass_source == ("orifice" if mode == "orifice" else "trace")
        m = pre_chamber_metrics(result, rpm=cycle.rpm)
        assert m.q_pc_gross_j == pytest.approx(cycle.q_pc_total_j, rel=1e-3)
        assert m.q_mc_gross_j == pytest.approx(cycle.q_mc_total_j, rel=1e-4)
        # The orifice mode integrates the mass from spline-interpolated pressures: drift < 0.1 %.
        rtol = 1e-3 if mode == "orifice" else 1e-6
        np.testing.assert_allclose(result.pre_chamber_mass_kg, cycle.pre_chamber_mass_kg, rtol=rtol)

    def test_discharge_coefficient_recovered_from_mass_trace(self, cycle) -> None:
        result = self.analyze(cycle, pre_chamber_mass_kg=cycle.pre_chamber_mass_kg,
                              main_chamber_mass_kg=cycle.main_chamber_mass_kg)
        assert estimate_discharge_coefficient(result) == pytest.approx(0.7, rel=1e-3)
        orifice = self.analyze(cycle, trapped_mass_kg=cycle.trapped_mass_kg)
        with pytest.raises(ValueError, match="traced zone masses"):
            estimate_discharge_coefficient(orifice)

    def test_jet_metrics_are_consistent(self, cycle) -> None:
        result = self.analyze(cycle, trapped_mass_kg=cycle.trapped_mass_kg)
        m = pre_chamber_metrics(result, rpm=cycle.rpm)
        dp = cycle.pre_chamber_pressure_pa - cycle.main_chamber_pressure_pa
        assert m.dp_max_pa == pytest.approx(float(dp.max()), rel=1e-12)
        assert m.dp_min_pa < 0.0  # flow back into the pre-chamber at some point of the cycle
        # The jet starts after pre-chamber combustion begins (-18 deg) and before the peak.
        assert -18.0 < m.jet_start_deg < m.theta_dp_max_deg < m.jet_end_deg
        assert 0.0 < m.jet_mass_fraction < 1.0
        assert m.pc_energy_share == pytest.approx(
            cycle.q_pc_total_j / (cycle.q_pc_total_j + cycle.q_mc_total_j), rel=1e-3)
        assert m.mc_ca50_deg > m.pc_ca50_deg

    def test_requires_a_mass_source(self, cycle) -> None:
        with pytest.raises(ValueError, match="trapped_mass_kg"):
            self.analyze(cycle)


@pytest.mark.strict_audit
class TestPreChamberCaseFiles:
    @pytest.mark.parametrize("with_temperatures", [True, False])
    def test_case_file_round_trip(self, cycle, tmp_path, with_temperatures) -> None:
        from combusta import load_case, write_pre_chamber_case_directory, write_pre_chamber_report

        path = write_pre_chamber_case_directory(
            cycle, GEOM, PC, tmp_path / "pc", name="synthetic_pc", fuel_key="hydrogen",
            lambda_air=2.0, include_temperatures=with_temperatures,
        )
        loaded = load_case(path)
        assert loaded.pre_chamber is not None
        assert loaded.pre_chamber.geometry == PC
        assert loaded.case.metadata["pre_chamber_sha256"] == loaded.pre_chamber.sha256
        report = write_pre_chamber_report(path, tmp_path / "out")
        assert report.result.mass_source == ("trace" if with_temperatures else "orifice")
        assert report.metrics.q_pc_gross_j == pytest.approx(cycle.q_pc_total_j, rel=2e-3)
        assert report.metrics.q_mc_gross_j == pytest.approx(cycle.q_mc_total_j, rel=1e-3)
        if with_temperatures:
            assert report.discharge_coefficient_estimate == pytest.approx(0.7, rel=5e-3)
        else:
            assert report.discharge_coefficient_estimate is None
        text = report.report_path.read_text(encoding="utf-8")
        assert "synthetic verification data" in text and "Jet start" in text
        assert report.figure_path.is_file()

    def test_same_file_requires_explicit_pressure_column(self, cycle, tmp_path) -> None:
        from combusta import CaseFileError, load_case, write_pre_chamber_case_directory

        path = write_pre_chamber_case_directory(cycle, GEOM, PC, tmp_path / "pc", name="x",
                                                fuel_key="hydrogen", lambda_air=2.0)
        text = path.read_text(encoding="utf-8").replace('file = "prechamber.csv"\n', "")
        text = text.replace('[pre_chamber.columns]\npressure = "PC Pressure"\ntemperature = "PC Temperature"\n',
                            "")
        path.write_text(text, encoding="utf-8")
        with pytest.raises(CaseFileError, match="must name the pre-chamber pressure column"):
            load_case(path)

    def test_cli(self, cycle, tmp_path, capsys) -> None:
        from combusta import write_pre_chamber_case_directory
        from combusta.cli import main

        path = write_pre_chamber_case_directory(cycle, GEOM, PC, tmp_path / "pc", name="x",
                                                fuel_key="hydrogen", lambda_air=2.0)
        assert main(["prechamber", str(path.parent), "--out", str(tmp_path / "o"), "--no-figure"]) == 0
        assert "Pre-chamber share of the heat release" in capsys.readouterr().out
        plain = tmp_path / "plain"
        from combusta import manufactured_cycle, write_case_directory

        write_case_directory(manufactured_cycle("p", GEOM, HYDROGEN, lambda_air=2.0, rpm=2000.0).case, plain,
                             fuel_key="hydrogen")
        assert main(["prechamber", str(plain), "--out", str(tmp_path / "o2")]) == 1
        assert "no [pre_chamber] table" in capsys.readouterr().err
