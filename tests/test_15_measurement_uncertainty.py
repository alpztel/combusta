"""Propagation of cylinder-pressure measurement errors into cycle metrics.

Exact properties of the measurement model are used as checks:

- a pegging offset does not change gross IMEP integrated from BDC to BDC
  (the volume returns to its starting value) and shifts peak pressure by the offset;
- a gain error scales peak pressure and IMEP by exactly (1 + epsilon);
- a TDC phasing error shifts the crank angle of peak pressure by exactly the error;
- the IMEP sensitivity to TDC phasing equals the integral of dP/dtheta dV/dtheta;
- for small uncertainties the Monte Carlo standard deviation agrees with linear
  propagation (JCGM 100 / JCGM 101).
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest
from scipy.integrate import trapezoid

from combusta import HYDROGEN, SliderCrankGeometry, manufactured_cycle
from combusta.measurement import (
    PressureMeasurementUncertainty,
    perturb_pressure,
    pressure_measurement_uncertainty,
)
from combusta.signal.filtering import compute_derivative

GEOM = SliderCrankGeometry(bore_m=0.086, stroke_m=0.086, conrod_m=0.145, compression_ratio=12.0)


@pytest.fixture(scope="module")
def bench_case():
    """Manufactured hydrogen cycle stripped to what a test bench measures (no temperature)."""
    case = manufactured_cycle("bench", GEOM, HYDROGEN, lambda_air=2.0, rpm=2000.0, wall_heat_loss_j=60.0).case
    return replace(case, trace=replace(case.trace, temperature_k=None), source="test_bench")


@pytest.mark.unit
class TestPerturbPressure:
    def test_offset_and_gain(self) -> None:
        theta = np.linspace(-10.0, 10.0, 41)
        p = np.full_like(theta, 2.0e6)
        out = perturb_pressure(theta, p, offset_pa=1.0e4, gain=0.02)
        np.testing.assert_allclose(out, 2.0e6 * 1.02 + 1.0e4)

    def test_phase_shift(self) -> None:
        theta = np.linspace(-90.0, 90.0, 721)
        p = 2.0e6 + 1.0e6 * np.cos(np.deg2rad(theta))
        out = perturb_pressure(theta, p, tdc_shift_deg=0.5)
        inner = slice(10, -10)
        np.testing.assert_allclose(out[inner], (2.0e6 + 1.0e6 * np.cos(np.deg2rad(theta + 0.5)))[inner], rtol=1e-9)

    def test_input_validation(self) -> None:
        from combusta import PhysicalBoundsError

        with pytest.raises(PhysicalBoundsError):
            PressureMeasurementUncertainty(tdc_offset_deg=-0.1)
        with pytest.raises(PhysicalBoundsError, match="not positive"):
            perturb_pressure(np.linspace(0, 1, 5), np.full(5, 1.0e5), offset_pa=-2.0e5)


@pytest.mark.strict_audit
class TestLinearBudget:
    def test_exact_properties(self, bench_case) -> None:
        inputs = PressureMeasurementUncertainty(tdc_offset_deg=0.3, pegging_offset_pa=8.0e3, gain_relative=0.015)
        r = pressure_measurement_uncertainty(bench_case, inputs, n_samples=0)
        imep, pmax = r.metrics["imep_gross"], r.metrics["p_max"]
        # Pegging offset: no effect on BDC-to-BDC IMEP, full effect on peak pressure.
        assert imep.contributions["offset"] < 1e-9 * imep.nominal
        assert pmax.contributions["offset"] == pytest.approx(8.0e3, rel=1e-9)
        # Gain: exact scaling.
        assert pmax.contributions["gain"] == pytest.approx(0.015 * pmax.nominal, rel=1e-9)
        assert imep.contributions["gain"] == pytest.approx(0.015 * imep.nominal, rel=1e-9)
        # TDC phasing shifts the location of peak pressure by the phasing error.
        assert r.metrics["theta_p_max"].contributions["tdc"] == pytest.approx(0.3, rel=1e-3)

    def test_imep_tdc_sensitivity_matches_integral(self, bench_case) -> None:
        """dIMEP/d(delta) = (1/V_d) * integral of dP/dtheta * dV/dtheta over the gross window."""
        theta = bench_case.trace.crank_angle_deg
        kin = bench_case.kinematics()
        dp = compute_derivative(theta, bench_case.trace.pressure_pa)
        expected = abs(trapezoid(dp * kin.dv_dtheta_m3_per_deg, theta)) / GEOM.displaced_volume_m3
        u = 0.25
        r = pressure_measurement_uncertainty(bench_case, PressureMeasurementUncertainty(tdc_offset_deg=u),
                                             n_samples=0)
        assert r.metrics["imep_gross"].contributions["tdc"] == pytest.approx(u * expected, rel=1e-2)
        # Of the order of several percent of IMEP per degree.
        per_degree = r.metrics["imep_gross"].contributions["tdc"] / u / r.metrics["imep_gross"].nominal
        assert 0.01 < per_degree < 0.15


@pytest.mark.strict_audit
class TestMonteCarlo:
    def test_agrees_with_linear_propagation(self, bench_case) -> None:
        inputs = PressureMeasurementUncertainty(tdc_offset_deg=0.1, pegging_offset_pa=3.0e3, gain_relative=0.005)
        r = pressure_measurement_uncertainty(bench_case, inputs, n_samples=400, seed=7)
        for key in ("p_max", "imep_gross", "ca50"):
            m = r.metrics[key]
            assert m.mc_valid == 400
            assert abs(m.mc_std - m.u_linear) < 3.0 * m.mc_std_se, key
            assert abs(m.mc_mean - m.nominal) < 4.0 * m.mc_std / np.sqrt(400), key

    def test_reproducible_and_noise_included(self, bench_case) -> None:
        from combusta.signal.filtering import savgol_zero_phase

        inputs = PressureMeasurementUncertainty(noise_pa=3.0e3)
        kw = dict(n_samples=40, seed=3, noise_filter=lambda th, p: savgol_zero_phase(p, 21, 3))
        a = pressure_measurement_uncertainty(bench_case, inputs, **kw)
        b = pressure_measurement_uncertainty(bench_case, inputs, **kw)
        assert a.metrics["ca50"].mc_std == b.metrics["ca50"].mc_std > 0.0
        assert np.isnan(a.metrics["ca50"].u_linear)  # noise has no linear contribution
        assert a.uncertainties()["ca50"] == a.metrics["ca50"].mc_std
        assert any("filter" in n for n in a.notes)

    def test_argument_validation(self, bench_case) -> None:
        with pytest.raises(ValueError, match="nothing to propagate"):
            pressure_measurement_uncertainty(bench_case, PressureMeasurementUncertainty())
        with pytest.raises(ValueError, match="Monte Carlo only"):
            pressure_measurement_uncertainty(bench_case, PressureMeasurementUncertainty(noise_pa=1e3), n_samples=0)
        with pytest.raises(ValueError, match="at least 2"):
            pressure_measurement_uncertainty(bench_case, PressureMeasurementUncertainty(gain_relative=0.01),
                                             n_samples=1)


@pytest.mark.strict_audit
class TestCaseFileAndCli:
    def write(self, bench_case, tmp_path, *, measurement: bool):
        from combusta import write_case_directory

        path = write_case_directory(bench_case, tmp_path / "bench", trace_format="csv", fuel_key="hydrogen")
        if measurement:
            text = path.read_text(encoding="utf-8")
            text += "\n[measurement]\ntdc_offset_deg = 0.2\npegging_offset_pa = 5000.0\ngain_relative = 0.01\n"
            path.write_text(text, encoding="utf-8")
        return path

    def test_case_file_table_and_report(self, bench_case, tmp_path) -> None:
        from combusta import load_case, write_measurement_report

        path = self.write(bench_case, tmp_path, measurement=True)
        loaded = load_case(path)
        assert loaded.measurement == PressureMeasurementUncertainty(
            tdc_offset_deg=0.2, pegging_offset_pa=5000.0, gain_relative=0.01)
        report = write_measurement_report(path, tmp_path / "out", n_samples=20)
        text = report.report_path.read_text(encoding="utf-8")
        assert "TDC phasing" in text and "Thermal shock" in text
        assert "imep_gross" in report.metrics_path.read_text(encoding="utf-8")

    def test_unknown_measurement_key_rejected(self, bench_case, tmp_path) -> None:
        from combusta import CaseFileError, load_case

        path = self.write(bench_case, tmp_path, measurement=True)
        path.write_text(path.read_text(encoding="utf-8").replace("gain_relative", "gain_percent"),
                        encoding="utf-8")
        with pytest.raises(CaseFileError, match="gain_percent"):
            load_case(path)

    def test_cli_options_override_and_errors(self, bench_case, tmp_path, capsys) -> None:
        from combusta.cli import main

        path = self.write(bench_case, tmp_path, measurement=False)
        assert main(["measurement", str(path), "--out", str(tmp_path / "o")]) == 1
        assert "no measurement uncertainties" in capsys.readouterr().err
        code = main(["measurement", str(path), "--tdc", "0.3", "--noise", "2000", "--samples", "10",
                     "--savgol", "21", "--out", str(tmp_path / "o")])
        assert code == 0
        out = capsys.readouterr().out
        assert "TDC phasing" in out and "M = 10" in out


@pytest.mark.strict_audit
def test_temperature_handling_is_reported() -> None:
    case = manufactured_cycle("cfd", GEOM, HYDROGEN, lambda_air=2.0, rpm=2000.0).case
    inputs = PressureMeasurementUncertainty(gain_relative=0.01)
    with_mass = pressure_measurement_uncertainty(case, inputs, n_samples=0)
    assert any("replaced" in n for n in with_mass.notes)
    no_mass = replace(case, operating_point=replace(case.operating_point, trapped_mass_kg=None))
    kept = pressure_measurement_uncertainty(no_mass, inputs, n_samples=0)
    assert any("kept" in n for n in kept.notes)
