"""Data intake: column mapping, generic trace tables, case files, reports, and the CLI.

The manufactured cycles of :mod:`combusta.synthetic` have a known heat release, so the
whole chain file -> reader -> case file -> metrics can be checked against the exact
answer (code verification by manufactured solutions; Roache, 1998).
"""

from __future__ import annotations

import hashlib
from io import StringIO
import math
from pathlib import Path

import numpy as np
import pytest

from combusta import (
    HYDROGEN,
    METHANE,
    CaseFileError,
    ColumnSpec,
    ForteParseError,
    SliderCrankGeometry,
    UnitDetectionError,
    compute_cycle_metrics,
    compute_heat_release,
    inspect_columns,
    load_case,
    make_product_mixture,
    make_reactant_mixture,
    manufactured_cycle,
    parse_thermo_csv,
    parse_trace_csv,
    write_case_directory,
    write_comparison_report,
)
from combusta.cli import main

GEOM = SliderCrankGeometry(bore_m=0.086, stroke_m=0.086, conrod_m=0.145, compression_ratio=12.0)


@pytest.fixture(scope="module")
def cycles():
    """Two manufactured hydrogen cycles with different burn rates and wall heat loss."""
    slow = manufactured_cycle(
        "slow", GEOM, HYDROGEN, lambda_air=2.0, rpm=2000.0,
        combustion_start_deg=-15.0, combustion_duration_deg=45.0, wall_heat_loss_j=60.0,
    )
    fast = manufactured_cycle(
        "fast", GEOM, HYDROGEN, lambda_air=2.0, rpm=2000.0,
        combustion_start_deg=-8.0, combustion_duration_deg=28.0, wall_heat_loss_j=60.0,
    )
    return slow, fast


@pytest.mark.unit
class TestColumnMapping:
    CSV = "Angle,Pcyl,Vcyl,Tgas\n" + "\n".join(
        f"{th},{2.0 + 0.001 * i},{400.0 + i},{350.0 + i}" for i, th in enumerate(range(-10, 10))
    )

    def test_explicit_columns_and_units(self) -> None:
        data = parse_thermo_csv(StringIO(self.CSV), columns={
            "crank_angle": ColumnSpec("Angle", "deg"),
            "pressure": ColumnSpec("Pcyl", "bar"),
            "volume": ColumnSpec("Vcyl", "cm^3"),
            "temperature": ColumnSpec("Tgas", "degC"),
        })
        assert data.pressure_pa[0] == pytest.approx(2.0e5)
        assert data.volume_m3[0] == pytest.approx(400.0e-6)
        assert data.temperature_k is not None
        assert data.temperature_k[0] == pytest.approx(350.0 + 273.15)
        assert data.source_units == {"crank_angle": "deg", "pressure": "bar", "volume": "cm^3", "temperature": "degC"}

    def test_unit_override_beats_header_unit(self) -> None:
        csv = "Crank Angle [deg],Pressure [kPa],Volume [m^3]\n-1,100,1e-4\n0,110,1e-4\n1,100,1e-4\n"
        data = parse_thermo_csv(StringIO(csv), columns={"pressure": ColumnSpec("Pressure", "bar")})
        assert data.pressure_pa[1] == pytest.approx(110.0e5)

    def test_missing_column_lists_available_headers(self) -> None:
        with pytest.raises(ForteParseError, match=r"Available columns: \['angle', 'pcyl', 'tgas', 'vcyl'\]"):
            parse_thermo_csv(StringIO(self.CSV))

    def test_mapped_column_not_in_file(self) -> None:
        with pytest.raises(ForteParseError, match="'Pressure' mapped to 'pressure' not found"):
            parse_thermo_csv(StringIO(self.CSV), columns={
                "crank_angle": ColumnSpec("Angle", "deg"), "pressure": "Pressure",
            })

    def test_unit_still_required(self) -> None:
        with pytest.raises(UnitDetectionError, match="missing or invalid unit"):
            parse_thermo_csv(StringIO(self.CSV), columns={
                "crank_angle": ColumnSpec("Angle", "deg"), "pressure": "Pcyl", "volume": ColumnSpec("Vcyl", "cc"),
            })

    def test_unknown_field_rejected(self) -> None:
        with pytest.raises(ValueError, match="Unknown field"):
            parse_thermo_csv(StringIO(self.CSV), columns={"presure": "Pcyl"})

    def test_duplicate_headers_are_suffixed_not_dropped(self) -> None:
        csv = "Crank Angle [deg],Pressure [bar],Pressure [bar],Volume [cc]\n0,10,20,50\n1,11,21,50\n"
        infos = inspect_columns(StringIO(csv))
        assert [i.name for i in infos] == ["crank_angle", "pressure", "pressure_2", "volume"]
        second = parse_thermo_csv(StringIO(csv), columns={"pressure": "pressure_2"})
        assert second.pressure_pa[0] == pytest.approx(20.0e5)

    def test_inspect_reports_matches(self) -> None:
        csv = "Crank Angle [deg],Pressure [MPa],Unknown thing,Total Cells\n0,1,2,1000\n1,1,2,1000\n"
        infos = {i.name: i for i in inspect_columns(StringIO(csv))}
        assert infos["crank_angle"].matched_fields == ("crank_angle",)
        assert infos["pressure"].unit == "MPa"
        assert infos["unknown_thing"].matched_fields == ()
        assert infos["total_cells"].matched_fields == ("cell_count",)


@pytest.mark.unit
class TestTraceCsv:
    def test_pressure_only_table_with_unit_row(self) -> None:
        csv = "Crank Angle,Cylinder Pressure\ndeg,bar\n-1,20\n0,25\n1,20\n"
        table = parse_trace_csv(StringIO(csv))
        assert table.volume_m3 is None and table.temperature_k is None
        np.testing.assert_allclose(table.pressure_pa, [20e5, 25e5, 20e5])

    def test_rejects_nonpositive_pressure(self) -> None:
        from combusta import PhysicalBoundsError

        with pytest.raises(PhysicalBoundsError):
            parse_trace_csv(StringIO("CA [deg],P [bar]\n0,1\n1,0\n"))


@pytest.mark.strict_audit
class TestCompositionDependentTemperature:
    """T = PV/(m R(x_b)) must be used when no temperature trace exists.

    Hydrogen combustion removes one third of the H2 + O2 moles (2 H2 + O2 -> 2 H2O), so a
    temperature from the unburned gas constant is biased high after combustion.
    """

    def test_recovers_heat_release_without_temperature_trace(self, cycles) -> None:
        slow, _ = cycles
        case = slow.case
        theta, p, v = case.trace.crank_angle_deg, case.trace.pressure_pa, case.trace.volume_m3
        unburned = make_reactant_mixture(HYDROGEN, 2.0)
        burned = make_product_mixture(HYDROGEN, 2.0)
        common = dict(
            dq_ht_j_per_deg=case.trace.wall_heat_loss_j_per_deg, mixture=unburned,
            burned_mixture=burned, combustion_window_deg=(-40.0, 90.0),
        )
        from_mass = compute_heat_release(theta, p, v, trapped_mass_kg=case.operating_point.trapped_mass_kg, **common)
        assert v is not None
        m = case.operating_point.trapped_mass_kg
        assert m is not None
        t_frozen_r = p * v / (m * unburned.specific_gas_constant)  # unburned R only
        frozen_r = compute_heat_release(theta, p, v, temperature_k=t_frozen_r, **common)
        q_true = slow.total_heat_release_j
        err_mass = abs(from_mass.cumulative_gross_j[-1] - q_true) / q_true
        err_frozen = abs(frozen_r.cumulative_gross_j[-1] - q_true) / q_true
        assert err_mass < 1.0e-4
        assert err_frozen > 20.0 * err_mass

    def test_trapped_mass_rejects_explicit_gamma(self) -> None:
        theta = np.linspace(-10.0, 10.0, 21)
        with pytest.raises(ValueError, match="explicit gamma"):
            compute_heat_release(theta, np.full(21, 1e5), np.full(21, 1e-4), gamma=1.35, trapped_mass_kg=1e-3)


@pytest.mark.strict_audit
class TestCaseFiles:
    def test_formats_units_and_signs_give_identical_metrics(self, cycles, tmp_path: Path) -> None:
        """Forte-style (MPa, cm^3, sign -1) and bench-style (bar, no volume) files agree."""
        slow, _ = cycles
        direct = compute_cycle_metrics(slow.case)
        forte = load_case(write_case_directory(slow.case, tmp_path / "forte", fuel_key="H2",
                                               wall_heat_loss_sign=-1.0))
        bench = load_case(write_case_directory(slow.case, tmp_path / "bench", trace_format="csv",
                                               fuel_key="hydrogen"))
        assert forte.case.trace.wall_heat_loss_j_per_deg is not None
        np.testing.assert_allclose(forte.case.trace.wall_heat_loss_j_per_deg,
                                   slow.case.trace.wall_heat_loss_j_per_deg, atol=1e-9)
        m_forte = compute_cycle_metrics(forte.case)
        m_bench = compute_cycle_metrics(bench.case)
        q_true = slow.total_heat_release_j
        for metrics in (direct, m_forte, m_bench):
            assert metrics["q_gross"] == pytest.approx(q_true, rel=2e-4)
        for key in ("p_max", "imep_gross", "ca50", "eta_indicated_gross"):
            assert m_forte[key] == pytest.approx(direct[key], rel=1e-6, abs=1e-6)
            assert m_bench[key] == pytest.approx(direct[key], rel=1e-4, abs=1e-3)

    def test_provenance_metadata(self, cycles, tmp_path: Path) -> None:
        slow, _ = cycles
        path = write_case_directory(slow.case, tmp_path / "c", fuel_key="hydrogen", uncertainty={"p_max": 3e4})
        loaded = load_case(tmp_path / "c")
        digest = hashlib.sha256((tmp_path / "c" / "thermo.csv").read_bytes()).hexdigest()
        assert loaded.trace_sha256 == digest == loaded.case.metadata["trace_sha256"]
        assert loaded.case_file == path.resolve()
        assert loaded.case.is_synthetic
        assert loaded.case.metadata["unit:pressure"] == "MPa"
        assert dict(loaded.uncertainties) == {"p_max": 3e4}

    @pytest.mark.parametrize(
        ("edit", "message"),
        [
            (lambda t: t.replace("bore_m", "bore_mm"), r"unknown key\(s\) \['bore_mm'\]"),
            (lambda t: t.replace('fuel = "hydrogen"', 'fuel = "unobtainium"'), "unknown fuel"),
            (lambda t: t.replace('source = "synthetic"', 'source = "lab"'), "'source' must be one of"),
            (lambda t: t.replace('format = "forte_thermo"', 'format = "xlsx"'), "'format' must be one of"),
            (lambda t: t.replace('file = "thermo.csv"', 'file = "missing.csv"'), "trace file not found"),
            (lambda t: t.replace("compression_ratio = 12.0", "compression_ratio = 0.9"), "Compression ratio"),
            (lambda t: t.replace("wall_heat_loss_sign = 1.0", "wall_heat_loss_sign = 2.0"), "must be \\+1 or -1"),
            (lambda t: t.replace("rpm = 2000.0", 'rpm = "fast"'), "must be a number"),
        ],
    )
    def test_invalid_case_files_are_rejected(self, cycles, tmp_path: Path, edit, message) -> None:
        slow, _ = cycles
        path = write_case_directory(slow.case, tmp_path / "bad", fuel_key="hydrogen")
        path.write_text(edit(path.read_text(encoding="utf-8")), encoding="utf-8")
        with pytest.raises(CaseFileError, match=message):
            load_case(path)

    def test_custom_fuel_table_and_column_map(self, tmp_path: Path) -> None:
        theta = np.arange(-180.0, 180.5, 0.5)
        pressure_bar = 1.0 + 0.0 * theta
        (tmp_path / "bench.csv").write_text(
            "angle,pcyl\n" + "\n".join(f"{a},{p}" for a, p in zip(theta, pressure_bar)) + "\n", encoding="utf-8"
        )
        (tmp_path / "case.toml").write_text(
            'name = "custom"\nsource = "test_bench"\n'
            "[geometry]\nbore_m = 0.08\nstroke_m = 0.09\nconrod_m = 0.15\ncompression_ratio = 10.0\n"
            "[operating_point]\nrpm = 1500.0\nlambda_air = 1.0\n"
            '[operating_point.fuel]\nname = "CH4 (custom)"\ncarbon_atoms = 1.0\nhydrogen_atoms = 4.0\n'
            "lhv_j_per_kg = 50.0e6\n"
            '[trace]\nfile = "bench.csv"\nformat = "csv"\n'
            '[trace.columns]\ncrank_angle = { column = "angle", unit = "deg" }\n'
            'pressure = { column = "pcyl", unit = "bar" }\n',
            encoding="utf-8",
        )
        loaded = load_case(tmp_path)
        fuel = loaded.case.operating_point.fuel
        assert fuel.name == "CH4 (custom)"
        assert fuel.stoichiometric_afr == pytest.approx(METHANE.stoichiometric_afr, rel=1e-12)
        assert loaded.case.trace.volume_m3 is None
        assert loaded.case.source == "test_bench"


@pytest.mark.strict_audit
class TestReportAndCli:
    def test_report_files_and_content(self, cycles, tmp_path: Path) -> None:
        slow, fast = cycles
        a = write_case_directory(slow.case, tmp_path / "a", fuel_key="hydrogen",
                                 uncertainty={"p_max": 3e4, "ca50": 0.3, "not_a_metric": 1.0})
        b = write_case_directory(fast.case, tmp_path / "b", trace_format="csv", fuel_key="hydrogen",
                                 uncertainty={"p_max": 3e4, "ca50": 0.3})
        report = write_comparison_report([a, b], tmp_path / "out", combustion_window_deg=(-40.0, 90.0))
        text = report.report_path.read_text(encoding="utf-8")
        assert "synthetic verification data" in text
        assert load_case(a).trace_sha256[:16] in text
        assert "| CA50 |" in text
        assert report.figure_path.is_file() and report.figure_path.stat().st_size > 10_000
        assert any("not_a_metric" in w for w in report.warnings)
        ca50 = report.comparison.row("ca50")
        assert ca50.significant[1] is True  # 0.9 deg shift vs. U = 0.42 deg
        q = report.comparison.row("q_gross")
        assert abs(q.deltas[1]) < 1e-3 * q.baseline_value  # same fuel energy, exact recovery
        import pandas as pd

        table = pd.read_csv(report.metrics_path)
        assert set(table["case"]) == {"slow", "fast"}
        assert math.isclose(
            float(table.query("metric == 'p_max' and case == 'fast'")["value_si"].iloc[0]),
            float(fast.case.trace.pressure_pa.max()), rel_tol=1e-9,
        )

    def test_cli_compare_and_inspect(self, cycles, tmp_path: Path, capsys) -> None:
        slow, fast = cycles
        write_case_directory(slow.case, tmp_path / "a", fuel_key="hydrogen")
        write_case_directory(fast.case, tmp_path / "b", fuel_key="hydrogen")
        code = main(["compare", str(tmp_path / "a"), str(tmp_path / "b"), "--out", str(tmp_path / "r"),
                     "--baseline", "fast", "--no-figure"])
        assert code == 0
        out = capsys.readouterr().out
        assert "fast (baseline)" in out
        assert (tmp_path / "r" / "report.md").is_file()
        assert not (tmp_path / "r" / "pressure_heat_release.png").exists()

        assert main(["inspect", str(tmp_path / "a" / "thermo.csv"), "--kind", "thermo"]) == 0
        out = capsys.readouterr().out
        assert "wall_heat_loss" in out and "MPa" in out

    def test_cli_reports_errors_without_traceback(self, tmp_path: Path, capsys) -> None:
        assert main(["compare", str(tmp_path / "x"), str(tmp_path / "y")]) == 1
        assert "Case file not found" in capsys.readouterr().err
