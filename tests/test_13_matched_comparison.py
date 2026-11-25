"""Matched-condition comparison: interpolation to a common CA50 and MBT timing.

Verification strategy:

- the quadratic vertex and its delta-method uncertainty are checked against exact
  parabolas and a Monte Carlo experiment (JCGM 100:2008, Sec. 5.1);
- PCHIP interpolation reproduces sweep members exactly and linear data exactly
  (Fritsch & Carlson, 1980);
- two manufactured sweeps with the same burn law but different ignition delay are
  identical at matched CA50 (to interpolation error) and clearly different at equal
  spark timing, which is the purpose of a matched comparison.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, cast

import numpy as np
import pytest

from combusta import (
    HYDROGEN,
    CycleMetrics,
    SliderCrankGeometry,
    Sweep,
    SweepError,
    build_sweep,
    compare_matched,
    find_mbt,
    manufactured_cycle,
    metrics_at,
    write_case_directory,
    write_matched_report,
)
from combusta.cli import main
from combusta.sweep import fit_quadratic_vertex

GEOM = SliderCrankGeometry(bore_m=0.086, stroke_m=0.086, conrod_m=0.145, compression_ratio=12.0)


def fake_sweep(name: str, rows: list[dict[str, float]], uncertainties: list[dict[str, float]] | None = None) -> Sweep:
    """Sweep from metric dictionaries only (heat-release results are not used by sweep functions)."""
    metrics = tuple(
        CycleMetrics(case_name=f"{name}_{i}", values=row, heat_release=cast(Any, None), burn=None)
        for i, row in enumerate(rows)
    )
    return Sweep(
        name=name,
        case_names=tuple(m.case_name for m in metrics),
        metrics=metrics,
        uncertainties=tuple(uncertainties or [{} for _ in rows]),
        contains_synthetic=True,
    )


@pytest.mark.unit
class TestQuadraticVertex:
    def test_exact_parabola(self) -> None:
        x = np.array([-32.0, -26.0, -20.0, -14.0, -8.0])
        y = 9.0e5 - 250.0 * (x + 21.5) ** 2
        vertex, u, y_vertex, (a, b, c), rms = fit_quadratic_vertex(x, y)
        assert vertex == pytest.approx(-21.5, abs=1e-9)
        assert y_vertex == pytest.approx(9.0e5, rel=1e-12)
        assert c == pytest.approx(-250.0, rel=1e-10)
        assert a + b * vertex + c * vertex**2 == pytest.approx(9.0e5, rel=1e-12)
        assert u == pytest.approx(0.0, abs=1e-6) and rms == pytest.approx(0.0, abs=1e-6)

    def test_three_points_have_no_uncertainty(self) -> None:
        vertex, u, *_ = fit_quadratic_vertex(np.array([-1.0, 0.0, 1.0]), np.array([0.0, 1.0, 0.0]))
        assert vertex == pytest.approx(0.0) and u is None

    def test_convex_data_rejected(self) -> None:
        with pytest.raises(SweepError, match="not concave"):
            fit_quadratic_vertex(np.array([-1.0, 0.0, 1.0, 2.0]), np.array([1.0, 0.0, 1.0, 4.0]))

    @pytest.mark.strict_audit
    def test_delta_method_matches_monte_carlo(self) -> None:
        """E[u^2] from the fit equals the sampling variance of the vertex (small noise)."""
        rng = np.random.default_rng(20261125)
        x = np.array([-32.0, -26.0, -20.0, -14.0, -8.0])
        truth = 10.0 - 0.01 * (x + 20.0) ** 2
        sigma = 0.02
        vertices, u2 = [], []
        for _ in range(4000):
            vertex, u, *_ = fit_quadratic_vertex(x, truth + rng.normal(0.0, sigma, x.size))
            assert u is not None
            vertices.append(vertex)
            u2.append(u**2)
        assert float(np.mean(vertices)) == pytest.approx(-20.0, abs=0.05)
        assert float(np.mean(u2)) == pytest.approx(float(np.var(vertices, ddof=1)), rel=0.15)


@pytest.mark.unit
class TestFindMbt:
    @staticmethod
    def sweep(timings: list[float], imep: list[float]) -> Sweep:
        return fake_sweep("s", [{"ignition_timing": t, "imep_gross": y} for t, y in zip(timings, imep)])

    def test_selects_points_around_maximum(self) -> None:
        t = [-40.0, -34.0, -28.0, -22.0, -16.0, -10.0, -4.0]
        y = [9.0e5 - 300.0 * (ti + 23.0) ** 2 for ti in t]  # sampled maximum at -22 deg (index 3)
        res = find_mbt(self.sweep(t, y), fit_points=5)
        assert res.theta_mbt_deg == pytest.approx(-23.0, abs=1e-9)
        assert res.fit_case_names == ("s_1", "s_2", "s_3", "s_4", "s_5")

    def test_unordered_input(self) -> None:
        t = [-10.0, -30.0, -20.0, -25.0, -15.0]
        y = [9.0e5 - 300.0 * (ti + 21.0) ** 2 for ti in t]
        assert find_mbt(self.sweep(t, y)).theta_mbt_deg == pytest.approx(-21.0, abs=1e-9)

    def test_unbracketed_maximum(self) -> None:
        with pytest.raises(SweepError, match="earliest spark timing"):
            find_mbt(self.sweep([-30.0, -20.0, -10.0, 0.0], [9.0, 8.0, 7.0, 6.0]))
        with pytest.raises(SweepError, match="latest spark timing"):
            find_mbt(self.sweep([-30.0, -20.0, -10.0, 0.0], [6.0, 7.0, 8.0, 9.0]))

    def test_missing_or_duplicate_timings(self) -> None:
        with pytest.raises(SweepError, match="ignition_timing_deg"):
            find_mbt(self.sweep([-30.0, float("nan"), -10.0], [1.0, 2.0, 1.0]))
        with pytest.raises(SweepError, match="distinct"):
            find_mbt(self.sweep([-30.0, -20.0, -20.0, -10.0], [1.0, 2.0, 2.0, 1.0]))


@pytest.mark.unit
class TestMetricsAt:
    ROWS = [
        {"ca50": 0.0, "ignition_timing": -30.0, "imep_gross": 8.0e5, "p_max": 6.0e6},
        {"ca50": 6.0, "ignition_timing": -24.0, "imep_gross": 8.6e5, "p_max": 5.4e6},
        {"ca50": 12.0, "ignition_timing": -18.0, "imep_gross": 9.2e5, "p_max": 4.8e6},
    ]

    def test_linear_data_reproduced_exactly(self) -> None:
        point = metrics_at(fake_sweep("s", self.ROWS), control="ca50", target=8.0)
        assert point.values["imep_gross"] == pytest.approx(8.8e5, rel=1e-12)
        assert point.values["ignition_timing"] == pytest.approx(-22.0, rel=1e-12)
        assert point.values["ca50"] == 8.0
        assert point.bracket == ("s_1", "s_2")
        assert point.interpolation_spread["imep_gross"] == pytest.approx(0.0, abs=1e-6)

    def test_member_value_reproduced(self) -> None:
        rows = [dict(r, eta=float(v)) for r, v in zip(self.ROWS, (0.40, 0.47, 0.41))]
        point = metrics_at(fake_sweep("s", rows), control="ca50", target=6.0)
        assert point.values["eta"] == pytest.approx(0.47, rel=1e-12)

    def test_uncertainty_interpolated_linearly(self) -> None:
        unc = [{"p_max": 1.0e4}, {"p_max": 2.0e4}, {"p_max": 4.0e4, "imep_gross": 1.0}]
        point = metrics_at(fake_sweep("s", self.ROWS, unc), control="ca50", target=9.0)
        assert point.uncertainties == {"p_max": pytest.approx(3.0e4)}

    def test_no_extrapolation(self) -> None:
        with pytest.raises(SweepError, match="outside the sweep range"):
            metrics_at(fake_sweep("s", self.ROWS), control="ca50", target=13.0)

    def test_duplicate_control_values(self) -> None:
        rows = [dict(self.ROWS[0]), dict(self.ROWS[1], ca50=0.0), dict(self.ROWS[2])]
        with pytest.raises(SweepError, match="distinct"):
            metrics_at(fake_sweep("s", rows), control="ca50", target=5.0)

    def test_non_monotonic_ca50_is_flagged(self) -> None:
        rows = [dict(self.ROWS[0]), dict(self.ROWS[1], ignition_timing=-35.0), dict(self.ROWS[2])]
        point = metrics_at(fake_sweep("s", rows), control="ca50", target=5.0)
        assert any("not monotonic" in n for n in point.notes)


@pytest.mark.unit
class TestCompareMatchedArguments:
    def test_exactly_one_mode(self) -> None:
        s = [fake_sweep("a", TestMetricsAt.ROWS), fake_sweep("b", TestMetricsAt.ROWS)]
        with pytest.raises(SweepError, match="exactly one"):
            compare_matched(s)
        with pytest.raises(SweepError, match="exactly one"):
            compare_matched(s, ca50_deg=5.0, at_mbt=True)

    def test_significance_uses_interpolated_uncertainty(self) -> None:
        unc = [{"imep_gross": 1.0e3} for _ in range(3)]
        shifted = [dict(r, imep_gross=r["imep_gross"] + 5.0e3) for r in TestMetricsAt.ROWS]
        cmp = compare_matched([fake_sweep("a", TestMetricsAt.ROWS, unc), fake_sweep("b", shifted, unc)],
                              ca50_deg=7.0)
        row = cmp.row("imep_gross")
        assert row.deltas[1] == pytest.approx(5.0e3)
        assert row.significant == (None, True)
        assert "CA50 = 7" in cmp.describe()


@pytest.fixture(scope="module")
def delay_sweeps():
    """Same Wiebe burn law, ignition delay 4 deg vs. 6 deg, spark steps of 4 deg."""
    sweeps = {}
    for name, delay in (("short_delay", 4.0), ("long_delay", 6.0)):
        cases = [
            manufactured_cycle(
                f"{name}_{s:+.0f}", GEOM, HYDROGEN, lambda_air=2.0, rpm=2000.0,
                combustion_start_deg=s + delay, combustion_duration_deg=36.0, wall_heat_loss_j=60.0,
                ignition_timing_deg=s,
            ).case
            for s in (-34.0, -30.0, -26.0, -22.0, -18.0, -14.0)
        ]
        sweeps[name] = cases
    return sweeps


@pytest.mark.strict_audit
class TestManufacturedSweeps:
    def test_matched_ca50_removes_phasing_difference(self, delay_sweeps) -> None:
        a = build_sweep("short_delay", delay_sweeps["short_delay"])
        b = build_sweep("long_delay", delay_sweeps["long_delay"])
        # Same burn law: at matched CA50 the cycles are the same up to interpolation error.
        matched = compare_matched([a, b], ca50_deg=3.0)
        for key, rel in (("imep_gross", 1e-3), ("eta_indicated_gross", 1e-3), ("p_max", 3e-3)):
            row = matched.row(key)
            assert abs(row.deltas[1]) < rel * abs(row.baseline_value), key
        assert matched.row("ignition_timing").deltas[1] == pytest.approx(-2.0, abs=0.05)
        # At equal spark timing the later ignition shows up as a large difference.
        same_spark = metrics_at(b, control="ignition_timing", target=-22.0)
        ref = metrics_at(a, control="ignition_timing", target=-22.0)
        assert abs(same_spark.values["p_max"] - ref.values["p_max"]) > 10.0 * abs(matched.row("p_max").deltas[1])

    def test_mbt_properties(self, delay_sweeps) -> None:
        a = build_sweep("short_delay", delay_sweeps["short_delay"])
        b = build_sweep("long_delay", delay_sweeps["long_delay"])
        res_a, res_b = find_mbt(a), find_mbt(b)
        # Shifting ignition by +2 deg shifts MBT spark timing by -2 deg (same CA50 at MBT).
        assert res_b.theta_mbt_deg - res_a.theta_mbt_deg == pytest.approx(-2.0, abs=0.3)
        assert res_a.metric_at_mbt >= float(np.max(a.values("imep_gross"))) - 1.0
        mbt = compare_matched([a, b], at_mbt=True)
        assert mbt.row("ignition_timing").values[0] == pytest.approx(res_a.theta_mbt_deg)
        assert abs(mbt.row("imep_gross").deltas[1]) < 2e-3 * mbt.row("imep_gross").baseline_value

    def test_report_and_cli(self, delay_sweeps, tmp_path: Path, capsys) -> None:
        dirs = []
        for name, cases in delay_sweeps.items():
            for case in cases:
                dirs.append(write_case_directory(case, tmp_path / name / case.name, fuel_key="hydrogen",
                                                 configuration=name, uncertainty={"imep_gross": 2.0e3}).parent)
        report = write_matched_report(dirs, tmp_path / "out", at_mbt=True)
        text = report.report_path.read_text(encoding="utf-8")
        assert "MBT timing [deg ATDC]" in text and "synthetic verification data" in text
        assert report.figure_path.is_file()
        assert "interpolation_spread_si" in report.metrics_path.read_text(encoding="utf-8")

        code = main(["matched", *map(str, dirs), "--ca50", "3", "--out", str(tmp_path / "cli"), "--no-figure"])
        assert code == 0
        assert "CA50 = 3" in capsys.readouterr().out

    def test_cases_need_configuration(self, delay_sweeps, tmp_path: Path, capsys) -> None:
        dirs = [write_case_directory(c, tmp_path / c.name, fuel_key="hydrogen").parent
                for c in delay_sweeps["short_delay"][:4]]
        assert main(["matched", *map(str, dirs), "--mbt", "--out", str(tmp_path / "x")]) == 1
        assert "without 'configuration'" in capsys.readouterr().err
