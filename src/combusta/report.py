"""Comparison reports from case files.

:func:`write_comparison_report` loads several ``case.toml`` files, compares them with
:func:`combusta.compare.compare_cases`, and writes

- ``report.md``: provenance (Combusta version, case files, trace SHA-256 digests),
  data checks, the comparison table, analysis notes, and the figure,
- ``metrics.csv``: one row per metric and case in SI and display units,
- ``pressure_heat_release.png``: pressure and gross heat-release overlay.

The report states what it does not know: missing uncertainties are reported as such,
and synthetic cases carry a banner.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from datetime import datetime, timezone
import math
from pathlib import Path
from typing import Final

from combusta.casefile import LoadedCase, load_case
from combusta.compare import CaseComparison, compare_cases
from combusta.plotting.publication_plots import (
    CaseOverlayTrace,
    SweepPlotSeries,
    plot_case_overlay,
    plot_sweep_overview,
)
from combusta.sweep import MatchedComparison, Sweep, build_sweep, compare_matched

VOLUME_CONSISTENCY_WARNING: Final[float] = 0.01
"""Relative traced-vs-geometric volume deviation above which the report warns [-]."""

OPERATING_POINT_SPREAD_WARNING: Final[float] = 0.01
"""Relative spread of speed, lambda, or fuel mass above which a matched report warns [-]."""

REPORT_FILE: Final[str] = "report.md"
METRICS_FILE: Final[str] = "metrics.csv"
FIGURE_FILE: Final[str] = "pressure_heat_release.png"
MATCHED_REPORT_FILE: Final[str] = "matched_report.md"
MATCHED_METRICS_FILE: Final[str] = "matched_metrics.csv"
SWEEP_FIGURE_FILE: Final[str] = "sweep_overview.png"


@dataclass(frozen=True, slots=True)
class ComparisonReport:
    """Result of :func:`write_comparison_report`.

    Attributes
    ----------
    comparison : CaseComparison
        The comparison that was reported.
    cases : tuple[LoadedCase, ...]
        Loaded cases in table order.
    markdown : str
        Report text.
    report_path, metrics_path, figure_path : Path
        Written files.
    warnings : tuple[str, ...]
        Data checks that failed or could not be made.
    """

    comparison: CaseComparison
    cases: tuple[LoadedCase, ...]
    markdown: str
    report_path: Path
    metrics_path: Path
    figure_path: Path
    warnings: tuple[str, ...]


def _display_path(path: Path) -> str:
    """Path relative to the working directory when possible (keeps reports portable)."""
    try:
        return path.relative_to(Path.cwd()).as_posix()
    except ValueError:
        return str(path)


def _data_checks(loaded: Sequence[LoadedCase], reported_keys: set[str]) -> list[str]:
    warnings: list[str] = []
    for item in loaded:
        case = item.case
        if case.trace.volume_m3 is not None:
            err = case.volume_consistency_error()
            if err > VOLUME_CONSISTENCY_WARNING:
                warnings.append(
                    f"{case.name}: traced volume deviates from the declared geometry by up to "
                    f"{100 * err:.2f} % (check compression ratio, connecting rod, crank-angle offset)."
                )
        if case.operating_point.fuel_mass_per_cycle_kg is None:
            warnings.append(f"{case.name}: no fuel mass per cycle; efficiency and ISFC are not computed.")
        if case.operating_point.lambda_air is None:
            warnings.append(f"{case.name}: no lambda_air; dry-air properties with frozen composition are used.")
    for item in loaded:
        unknown = sorted(set(item.uncertainties) - reported_keys)
        if unknown:
            warnings.append(
                f"{item.case.name}: [uncertainty] keys {unknown} match no reported metric and were ignored."
            )
    if not any(item.uncertainties for item in loaded):
        warnings.append(
            "No numerical uncertainties given: differences cannot be tested for significance."
        )
    return warnings


def _provenance_table(loaded: Sequence[LoadedCase]) -> list[str]:
    lines = [
        "| Case | Source | Fuel | λ | Speed [rpm] | Trace file | SHA-256 |",
        "|---|---|---|---|---|---|---|",
    ]
    for item in loaded:
        case, op = item.case, item.case.operating_point
        lam = "–" if op.lambda_air is None else f"{op.lambda_air:g}"
        lines.append(
            f"| {case.name} | {case.source} | {op.fuel.name} | {lam} | {op.rpm:g} | "
            f"`{item.trace_file.name}` ({case.metadata.get('trace_format', '')}) | `{item.trace_sha256[:16]}` |"
        )
    return lines


def render_report(
    loaded: Sequence[LoadedCase],
    comparison: CaseComparison,
    *,
    warnings: Sequence[str],
    figure_name: str | None,
    analysis_settings: Mapping[str, object],
    generated_utc: datetime | None = None,
) -> str:
    """Render the Markdown text of a comparison report."""
    from combusta import __version__

    stamp = (generated_utc or datetime.now(timezone.utc)).strftime("%Y-%m-%d %H:%M UTC")
    lines = [
        "# Engine case comparison",
        "",
        f"Generated by Combusta {__version__} on {stamp}.",
        "",
    ]
    if comparison.contains_synthetic:
        lines += [
            "> [!WARNING]",
            "> This comparison contains **synthetic verification data**. The numbers are not engineering results.",
            "",
        ]
    lines += ["## Cases", ""] + _provenance_table(loaded) + [""]
    lines += ["Case files:", ""] + [f"- `{item.case.name}`: `{_display_path(item.case_file)}`" for item in loaded] + [""]

    lines += ["## Data checks", ""]
    lines += [f"- {w}" for w in warnings] if warnings else ["- All checks passed."]
    lines += [""]

    lines += ["## Results", "", comparison.to_markdown(), ""]

    lines += ["## Analysis settings", ""]
    lines += [f"- `{key}` = `{value}`" for key, value in analysis_settings.items()]
    lines += [""]
    lines += ["## Analysis notes", ""]
    for metrics in comparison.metrics:
        lines.append(f"**{metrics.case_name}**")
        lines.append("")
        lines += [f"- {note}" for note in metrics.notes] or ["- none"]
        lines.append("")

    if figure_name:
        lines += ["## Figure", "", f"![Pressure and gross heat-release rate]({figure_name})", ""]
    return "\n".join(lines)


def write_comparison_report(
    case_paths: Sequence[str | Path],
    out_dir: str | Path,
    *,
    baseline: int | str = 0,
    figure: bool = True,
    combustion_window_deg: tuple[float, float] = (-40.0, 90.0),
    gross_window_deg: tuple[float, float] = (-180.0, 180.0),
    revolutions_per_cycle: int = 2,
    variable_composition: bool = True,
) -> ComparisonReport:
    """Load case files, compare them, and write the report files to ``out_dir``.

    Parameters
    ----------
    case_paths : Sequence[str | Path]
        Two or more ``case.toml`` files or directories containing one.
    out_dir : str | Path
        Output directory (created if needed).
    baseline : int | str, default=0
        Baseline case index or name.
    figure : bool, default=True
        Write the pressure / heat-release overlay figure.
    combustion_window_deg, gross_window_deg, revolutions_per_cycle, variable_composition
        Analysis settings forwarded to :func:`combusta.compare.compute_cycle_metrics`.
    """
    settings: dict[str, object] = {
        "combustion_window_deg": combustion_window_deg,
        "gross_window_deg": gross_window_deg,
        "revolutions_per_cycle": revolutions_per_cycle,
        "variable_composition": variable_composition,
    }
    loaded = tuple(load_case(p) for p in case_paths)
    comparison = compare_cases(
        [item.case for item in loaded],
        baseline=baseline,
        uncertainties={item.case.name: dict(item.uncertainties) for item in loaded if item.uncertainties},
        combustion_window_deg=combustion_window_deg,
        gross_window_deg=gross_window_deg,
        revolutions_per_cycle=revolutions_per_cycle,
        variable_composition=variable_composition,
    )
    warnings = _data_checks(loaded, {row.key for row in comparison.rows})

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    figure_path = out / FIGURE_FILE
    if figure:
        overlay = [
            CaseOverlayTrace(
                label=m.case_name,
                theta_deg=item.case.trace.crank_angle_deg,
                pressure_pa=item.case.trace.pressure_pa,
                hrr_result=m.heat_release,
                ca50_deg=float(m.values.get("ca50", math.nan)),
            )
            for item, m in zip(loaded, comparison.metrics)
        ]
        plot_case_overlay(overlay, theta_limits_deg=combustion_window_deg, save_path=figure_path)

    settings = {"baseline": comparison.case_names[comparison.baseline_index], **settings}
    markdown = render_report(
        loaded,
        comparison,
        warnings=warnings,
        figure_name=FIGURE_FILE if figure else None,
        analysis_settings=settings,
    )
    report_path = out / REPORT_FILE
    report_path.write_text(markdown + "\n", encoding="utf-8")
    metrics_path = out / METRICS_FILE
    comparison.to_dataframe().to_csv(metrics_path, index=False)

    return ComparisonReport(
        comparison=comparison,
        cases=loaded,
        markdown=markdown,
        report_path=report_path,
        metrics_path=metrics_path,
        figure_path=figure_path,
        warnings=tuple(warnings),
    )


@dataclass(frozen=True, slots=True)
class MatchedReport:
    """Result of :func:`write_matched_report`.

    Attributes
    ----------
    comparison : MatchedComparison
        The matched comparison.
    sweeps : tuple[Sweep, ...]
        Sweeps per configuration.
    cases : tuple[LoadedCase, ...]
        All loaded cases.
    markdown : str
        Report text.
    report_path, metrics_path, figure_path : Path
        Written files.
    warnings : tuple[str, ...]
        Data checks that failed.
    """

    comparison: MatchedComparison
    sweeps: tuple[Sweep, ...]
    cases: tuple[LoadedCase, ...]
    markdown: str
    report_path: Path
    metrics_path: Path
    figure_path: Path
    warnings: tuple[str, ...]


def _operating_point_checks(loaded: Sequence[LoadedCase]) -> list[str]:
    warnings = []
    for label, getter in (
        ("speed", lambda op: op.rpm),
        ("lambda", lambda op: op.lambda_air),
        ("fuel mass per cycle", lambda op: op.fuel_mass_per_cycle_kg),
    ):
        vals = [getter(item.case.operating_point) for item in loaded]
        if any(v is None for v in vals):
            continue
        lo, hi = min(vals), max(vals)
        if hi > 0.0 and (hi - lo) / hi > OPERATING_POINT_SPREAD_WARNING:
            warnings.append(
                f"Cases differ in {label} ({lo:.4g} to {hi:.4g}); a matched comparison assumes equal "
                "speed, lambda, and fuel mass, so part of the differences is not due to the configuration."
            )
    return warnings


def write_matched_report(
    case_paths: Sequence[str | Path],
    out_dir: str | Path,
    *,
    ca50_deg: float | None = None,
    at_mbt: bool = False,
    baseline: int | str = 0,
    fit_points: int = 5,
    figure: bool = True,
    combustion_window_deg: tuple[float, float] = (-40.0, 90.0),
    gross_window_deg: tuple[float, float] = (-180.0, 180.0),
    revolutions_per_cycle: int = 2,
) -> MatchedReport:
    """Group case files by ``configuration``, compare at matched CA50 or MBT, write a report.

    Every case file must set ``configuration``; each configuration needs at least two
    cases (a spark-timing sweep). Writes ``matched_report.md``, ``matched_metrics.csv``,
    and ``sweep_overview.png`` to ``out_dir``.
    """
    loaded = tuple(load_case(p) for p in case_paths)
    missing = [item.case.name for item in loaded if item.configuration is None]
    if missing:
        raise ValueError(f"Cases without 'configuration' in their case file: {missing}.")
    groups: dict[str, list[LoadedCase]] = {}
    for item in loaded:
        assert item.configuration is not None
        groups.setdefault(item.configuration, []).append(item)
    kwargs: dict[str, object] = {
        "combustion_window_deg": combustion_window_deg,
        "gross_window_deg": gross_window_deg,
        "revolutions_per_cycle": revolutions_per_cycle,
    }
    sweeps = tuple(
        build_sweep(
            name,
            [item.case for item in items],
            uncertainties={item.case.name: dict(item.uncertainties) for item in items if item.uncertainties},
            **kwargs,
        )
        for name, items in groups.items()
    )
    comparison = compare_matched(sweeps, ca50_deg=ca50_deg, at_mbt=at_mbt, baseline=baseline,
                                 fit_points=fit_points)

    reported = {row.key for row in comparison.rows}
    warnings = _data_checks(loaded, reported) + _operating_point_checks(loaded)
    warnings += [f"{p.sweep_name}: {note}" for p in comparison.points for note in p.notes]

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    figure_path = out / SWEEP_FIGURE_FILE
    if figure:
        series = []
        for sweep, mbt in zip(sweeps, comparison.mbt):
            fit_range = None
            if mbt is not None:
                fitted = [sweep.values("ignition_timing")[sweep.case_names.index(n)] for n in mbt.fit_case_names]
                fit_range = (float(min(fitted)), float(max(fitted)))
            series.append(SweepPlotSeries(
                label=sweep.name,
                ignition_timing_deg=sweep.values("ignition_timing"),
                imep_pa=sweep.values("imep_gross"),
                ca50_deg=sweep.values("ca50"),
                efficiency=sweep.values("eta_indicated_gross"),
                mbt_deg=math.nan if mbt is None else mbt.theta_mbt_deg,
                mbt_coefficients=None if mbt is None else mbt.coefficients,
                fit_range_deg=fit_range,
            ))
        plot_sweep_overview(series, ca50_target_deg=ca50_deg, save_path=figure_path)

    from combusta import __version__

    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    lines = ["# Matched engine configuration comparison", "",
             f"Generated by Combusta {__version__} on {stamp}.", ""]
    if comparison.contains_synthetic:
        lines += ["> [!WARNING]",
                  "> This comparison contains **synthetic verification data**. The numbers are not engineering results.",
                  ""]
    lines += ["## Sweeps", ""]
    for sweep in sweeps:
        timings = ", ".join(f"{t:g}" for t in sorted(sweep.values("ignition_timing")) if math.isfinite(t))
        lines.append(f"- **{sweep.name}**: {len(sweep.case_names)} cases"
                     + (f"; spark timings {timings} deg ATDC" if timings else ""))
    lines += ["", "## Cases", ""] + _provenance_table(loaded) + [""]
    lines += ["## Data checks", ""] + ([f"- {w}" for w in warnings] or ["- All checks passed."]) + [""]
    lines += ["## Results", "", comparison.to_markdown(), ""]
    lines += ["Values are interpolated with a monotone cubic (PCHIP) between the bracketing cases; "
              "numerical uncertainties are interpolated linearly between them. "
              + ("The MBT timing uncertainty is not propagated into the metric uncertainties. " if at_mbt else "")
              + "`interpolation_spread_si` in the CSV file shows how much linear and cubic interpolation "
              "differ; large values mean the sweep is too coarse.", ""]
    lines += ["## Analysis settings", ""]
    settings = {"mode": comparison.mode, "ca50_deg": ca50_deg, "fit_points": fit_points,
                "baseline": comparison.case_names[comparison.baseline_index], **kwargs}
    lines += [f"- `{k}` = `{v}`" for k, v in settings.items()] + [""]
    if figure:
        lines += ["## Figure", "", f"![Sweep overview]({SWEEP_FIGURE_FILE})", ""]
    markdown = "\n".join(lines)
    report_path = out / MATCHED_REPORT_FILE
    report_path.write_text(markdown + "\n", encoding="utf-8")
    metrics_path = out / MATCHED_METRICS_FILE
    comparison.to_dataframe().to_csv(metrics_path, index=False)
    return MatchedReport(
        comparison=comparison, sweeps=sweeps, cases=loaded, markdown=markdown,
        report_path=report_path, metrics_path=metrics_path, figure_path=figure_path,
        warnings=tuple(warnings),
    )
