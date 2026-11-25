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
    plot_pre_chamber,
    plot_sweep_overview,
)
from combusta.measurement import (
    FloatArray,
    MeasurementUncertaintyResult,
    NoiseFilter,
    PressureMeasurementUncertainty,
    pressure_measurement_uncertainty,
)
from combusta.prechamber import (
    PreChamberMetrics,
    PreChamberResult,
    analyze_pre_chamber,
    estimate_discharge_coefficient,
    pre_chamber_metrics,
)
from combusta.signal.filtering import savgol_zero_phase
from combusta.sweep import MatchedComparison, Sweep, build_sweep, compare_matched
from combusta.thermo.fuels import make_reactant_mixture

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


PRE_CHAMBER_REPORT_FILE: Final[str] = "prechamber_report.md"
PRE_CHAMBER_METRICS_FILE: Final[str] = "prechamber_metrics.csv"
PRE_CHAMBER_FIGURE_FILE: Final[str] = "prechamber.png"

_PRE_CHAMBER_DISPLAY: Final[tuple[tuple[str, str, str, float], ...]] = (
    ("dp_max_pa", "Peak pressure difference p_pc - p_mc", "bar", 1e-5),
    ("theta_dp_max_deg", "Crank angle of peak pressure difference", "deg ATDC", 1.0),
    ("dp_min_pa", "Most negative pressure difference (flow into the pre-chamber)", "bar", 1e-5),
    ("theta_dp_min_deg", "Crank angle of most negative pressure difference", "deg ATDC", 1.0),
    ("jet_start_deg", "Jet start (Δp > 0)", "deg ATDC", 1.0),
    ("jet_end_deg", "Jet end", "deg ATDC", 1.0),
    ("jet_mass_kg", "Mass leaving the pre-chamber during the jet", "mg", 1e6),
    ("jet_mass_fraction", "Jet mass / pre-chamber mass at jet start", "%", 100.0),
    ("peak_jet_rate_kg_per_s", "Peak jet mass flow", "g/s", 1e3),
    ("choked_fraction", "Share of the jet interval with choked flow", "%", 100.0),
    ("q_pc_gross_j", "Pre-chamber gross heat release", "J", 1.0),
    ("q_mc_gross_j", "Main-chamber gross heat release", "J", 1.0),
    ("pc_energy_share", "Pre-chamber share of the heat release", "%", 100.0),
    ("pc_ca50_deg", "Pre-chamber CA50", "deg ATDC", 1.0),
    ("mc_ca10_deg", "Main-chamber CA10", "deg ATDC", 1.0),
    ("mc_ca50_deg", "Main-chamber CA50", "deg ATDC", 1.0),
    ("mc_ca90_deg", "Main-chamber CA90", "deg ATDC", 1.0),
)


def analyze_case_pre_chamber(loaded: LoadedCase) -> PreChamberResult:
    """Two-volume analysis of a case file with a ``[pre_chamber]`` table.

    The case trace is the main chamber. Zone masses come from the two temperature
    traces when both exist, otherwise from the orifice model with the case's
    ``trapped_mass_kg`` (total mass of both volumes). Both volumes use the unburned
    charge of the case's fuel and λ (frozen composition).
    """
    pc = loaded.pre_chamber
    if pc is None:
        raise ValueError(f"{loaded.case.name}: the case file has no [pre_chamber] table.")
    case, op = loaded.case, loaded.case.operating_point
    if op.lambda_air is None:
        raise ValueError(f"{case.name}: lambda_air is required for the pre-chamber gas model.")
    mixture = make_reactant_mixture(op.fuel, op.lambda_air)
    kinematics = case.kinematics()
    temperatures = pc.temperature_k is not None and case.trace.temperature_k is not None
    return analyze_pre_chamber(
        case.trace.crank_angle_deg,
        pc.pressure_pa,
        case.trace.pressure_pa,
        case.volume_m3(),
        pc.geometry,
        rpm=op.rpm,
        mixture_pc=mixture,
        mixture_mc=mixture,
        pre_chamber_temperature_k=pc.temperature_k if temperatures else None,
        main_chamber_temperature_k=case.trace.temperature_k if temperatures else None,
        trapped_mass_kg=None if temperatures else op.trapped_mass_kg,
        wall_heat_loss_pc_j_per_deg=pc.wall_heat_loss_j_per_deg,
        wall_heat_loss_mc_j_per_deg=case.trace.wall_heat_loss_j_per_deg,
        main_chamber_dv_dtheta_m3_per_deg=(
            kinematics.dv_dtheta_m3_per_deg if case.trace.volume_m3 is None else None
        ),
    )


@dataclass(frozen=True, slots=True)
class PreChamberReport:
    """Result of :func:`write_pre_chamber_report`.

    Attributes
    ----------
    case : LoadedCase
        The analyzed case.
    result : PreChamberResult
        Crank-angle-resolved two-volume analysis.
    metrics : PreChamberMetrics
        Scalar metrics.
    discharge_coefficient_estimate : float | None
        $C_D$ fitted to the traced mass transfer (``None`` with the orifice model).
    markdown : str
        Report text.
    report_path, metrics_path, figure_path : Path
        Written files.
    """

    case: LoadedCase
    result: PreChamberResult
    metrics: PreChamberMetrics
    discharge_coefficient_estimate: float | None
    markdown: str
    report_path: Path
    metrics_path: Path
    figure_path: Path


def write_pre_chamber_report(
    case_path: str | Path,
    out_dir: str | Path,
    *,
    window_deg: tuple[float, float] = (-40.0, 90.0),
    figure: bool = True,
) -> PreChamberReport:
    """Analyze a pre-chamber case file and write ``prechamber_report.md``, metrics CSV, and figure."""
    loaded = load_case(case_path)
    result = analyze_case_pre_chamber(loaded)
    rpm = loaded.case.operating_point.rpm
    metrics = pre_chamber_metrics(result, rpm=rpm, window_deg=window_deg)
    cd_fit = estimate_discharge_coefficient(result) if result.mass_source == "trace" else None
    pc = loaded.pre_chamber
    assert pc is not None

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    figure_path = out / PRE_CHAMBER_FIGURE_FILE
    if figure:
        plot_pre_chamber(
            result.crank_angle_deg, pc.pressure_pa, loaded.case.trace.pressure_pa,
            result.transfer_rate_kg_per_deg * 6.0 * rpm,
            result.dq_pc_gross_j_per_deg, result.dq_mc_gross_j_per_deg,
            jet_window_deg=(metrics.jet_start_deg, metrics.jet_end_deg),
            theta_limits_deg=window_deg, title=loaded.case.name, save_path=figure_path,
        )

    from combusta import __version__

    values = metrics.as_dict()
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    geo = pc.geometry
    lines = ["# Pre-chamber analysis", "", f"Generated by Combusta {__version__} on {stamp}.", ""]
    if loaded.case.is_synthetic:
        lines += ["> [!WARNING]",
                  "> This case contains **synthetic verification data**. The numbers are not engineering results.", ""]
    lines += ["## Case", "", *_provenance_table([loaded]), "",
              f"- Pre-chamber signals: `{pc.file.name}` (SHA-256 `{pc.sha256[:16]}`)",
              f"- Pre-chamber volume {geo.volume_m3 * 1e6:.4g} cm³; {geo.orifice_count} orifices of "
              f"{geo.orifice_diameter_m * 1e3:.4g} mm (total area {geo.flow_area_m2 * 1e6:.4g} mm²); "
              f"C_D = {geo.discharge_coefficient:g}",
              f"- Mass exchange from: **{'zone temperature traces' if result.mass_source == 'trace' else 'orifice model with the trapped mass'}**",
              ""]
    lines += ["## Results", "", "| Metric | Unit | Value |", "|---|---|---|"]
    for key, label, unit, factor in _PRE_CHAMBER_DISPLAY:
        v = values[key]
        lines.append(f"| {label} | {unit} | {'–' if not math.isfinite(v) else f'{v * factor:.4g}'} |")
    lines.append("")
    if cd_fit is not None:
        lines += [f"Effective discharge coefficient fitted to the traced mass transfer: **{cd_fit:.3f}** "
                  f"(case file: {geo.discharge_coefficient:g}). A large difference means the orifice model "
                  "with the declared C_D does not reproduce the exported mass exchange.", ""]
    lines += ["## Assumptions", "",
              "- The case trace describes the main chamber only; the pre-chamber volume is not part of it.",
              "- Each volume is one zone of uniform state with the frozen unburned-charge gas model; gas "
              "entering a volume takes on that volume's gas model. With hydrogen this biases the absolute "
              "heat release; timing quantities are less sensitive.",
              "- Orifice flow is quasi-steady and isentropic with a constant discharge coefficient.", ""]
    if figure:
        lines += ["## Figure", "", f"![Pre-chamber analysis]({PRE_CHAMBER_FIGURE_FILE})", ""]
    markdown = "\n".join(lines)
    report_path = out / PRE_CHAMBER_REPORT_FILE
    report_path.write_text(markdown + "\n", encoding="utf-8")
    metrics_path = out / PRE_CHAMBER_METRICS_FILE
    import pandas as pd

    pd.DataFrame([{"metric": k, "label": label, "value_si": values[k], "display_unit": unit,
                   "value_display": values[k] * f} for k, label, unit, f in _PRE_CHAMBER_DISPLAY]
                 ).to_csv(metrics_path, index=False)
    return PreChamberReport(
        case=loaded, result=result, metrics=metrics, discharge_coefficient_estimate=cd_fit,
        markdown=markdown, report_path=report_path, metrics_path=metrics_path, figure_path=figure_path,
    )


MEASUREMENT_REPORT_FILE: Final[str] = "measurement_uncertainty.md"
MEASUREMENT_METRICS_FILE: Final[str] = "measurement_uncertainty.csv"


@dataclass(frozen=True, slots=True)
class MeasurementReport:
    """Result of :func:`write_measurement_report`.

    Attributes
    ----------
    case : LoadedCase
        The analyzed case.
    result : MeasurementUncertaintyResult
        Uncertainty budget.
    markdown : str
        Report text.
    report_path, metrics_path : Path
        Written files.
    """

    case: LoadedCase
    result: MeasurementUncertaintyResult
    markdown: str
    report_path: Path
    metrics_path: Path


def write_measurement_report(
    case_path: str | Path,
    out_dir: str | Path,
    *,
    inputs: PressureMeasurementUncertainty | None = None,
    n_samples: int = 500,
    seed: int | None = 0,
    savgol_window: int | None = None,
    combustion_window_deg: tuple[float, float] = (-40.0, 90.0),
    gross_window_deg: tuple[float, float] = (-180.0, 180.0),
    revolutions_per_cycle: int = 2,
) -> MeasurementReport:
    """Propagate pressure-measurement errors for a case file and write the uncertainty budget.

    ``inputs`` overrides the case file's ``[measurement]`` table. ``savgol_window``
    applies a zero-phase Savitzky-Golay filter (cubic) of that odd window length to
    every trace, as in the processing of real data.
    """
    loaded = load_case(case_path)
    meas = inputs or loaded.measurement
    if meas is None:
        raise ValueError(f"{loaded.case.name}: no measurement uncertainties given ([measurement] or options).")
    noise_filter: NoiseFilter | None = None
    if savgol_window is not None:
        window = savgol_window

        def noise_filter(theta: FloatArray, p: FloatArray) -> FloatArray:  # noqa: ARG001
            return savgol_zero_phase(p, window, 3)

    result = pressure_measurement_uncertainty(
        loaded.case, meas, n_samples=n_samples, seed=seed, noise_filter=noise_filter,
        combustion_window_deg=combustion_window_deg, gross_window_deg=gross_window_deg,
        revolutions_per_cycle=revolutions_per_cycle,
    )

    from combusta import __version__

    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    lines = ["# Pressure-measurement uncertainty", "", f"Generated by Combusta {__version__} on {stamp}.", ""]
    if loaded.case.is_synthetic:
        lines += ["> [!WARNING]",
                  "> This case contains **synthetic verification data**. The numbers are not engineering results.", ""]
    lines += ["## Case", "", *_provenance_table([loaded]), ""]
    lines += ["## Inputs (standard uncertainties, normal, independent)", "",
              f"- TDC phasing: {meas.tdc_offset_deg:g} deg",
              f"- Pegging offset: {meas.pegging_offset_pa * 1e-5:g} bar",
              f"- Transducer gain: {meas.gain_relative * 100:g} %",
              f"- Sample noise: {meas.noise_pa * 1e-5:g} bar"
              + (f" (filtered with a {savgol_window}-point zero-phase Savitzky-Golay filter)" if savgol_window else ""),
              ""]
    lines += ["## Uncertainty budget", "", result.to_markdown(), ""]
    lines += ["Contribution columns are |c_i u_i| from central differences (linear propagation, JCGM 100:2008). "
              "The Monte Carlo column combines all sources including noise; its ± value is the standard error "
              "of the standard deviation. The number of trials is far below what JCGM 101:2008 recommends for "
              "coverage intervals, so none are reported.", ""]
    if result.notes:
        lines += ["## Notes", ""] + [f"- {n}" for n in result.notes] + [""]
    lines += ["## Not included", "",
              "- Thermal shock of the transducer (an intra-cycle, load-dependent error).",
              "- Crank-angle encoder resolution and errors of the declared geometry or volume.",
              "- Cycle-to-cycle variation: combine with the Type A uncertainty of the cycle mean "
              "(`combusta.vv.experimental_uncertainty`).", ""]
    markdown = "\n".join(lines)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    report_path = out / MEASUREMENT_REPORT_FILE
    report_path.write_text(markdown + "\n", encoding="utf-8")
    metrics_path = out / MEASUREMENT_METRICS_FILE
    result.to_dataframe().to_csv(metrics_path, index=False)
    return MeasurementReport(case=loaded, result=result, markdown=markdown, report_path=report_path,
                             metrics_path=metrics_path)
