"""Cycle metrics for one engine case and side-by-side comparison of several cases.

:func:`compute_cycle_metrics` reduces an :class:`~combusta.case.EngineCase` to the
standard indicated-cycle quantities (Heywood, 2018, Ch. 2 and 9): peak pressure and
its phasing, maximum pressure-rise rate, gross/net IMEP, gross indicated power,
efficiency and specific fuel consumption, CA10/CA50/CA90, and cumulative heat release.

:func:`compare_cases` tabulates these metrics for several cases against a baseline.
When numerical uncertainties are supplied (for example from a three-grid GCI study,
see :func:`uncertainty_from_gci`), each difference is tested against the combined
uncertainty of the two cases,

$$U_\\Delta = \\sqrt{U_a^2 + U_b^2},$$

and flagged significant when $|\\Delta| > U_\\Delta$. Combining in quadrature assumes
independent numerical errors; when both designs share the same meshing strategy their
errors are usually positively correlated, which makes this test conservative.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
import math
from typing import Final, TypeAlias

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import trapezoid

from combusta.case import EngineCase
from combusta.emissions import (
    combustion_efficiency,
    emission_index,
    has_heating_value,
    specific_emission,
)
from combusta.exceptions import PhysicalBoundsError
from combusta.signal.filtering import compute_derivative
from combusta.thermo.fuels import make_product_mixture, make_reactant_mixture
from combusta.thermo.heat_release import (
    BurnFractionMetrics,
    HeatReleaseResult,
    compute_heat_release,
    compute_imep,
    extract_burn_fractions,
)
from combusta.thermo.nasa_polynomials import DRY_AIR, NASAGasMixture
from combusta.vv.gci import GCIResult

FloatArray: TypeAlias = NDArray[np.float64]

J_PER_KWH: Final[float] = 3.6e6


@dataclass(frozen=True, slots=True)
class MetricSpec:
    """Display metadata of a cycle metric; values are always stored in SI units."""

    label: str
    display_unit: str
    si_to_display: float


METRIC_SPECS: Final[Mapping[str, MetricSpec]] = {
    "p_max": MetricSpec("Peak cylinder pressure", "MPa", 1.0e-6),
    "theta_p_max": MetricSpec("Crank angle of peak pressure", "deg ATDC", 1.0),
    "max_prr": MetricSpec("Maximum pressure-rise rate", "bar/deg", 1.0e-5),
    "imep_gross": MetricSpec("Gross IMEP", "bar", 1.0e-5),
    "imep_net": MetricSpec("Net IMEP (full trace)", "bar", 1.0e-5),
    "indicated_power_gross": MetricSpec("Gross indicated power per cylinder", "kW", 1.0e-3),
    "eta_indicated_gross": MetricSpec("Gross indicated fuel conversion efficiency", "%", 100.0),
    "isfc_gross": MetricSpec("Gross indicated specific fuel consumption", "g/kWh", 1.0e3 * J_PER_KWH),
    "ca10": MetricSpec("CA10", "deg ATDC", 1.0),
    "ca50": MetricSpec("CA50", "deg ATDC", 1.0),
    "ca90": MetricSpec("CA90", "deg ATDC", 1.0),
    "burn_duration_10_90": MetricSpec("Combustion duration CA10-CA90", "deg", 1.0),
    "q_gross": MetricSpec("Cumulative gross heat release", "J", 1.0),
    "heat_release_fraction": MetricSpec("Gross heat release / fuel energy", "%", 100.0),
    "wall_heat_loss": MetricSpec("Cumulative wall heat loss", "J", 1.0),
    "combustion_efficiency": MetricSpec("Combustion efficiency (exhaust chemical energy)", "%", 100.0),
}

SPECIFIC_EMISSION_PREFIX: Final[str] = "is_"
EMISSION_INDEX_PREFIX: Final[str] = "ei_"


@dataclass(frozen=True, slots=True)
class CycleMetrics:
    """Indicated-cycle metrics of one case (SI units; ``nan`` when not computable).

    Attributes
    ----------
    case_name : str
        Name of the analyzed case.
    values : Mapping[str, float]
        Metric values keyed as in ``METRIC_SPECS`` plus any case scalars.
    heat_release : HeatReleaseResult
        Underlying heat-release analysis.
    burn : BurnFractionMetrics | None
        Burn-fraction analysis, or ``None`` when no positive heat release was found.
    notes : tuple[str, ...]
        Assumptions and fallbacks applied during the analysis.
    """

    case_name: str
    values: Mapping[str, float]
    heat_release: HeatReleaseResult
    burn: BurnFractionMetrics | None
    notes: tuple[str, ...] = ()

    def __getitem__(self, key: str) -> float:
        return self.values[key]


def _mixtures(case: EngineCase) -> tuple[NASAGasMixture, NASAGasMixture | None, str]:
    op = case.operating_point
    if op.lambda_air is None:
        return DRY_AIR, None, "lambda_air not given: dry-air properties, frozen composition."
    unburned = make_reactant_mixture(op.fuel, op.lambda_air)
    burned = make_product_mixture(op.fuel, op.lambda_air)
    return unburned, burned, (
        f"{op.fuel.name} at lambda={op.lambda_air:g}: x_b-weighted unburned/burned properties."
    )


def compute_cycle_metrics(
    case: EngineCase,
    *,
    combustion_window_deg: tuple[float, float] = (-40.0, 90.0),
    gross_window_deg: tuple[float, float] = (-180.0, 180.0),
    revolutions_per_cycle: int = 2,
    variable_composition: bool = True,
    min_net_to_absolute: float = 0.5,
) -> CycleMetrics:
    """Reduce one :class:`EngineCase` to standard indicated-cycle metrics.

    Parameters
    ----------
    case : EngineCase
        Case to analyze. Needs a temperature trace, or trapped mass and ``lambda_air``.
    combustion_window_deg : tuple[float, float], default=(-40, 90)
        Window [deg ATDC] for the mass fraction burned and CA10/50/90; clipped to the trace.
    gross_window_deg : tuple[float, float], default=(-180, 180)
        Compression-expansion window for gross indicated work (Heywood, 2018, Sec. 2.4).
    revolutions_per_cycle : int, default=2
        Crank revolutions per cycle $n_R$ (2 for four-stroke, 1 for two-stroke).
    variable_composition : bool, default=True
        Use $x_b$-weighted unburned/burned gas properties when ``lambda_air`` is known.
    min_net_to_absolute : float, default=0.5
        Plausibility check (a Combusta heuristic, not a literature criterion): burn
        fractions are reported only if $\\int dQ \\ge$ ``min_net_to_absolute`` $\\cdot
        \\int |dQ|$ over the window. Real combustion gives a ratio near 1; a motored or
        misfiring trace gives a small ratio because its apparent heat release is
        dominated by sign-alternating property and noise errors.

    Notes
    -----
    Gross indicated power per cylinder is $P_{ig} = W_{c,ig} N / (60 n_R)$, the gross
    indicated fuel conversion efficiency $\\eta_{f,ig} = W_{c,ig} / (m_f Q_{LHV})$, and
    the indicated specific fuel consumption $m_f / W_{c,ig}$ (Heywood, 2018, Eqs. 2.13,
    2.22, 2.24). The trace wall heat-loss channel, when present, must be positive for
    heat leaving the gas.
    """
    if revolutions_per_cycle not in (1, 2):
        raise ValueError("revolutions_per_cycle must be 1 (two-stroke) or 2 (four-stroke).")
    op = case.operating_point
    theta = case.trace.crank_angle_deg
    pressure = case.trace.pressure_pa
    notes: list[str] = []

    if case.trace.volume_m3 is not None:
        volume = case.trace.volume_m3
        dv = compute_derivative(theta, volume)
        notes.append("Volume and dV/dtheta from the trace.")
    else:
        kin = case.kinematics()
        volume, dv = kin.volume_m3, kin.dv_dtheta_m3_per_deg
        notes.append("Volume and dV/dtheta from slider-crank geometry.")

    unburned, burned, mixture_note = _mixtures(case)
    if not variable_composition:
        burned = None
        mixture_note = "Frozen unburned-gas composition (variable_composition=False)."
    notes.append(mixture_note)

    w_start = max(combustion_window_deg[0], float(theta[0]))
    w_end = min(combustion_window_deg[1], float(theta[-1]))
    if w_end <= w_start:
        raise ValueError("combustion_window_deg does not overlap the trace.")

    temperature = case.temperature_k()
    dq_ht = case.trace.wall_heat_loss_j_per_deg

    def heat_release(burned_gas: NASAGasMixture | None) -> HeatReleaseResult:
        return compute_heat_release(
            theta,
            pressure,
            volume,
            temperature_k=temperature,
            dv_dtheta_m3_per_deg=dv,
            dq_ht_j_per_deg=dq_ht,
            mixture=unburned,
            burned_mixture=burned_gas,
            combustion_window_deg=(w_start, w_end) if burned_gas is not None else None,
        )

    try:
        hrr = heat_release(burned)
    except PhysicalBoundsError:
        if burned is None:
            raise
        hrr = heat_release(None)
        notes.append("No positive heat release in the window: frozen composition used.")

    burn: BurnFractionMetrics | None = None
    if _combustion_detected(theta, hrr.dq_gross_j_per_deg, w_start, w_end, min_net_to_absolute):
        try:
            burn = extract_burn_fractions(theta, hrr.dq_gross_j_per_deg, start_deg=w_start, end_deg=w_end)
        except (PhysicalBoundsError, ValueError):
            burn = None
    if burn is None:
        notes.append("No combustion detected in the window: burn fractions undefined.")

    imep = compute_imep(
        theta, pressure, volume,
        displaced_volume_m3=case.geometry.displaced_volume_m3,
        dv_dtheta_m3_per_deg=dv,
        gross_window_deg=gross_window_deg,
    )
    dp = hrr.dp_dtheta_pa_per_deg
    i_max = int(np.argmax(pressure))

    nan = float("nan")
    values: dict[str, float] = {
        "p_max": float(pressure[i_max]),
        "theta_p_max": _refined_peak_angle(theta, pressure, i_max),
        "max_prr": float(np.max(dp)),
        "imep_gross": imep.imep_gross_pa,
        "imep_net": imep.imep_net_pa,
        "indicated_power_gross": imep.indicated_work_gross_j * op.rpm / (60.0 * revolutions_per_cycle),
        "eta_indicated_gross": nan,
        "isfc_gross": nan,
        "ca10": burn.ca10_deg if burn else nan,
        "ca50": burn.ca50_deg if burn else nan,
        "ca90": burn.ca90_deg if burn else nan,
        "burn_duration_10_90": burn.burn_duration_10_90_deg if burn else nan,
        "q_gross": _window_integral(hrr.cumulative_gross_j, theta, w_start, w_end),
        "heat_release_fraction": nan,
        "wall_heat_loss": nan,
    }
    if op.fuel_mass_per_cycle_kg is not None:
        q_fuel = op.fuel_energy_per_cycle_j
        values["eta_indicated_gross"] = imep.indicated_work_gross_j / q_fuel
        values["heat_release_fraction"] = values["q_gross"] / q_fuel
        if imep.indicated_work_gross_j > 0.0:
            values["isfc_gross"] = op.fuel_mass_per_cycle_kg / imep.indicated_work_gross_j
    if case.trace.wall_heat_loss_j_per_deg is not None:
        values["wall_heat_loss"] = float(trapezoid(case.trace.wall_heat_loss_j_per_deg, theta))
    values.update(_emission_metrics(case, imep.indicated_work_gross_j))
    values.update(case.scalars)

    return CycleMetrics(
        case_name=case.name, values=values, heat_release=hrr, burn=burn, notes=tuple(notes)
    )


def _emission_metrics(case: EngineCase, work_gross_j: float) -> dict[str, float]:
    """Gross indicated specific emissions, emission indices, and combustion efficiency."""
    emissions = case.emissions_kg_per_cycle
    m_fuel = case.operating_point.fuel_mass_per_cycle_kg
    out: dict[str, float] = {}
    for species, mass in emissions.items():
        out[f"{SPECIFIC_EMISSION_PREFIX}{species}"] = (
            specific_emission(mass, work_gross_j) if work_gross_j > 0.0 else float("nan")
        )
        out[f"{EMISSION_INDEX_PREFIX}{species}"] = (
            emission_index(mass, m_fuel) if m_fuel is not None else float("nan")
        )
    if m_fuel is not None and any(has_heating_value(sp) for sp in emissions):
        out["combustion_efficiency"] = combustion_efficiency(
            case.operating_point.fuel, m_fuel, emissions
        )
    return out


def _dynamic_spec(key: str) -> MetricSpec | None:
    if key.startswith(SPECIFIC_EMISSION_PREFIX):
        species = key.removeprefix(SPECIFIC_EMISSION_PREFIX)
        return MetricSpec(f"Gross indicated specific {species}", "g/kWh", 1.0e3 * J_PER_KWH)
    if key.startswith(EMISSION_INDEX_PREFIX):
        species = key.removeprefix(EMISSION_INDEX_PREFIX)
        return MetricSpec(f"Emission index {species}", "g/kg fuel", 1.0e3)
    return None


def _combustion_detected(
    theta: FloatArray, dq: FloatArray, start: float, end: float, min_ratio: float
) -> bool:
    mask = (theta >= start) & (theta <= end)
    if int(np.count_nonzero(mask)) < 5:
        return False
    net = float(trapezoid(dq[mask], theta[mask]))
    absolute = float(trapezoid(np.abs(dq[mask]), theta[mask]))
    return absolute > 0.0 and net >= min_ratio * absolute


def _refined_peak_angle(theta: FloatArray, values: FloatArray, idx: int) -> float:
    """Sub-grid peak location from a parabola through the three samples around ``idx``."""
    if idx == 0 or idx == theta.size - 1:
        return float(theta[idx])
    x = theta[idx - 1 : idx + 2]
    y = values[idx - 1 : idx + 2]
    a, b, _ = np.polyfit(x - x[1], y, 2)
    if a >= 0.0:
        return float(theta[idx])
    return float(x[1] - b / (2.0 * a))


def _window_integral(cumulative: FloatArray, theta: FloatArray, start: float, end: float) -> float:
    return float(np.interp(end, theta, cumulative) - np.interp(start, theta, cumulative))


def uncertainty_from_gci(gci: GCIResult) -> float:
    """Absolute numerical uncertainty $U = \\mathrm{GCI}_{fine}^{21} |\\phi_1|$ of the fine-grid value.

    Celik et al. (2008) report the GCI as a fraction of the fine-grid solution
    $\\phi_1$, so the case value to compare must be the fine-grid result.
    """
    if not math.isfinite(gci.gci_fine_21):
        raise PhysicalBoundsError("GCI is not finite; no numerical uncertainty can be derived.")
    if gci.convergence_regime == "DIVERGENT":
        raise PhysicalBoundsError("Divergent grid triplet: the GCI is not a valid uncertainty.")
    return float(gci.gci_fine_21 * abs(gci.phi1))


@dataclass(frozen=True, slots=True)
class ComparisonRow:
    """One metric across all cases, with differences relative to the baseline (SI units).

    Attributes
    ----------
    key, label, display_unit : str
        Metric identifier and display metadata.
    si_to_display : float
        Factor converting SI values to ``display_unit``.
    baseline_value : float
        Value of the baseline case.
    values : tuple[float, ...]
        Metric value per case (case order of the comparison).
    deltas : tuple[float, ...]
        Difference to the baseline value per case.
    uncertainties : tuple[float | None, ...]
        Absolute numerical uncertainty $U$ per case, if supplied.
    significant : tuple[bool | None, ...]
        ``True``/``False`` when $|\\Delta|$ exceeds / stays within $U_\\Delta$; ``None``
        when an uncertainty is missing (and for the baseline itself).
    """

    key: str
    label: str
    display_unit: str
    si_to_display: float
    baseline_value: float
    values: tuple[float, ...]
    deltas: tuple[float, ...]
    uncertainties: tuple[float | None, ...]
    significant: tuple[bool | None, ...]

    def relative_deltas(self) -> tuple[float, ...]:
        """Differences relative to the baseline value [-] (``nan`` for a zero baseline)."""
        base = self.baseline_value
        if base == 0.0 or not math.isfinite(base):
            return tuple(float("nan") for _ in self.deltas)
        return tuple(d / base for d in self.deltas)


@dataclass(frozen=True, slots=True)
class CaseComparison:
    """Side-by-side comparison of several engine cases.

    Attributes
    ----------
    case_names : tuple[str, ...]
        Case labels in table order.
    baseline_index : int
        Index of the reference case for differences.
    metrics : tuple[CycleMetrics, ...]
        Per-case metrics.
    rows : tuple[ComparisonRow, ...]
        Metric rows (only metrics finite for at least one case).
    contains_synthetic : bool
        ``True`` if any case is a synthetic verification fixture.
    """

    case_names: tuple[str, ...]
    baseline_index: int
    metrics: tuple[CycleMetrics, ...]
    rows: tuple[ComparisonRow, ...]
    contains_synthetic: bool

    def row(self, key: str) -> ComparisonRow:
        """Return the row of metric ``key``."""
        for r in self.rows:
            if r.key == key:
                return r
        raise KeyError(f"Metric '{key}' not in comparison.")

    def to_markdown(self, *, digits: int = 4) -> str:
        """Render the comparison as a Markdown table in display units.

        Each non-baseline cell shows the value, the difference to the baseline, and a
        marker: ``*`` significant (|Δ| > U_Δ), ``~`` within numerical uncertainty.
        """
        base = self.case_names[self.baseline_index]
        header = ["Metric", "Unit"] + [
            f"{n} (baseline)" if i == self.baseline_index else n for i, n in enumerate(self.case_names)
        ]
        lines = []
        if self.contains_synthetic:
            lines.append(
                "> **Synthetic verification data** — not engineering results.\n"
            )
        lines.append("| " + " | ".join(header) + " |")
        lines.append("|" + "---|" * len(header))
        for r in self.rows:
            cells = [r.label, r.display_unit]
            for i, (v, d, s) in enumerate(zip(r.values, r.deltas, r.significant)):
                if not math.isfinite(v):
                    cells.append("–")
                    continue
                text = f"{v * r.si_to_display:.{digits}g}"
                if i != self.baseline_index:
                    mark = {True: " *", False: " ~", None: ""}[s]
                    text += f" (Δ {d * r.si_to_display:+.{digits}g}{mark})"
                cells.append(text)
            lines.append("| " + " | ".join(cells) + " |")
        lines.append("")
        lines.append(f"Δ relative to **{base}**. `*` exceeds the combined numerical uncertainty, "
                     "`~` lies within it; no marker means no uncertainty was supplied.")
        return "\n".join(lines)

    def to_dataframe(self):  # noqa: ANN201 - pandas is an optional-at-call-time import
        """Return a tidy ``pandas.DataFrame`` (one row per metric and case, SI and display units)."""
        import pandas as pd

        records = []
        for r in self.rows:
            for name, v, d, u, s in zip(self.case_names, r.values, r.deltas, r.uncertainties, r.significant):
                records.append({
                    "metric": r.key,
                    "label": r.label,
                    "case": name,
                    "value_si": v,
                    "delta_si": d,
                    "uncertainty_si": u,
                    "significant": s,
                    "value_display": v * r.si_to_display,
                    "display_unit": r.display_unit,
                })
        return pd.DataFrame.from_records(records)


def compare_cases(
    cases: Sequence[EngineCase],
    *,
    baseline: int | str = 0,
    uncertainties: Mapping[str, Mapping[str, float]] | None = None,
    scalar_units: Mapping[str, tuple[str, float]] | None = None,
    **metric_kwargs: object,
) -> CaseComparison:
    """Compare engine cases metric by metric against a baseline case.

    Parameters
    ----------
    cases : Sequence[EngineCase]
        Two or more cases with unique names.
    baseline : int | str, default=0
        Baseline case index or name.
    uncertainties : Mapping[str, Mapping[str, float]] | None
        Absolute numerical uncertainty $U$ [SI] per case name and metric key, e.g.
        ``{"open": {"p_max": uncertainty_from_gci(gci_open)}}``.
    scalar_units : Mapping[str, tuple[str, float]] | None
        Display unit and SI-to-display factor for case scalars, e.g.
        ``{"NOx_g_per_kWh": ("g/kWh", 1.0)}``. Unlisted scalars are shown as stored.
    **metric_kwargs
        Forwarded to :func:`compute_cycle_metrics`.
    """
    if len(cases) < 2:
        raise ValueError("Provide at least two cases to compare.")
    names = tuple(c.name for c in cases)
    if len(set(names)) != len(names):
        raise ValueError(f"Case names must be unique, got {names}.")
    base_idx = names.index(baseline) if isinstance(baseline, str) else baseline
    if not 0 <= base_idx < len(cases):
        raise IndexError("baseline index out of range.")
    unc = uncertainties or {}
    unknown = set(unc) - set(names)
    if unknown:
        raise KeyError(f"Uncertainties given for unknown cases: {sorted(unknown)}.")

    metrics = tuple(compute_cycle_metrics(c, **metric_kwargs) for c in cases)  # type: ignore[arg-type]
    keys = list(METRIC_SPECS)
    keys += sorted({k for m in metrics for k in m.values} - set(METRIC_SPECS))
    units = dict(scalar_units or {})

    rows = []
    for key in keys:
        vals = tuple(float(m.values.get(key, float("nan"))) for m in metrics)
        if not any(math.isfinite(v) for v in vals):
            continue
        spec = METRIC_SPECS.get(key) or _dynamic_spec(key)
        if spec is None:
            unit, factor = units.get(key, ("", 1.0))
            spec = MetricSpec(key, unit, factor)
        base_val = vals[base_idx]
        deltas = tuple(v - base_val for v in vals)
        u_vals = tuple(unc.get(n, {}).get(key) for n in names)
        u_base = u_vals[base_idx]
        significant: list[bool | None] = []
        for i, (d, u) in enumerate(zip(deltas, u_vals)):
            if i == base_idx or u is None or u_base is None or not math.isfinite(d):
                significant.append(None)
            else:
                significant.append(abs(d) > math.hypot(u, u_base))
        rows.append(
            ComparisonRow(
                key=key,
                label=spec.label,
                display_unit=spec.display_unit,
                si_to_display=spec.si_to_display,
                baseline_value=base_val,
                values=vals,
                deltas=deltas,
                uncertainties=u_vals,
                significant=tuple(significant),
            )
        )

    return CaseComparison(
        case_names=names,
        baseline_index=base_idx,
        metrics=metrics,
        rows=tuple(rows),
        contains_synthetic=any(c.is_synthetic for c in cases),
    )
