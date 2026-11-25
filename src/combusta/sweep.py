"""Comparison of engine configurations at matched combustion phasing or at MBT timing.

Two combustion systems with different burn rates (for example an open chamber and a
pre-chamber) reach different combustion phasing at the same spark timing, so a
comparison at equal spark timing mixes the effect of the combustion system with the
effect of phasing. The usual practice is to run a spark-timing sweep per configuration
and compare either

- at **matched CA50**: every metric is interpolated to the same crank angle of 50 %
  mass fraction burned, or
- at **MBT timing**: each configuration at the spark timing that maximizes its torque
  (maximum brake torque, MBT; Heywood, 2018, Sec. 9.2).

A :class:`Sweep` holds the cycle metrics of one configuration at several spark
timings. :func:`metrics_at` interpolates all metrics of a sweep to a target CA50 or
spark timing with a monotone piecewise-cubic Hermite interpolant (PCHIP; Fritsch &
Carlson, 1980); targets outside the sweep are rejected rather than extrapolated.
:func:`find_mbt` fits a quadratic to gross IMEP against spark timing around the
sweep maximum and returns its vertex,

$$\\theta_\\mathrm{MBT} = -\\frac{b}{2c}, \\qquad
  u^2(\\theta_\\mathrm{MBT}) = \\mathbf{g}^\\mathsf{T}\\,\\boldsymbol{\\Sigma}\\,\\mathbf{g},
  \\qquad \\mathbf{g} = \\left(0,\\; -\\frac{1}{2c},\\; \\frac{b}{2c^2}\\right),$$

with $\\boldsymbol{\\Sigma} = s^2 (\\mathbf{X}^\\mathsf{T}\\mathbf{X})^{-1}$ the covariance of
the least-squares coefficients $(a, b, c)$ and $s^2$ the residual variance (first-order
propagation, JCGM 100:2008, Sec. 5.1). MBT is defined on gross IMEP, which assumes that
friction and pumping work do not change with spark timing.
:func:`compare_matched` tabulates several configurations at the matched condition with
the same significance test as :func:`combusta.compare.compare_cases`.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
import math
from typing import Final, Literal, TypeAlias

import numpy as np
from numpy.typing import NDArray
from scipy.interpolate import PchipInterpolator

from combusta.case import EngineCase
from combusta.compare import (
    ComparisonRow,
    CycleMetrics,
    build_comparison_rows,
    compute_cycle_metrics,
    rows_to_dataframe,
    rows_to_markdown,
)
from combusta.exceptions import CombustaError

FloatArray: TypeAlias = NDArray[np.float64]
Control: TypeAlias = Literal["ca50", "ignition_timing"]

DEFAULT_MBT_FIT_POINTS: Final[int] = 5
"""Number of sweep points around the IMEP maximum used for the quadratic MBT fit."""


class SweepError(CombustaError, ValueError):
    """Raised when a sweep cannot support the requested matching or MBT search."""


@dataclass(frozen=True, slots=True)
class Sweep:
    """Cycle metrics of one engine configuration at several operating points.

    Attributes
    ----------
    name : str
        Configuration name (column label in comparisons).
    case_names : tuple[str, ...]
        Member case names.
    metrics : tuple[CycleMetrics, ...]
        Cycle metrics per member.
    uncertainties : tuple[Mapping[str, float], ...]
        Absolute numerical uncertainty $U$ [SI] per member and metric.
    contains_synthetic : bool
        ``True`` if any member is a synthetic verification case.
    """

    name: str
    case_names: tuple[str, ...]
    metrics: tuple[CycleMetrics, ...]
    uncertainties: tuple[Mapping[str, float], ...]
    contains_synthetic: bool

    def values(self, key: str) -> FloatArray:
        """Metric ``key`` [SI] of every member, in member order (``nan`` if missing)."""
        return np.array([float(m.values.get(key, math.nan)) for m in self.metrics])


def build_sweep(
    name: str,
    cases: Sequence[EngineCase],
    *,
    uncertainties: Mapping[str, Mapping[str, float]] | None = None,
    **metric_kwargs: object,
) -> Sweep:
    """Compute the cycle metrics of every case of one configuration.

    Parameters
    ----------
    name : str
        Configuration name.
    cases : Sequence[EngineCase]
        Two or more cases of the same configuration (typically a spark-timing sweep).
    uncertainties : Mapping[str, Mapping[str, float]] | None
        Absolute numerical uncertainty $U$ [SI] per case name and metric key.
    **metric_kwargs
        Forwarded to :func:`combusta.compare.compute_cycle_metrics`.
    """
    if len(cases) < 2:
        raise SweepError(f"Sweep '{name}' needs at least two cases.")
    names = tuple(c.name for c in cases)
    if len(set(names)) != len(names):
        raise SweepError(f"Case names in sweep '{name}' must be unique, got {names}.")
    unc = uncertainties or {}
    unknown = set(unc) - set(names)
    if unknown:
        raise SweepError(f"Uncertainties given for cases not in sweep '{name}': {sorted(unknown)}.")
    metrics = tuple(compute_cycle_metrics(c, **metric_kwargs) for c in cases)  # type: ignore[arg-type]
    return Sweep(
        name=name,
        case_names=names,
        metrics=metrics,
        uncertainties=tuple(dict(unc.get(n, {})) for n in names),
        contains_synthetic=any(c.is_synthetic for c in cases),
    )


@dataclass(frozen=True, slots=True)
class MBTResult:
    """Spark timing of maximum gross IMEP from a quadratic least-squares fit.

    Attributes
    ----------
    theta_mbt_deg : float
        Vertex of the fitted parabola [deg ATDC].
    u_theta_mbt_deg : float | None
        Standard uncertainty of the vertex from the fit residuals [deg]; ``None`` when
        the fit has no residual degrees of freedom (three points).
    metric_at_mbt : float
        Fitted metric value at the vertex [SI].
    coefficients : tuple[float, float, float]
        $(a, b, c)$ of $y = a + b\\,\\theta + c\\,\\theta^2$ [SI, deg].
    fit_case_names : tuple[str, ...]
        Cases used in the fit.
    residual_rms : float
        Root-mean-square fit residual [SI].
    metric : str
        Fitted metric key (``"imep_gross"`` by default).
    """

    theta_mbt_deg: float
    u_theta_mbt_deg: float | None
    metric_at_mbt: float
    coefficients: tuple[float, float, float]
    fit_case_names: tuple[str, ...]
    residual_rms: float
    metric: str


def fit_quadratic_vertex(
    x: FloatArray, y: FloatArray
) -> tuple[float, float | None, float, tuple[float, float, float], float]:
    """Least-squares parabola through ``(x, y)``; return vertex, its uncertainty, and fit data.

    Returns
    -------
    tuple
        ``(x_vertex, u_x_vertex or None, y_vertex, (a, b, c), residual_rms)``.

    Raises
    ------
    SweepError
        Fewer than three points, or a parabola that is not concave.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    n = x.size
    if n < 3:
        raise SweepError("A quadratic fit needs at least three points.")
    x0 = float(np.mean(x))
    xc = x - x0
    design = np.column_stack([np.ones(n), xc, xc**2])
    coef, *_ = np.linalg.lstsq(design, y, rcond=None)
    a_c, b_c, c_c = (float(v) for v in coef)
    if not c_c < 0.0:
        raise SweepError("The fitted parabola is not concave; the sweep shows no maximum.")
    vertex_c = -b_c / (2.0 * c_c)
    residual = y - design @ coef
    rms = float(math.sqrt(float(np.mean(residual**2))))
    u_vertex: float | None = None
    if n > 3:
        s2 = float(residual @ residual) / (n - 3)
        cov = s2 * np.linalg.inv(design.T @ design)
        grad = np.array([0.0, -1.0 / (2.0 * c_c), b_c / (2.0 * c_c**2)])
        u_vertex = float(math.sqrt(max(float(grad @ cov @ grad), 0.0)))
    y_vertex = a_c + b_c * vertex_c + c_c * vertex_c**2
    # Coefficients in the original (uncentred) variable.
    a = a_c - b_c * x0 + c_c * x0**2
    b = b_c - 2.0 * c_c * x0
    return vertex_c + x0, u_vertex, y_vertex, (a, b, c_c), rms


def find_mbt(
    sweep: Sweep,
    *,
    metric: str = "imep_gross",
    fit_points: int = DEFAULT_MBT_FIT_POINTS,
) -> MBTResult:
    """Spark timing of maximum ``metric`` (gross IMEP) from a quadratic fit around the maximum.

    The ``fit_points`` sweep points nearest to the sampled maximum (in spark-timing
    order) are fitted. The maximum must be bracketed: if the largest value lies at the
    earliest or latest spark timing, the MBT timing is outside the sweep and no value is
    returned. Knock- or pre-ignition-limited configurations often cannot reach MBT; use
    a matched-CA50 comparison for them.

    Raises
    ------
    SweepError
        Missing spark timings or metric values, duplicate timings, an unbracketed
        maximum, a non-concave fit, or a vertex outside the fitted range.
    """
    if fit_points < 3:
        raise SweepError("fit_points must be at least 3.")
    timing = sweep.values("ignition_timing")
    y = sweep.values(metric)
    if not np.all(np.isfinite(timing)):
        raise SweepError(f"Sweep '{sweep.name}': every case needs ignition_timing_deg for an MBT search.")
    if not np.all(np.isfinite(y)):
        raise SweepError(f"Sweep '{sweep.name}': metric '{metric}' is missing for some cases.")
    order = np.argsort(timing)
    t_sorted, y_sorted = timing[order], y[order]
    if np.any(np.diff(t_sorted) <= 0.0):
        raise SweepError(f"Sweep '{sweep.name}': spark timings must be distinct.")
    i_max = int(np.argmax(y_sorted))
    if i_max in (0, t_sorted.size - 1):
        edge = "earliest" if i_max == 0 else "latest"
        raise SweepError(
            f"Sweep '{sweep.name}': the maximum of '{metric}' is at the {edge} spark timing "
            f"({t_sorted[i_max]:g} deg ATDC); MBT is not bracketed by the sweep."
        )
    k = min(fit_points, t_sorted.size)
    lo = min(max(i_max - k // 2, 0), t_sorted.size - k)
    sel = slice(lo, lo + k)
    vertex, u_vertex, y_vertex, coef, rms = fit_quadratic_vertex(t_sorted[sel], y_sorted[sel])
    if not t_sorted[sel][0] <= vertex <= t_sorted[sel][-1]:
        raise SweepError(
            f"Sweep '{sweep.name}': the fitted MBT timing {vertex:.2f} deg ATDC lies outside the "
            f"fitted range [{t_sorted[sel][0]:g}, {t_sorted[sel][-1]:g}]."
        )
    names = tuple(sweep.case_names[int(i)] for i in order[sel])
    return MBTResult(
        theta_mbt_deg=vertex,
        u_theta_mbt_deg=u_vertex,
        metric_at_mbt=y_vertex,
        coefficients=coef,
        fit_case_names=names,
        residual_rms=rms,
        metric=metric,
    )


@dataclass(frozen=True, slots=True)
class MatchedPoint:
    """All metrics of one sweep interpolated to a target value of a control variable.

    Attributes
    ----------
    sweep_name : str
        Configuration name.
    control : {"ca50", "ignition_timing"}
        Control variable.
    target : float
        Target value of the control variable [deg ATDC].
    values : Mapping[str, float]
        Interpolated metric values [SI] (PCHIP).
    uncertainties : Mapping[str, float]
        Numerical uncertainty [SI] linearly interpolated between the two bracketing
        cases, for metrics where both give one.
    interpolation_spread : Mapping[str, float]
        $|y_\\mathrm{PCHIP} - y_\\mathrm{linear}|$ per metric [SI]. A heuristic indicator
        of how strongly the result depends on the interpolant (large values mean the
        sweep is too coarse); it is not an uncertainty estimate.
    bracket : tuple[str, str]
        Names of the cases on either side of the target.
    notes : tuple[str, ...]
        Warnings, e.g. CA50 not monotonic in spark timing.
    """

    sweep_name: str
    control: Control
    target: float
    values: Mapping[str, float]
    uncertainties: Mapping[str, float]
    interpolation_spread: Mapping[str, float]
    bracket: tuple[str, str]
    notes: tuple[str, ...] = ()


def metrics_at(sweep: Sweep, *, control: Control, target: float) -> MatchedPoint:
    """Interpolate every metric of ``sweep`` to ``control = target``.

    Parameters
    ----------
    sweep : Sweep
        Configuration sweep.
    control : {"ca50", "ignition_timing"}
        Control variable; its values must be finite and distinct across the sweep.
    target : float
        Target [deg ATDC]; must lie within the sweep range (no extrapolation).

    Notes
    -----
    Only metrics that are finite for every member are interpolated. At a member's own
    control value the result equals that member's metrics.
    """
    if control not in ("ca50", "ignition_timing"):
        raise SweepError(f"Unknown control variable '{control}'.")
    x = sweep.values(control)
    if not np.all(np.isfinite(x)):
        raise SweepError(f"Sweep '{sweep.name}': '{control}' is undefined for some cases.")
    order = np.argsort(x)
    xs = x[order]
    if np.any(np.diff(xs) <= 0.0):
        raise SweepError(f"Sweep '{sweep.name}': '{control}' values must be distinct.")
    if not xs[0] <= target <= xs[-1]:
        raise SweepError(
            f"Sweep '{sweep.name}': target {control} = {target:g} lies outside the sweep range "
            f"[{xs[0]:g}, {xs[-1]:g}]; extend the sweep instead of extrapolating."
        )
    notes: list[str] = []
    if control == "ca50":
        timing = sweep.values("ignition_timing")
        if np.all(np.isfinite(timing)):
            t_by_ca50 = timing[order]
            d = np.diff(t_by_ca50)
            if not (np.all(d < 0.0) or np.all(d > 0.0)):
                notes.append("CA50 is not monotonic in spark timing across the sweep; check for "
                             "partial burns or misfires before trusting the interpolation.")

    keys = sorted({k for m in sweep.metrics for k in m.values})
    values: dict[str, float] = {}
    spread: dict[str, float] = {}
    for key in keys:
        y = sweep.values(key)[order]
        if not np.all(np.isfinite(y)):
            continue
        pchip = float(PchipInterpolator(xs, y)(target))
        linear = float(np.interp(target, xs, y))
        values[key] = pchip
        spread[key] = abs(pchip - linear)
    values[control] = float(target)
    spread[control] = 0.0

    hi = int(np.searchsorted(xs, target, side="left"))
    hi = min(max(hi, 1), xs.size - 1)
    lo = hi - 1
    i_lo, i_hi = int(order[lo]), int(order[hi])
    w = (target - xs[lo]) / (xs[hi] - xs[lo])
    unc_lo, unc_hi = sweep.uncertainties[i_lo], sweep.uncertainties[i_hi]
    uncertainties = {
        key: (1.0 - w) * unc_lo[key] + w * unc_hi[key]
        for key in set(unc_lo) & set(unc_hi)
        if key in values
    }
    return MatchedPoint(
        sweep_name=sweep.name,
        control=control,
        target=float(target),
        values=values,
        uncertainties=uncertainties,
        interpolation_spread=spread,
        bracket=(sweep.case_names[i_lo], sweep.case_names[i_hi]),
        notes=tuple(notes),
    )


@dataclass(frozen=True, slots=True)
class MatchedComparison:
    """Configurations compared at matched CA50 or at their own MBT timing.

    Attributes
    ----------
    mode : {"ca50", "mbt"}
        Matching mode.
    case_names : tuple[str, ...]
        Configuration names (table columns).
    baseline_index : int
        Reference configuration.
    points : tuple[MatchedPoint, ...]
        Interpolated metrics per configuration.
    mbt : tuple[MBTResult | None, ...]
        MBT fit per configuration (``None`` in CA50 mode).
    rows : tuple[ComparisonRow, ...]
        Metric rows with differences and significance.
    contains_synthetic : bool
        ``True`` if any sweep contains synthetic cases.
    """

    mode: Literal["ca50", "mbt"]
    case_names: tuple[str, ...]
    baseline_index: int
    points: tuple[MatchedPoint, ...]
    mbt: tuple[MBTResult | None, ...]
    rows: tuple[ComparisonRow, ...]
    contains_synthetic: bool

    def row(self, key: str) -> ComparisonRow:
        """Return the row of metric ``key``."""
        for r in self.rows:
            if r.key == key:
                return r
        raise KeyError(f"Metric '{key}' not in comparison.")

    def describe(self) -> str:
        """One sentence stating the matching condition."""
        if self.mode == "ca50":
            return f"All configurations interpolated to CA50 = {self.points[0].target:g} deg ATDC."
        return "Each configuration at its own MBT spark timing (maximum gross IMEP)."

    def to_markdown(self, *, digits: int = 4) -> str:
        """Matching condition, comparison table, bracketing cases, and MBT fits."""
        lines = [self.describe(), ""]
        lines.append(rows_to_markdown(self.case_names, self.baseline_index, self.rows,
                                      self.contains_synthetic, digits=digits))
        lines += ["", "| Configuration | Interpolated between | Notes |", "|---|---|---|"]
        for p in self.points:
            lines.append(f"| {p.sweep_name} | {p.bracket[0]}, {p.bracket[1]} | {' '.join(p.notes) or '–'} |")
        if self.mode == "mbt":
            lines += ["", "| Configuration | MBT timing [deg ATDC] | u [deg] | Fit points | RMS residual [bar] |",
                      "|---|---|---|---|---|"]
            for name, res in zip(self.case_names, self.mbt):
                if res is None:
                    continue
                u = "–" if res.u_theta_mbt_deg is None else f"{res.u_theta_mbt_deg:.2f}"
                lines.append(f"| {name} | {res.theta_mbt_deg:.2f} | {u} | {len(res.fit_case_names)} | "
                             f"{res.residual_rms * 1e-5:.3g} |")
        return "\n".join(lines)

    def to_dataframe(self):  # noqa: ANN201 - pandas imported at call time
        """Tidy ``pandas.DataFrame`` with one row per metric and configuration."""
        frame = rows_to_dataframe(self.case_names, self.rows)
        spread = {(p.sweep_name, k): v for p in self.points for k, v in p.interpolation_spread.items()}
        frame["interpolation_spread_si"] = [spread.get((c, m), math.nan) for c, m in zip(frame["case"], frame["metric"])]
        return frame


def compare_matched(
    sweeps: Sequence[Sweep],
    *,
    ca50_deg: float | None = None,
    at_mbt: bool = False,
    baseline: int | str = 0,
    mbt_metric: str = "imep_gross",
    fit_points: int = DEFAULT_MBT_FIT_POINTS,
    scalar_units: Mapping[str, tuple[str, float]] | None = None,
) -> MatchedComparison:
    """Compare configurations at matched CA50 (``ca50_deg``) or at MBT timing (``at_mbt``).

    Exactly one of ``ca50_deg`` and ``at_mbt`` must be given. In MBT mode each
    configuration is evaluated at its own fitted MBT timing (:func:`find_mbt`) by
    interpolating all metrics to that spark timing.
    """
    if (ca50_deg is None) == (not at_mbt):
        raise SweepError("Give exactly one of ca50_deg and at_mbt=True.")
    if len(sweeps) < 2:
        raise SweepError("Provide at least two sweeps to compare.")
    names = tuple(s.name for s in sweeps)
    if len(set(names)) != len(names):
        raise SweepError(f"Sweep names must be unique, got {names}.")
    base_idx = names.index(baseline) if isinstance(baseline, str) else baseline
    if not 0 <= base_idx < len(sweeps):
        raise SweepError("baseline index out of range.")

    if at_mbt:
        mbt = tuple(find_mbt(s, metric=mbt_metric, fit_points=fit_points) for s in sweeps)
        points = tuple(
            metrics_at(s, control="ignition_timing", target=r.theta_mbt_deg) for s, r in zip(sweeps, mbt)
        )
        mbt_out: tuple[MBTResult | None, ...] = mbt
        mode: Literal["ca50", "mbt"] = "mbt"
    else:
        assert ca50_deg is not None
        points = tuple(metrics_at(s, control="ca50", target=ca50_deg) for s in sweeps)
        mbt_out = tuple(None for _ in sweeps)
        mode = "ca50"

    rows = build_comparison_rows(
        names, [p.values for p in points], base_idx, [p.uncertainties for p in points],
        scalar_units=scalar_units,
    )
    return MatchedComparison(
        mode=mode,
        case_names=names,
        baseline_index=base_idx,
        points=points,
        mbt=mbt_out,
        rows=rows,
        contains_synthetic=any(s.contains_synthetic for s in sweeps),
    )
