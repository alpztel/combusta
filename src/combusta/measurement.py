"""Uncertainty of cycle metrics caused by cylinder-pressure measurement errors.

Pressure-derived metrics of a test-bench case inherit the errors of the pressure
measurement. Brunt & Emtage (1996) identify thermal shock, crank-angle phasing, and
transducer sensitivity as the main sources of IMEP error; absolute referencing
("pegging") of piezoelectric transducers is discussed by Randolph (1990). This module
propagates four of these error sources into every metric of
:func:`combusta.compare.compute_cycle_metrics`. The measured trace is modelled as

$$p_\\mathrm{meas}(\\theta) = (1 + \\varepsilon_g)\\,p(\\theta + \\delta_\\theta) + \\Delta p_0
  + n(\\theta),$$

with crank-angle phasing (TDC) error $\\delta_\\theta$, transducer gain error
$\\varepsilon_g$, pegging offset $\\Delta p_0$, and uncorrelated sample noise $n$. The
inputs are standard uncertainties of normally distributed, mutually independent errors.

Two estimates are reported:

- **Linear propagation** (JCGM 100:2008, Sec. 5.1) of the three systematic sources.
  Each contribution $c_i u_i$ is a central difference of the metric at
  $x_i = \\pm u_i$; their root-sum-square is $u_\\mathrm{lin}$. This is the uncertainty
  budget: it shows which source dominates.
- **Monte Carlo** propagation of all four sources (in the manner of JCGM 101:2008)
  with $M$ trials. Its standard deviation $u_\\mathrm{MC}$ includes noise and
  nonlinearity. Because each trial is a full heat-release analysis, $M$ is far below
  the number JCGM 101 recommends for coverage intervals; only the mean and standard
  deviation are reported, the latter with its standard error
  $u_\\mathrm{MC} / \\sqrt{2(M - 1)}$.

Not modelled: thermal shock (an intra-cycle, load-dependent error), crank-angle
encoder resolution, and errors of the volume or geometry. Noise is added to the
pressure before the analysis; if real data are filtered, pass the same filter as
``noise_filter`` so that its effect on the noise is reproduced.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, replace
import math
from typing import Final, Literal, TypeAlias

import numpy as np
from numpy.typing import NDArray
from scipy.interpolate import CubicSpline

from combusta.case import EngineCase
from combusta.compare import METRIC_SPECS, compute_cycle_metrics
from combusta.exceptions import PhysicalBoundsError

FloatArray: TypeAlias = NDArray[np.float64]
Source: TypeAlias = Literal["tdc", "offset", "gain"]
NoiseFilter: TypeAlias = Callable[[FloatArray, FloatArray], FloatArray]

DEFAULT_METRICS: Final[tuple[str, ...]] = (
    "p_max", "theta_p_max", "max_prr", "imep_gross", "imep_net", "eta_indicated_gross",
    "ca10", "ca50", "ca90", "burn_duration_10_90", "q_gross",
)
"""Metrics evaluated by default."""

SOURCE_LABELS: Final[Mapping[str, str]] = {
    "tdc": "TDC phasing",
    "offset": "Pegging offset",
    "gain": "Transducer gain",
}


@dataclass(frozen=True, slots=True)
class PressureMeasurementUncertainty:
    """Standard uncertainties of the cylinder-pressure measurement.

    Parameters
    ----------
    tdc_offset_deg : float, default=0
        Crank-angle phasing (TDC) error [deg].
    pegging_offset_pa : float, default=0
        Absolute pressure referencing error [Pa].
    gain_relative : float, default=0
        Relative transducer sensitivity error [-] (0.01 = 1 %).
    noise_pa : float, default=0
        Standard deviation of uncorrelated sample-to-sample noise [Pa].
    """

    tdc_offset_deg: float = 0.0
    pegging_offset_pa: float = 0.0
    gain_relative: float = 0.0
    noise_pa: float = 0.0

    def __post_init__(self) -> None:
        for name in ("tdc_offset_deg", "pegging_offset_pa", "gain_relative", "noise_pa"):
            value = getattr(self, name)
            if not (math.isfinite(value) and value >= 0.0):
                raise PhysicalBoundsError(f"{name} must be finite and >= 0, got {value}.")
        if self.gain_relative >= 0.5:
            raise PhysicalBoundsError("gain_relative must be below 0.5.")

    def systematic(self) -> dict[Source, float]:
        """Non-zero systematic standard uncertainties by source."""
        values: dict[Source, float] = {
            "tdc": self.tdc_offset_deg, "offset": self.pegging_offset_pa, "gain": self.gain_relative,
        }
        return {k: v for k, v in values.items() if v > 0.0}


def perturb_pressure(
    crank_angle_deg: FloatArray,
    pressure_pa: FloatArray,
    *,
    tdc_shift_deg: float = 0.0,
    offset_pa: float = 0.0,
    gain: float = 0.0,
    noise_pa: FloatArray | None = None,
) -> FloatArray:
    """Measured pressure $(1 + \\varepsilon_g) p(\\theta + \\delta_\\theta) + \\Delta p_0 + n$ [Pa].

    The phase shift interpolates the trace with a cubic spline (extrapolating by at most
    ``tdc_shift_deg`` at the ends).
    """
    theta = np.asarray(crank_angle_deg, dtype=np.float64)
    p = np.asarray(pressure_pa, dtype=np.float64)
    shifted = p if tdc_shift_deg == 0.0 else np.asarray(CubicSpline(theta, p)(theta + tdc_shift_deg))
    out = (1.0 + gain) * shifted + offset_pa
    if noise_pa is not None:
        out = out + np.asarray(noise_pa, dtype=np.float64)
    if np.any(out <= 0.0):
        raise PhysicalBoundsError("Perturbed pressure is not positive; the uncertainties are too large.")
    return np.asarray(out, dtype=np.float64)


@dataclass(frozen=True, slots=True)
class MetricUncertainty:
    """Measurement uncertainty of one metric (SI units).

    Attributes
    ----------
    key : str
        Metric key.
    nominal : float
        Metric of the unperturbed case.
    contributions : Mapping[str, float]
        $|c_i u_i|$ per systematic source.
    u_linear : float
        Root-sum-square of the contributions.
    mc_mean, mc_std, mc_std_se : float
        Monte Carlo mean, standard deviation, and the standard error of the standard
        deviation (``nan`` without Monte Carlo trials).
    mc_valid : int
        Number of trials in which the metric was defined.
    """

    key: str
    nominal: float
    contributions: Mapping[str, float]
    u_linear: float
    mc_mean: float
    mc_std: float
    mc_std_se: float
    mc_valid: int

    @property
    def standard_uncertainty(self) -> float:
        """Monte Carlo standard deviation if available, else the linear estimate."""
        return self.mc_std if math.isfinite(self.mc_std) else self.u_linear


@dataclass(frozen=True, slots=True)
class MeasurementUncertaintyResult:
    """Uncertainty of the cycle metrics of one case due to pressure-measurement errors.

    Attributes
    ----------
    case_name : str
        Analyzed case.
    inputs : PressureMeasurementUncertainty
        Standard uncertainties of the error sources.
    metrics : Mapping[str, MetricUncertainty]
        Result per metric key.
    n_samples : int
        Monte Carlo trials (0 = linear propagation only).
    seed : int | None
        Random seed of the Monte Carlo trials.
    notes : tuple[str, ...]
        Assumptions applied (e.g. temperature reconstruction).
    """

    case_name: str
    inputs: PressureMeasurementUncertainty
    metrics: Mapping[str, MetricUncertainty]
    n_samples: int
    seed: int | None
    notes: tuple[str, ...] = ()

    def uncertainties(self) -> dict[str, float]:
        """Standard uncertainty per metric [SI], e.g. as ``u_D`` for :func:`combusta.vv.validate`."""
        return {k: m.standard_uncertainty for k, m in self.metrics.items()
                if math.isfinite(m.standard_uncertainty)}

    def to_markdown(self, *, digits: int = 3) -> str:
        """Uncertainty budget in display units."""
        sources = list(self.inputs.systematic())
        header = ["Metric", "Unit", "Nominal"] + [SOURCE_LABELS[s] for s in sources] + ["u (linear)"]
        if self.n_samples:
            header += [f"u (Monte Carlo, M = {self.n_samples})"]
        lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
        for key, m in self.metrics.items():
            spec = METRIC_SPECS.get(key)
            unit, factor = (spec.display_unit, spec.si_to_display) if spec else ("", 1.0)

            def fmt(v: float) -> str:
                return "–" if not math.isfinite(v) else f"{v * factor:.{digits}g}"

            cells = [spec.label if spec else key, unit, fmt(m.nominal)]
            cells += [fmt(m.contributions.get(s, math.nan)) for s in sources] + [fmt(m.u_linear)]
            if self.n_samples:
                cells.append("–" if not math.isfinite(m.mc_std)
                             else f"{m.mc_std * factor:.{digits}g} ± {m.mc_std_se * factor:.1g}")
            lines.append("| " + " | ".join(cells) + " |")
        return "\n".join(lines)

    def to_dataframe(self):  # noqa: ANN201 - pandas imported at call time
        """Tidy ``pandas.DataFrame`` with one row per metric."""
        import pandas as pd

        rows = []
        for key, m in self.metrics.items():
            row: dict[str, object] = {"metric": key, "nominal_si": m.nominal, "u_linear_si": m.u_linear,
                                      "mc_mean_si": m.mc_mean, "mc_std_si": m.mc_std,
                                      "mc_std_se_si": m.mc_std_se, "mc_valid": m.mc_valid}
            row.update({f"contribution_{s}_si": v for s, v in m.contributions.items()})
            rows.append(row)
        return pd.DataFrame(rows)


def _perturbed_case(case: EngineCase, pressure: FloatArray, drop_temperature: bool) -> EngineCase:
    trace = replace(case.trace, pressure_pa=pressure,
                    temperature_k=None if drop_temperature else case.trace.temperature_k)
    return replace(case, trace=trace)


def _evaluate(case: EngineCase, keys: Sequence[str], metric_kwargs: Mapping[str, object]) -> dict[str, float]:
    values = compute_cycle_metrics(case, **metric_kwargs).values  # type: ignore[arg-type]
    return {k: float(values.get(k, math.nan)) for k in keys}


def pressure_measurement_uncertainty(
    case: EngineCase,
    inputs: PressureMeasurementUncertainty,
    *,
    n_samples: int = 500,
    seed: int | None = 0,
    metrics: Sequence[str] = DEFAULT_METRICS,
    noise_filter: NoiseFilter | None = None,
    **metric_kwargs: object,
) -> MeasurementUncertaintyResult:
    """Propagate pressure-measurement errors into the cycle metrics of ``case``.

    Parameters
    ----------
    case : EngineCase
        Case whose trace is the measured pressure.
    inputs : PressureMeasurementUncertainty
        Standard uncertainties of the error sources.
    n_samples : int, default=500
        Monte Carlo trials; 0 skips the Monte Carlo propagation.
    seed : int | None, default=0
        Random seed (``None`` for a non-reproducible run).
    metrics : Sequence[str]
        Metric keys to evaluate.
    noise_filter : callable, optional
        ``f(theta, p) -> p`` applied to every perturbed trace (and to the nominal one),
        e.g. the zero-phase filter used for the real data.
    **metric_kwargs
        Forwarded to :func:`combusta.compare.compute_cycle_metrics`.

    Notes
    -----
    A traced gas temperature cannot follow a perturbed pressure. When the operating
    point has a trapped mass, the temperature is dropped and reconstructed from the
    perturbed pressure; otherwise the traced temperature is kept, which understates the
    effect of the pressure errors on the gas properties (stated in ``notes``).
    """
    if n_samples < 0:
        raise ValueError("n_samples must be >= 0.")
    if n_samples == 1:
        raise ValueError("n_samples must be 0 or at least 2.")
    systematic = inputs.systematic()
    if not systematic and inputs.noise_pa == 0.0:
        raise ValueError("All standard uncertainties are zero; nothing to propagate.")
    if inputs.noise_pa > 0.0 and n_samples == 0:
        raise ValueError("Noise is propagated by Monte Carlo only; set n_samples >= 2.")

    theta = case.trace.crank_angle_deg
    p_nominal = case.trace.pressure_pa
    notes: list[str] = []
    drop_temperature = case.trace.temperature_k is not None and case.operating_point.trapped_mass_kg is not None
    if drop_temperature:
        notes.append("Traced temperature replaced by the ideal-gas reconstruction from the perturbed pressure.")
    elif case.trace.temperature_k is not None:
        notes.append("Traced temperature kept (no trapped mass); gas-property effects of pressure errors are "
                     "not propagated.")
    if noise_filter is not None:
        notes.append("Noise filter applied to every trace, including the nominal one.")

    def run(tdc: float = 0.0, offset: float = 0.0, gain: float = 0.0,
            noise: FloatArray | None = None) -> dict[str, float]:
        p = perturb_pressure(theta, p_nominal, tdc_shift_deg=tdc, offset_pa=offset, gain=gain, noise_pa=noise)
        if noise_filter is not None:
            p = np.asarray(noise_filter(theta, p), dtype=np.float64)
        return _evaluate(_perturbed_case(case, p, drop_temperature), metrics, metric_kwargs)

    nominal = run()
    contributions: dict[str, dict[str, float]] = {k: {} for k in metrics}
    def run_source(source: Source, value: float) -> dict[str, float]:
        if source == "tdc":
            return run(tdc=value)
        if source == "offset":
            return run(offset=value)
        return run(gain=value)

    for source, u in systematic.items():
        plus = run_source(source, u)
        minus = run_source(source, -u)
        for k in metrics:
            contributions[k][source] = abs(plus[k] - minus[k]) / 2.0

    samples: dict[str, list[float]] = {k: [] for k in metrics}
    if n_samples:
        rng = np.random.default_rng(seed)
        for _ in range(n_samples):
            noise = rng.normal(0.0, inputs.noise_pa, theta.size) if inputs.noise_pa > 0.0 else None
            values = run(
                tdc=float(rng.normal(0.0, inputs.tdc_offset_deg)) if inputs.tdc_offset_deg else 0.0,
                offset=float(rng.normal(0.0, inputs.pegging_offset_pa)) if inputs.pegging_offset_pa else 0.0,
                gain=float(rng.normal(0.0, inputs.gain_relative)) if inputs.gain_relative else 0.0,
                noise=noise,
            )
            for k in metrics:
                samples[k].append(values[k])

    result: dict[str, MetricUncertainty] = {}
    for k in metrics:
        contrib = contributions[k]
        u_lin = math.sqrt(sum(v**2 for v in contrib.values())) if contrib else math.nan
        arr = np.array(samples[k], dtype=np.float64)
        valid = arr[np.isfinite(arr)]
        if valid.size >= 2:
            std = float(np.std(valid, ddof=1))
            mean = float(np.mean(valid))
            se = std / math.sqrt(2.0 * (valid.size - 1))
        else:
            std = mean = se = math.nan
        result[k] = MetricUncertainty(
            key=k, nominal=nominal[k], contributions=contrib, u_linear=u_lin,
            mc_mean=mean, mc_std=std, mc_std_se=se, mc_valid=valid.size,
        )
    return MeasurementUncertaintyResult(
        case_name=case.name, inputs=inputs, metrics=result, n_samples=n_samples, seed=seed, notes=tuple(notes),
    )
