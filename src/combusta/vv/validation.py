"""Validation of a simulation result against experimental data (ASME V&V 20-2009).

ASME V&V 20 compares a simulation value $S$ with an experimental value $D$ through the
comparison error and the validation uncertainty

$$E = S - D, \\qquad u_{val} = \\sqrt{u_{num}^2 + u_{input}^2 + u_D^2},$$

where $u_{num}$ is the numerical (discretization) uncertainty, $u_{input}$ the
uncertainty propagated from simulation inputs, and $u_D$ the experimental uncertainty,
all as *standard* uncertainties. The model error $\\delta_{model}$ then lies in
$E \\pm k\\,u_{val}$ for a chosen coverage factor $k$. If $|E| \\le k\\,u_{val}$ the model
error cannot be distinguished from the validation uncertainty; otherwise $E$ estimates
the model error. Combining in quadrature assumes the three contributions are
independent.

Helper functions build the three contributions:

* :func:`numerical_uncertainty_from_gci` converts a GCI (an expanded uncertainty) to a
  standard uncertainty by dividing by a coverage factor.
* :func:`input_uncertainty` propagates input standard uncertainties with sensitivity
  coefficients (first-order Taylor series, the sensitivity-coefficient method of
  V&V 20, Sec. 4); :func:`sensitivity_coefficient` estimates them by central differences.
* :func:`experimental_uncertainty` combines the Type A standard uncertainty of a mean,
  $s/\\sqrt{n}$, with a Type B (systematic) contribution (JCGM 100:2008, GUM).

All functions accept scalars or arrays, so crank-angle-resolved traces (for example a
pressure curve against measured pressure with a per-sample $u_D$) can be validated
point by point.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Final, TypeAlias

import numpy as np
from numpy.typing import ArrayLike, NDArray

from combusta.exceptions import PhysicalBoundsError
from combusta.vv.gci import GCIResult

FloatArray: TypeAlias = NDArray[np.float64]

#: Coverage factor relating a three-grid GCI (Fs = 1.25) to a standard uncertainty.
#: The GCI is treated as an expanded uncertainty at roughly 95 % coverage. NIST IR 8298
#: (2020), summarizing V&V 20, gives k ≈ 2 for a Gaussian error distribution about the
#: fine-grid solution and k ≈ 1.1-1.15 for a shifted one; 1.15 is the more conservative
#: (larger u_num) end of the latter. The value is a user choice and must be reported.
DEFAULT_GCI_COVERAGE_FACTOR: Final[float] = 1.15

#: Coverage factor for a ~95 % interval of an approximately normal distribution.
DEFAULT_VALIDATION_COVERAGE_FACTOR: Final[float] = 2.0


def _as_nonnegative(values: ArrayLike, name: str) -> FloatArray:
    arr = np.asarray(values, dtype=np.float64)
    if not np.all(np.isfinite(arr)) or np.any(arr < 0.0):
        raise PhysicalBoundsError(f"{name} must be finite and >= 0.")
    return arr


def numerical_uncertainty_from_gci(
    gci: GCIResult,
    *,
    coverage_factor: float = DEFAULT_GCI_COVERAGE_FACTOR,
) -> float:
    """Standard numerical uncertainty $u_{num} = \\mathrm{GCI}_{fine}^{21} |\\phi_1| / k$.

    The fine-grid value $\\phi_1$ is the simulation result $S$ to validate.
    """
    if coverage_factor <= 0.0:
        raise PhysicalBoundsError("coverage_factor must be > 0.")
    if gci.convergence_regime == "DIVERGENT" or not np.isfinite(gci.gci_fine_21):
        raise PhysicalBoundsError("A divergent or non-finite GCI is not a valid uncertainty.")
    return float(gci.gci_fine_21 * abs(gci.phi1) / coverage_factor)


def sensitivity_coefficient(
    result_plus: ArrayLike,
    result_minus: ArrayLike,
    input_step: float,
) -> FloatArray:
    """Central-difference sensitivity $\\partial S / \\partial X \\approx (S_+ - S_-) / (2 \\Delta X)$.

    ``result_plus`` and ``result_minus`` are simulation results at $X \\pm \\Delta X$.
    The estimate is second-order accurate and exact for results quadratic in $X$.
    """
    if not np.isfinite(input_step) or input_step <= 0.0:
        raise PhysicalBoundsError("input_step must be finite and > 0.")
    plus = np.asarray(result_plus, dtype=np.float64)
    minus = np.asarray(result_minus, dtype=np.float64)
    return np.asarray((plus - minus) / (2.0 * input_step), dtype=np.float64)


def input_uncertainty(
    sensitivities: Sequence[ArrayLike],
    input_standard_uncertainties: Sequence[float],
) -> FloatArray:
    """Propagated input uncertainty $u_{input} = \\sqrt{\\sum_i (\\partial S/\\partial X_i\\, u_{X_i})^2}$.

    First-order propagation for independent inputs (JCGM 100:2008, Eq. 10).
    """
    if len(sensitivities) != len(input_standard_uncertainties) or len(sensitivities) == 0:
        raise ValueError("Provide one input uncertainty per sensitivity coefficient.")
    u_x = _as_nonnegative(np.asarray(input_standard_uncertainties), "input_standard_uncertainties")
    terms = [np.asarray(c, dtype=np.float64) * u for c, u in zip(sensitivities, u_x)]
    return np.asarray(np.sqrt(sum(t**2 for t in terms)), dtype=np.float64)


def experimental_uncertainty(
    *,
    samples: ArrayLike | None = None,
    standard_deviation: ArrayLike | None = None,
    n_samples: int | None = None,
    systematic_standard_uncertainty: ArrayLike = 0.0,
) -> FloatArray:
    """Experimental standard uncertainty $u_D = \\sqrt{(s/\\sqrt{n})^2 + u_B^2}$ of a mean value.

    Give either repeated ``samples`` (axis 0 = repetition) or their
    ``standard_deviation`` and ``n_samples``. ``systematic_standard_uncertainty`` is the
    Type B contribution (sensor calibration, reference uncertainty).
    """
    if (samples is None) == (standard_deviation is None):
        raise ValueError("Provide either samples or standard_deviation with n_samples.")
    if samples is not None:
        arr = np.asarray(samples, dtype=np.float64)
        n = arr.shape[0]
        if n < 2:
            raise ValueError("At least two repeated samples are needed for a Type A estimate.")
        s = np.std(arr, axis=0, ddof=1)
    else:
        if n_samples is None or n_samples < 1:
            raise ValueError("n_samples must be >= 1 with standard_deviation.")
        n = n_samples
        s = _as_nonnegative(np.asarray(standard_deviation), "standard_deviation")
    u_a = s / np.sqrt(n)
    u_b = _as_nonnegative(systematic_standard_uncertainty, "systematic_standard_uncertainty")
    return np.asarray(np.sqrt(u_a**2 + u_b**2), dtype=np.float64)


@dataclass(frozen=True, slots=True)
class ValidationResult:
    """ASME V&V 20 validation comparison (scalars or arrays, SI units).

    Attributes
    ----------
    comparison_error : FloatArray
        $E = S - D$.
    u_num, u_input, u_d : FloatArray
        Standard numerical, input, and experimental uncertainties.
    u_val : FloatArray
        Validation standard uncertainty $\\sqrt{u_{num}^2 + u_{input}^2 + u_D^2}$.
    coverage_factor : float
        $k$ used for the model-error interval.
    model_error_lower, model_error_upper : FloatArray
        Interval $E \\mp k\\,u_{val}$ containing $\\delta_{model}$.
    error_exceeds_uncertainty : NDArray[np.bool_]
        ``True`` where $|E| > k\\,u_{val}$ (model error resolved above the noise level).
    """

    comparison_error: FloatArray
    u_num: FloatArray
    u_input: FloatArray
    u_d: FloatArray
    u_val: FloatArray
    coverage_factor: float
    model_error_lower: FloatArray
    model_error_upper: FloatArray
    error_exceeds_uncertainty: NDArray[np.bool_]

    def summary(self) -> str:
        """One-line interpretation for scalar comparisons."""
        if self.comparison_error.ndim != 0:
            share = float(np.mean(self.error_exceeds_uncertainty)) * 100.0
            return f"|E| exceeds k*u_val at {share:.1f} % of the samples (k = {self.coverage_factor:g})."
        e = float(self.comparison_error)
        band = self.coverage_factor * float(self.u_val)
        if bool(self.error_exceeds_uncertainty):
            return (
                f"E = {e:.4g} exceeds k*u_val = {band:.4g}: model error estimated as "
                f"{e:.4g}, within [{float(self.model_error_lower):.4g}, {float(self.model_error_upper):.4g}]."
            )
        return (
            f"|E| = {abs(e):.4g} <= k*u_val = {band:.4g}: model error is not distinguishable "
            "from the validation uncertainty."
        )


def validate(
    simulation: ArrayLike,
    experiment: ArrayLike,
    *,
    u_num: ArrayLike = 0.0,
    u_input: ArrayLike = 0.0,
    u_d: ArrayLike,
    coverage_factor: float = DEFAULT_VALIDATION_COVERAGE_FACTOR,
) -> ValidationResult:
    """Compare a simulation result with experimental data following ASME V&V 20-2009.

    Parameters
    ----------
    simulation, experiment : ArrayLike
        Simulation value(s) $S$ and experimental value(s) $D$ on the same support.
    u_num, u_input, u_d : ArrayLike
        Standard uncertainties (broadcast against ``simulation``). ``u_d`` is required;
        leaving ``u_num`` or ``u_input`` at zero must be justified in the report.
    coverage_factor : float, default=2.0
        $k$ for the model-error interval.
    """
    s = np.asarray(simulation, dtype=np.float64)
    d = np.asarray(experiment, dtype=np.float64)
    if s.shape != d.shape:
        raise ValueError("simulation and experiment must have the same shape.")
    if not (np.all(np.isfinite(s)) and np.all(np.isfinite(d))):
        raise PhysicalBoundsError("simulation and experiment must be finite.")
    if coverage_factor <= 0.0:
        raise PhysicalBoundsError("coverage_factor must be > 0.")
    un = np.broadcast_to(_as_nonnegative(u_num, "u_num"), s.shape)
    ui = np.broadcast_to(_as_nonnegative(u_input, "u_input"), s.shape)
    ud = np.broadcast_to(_as_nonnegative(u_d, "u_d"), s.shape)

    e = s - d
    u_val = np.sqrt(un**2 + ui**2 + ud**2)
    band = coverage_factor * u_val
    return ValidationResult(
        comparison_error=np.asarray(e, dtype=np.float64),
        u_num=np.asarray(un, dtype=np.float64),
        u_input=np.asarray(ui, dtype=np.float64),
        u_d=np.asarray(ud, dtype=np.float64),
        u_val=np.asarray(u_val, dtype=np.float64),
        coverage_factor=float(coverage_factor),
        model_error_lower=np.asarray(e - band, dtype=np.float64),
        model_error_upper=np.asarray(e + band, dtype=np.float64),
        error_exceeds_uncertainty=np.asarray(np.abs(e) > band),
    )
