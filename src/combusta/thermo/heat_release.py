"""Single-zone First-Law AHRR, burn-fraction (CA10/50/90), IMEP, and knock metrics.

Implements Heywood (2018) single-zone Apparent & Gross Heat Release Rate with
variable $\\gamma(T)$, sub-degree cubic spline Mass Fraction Burned (MFB)
extraction, cycle-integrated IMEP, and MAPO / Eng (2002) Ringing Intensity.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.interpolate import CubicSpline

from combusta.exceptions import PhysicalBoundsError
from combusta.signal.filtering import (
    DerivativeMethod,
    butterworth_zero_phase,
    compute_derivative,
)
from combusta.thermo.fuels import blend_burned_unburned
from combusta.thermo.nasa_polynomials import (
    DIESEL_BURNED_PRODUCTS,
    DRY_AIR,
    NASAGasMixture,
)

FloatArray = NDArray[np.float64]


@dataclass(frozen=True, slots=True)
class HeatReleaseResult:
    """First-Law Apparent and Gross Heat Release Rate results in strict SI units.

    Attributes
    ----------
    crank_angle_deg : NDArray[np.float64]
        Crank angle $\\theta$ [deg].
    gamma : NDArray[np.float64]
        Instantaneous specific heat ratio $\\gamma(T)$ [-].
    dp_dtheta_pa_per_deg : NDArray[np.float64]
        Pressure derivative $dP/d\\theta$ [$\\text{Pa/deg}$].
    dv_dtheta_m3_per_deg : NDArray[np.float64]
        Volume derivative $dV/d\\theta$ [$\\text{m}^3/\\text{deg}$].
    dq_net_j_per_deg : NDArray[np.float64]
        Apparent net heat-release rate $dQ_n/d\\theta$ [$\\text{J/deg}$].
    dq_ht_j_per_deg : NDArray[np.float64]
        Wall heat-transfer loss rate $dQ_{ht}/d\\theta$ [$\\text{J/deg}$].
    dq_gross_j_per_deg : NDArray[np.float64]
        Gross chemical heat-release rate $dQ_{\\text{gross}}/d\\theta$ [$\\text{J/deg}$].
    cumulative_net_j : NDArray[np.float64]
        Cumulative net heat release $Q_n(\\theta)$ [J].
    cumulative_gross_j : NDArray[np.float64]
        Cumulative gross heat release $Q_{\\text{gross}}(\\theta)$ [J].
    mass_fraction_burned : NDArray[np.float64] | None
        Mass fraction burned $x_b(\\theta)$ [-] used to blend unburned/burned gas
        properties, or ``None`` when a single mixture was used.
    property_iterations : int
        Fixed-point iterations needed for the $x_b$-dependent $\\gamma$ to converge
        (0 for a single mixture).
    """

    crank_angle_deg: FloatArray
    gamma: FloatArray
    dp_dtheta_pa_per_deg: FloatArray
    dv_dtheta_m3_per_deg: FloatArray
    dq_net_j_per_deg: FloatArray
    dq_ht_j_per_deg: FloatArray
    dq_gross_j_per_deg: FloatArray
    cumulative_net_j: FloatArray
    cumulative_gross_j: FloatArray
    mass_fraction_burned: FloatArray | None = None
    property_iterations: int = 0


@dataclass(frozen=True, slots=True)
class BurnFractionMetrics:
    """Sub-degree Mass Fraction Burned (MFB) combustion phasing metrics.

    Attributes
    ----------
    ca10_deg : float
        Crank angle at 10% cumulative heat release [deg].
    ca50_deg : float
        Crank angle at 50% cumulative heat release [deg].
    ca90_deg : float
        Crank angle at 90% cumulative heat release [deg].
    burn_duration_10_90_deg : float
        10%–90% combustion duration $\\Delta\\theta_{10\\text{--}90}$ [deg].
    total_heat_j : float
        Total integrated heat release across combustion window [J].
    crank_angle_deg : NDArray[np.float64]
        Window crank angle array [deg].
    mfb : NDArray[np.float64]
        Normalized cumulative burn fraction $x_b(\\theta)$ [-].
    """

    ca10_deg: float
    ca50_deg: float
    ca90_deg: float
    burn_duration_10_90_deg: float
    total_heat_j: float
    crank_angle_deg: FloatArray
    mfb: FloatArray


@dataclass(frozen=True, slots=True)
class IMEPResult:
    """Indicated Mean Effective Pressure and cycle work summary in SI units.

    Attributes
    ----------
    imep_gross_pa : float
        Closed compression/expansion ($[-180^\\circ, +180^\\circ]$) IMEP [Pa].
    imep_net_pa : float
        Full-domain net IMEP [Pa].
    pmep_pa : float
        Pumping Mean Effective Pressure $\\text{IMEP}_g - \\text{IMEP}_n$ [Pa].
    indicated_work_gross_j : float
        Closed compression/expansion indicated boundary work [J].
    indicated_work_net_j : float
        Full-domain net indicated boundary work $\\oint P\\,dV$ [J].
    """

    imep_gross_pa: float
    imep_net_pa: float
    pmep_pa: float
    indicated_work_gross_j: float
    indicated_work_net_j: float


@dataclass(frozen=True, slots=True)
class KnockMetrics:
    """Acoustic knock amplitude (MAPO) and Eng (2002) Ringing Intensity metrics.

    Attributes
    ----------
    mapo_pa : float
        Maximum Amplitude of Pressure Oscillations $\\max|\\Delta P_{\\text{knock}}|$ [Pa].
    ringing_intensity_w_m2 : float
        Eng (2002) acoustic Ringing Intensity $RI$ [$\\text{W/m}^2$].
    max_dp_dtheta_pa_per_deg : float
        Peak cylinder pressure-rise rate $(dP/d\\theta)_{\\max}$ [$\\text{Pa/deg}$].
    knock_pressure_pa : NDArray[np.float64] | None
        Zero-phase bandpass-isolated acoustic pressure oscillation trace [Pa],
        or ``None`` when ``use_prr_correlation=True`` or Nyquist-limited.
    """

    mapo_pa: float
    ringing_intensity_w_m2: float
    max_dp_dtheta_pa_per_deg: float
    knock_pressure_pa: FloatArray | None = None


def compute_heat_release(
    crank_angle_deg: FloatArray,
    pressure_pa: FloatArray,
    volume_m3: FloatArray,
    *,
    temperature_k: FloatArray | None = None,
    gamma: FloatArray | float | None = None,
    dv_dtheta_m3_per_deg: FloatArray | None = None,
    dq_ht_j_per_deg: FloatArray | None = None,
    mixture: NASAGasMixture = DRY_AIR,
    burned_mixture: NASAGasMixture | None = None,
    combustion_window_deg: tuple[float, float] | None = None,
    max_property_iterations: int = 50,
    property_tolerance: float = 1.0e-9,
    derivative_method: DerivativeMethod = "auto",
) -> HeatReleaseResult:
    """Compute single-zone Apparent & Gross Heat Release Rate (Heywood, 2018).

    Evaluates the First-Law closed-system energy balance (Heywood, 2018, Eq. 9.27):

    $$\\frac{dQ_n}{d\\theta} = \\frac{\\gamma(T)}{\\gamma(T) - 1} P \\frac{dV}{d\\theta}
      + \\frac{1}{\\gamma(T) - 1} V \\frac{dP}{d\\theta}$$

    $$\\frac{dQ_{\\text{gross}}}{d\\theta} = \\frac{dQ_n}{d\\theta} + \\frac{dQ_{ht}}{d\\theta}$$

    Parameters
    ----------
    crank_angle_deg : NDArray[np.float64]
        Strictly increasing crank angle $\\theta$ [deg].
    pressure_pa : NDArray[np.float64]
        In-cylinder absolute pressure $P(\\theta)$ [Pa] (optionally pre-filtered).
    volume_m3 : NDArray[np.float64]
        Instantaneous cylinder volume $V(\\theta)$ [$\\text{m}^3$].
    temperature_k : NDArray[np.float64] | None, optional
        Mass-averaged gas temperature $T(\\theta)$ [K] for NASA $\\gamma(T)$ evaluation;
        mandatory unless explicit ``gamma`` is provided.
    gamma : NDArray[np.float64] | float | None, optional
        Explicit specific heat ratio $\\gamma$ [-] (used for analytical polytropic
        verification where $\\gamma = \\kappa$). When ``None``, evaluated via NASA
        7-term polynomials from ``temperature_k``.
    dv_dtheta_m3_per_deg : NDArray[np.float64] | None, optional
        Analytical volume derivative $dV/d\\theta$ [$\\text{m}^3/\\text{deg}$]; computed
        numerically via 4th-order differentiation if omitted.
    dq_ht_j_per_deg : NDArray[np.float64] | None, optional
        Wall heat-loss rate $dQ_{ht}/d\\theta$ [$\\text{J/deg}$] from Woschni/Hohenberg.
    mixture : NASAGasMixture, default=DRY_AIR
        NASA polynomial gas model used when ``gamma is None``. When ``burned_mixture``
        is given, this is the unburned charge (e.g. from ``make_reactant_mixture``).
    burned_mixture : NASAGasMixture | None, optional
        Burned-gas model (e.g. from ``make_product_mixture``). When provided,
        $\\gamma(T, x_b)$ blends unburned and burned properties by mass fraction burned
        (``blend_burned_unburned``), with $x_b(\\theta)$ taken as the normalized
        cumulative gross heat release over ``combustion_window_deg``. Because $x_b$
        depends on the heat release it is resolved by fixed-point iteration
        (Gatowski et al., 1984, SAE 841359; Brunt et al., 1998, SAE 981052).
    combustion_window_deg : tuple[float, float] | None, optional
        Crank-angle window [deg] over which $x_b$ rises from 0 to 1; defaults to the
        full trace. $x_b = 0$ before and $x_b = 1$ after the window.
    max_property_iterations : int, default=50
        Maximum fixed-point iterations for $x_b$.
    property_tolerance : float, default=1e-9
        Convergence threshold on $\\max|\\Delta x_b|$ [-].
    derivative_method : {"auto", "central4", "spline"}, default="auto"
        Numerical differentiation scheme for $dP/d\\theta$ (and $dV/d\\theta$ if needed).

    Returns
    -------
    HeatReleaseResult
        Complete First-Law net/gross HRR and cumulative energy arrays in SI units.
    """
    theta = np.asarray(crank_angle_deg, dtype=np.float64)
    p_pa = np.asarray(pressure_pa, dtype=np.float64)
    v_m3 = np.asarray(volume_m3, dtype=np.float64)

    if theta.ndim != 1 or not (theta.shape == p_pa.shape == v_m3.shape):
        raise ValueError("crank_angle_deg, pressure_pa, and volume_m3 must be 1D arrays of identical length.")
    if not (np.all(p_pa > 0.0) and np.all(v_m3 > 0.0)):
        raise PhysicalBoundsError("pressure_pa and volume_m3 must be finite and strictly positive.")

    if gamma is None:
        if temperature_k is None:
            raise ValueError("Provide temperature_k for NASA gamma(T) or explicit gamma.")
        t_arr = np.asarray(temperature_k, dtype=np.float64)
        if t_arr.shape != theta.shape:
            raise ValueError("temperature_k must share the same shape as crank_angle_deg.")
        gamma_arr = mixture.gamma(t_arr)
    elif isinstance(gamma, (int, float)):
        if not np.isfinite(gamma) or gamma <= 1.0:
            raise PhysicalBoundsError("Specific heat ratio gamma must be finite and > 1.0.")
        gamma_arr = np.full_like(theta, float(gamma))
    else:
        gamma_arr = np.asarray(gamma, dtype=np.float64)
        if gamma_arr.shape != theta.shape:
            raise ValueError("gamma array must share the same shape as crank_angle_deg.")
        if not np.all(gamma_arr > 1.0):
            raise PhysicalBoundsError("Specific heat ratio gamma must be finite and > 1.0 everywhere.")

    dp_dtheta = compute_derivative(theta, p_pa, method=derivative_method)
    if dv_dtheta_m3_per_deg is not None:
        dv_dtheta = np.asarray(dv_dtheta_m3_per_deg, dtype=np.float64)
        if dv_dtheta.shape != theta.shape:
            raise ValueError("dv_dtheta_m3_per_deg must match crank_angle_deg shape.")
    else:
        dv_dtheta = compute_derivative(theta, v_m3, method=derivative_method)

    if dq_ht_j_per_deg is not None:
        dq_ht = np.asarray(dq_ht_j_per_deg, dtype=np.float64)
        if dq_ht.shape != theta.shape:
            raise ValueError("dq_ht_j_per_deg must match crank_angle_deg shape.")
    else:
        dq_ht = np.zeros_like(theta)

    def first_law(gamma_vals: FloatArray) -> tuple[FloatArray, FloatArray]:
        inv_gm1 = 1.0 / (gamma_vals - 1.0)
        net = gamma_vals * inv_gm1 * p_pa * dv_dtheta + inv_gm1 * v_m3 * dp_dtheta
        return net, net + dq_ht

    dq_net, dq_gross = first_law(gamma_arr)
    x_b_result: FloatArray | None = None
    iterations = 0

    if burned_mixture is not None:
        if gamma is not None or temperature_k is None:
            raise ValueError("burned_mixture requires temperature_k and no explicit gamma.")
        if max_property_iterations < 1 or property_tolerance <= 0.0:
            raise ValueError("max_property_iterations must be >= 1 and property_tolerance > 0.")
        t_arr = np.asarray(temperature_k, dtype=np.float64)
        w_start, w_end = combustion_window_deg or (float(theta[0]), float(theta[-1]))
        if not (theta[0] <= w_start < w_end <= theta[-1]):
            raise ValueError("combustion_window_deg must be an increasing pair inside the trace.")

        x_b = np.zeros_like(theta)
        converged = False
        change = float("inf")
        while iterations < max_property_iterations:
            iterations += 1
            gamma_arr = blend_burned_unburned(t_arr, x_b, mixture, burned_mixture).gamma
            dq_net, dq_gross = first_law(gamma_arr)
            x_b_new = _mass_fraction_burned(theta, dq_gross, w_start, w_end)
            change = float(np.max(np.abs(x_b_new - x_b)))
            x_b = x_b_new
            if change <= property_tolerance:
                converged = True
                break
        if not converged:
            raise PhysicalBoundsError(
                f"x_b-dependent gamma did not converge within {max_property_iterations} "
                f"iterations (last max|dx_b| = {change:.3e})."
            )
        gamma_arr = blend_burned_unburned(t_arr, x_b, mixture, burned_mixture).gamma
        dq_net, dq_gross = first_law(gamma_arr)
        x_b_result = x_b

    net_anti = CubicSpline(theta, dq_net, bc_type="not-a-knot").antiderivative()
    gross_anti = CubicSpline(theta, dq_gross, bc_type="not-a-knot").antiderivative()
    cum_net = np.asarray(net_anti(theta) - net_anti(theta[0]), dtype=np.float64)
    cum_gross = np.asarray(gross_anti(theta) - gross_anti(theta[0]), dtype=np.float64)

    return HeatReleaseResult(
        crank_angle_deg=theta,
        gamma=gamma_arr,
        dp_dtheta_pa_per_deg=dp_dtheta,
        dv_dtheta_m3_per_deg=dv_dtheta,
        dq_net_j_per_deg=dq_net,
        dq_ht_j_per_deg=dq_ht,
        dq_gross_j_per_deg=dq_gross,
        cumulative_net_j=cum_net,
        cumulative_gross_j=cum_gross,
        mass_fraction_burned=x_b_result,
        property_iterations=iterations,
    )


def _mass_fraction_burned(
    theta: FloatArray, dq_gross: FloatArray, start_deg: float, end_deg: float
) -> FloatArray:
    """Normalized cumulative gross heat release, clipped to a monotone [0, 1] profile."""
    anti = CubicSpline(theta, dq_gross, bc_type="not-a-knot").antiderivative()
    lo = float(anti(start_deg))
    total = float(anti(end_deg)) - lo
    if total <= 0.0:
        raise PhysicalBoundsError(
            "Gross heat release over combustion_window_deg is <= 0 J; cannot define x_b."
        )
    q = (np.asarray(anti(np.clip(theta, start_deg, end_deg)), dtype=np.float64) - lo) / total
    return np.asarray(np.maximum.accumulate(np.clip(q, 0.0, 1.0)), dtype=np.float64)


def _solve_mfb_crossing(theta_w: FloatArray, mfb: FloatArray, target: float) -> float:
    spline = CubicSpline(theta_w, mfb - target, bc_type="not-a-knot")
    roots = spline.roots(extrapolate=False)
    if roots.size == 0:
        raise ValueError(f"MFB profile never crosses target fraction {target:.2f}.")
    slopes = np.asarray(spline.derivative(1)(roots), dtype=np.float64)
    rising = roots[slopes > 0.0]
    candidates = rising if rising.size > 0 else roots
    idx_approx = int(np.argmin(np.abs(mfb - target)))
    best = int(np.argmin(np.abs(candidates - theta_w[idx_approx])))
    return float(candidates[best])


def extract_burn_fractions(
    crank_angle_deg: FloatArray,
    dq_dtheta_j_per_deg: FloatArray,
    *,
    start_deg: float = -40.0,
    end_deg: float = 90.0,
) -> BurnFractionMetrics:
    """Extract sub-degree CA10, CA50, CA90, and burn duration via cubic splines.

    Integrates $dQ/d\\theta$ across ``[start_deg, end_deg]`` using exact cubic-spline
    antiderivatives and resolves continuous roots of $x_b(\\theta) - f = 0$ for
    $f \\in \\{0.10, 0.50, 0.90\\}$ (Heywood, 2018, Ch. 9).

    Parameters
    ----------
    crank_angle_deg : NDArray[np.float64]
        Crank angle $\\theta$ [deg].
    dq_dtheta_j_per_deg : NDArray[np.float64]
        Net or gross heat-release rate $dQ/d\\theta$ [$\\text{J/deg}$].
    start_deg : float, default=-40.0
        Start of combustion integration window [deg].
    end_deg : float, default=90.0
        End of combustion integration window [deg].

    Returns
    -------
    BurnFractionMetrics
        Sub-degree CA10, CA50, CA90, $\\Delta\\theta_{10\\text{--}90}$, and MFB curve.
    """
    theta = np.asarray(crank_angle_deg, dtype=np.float64)
    dq = np.asarray(dq_dtheta_j_per_deg, dtype=np.float64)
    if theta.ndim != 1 or theta.shape != dq.shape:
        raise ValueError("crank_angle_deg and dq_dtheta_j_per_deg must be 1D arrays of identical length.")
    if start_deg >= end_deg:
        raise ValueError("start_deg must be strictly less than end_deg.")

    mask = (theta >= start_deg) & (theta <= end_deg)
    theta_w = theta[mask]
    dq_w = dq[mask]
    if theta_w.size < 5:
        raise ValueError("Combustion window [start_deg, end_deg] contains fewer than 5 points.")

    anti = CubicSpline(theta_w, dq_w, bc_type="not-a-knot").antiderivative()
    q_cum = np.asarray(anti(theta_w) - anti(theta_w[0]), dtype=np.float64)
    total_q = float(q_cum[-1])
    if total_q <= 0.0:
        raise PhysicalBoundsError("Total integrated heat release over combustion window is <= 0 J.")

    mfb = q_cum / total_q
    ca10 = _solve_mfb_crossing(theta_w, mfb, 0.10)
    ca50 = _solve_mfb_crossing(theta_w, mfb, 0.50)
    ca90 = _solve_mfb_crossing(theta_w, mfb, 0.90)
    if not (ca10 < ca50 < ca90):
        raise PhysicalBoundsError(
            f"Non-monotonic burn-fraction crossings detected: CA10={ca10:.3f}, "
            f"CA50={ca50:.3f}, CA90={ca90:.3f} deg."
        )

    return BurnFractionMetrics(
        ca10_deg=ca10,
        ca50_deg=ca50,
        ca90_deg=ca90,
        burn_duration_10_90_deg=ca90 - ca10,
        total_heat_j=total_q,
        crank_angle_deg=theta_w,
        mfb=mfb,
    )


def compute_imep(
    crank_angle_deg: FloatArray,
    pressure_pa: FloatArray,
    volume_m3: FloatArray,
    *,
    displaced_volume_m3: float | None = None,
    dv_dtheta_m3_per_deg: FloatArray | None = None,
    gross_window_deg: tuple[float, float] = (-180.0, 180.0),
) -> IMEPResult:
    """Compute Gross & Net Indicated Mean Effective Pressure (Heywood, 2018, Eq. 2.14).

    $$\\text{IMEP} = \\frac{1}{V_d} \\int_{\\theta_1}^{\\theta_2} P(\\theta) \\frac{dV}{d\\theta}\\,d\\theta$$
    """
    theta = np.asarray(crank_angle_deg, dtype=np.float64)
    p_pa = np.asarray(pressure_pa, dtype=np.float64)
    v_m3 = np.asarray(volume_m3, dtype=np.float64)
    if theta.ndim != 1 or not (theta.shape == p_pa.shape == v_m3.shape):
        raise ValueError("crank_angle_deg, pressure_pa, and volume_m3 must be 1D arrays of identical length.")

    v_d = (
        displaced_volume_m3
        if displaced_volume_m3 is not None
        else float(np.max(v_m3) - np.min(v_m3))
    )
    if v_d <= 0.0:
        raise PhysicalBoundsError("Displaced volume V_d must be > 0 m^3.")

    dv_dtheta = (
        np.asarray(dv_dtheta_m3_per_deg, dtype=np.float64)
        if dv_dtheta_m3_per_deg is not None
        else compute_derivative(theta, v_m3)
    )
    integrand_spline = CubicSpline(theta, p_pa * dv_dtheta, bc_type="not-a-knot")

    work_net_j = float(integrand_spline.integrate(float(theta[0]), float(theta[-1])))
    g_start = max(float(theta[0]), gross_window_deg[0])
    g_end = min(float(theta[-1]), gross_window_deg[1])
    work_gross_j = (
        float(integrand_spline.integrate(g_start, g_end)) if g_end > g_start else work_net_j
    )

    imep_gross = work_gross_j / v_d
    imep_net = work_net_j / v_d
    return IMEPResult(
        imep_gross_pa=imep_gross,
        imep_net_pa=imep_net,
        pmep_pa=imep_gross - imep_net,
        indicated_work_gross_j=work_gross_j,
        indicated_work_net_j=work_net_j,
    )


def compute_knock_metrics(
    crank_angle_deg: FloatArray,
    pressure_pa: FloatArray,
    temperature_k: FloatArray,
    rpm: float,
    *,
    band_hz: tuple[float, float] = (4_000.0, 20_000.0),
    mixture: NASAGasMixture = DIESEL_BURNED_PRODUCTS,
    use_prr_correlation: bool = False,
    beta_s: float = 5.0e-5,
) -> KnockMetrics:
    """Compute MAPO [Pa] and Eng (2002) Ringing Intensity $RI$ [$\\text{W/m}^2$].

    Implements Eng (2002, SAE Paper 2002-01-2859, Eq. 5 & Eq. 7):

    $$RI = \\frac{\\Delta P_{\\max}^2}{2 \\gamma_{\\max} P_{\\max}} \\sqrt{\\gamma_{\\max} R_{\\text{spec}} T_{\\max}}$$

    Isolates acoustic oscillations $\\Delta P_{\\text{knock}}$ via zero-phase Butterworth
    bandpass filtering when Nyquist sampling permits ($f_s > 2 f_{\\text{high}}$),
    or applies Eng's pressure-rise correlation $\\Delta P_{\\max} = \\beta (dP/dt)_{\\max}$
    ($\\beta = 0.05\\,\\text{ms}$) when ``use_prr_correlation=True`` or sampling is coarse.

    Scope: MAPO needs a *local* pressure signal (a transducer or a CFD monitor point).
    A spatially averaged CFD cylinder pressure, such as the one in ANSYS Forte
    ``thermo.csv``, averages out the acoustic modes, so its MAPO and bandpass-based
    ringing intensity are not meaningful.
    """
    if rpm <= 0.0:
        raise PhysicalBoundsError("Engine speed rpm must be > 0.")
    if beta_s < 0.0:
        raise PhysicalBoundsError("Correlation constant beta_s must be >= 0 s.")

    theta = np.asarray(crank_angle_deg, dtype=np.float64)
    p_pa = np.asarray(pressure_pa, dtype=np.float64)
    t_k = np.asarray(temperature_k, dtype=np.float64)
    if theta.ndim != 1 or not (theta.shape == p_pa.shape == t_k.shape):
        raise ValueError("crank_angle_deg, pressure_pa, and temperature_k must match in shape.")
    if not (np.all(p_pa > 0.0) and np.all(t_k > 0.0)):
        raise PhysicalBoundsError("pressure_pa and temperature_k must be strictly positive.")

    dp_dtheta = compute_derivative(theta, p_pa)
    max_dp_dtheta = float(np.max(dp_dtheta))

    dtheta = np.diff(theta)
    is_uniform = theta.size >= 5 and np.allclose(dtheta, dtheta[0], rtol=1e-4, atol=1e-9)
    dtheta_mean = float(np.mean(dtheta))
    fs_hz = 6.0 * rpm / dtheta_mean

    knock_trace: FloatArray | None = None
    if not use_prr_correlation and is_uniform and fs_hz > 2.0 * band_hz[1]:
        knock_trace = butterworth_zero_phase(p_pa, fs=fs_hz, cutoff=band_hz, btype="bandpass")
        mapo_pa = float(np.max(np.abs(knock_trace)))
    else:
        max_dp_dt = 6.0 * rpm * max_dp_dtheta
        mapo_pa = beta_s * max(max_dp_dt, 0.0)

    p_max = float(np.max(p_pa))
    t_max = float(np.max(t_k))
    gamma_max = float(mixture.gamma(np.array([t_max], dtype=np.float64))[0])
    speed_of_sound = float(np.sqrt(gamma_max * mixture.specific_gas_constant * t_max))
    ri_w_m2 = ((mapo_pa**2) / (2.0 * gamma_max * p_max)) * speed_of_sound

    return KnockMetrics(
        mapo_pa=mapo_pa,
        ringing_intensity_w_m2=ri_w_m2,
        max_dp_dtheta_pa_per_deg=max_dp_dtheta,
        knock_pressure_pa=knock_trace,
    )
