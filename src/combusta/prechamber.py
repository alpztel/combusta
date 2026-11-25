"""Two-volume analysis of pre-chamber (turbulent jet ignition) engines.

The pre-chamber (index ``pc``, constant volume $V_{pc}$) and the main chamber (``mc``,
volume $V_{mc}(\\theta)$) are treated as two single-zone open systems connected by the
orifices. Pre-chamber jet ignition systems are reviewed by Toulson, Schock and Attard
(2010).

**Orifice flow.** Quasi-steady isentropic flow of an ideal gas through a restriction
with discharge coefficient $C_D$ and total flow area $A$ (Heywood, 2018, App. C). With
upstream state $(p_0, T_0)$ and downstream pressure $p_T$,

$$\\dot m = \\frac{C_D A p_0}{\\sqrt{R T_0}} \\left(\\frac{p_T}{p_0}\\right)^{1/\\gamma}
  \\left\\{\\frac{2\\gamma}{\\gamma - 1}\\left[1 - \\left(\\frac{p_T}{p_0}\\right)^{(\\gamma-1)/\\gamma}\\right]\\right\\}^{1/2}$$

and, for $p_T/p_0 \\le [2/(\\gamma+1)]^{\\gamma/(\\gamma-1)}$ (choked flow),

$$\\dot m = \\frac{C_D A p_0}{\\sqrt{R T_0}}\\,\\gamma^{1/2}
  \\left(\\frac{2}{\\gamma + 1}\\right)^{(\\gamma+1)/(2(\\gamma-1))}.$$

**Energy balance.** For a zone with mass $m$, frozen composition, $u(T)$, and a net
inflow $dm$ carrying specific enthalpy $h_f$ (the upstream zone's enthalpy), the
First Law with $pV = mRT$ gives the apparent heat release

$$dQ = \\frac{V\\,dp}{\\gamma - 1} + \\frac{\\gamma\\,p\\,dV}{\\gamma - 1}
  + \\left(u - c_v T - h_f\\right) dm .$$

It follows from $dU = m c_v\\,dT + u\\,dm$ and $m\\,dT = (V dp + p dV - RT\\,dm)/R$. For
outflow ($h_f = h = u + RT$) the bracket is $-c_p T$; in a rigid adiabatic zone this
reproduces isentropic blowdown, $dp/p = \\gamma\\,dm/m$. The gross heat release adds
the wall heat loss of the zone. Mass exchange is taken from traced zone masses (or
temperatures) when available, otherwise from the orifice model integrated over the
cycle.

Assumptions: each zone has uniform state and frozen composition (no burned/unburned
blending), the gas entering a zone takes on that zone's gas model, and the orifice
flow is quasi-steady. With hydrogen the frozen composition biases the absolute heat
release (see README, Section 3); timing quantities are less sensitive.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Literal, TypeAlias

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import solve_ivp, trapezoid
from scipy.interpolate import CubicSpline

from combusta.exceptions import PhysicalBoundsError
from combusta.signal.filtering import compute_derivative
from combusta.thermo.heat_release import extract_burn_fractions
from combusta.thermo.nasa_polynomials import NASAGasMixture

FloatArray: TypeAlias = NDArray[np.float64]
MassSource: TypeAlias = Literal["trace", "orifice"]


def critical_pressure_ratio(gamma: ArrayLike) -> FloatArray:
    """Critical downstream-to-upstream pressure ratio $[2/(\\gamma+1)]^{\\gamma/(\\gamma-1)}$ [-]."""
    g = np.asarray(gamma, dtype=np.float64)
    return np.asarray((2.0 / (g + 1.0)) ** (g / (g - 1.0)), dtype=np.float64)


def orifice_mass_flow(
    p_upstream_pa: ArrayLike,
    t_upstream_k: ArrayLike,
    p_downstream_pa: ArrayLike,
    *,
    gamma: ArrayLike,
    r_specific: float,
    flow_area_m2: float,
    discharge_coefficient: float,
) -> FloatArray:
    """Quasi-steady isentropic mass flow through a restriction [kg/s] (Heywood, 2018, App. C).

    Requires $p_\\mathrm{downstream} \\le p_\\mathrm{upstream}$; the flow is choked when the
    pressure ratio is at or below :func:`critical_pressure_ratio`.
    """
    p0 = np.asarray(p_upstream_pa, dtype=np.float64)
    t0 = np.asarray(t_upstream_k, dtype=np.float64)
    pt = np.asarray(p_downstream_pa, dtype=np.float64)
    g = np.asarray(gamma, dtype=np.float64)
    if np.any(p0 <= 0.0) or np.any(t0 <= 0.0) or np.any(pt <= 0.0):
        raise PhysicalBoundsError("Pressures and temperature must be > 0.")
    if np.any(pt > p0 * (1.0 + 1e-12)):
        raise ValueError("Downstream pressure exceeds upstream pressure; swap the states.")
    if np.any(g <= 1.0) or r_specific <= 0.0 or flow_area_m2 <= 0.0 or discharge_coefficient <= 0.0:
        raise PhysicalBoundsError("gamma > 1, R > 0, area > 0, and C_D > 0 are required.")
    prefactor = discharge_coefficient * flow_area_m2 * p0 / np.sqrt(r_specific * t0)
    return np.asarray(prefactor * _flow_function(np.minimum(pt / p0, 1.0), g), dtype=np.float64)


@dataclass(frozen=True, slots=True)
class PreChamberGeometry:
    """Pre-chamber volume and orifices.

    Parameters
    ----------
    volume_m3 : float
        Pre-chamber volume $V_{pc}$ [m^3], including the orifices.
    orifice_count : int
        Number of orifices.
    orifice_diameter_m : float
        Orifice diameter [m].
    discharge_coefficient : float
        Discharge coefficient $C_D$ [-] of one orifice, in (0, 1]. It depends on the
        orifice geometry and flow direction and should come from measurement, CFD, or
        :func:`estimate_discharge_coefficient`; Combusta provides no default.
    """

    volume_m3: float
    orifice_count: int
    orifice_diameter_m: float
    discharge_coefficient: float

    def __post_init__(self) -> None:
        if not (self.volume_m3 > 0.0 and self.orifice_diameter_m > 0.0):
            raise PhysicalBoundsError("Pre-chamber volume and orifice diameter must be > 0.")
        if self.orifice_count < 1:
            raise PhysicalBoundsError("orifice_count must be >= 1.")
        if not 0.0 < self.discharge_coefficient <= 1.0:
            raise PhysicalBoundsError("discharge_coefficient must lie in (0, 1].")

    @property
    def flow_area_m2(self) -> float:
        """Total geometric orifice area $A = n \\pi d^2 / 4$ [m^2]."""
        return self.orifice_count * math.pi * self.orifice_diameter_m**2 / 4.0


def transfer_rate(
    p_pc: ArrayLike,
    t_pc: ArrayLike,
    p_mc: ArrayLike,
    t_mc: ArrayLike,
    *,
    mixture_pc: NASAGasMixture,
    mixture_mc: NASAGasMixture,
    geometry: PreChamberGeometry,
    rpm: float,
    discharge_coefficient: float | None = None,
) -> FloatArray:
    """Orifice mass transfer per crank angle [kg/deg], positive from pre-chamber to main chamber.

    The upstream zone is the one with the higher pressure; its temperature, $\\gamma(T)$,
    and gas constant enter :func:`orifice_mass_flow`. $d\\theta/dt = 6N$ deg/s.
    """
    p_a, t_a = np.asarray(p_pc, dtype=np.float64), np.asarray(t_pc, dtype=np.float64)
    p_b, t_b = np.asarray(p_mc, dtype=np.float64), np.asarray(t_mc, dtype=np.float64)
    cd = geometry.discharge_coefficient if discharge_coefficient is None else discharge_coefficient
    out_pc = p_a >= p_b
    p_up = np.where(out_pc, p_a, p_b)
    p_dn = np.where(out_pc, p_b, p_a)
    t_up = np.where(out_pc, t_a, t_b)
    gamma = np.where(out_pc, mixture_pc.gamma(t_up), mixture_mc.gamma(t_up))
    r_up = np.where(out_pc, mixture_pc.specific_gas_constant, mixture_mc.specific_gas_constant)
    mdot = (
        cd * geometry.flow_area_m2 * p_up / np.sqrt(r_up * t_up)
        * _flow_function(np.minimum(p_dn / p_up, 1.0), gamma)
    )
    return np.asarray(np.where(out_pc, mdot, -mdot) / (6.0 * rpm), dtype=np.float64)


def _flow_function(ratio: FloatArray, gamma: FloatArray) -> FloatArray:
    """Dimensionless isentropic flow function $\\dot m \\sqrt{RT_0} / (C_D A p_0)$."""
    r = np.maximum(ratio, critical_pressure_ratio(gamma))
    return np.asarray(
        r ** (1.0 / gamma) * np.sqrt(np.maximum(2.0 * gamma / (gamma - 1.0) * (1.0 - r ** ((gamma - 1.0) / gamma)), 0.0)),
        dtype=np.float64,
    )


def zone_apparent_heat_release(
    pressure_pa: FloatArray,
    dp_dtheta: FloatArray,
    volume_m3: FloatArray,
    dv_dtheta: FloatArray,
    temperature_k: FloatArray,
    inflow_kg_per_deg: FloatArray,
    inflow_enthalpy_j_per_kg: FloatArray,
    mixture: NASAGasMixture,
) -> FloatArray:
    """Apparent heat-release rate of an open single zone [J/deg] (module docstring)."""
    gamma = mixture.gamma(temperature_k)
    u = mixture.internal_energy(temperature_k)
    cv = mixture.cv(temperature_k)
    return np.asarray(
        volume_m3 * dp_dtheta / (gamma - 1.0)
        + gamma * pressure_pa * dv_dtheta / (gamma - 1.0)
        + (u - cv * temperature_k - inflow_enthalpy_j_per_kg) * inflow_kg_per_deg,
        dtype=np.float64,
    )


@dataclass(frozen=True, slots=True)
class PreChamberResult:
    """Two-volume pre-chamber analysis on a crank-angle grid (SI units).

    Attributes
    ----------
    crank_angle_deg : NDArray[np.float64]
        Crank angle [deg ATDC].
    pressure_difference_pa : NDArray[np.float64]
        $\\Delta p = p_{pc} - p_{mc}$ [Pa].
    transfer_rate_kg_per_deg : NDArray[np.float64]
        Mass transfer from pre-chamber to main chamber [kg/deg] used in the energy
        balance (from the mass trace or the orifice model, see ``mass_source``).
    model_transfer_rate_kg_per_deg : NDArray[np.float64]
        Orifice-model transfer with the geometry's $C_D$ [kg/deg].
    choked : NDArray[np.bool_]
        Orifice flow choked according to the critical pressure ratio.
    pre_chamber_mass_kg, main_chamber_mass_kg : NDArray[np.float64]
        Zone masses [kg].
    pre_chamber_temperature_k, main_chamber_temperature_k : NDArray[np.float64]
        Zone temperatures [K].
    dq_pc_gross_j_per_deg, dq_mc_gross_j_per_deg : NDArray[np.float64]
        Gross heat-release rates of the zones [J/deg].
    dq_pc_apparent_j_per_deg, dq_mc_apparent_j_per_deg : NDArray[np.float64]
        Apparent (net) heat-release rates [J/deg].
    mass_source : {"trace", "orifice"}
        Origin of the mass exchange.
    discharge_coefficient : float
        $C_D$ used for ``model_transfer_rate_kg_per_deg`` [-].
    """

    crank_angle_deg: FloatArray
    pressure_difference_pa: FloatArray
    transfer_rate_kg_per_deg: FloatArray
    model_transfer_rate_kg_per_deg: FloatArray
    choked: NDArray[np.bool_]
    pre_chamber_mass_kg: FloatArray
    main_chamber_mass_kg: FloatArray
    pre_chamber_temperature_k: FloatArray
    main_chamber_temperature_k: FloatArray
    dq_pc_apparent_j_per_deg: FloatArray
    dq_mc_apparent_j_per_deg: FloatArray
    dq_pc_gross_j_per_deg: FloatArray
    dq_mc_gross_j_per_deg: FloatArray
    mass_source: MassSource
    discharge_coefficient: float


def _as_trace(values: ArrayLike, shape: tuple[int, ...], name: str, *, positive: bool) -> FloatArray:
    arr = np.asarray(values, dtype=np.float64)
    if arr.shape != shape:
        raise ValueError(f"{name} must match the crank-angle grid {shape}.")
    if not np.all(np.isfinite(arr)) or (positive and np.any(arr <= 0.0)):
        raise PhysicalBoundsError(f"{name} must be finite" + (" and > 0." if positive else "."))
    return arr


def _integrate_pc_mass(
    theta: FloatArray,
    p_pc: FloatArray,
    p_mc: FloatArray,
    v_mc: FloatArray,
    geometry: PreChamberGeometry,
    rpm: float,
    mixture_pc: NASAGasMixture,
    mixture_mc: NASAGasMixture,
    trapped_mass_kg: float,
) -> FloatArray:
    """Integrate the pre-chamber mass with the orifice model from a uniform initial temperature."""
    r_pc, r_mc = mixture_pc.specific_gas_constant, mixture_mc.specific_gas_constant
    # Uniform initial temperature: m_i = p_i V_i / (R_i T0) with m_pc + m_mc = m_total.
    t0 = (p_pc[0] * geometry.volume_m3 / r_pc + p_mc[0] * v_mc[0] / r_mc) / trapped_mass_kg
    m_pc0 = p_pc[0] * geometry.volume_m3 / (r_pc * t0)
    s_pc, s_mc, s_v = CubicSpline(theta, p_pc), CubicSpline(theta, p_mc), CubicSpline(theta, v_mc)

    def rhs(th: float, y: FloatArray) -> list[float]:
        m_pc = float(y[0])
        m_mc = trapped_mass_kg - m_pc
        if m_pc <= 0.0 or m_mc <= 0.0:
            raise PhysicalBoundsError("Orifice integration produced a non-positive zone mass.")
        pp, pm = float(s_pc(th)), float(s_mc(th))
        tp = pp * geometry.volume_m3 / (m_pc * r_pc)
        tm = pm * float(s_v(th)) / (m_mc * r_mc)
        x = transfer_rate(pp, tp, pm, tm, mixture_pc=mixture_pc, mixture_mc=mixture_mc,
                          geometry=geometry, rpm=rpm)
        return [-float(x)]

    sol = solve_ivp(rhs, (float(theta[0]), float(theta[-1])), [m_pc0], t_eval=theta,
                    method="LSODA", rtol=1e-9, atol=1e-14)
    if not sol.success:
        raise PhysicalBoundsError(f"Orifice mass integration failed: {sol.message}")
    return np.asarray(sol.y[0], dtype=np.float64)


def analyze_pre_chamber(
    crank_angle_deg: ArrayLike,
    pre_chamber_pressure_pa: ArrayLike,
    main_chamber_pressure_pa: ArrayLike,
    main_chamber_volume_m3: ArrayLike,
    geometry: PreChamberGeometry,
    *,
    rpm: float,
    mixture_pc: NASAGasMixture,
    mixture_mc: NASAGasMixture,
    pre_chamber_mass_kg: ArrayLike | None = None,
    main_chamber_mass_kg: ArrayLike | None = None,
    pre_chamber_temperature_k: ArrayLike | None = None,
    main_chamber_temperature_k: ArrayLike | None = None,
    trapped_mass_kg: float | None = None,
    wall_heat_loss_pc_j_per_deg: ArrayLike | None = None,
    wall_heat_loss_mc_j_per_deg: ArrayLike | None = None,
    main_chamber_dv_dtheta_m3_per_deg: ArrayLike | None = None,
) -> PreChamberResult:
    """Heat release of pre-chamber and main chamber with orifice mass exchange.

    Zone masses are taken, in this order of preference, from

    1. ``pre_chamber_mass_kg`` and ``main_chamber_mass_kg`` traces;
    2. zone temperature traces, via $m = pV/(RT)$;
    3. the orifice model integrated from a uniform initial temperature, which needs
       ``trapped_mass_kg`` (total mass of both zones).

    Wall heat losses [J/deg] are positive when heat leaves the gas.

    Raises
    ------
    ValueError
        Missing or inconsistent inputs.
    PhysicalBoundsError
        Non-physical inputs or a failed mass integration.
    """
    theta = np.asarray(crank_angle_deg, dtype=np.float64)
    if theta.ndim != 1 or theta.size < 5 or np.any(np.diff(theta) <= 0.0):
        raise ValueError("crank_angle_deg must be strictly increasing with at least 5 samples.")
    if not (np.isfinite(rpm) and rpm > 0.0):
        raise PhysicalBoundsError("rpm must be finite and > 0.")
    shape = theta.shape
    p_pc = _as_trace(pre_chamber_pressure_pa, shape, "pre_chamber_pressure_pa", positive=True)
    p_mc = _as_trace(main_chamber_pressure_pa, shape, "main_chamber_pressure_pa", positive=True)
    v_mc = _as_trace(main_chamber_volume_m3, shape, "main_chamber_volume_m3", positive=True)
    v_pc = geometry.volume_m3
    r_pc, r_mc = mixture_pc.specific_gas_constant, mixture_mc.specific_gas_constant

    source: MassSource
    if pre_chamber_mass_kg is not None and main_chamber_mass_kg is not None:
        m_pc = _as_trace(pre_chamber_mass_kg, shape, "pre_chamber_mass_kg", positive=True)
        m_mc = _as_trace(main_chamber_mass_kg, shape, "main_chamber_mass_kg", positive=True)
        source = "trace"
    elif pre_chamber_temperature_k is not None and main_chamber_temperature_k is not None:
        t_pc_in = _as_trace(pre_chamber_temperature_k, shape, "pre_chamber_temperature_k", positive=True)
        t_mc_in = _as_trace(main_chamber_temperature_k, shape, "main_chamber_temperature_k", positive=True)
        m_pc = p_pc * v_pc / (r_pc * t_pc_in)
        m_mc = p_mc * v_mc / (r_mc * t_mc_in)
        source = "trace"
    elif trapped_mass_kg is not None:
        if not (np.isfinite(trapped_mass_kg) and trapped_mass_kg > 0.0):
            raise PhysicalBoundsError("trapped_mass_kg must be finite and > 0.")
        m_pc = _integrate_pc_mass(theta, p_pc, p_mc, v_mc, geometry, rpm, mixture_pc, mixture_mc,
                                  trapped_mass_kg)
        m_mc = trapped_mass_kg - m_pc
        source = "orifice"
    else:
        raise ValueError(
            "Provide zone mass traces, zone temperature traces, or trapped_mass_kg for the orifice model."
        )

    t_pc = p_pc * v_pc / (m_pc * r_pc)
    t_mc = p_mc * v_mc / (m_mc * r_mc)
    model_x = transfer_rate(p_pc, t_pc, p_mc, t_mc, mixture_pc=mixture_pc, mixture_mc=mixture_mc,
                            geometry=geometry, rpm=rpm)
    x = model_x if source == "orifice" else -compute_derivative(theta, m_pc)

    gamma_up = np.where(p_pc >= p_mc, mixture_pc.gamma(t_pc), mixture_mc.gamma(t_mc))
    ratio = np.minimum(p_pc, p_mc) / np.maximum(p_pc, p_mc)
    choked = np.asarray(ratio <= critical_pressure_ratio(gamma_up))

    # Enthalpy carried by the flow: that of the zone it leaves.
    h_pc, h_mc = mixture_pc.enthalpy(t_pc), mixture_mc.enthalpy(t_mc)
    h_flow = np.where(x >= 0.0, h_pc, h_mc)
    dp_pc = compute_derivative(theta, p_pc)
    dp_mc = compute_derivative(theta, p_mc)
    dv_mc = (
        compute_derivative(theta, v_mc)
        if main_chamber_dv_dtheta_m3_per_deg is None
        else _as_trace(main_chamber_dv_dtheta_m3_per_deg, shape, "main_chamber_dv_dtheta_m3_per_deg", positive=False)
    )
    zeros = np.zeros_like(theta)
    dq_pc = zone_apparent_heat_release(p_pc, dp_pc, np.full_like(theta, v_pc), zeros, t_pc, -x, h_flow, mixture_pc)
    dq_mc = zone_apparent_heat_release(p_mc, dp_mc, v_mc, dv_mc, t_mc, x, h_flow, mixture_mc)
    q_w_pc = zeros if wall_heat_loss_pc_j_per_deg is None else _as_trace(
        wall_heat_loss_pc_j_per_deg, shape, "wall_heat_loss_pc_j_per_deg", positive=False)
    q_w_mc = zeros if wall_heat_loss_mc_j_per_deg is None else _as_trace(
        wall_heat_loss_mc_j_per_deg, shape, "wall_heat_loss_mc_j_per_deg", positive=False)

    return PreChamberResult(
        crank_angle_deg=theta,
        pressure_difference_pa=p_pc - p_mc,
        transfer_rate_kg_per_deg=np.asarray(x, dtype=np.float64),
        model_transfer_rate_kg_per_deg=model_x,
        choked=choked,
        pre_chamber_mass_kg=np.asarray(m_pc, dtype=np.float64),
        main_chamber_mass_kg=np.asarray(m_mc, dtype=np.float64),
        pre_chamber_temperature_k=np.asarray(t_pc, dtype=np.float64),
        main_chamber_temperature_k=np.asarray(t_mc, dtype=np.float64),
        dq_pc_apparent_j_per_deg=dq_pc,
        dq_mc_apparent_j_per_deg=dq_mc,
        dq_pc_gross_j_per_deg=dq_pc + q_w_pc,
        dq_mc_gross_j_per_deg=dq_mc + q_w_mc,
        mass_source=source,
        discharge_coefficient=geometry.discharge_coefficient,
    )


def estimate_discharge_coefficient(
    result: PreChamberResult,
    *,
    min_relative_pressure_difference: float = 0.05,
) -> float:
    """Least-squares $C_D$ that best maps the orifice model onto the traced mass transfer.

    $C_D = \\sum x_t x_1 / \\sum x_1^2$, where $x_t$ is the transfer from the zone-mass
    trace and $x_1$ the orifice model with $C_D = 1$, over samples with
    $|\\Delta p| \\ge$ ``min_relative_pressure_difference`` $\\cdot \\max|\\Delta p|$ (near
    $\\Delta p = 0$ the traced transfer is dominated by differentiation noise). It is an
    effective, cycle-averaged value for both flow directions.
    """
    if result.mass_source != "trace":
        raise ValueError("A discharge coefficient can only be estimated from traced zone masses.")
    dp = np.abs(result.pressure_difference_pa)
    mask = dp >= min_relative_pressure_difference * float(np.max(dp))
    # The model rate is proportional to C_D: fit x_t = k * x_model, then C_D = k * C_D(model).
    x_t = result.transfer_rate_kg_per_deg[mask]
    x_m = result.model_transfer_rate_kg_per_deg[mask]
    denom = float(np.dot(x_m, x_m))
    if denom == 0.0:
        raise ValueError("The orifice model predicts no flow.")
    return float(np.dot(x_t, x_m) / denom) * result.discharge_coefficient


@dataclass(frozen=True, slots=True)
class PreChamberMetrics:
    """Scalar pre-chamber metrics (SI units; ``nan`` when not defined).

    Attributes
    ----------
    dp_max_pa, theta_dp_max_deg : float
        Largest $p_{pc} - p_{mc}$ [Pa] and its crank angle [deg ATDC] (jet driving pressure).
    dp_min_pa, theta_dp_min_deg : float
        Most negative $\\Delta p$ [Pa] and its crank angle (flow into the pre-chamber).
    jet_start_deg, jet_end_deg : float
        Zero crossings of $\\Delta p$ around its maximum: start and end of outflow from
        the pre-chamber [deg ATDC] (linear interpolation between samples).
    jet_mass_kg : float
        Mass leaving the pre-chamber between jet start and end [kg].
    jet_mass_fraction : float
        ``jet_mass_kg`` relative to the pre-chamber mass at jet start [-].
    peak_jet_rate_kg_per_s : float
        Largest outflow rate [kg/s].
    choked_fraction : float
        Fraction of the jet interval with choked orifice flow [-].
    q_pc_gross_j, q_mc_gross_j : float
        Cumulative gross heat release of the zones over the window [J].
    pc_energy_share : float
        $Q_{pc} / (Q_{pc} + Q_{mc})$ [-].
    pc_ca50_deg, mc_ca10_deg, mc_ca50_deg, mc_ca90_deg : float
        Burn fractions of the pre-chamber and main-chamber heat release [deg ATDC].
    """

    dp_max_pa: float
    theta_dp_max_deg: float
    dp_min_pa: float
    theta_dp_min_deg: float
    jet_start_deg: float
    jet_end_deg: float
    jet_mass_kg: float
    jet_mass_fraction: float
    peak_jet_rate_kg_per_s: float
    choked_fraction: float
    q_pc_gross_j: float
    q_mc_gross_j: float
    pc_energy_share: float
    pc_ca50_deg: float
    mc_ca10_deg: float
    mc_ca50_deg: float
    mc_ca90_deg: float

    def as_dict(self) -> dict[str, float]:
        """Metrics as a plain dictionary (keys as attribute names)."""
        return {name: float(getattr(self, name)) for name in self.__slots__}  # type: ignore[attr-defined]


def _crossing(theta: FloatArray, y: FloatArray, i: int) -> float:
    """Linear zero crossing between samples ``i`` and ``i + 1``."""
    return float(theta[i] - y[i] * (theta[i + 1] - theta[i]) / (y[i + 1] - y[i]))


def pre_chamber_metrics(
    result: PreChamberResult,
    *,
    rpm: float,
    window_deg: tuple[float, float] = (-40.0, 90.0),
) -> PreChamberMetrics:
    """Jet timing, transferred mass, and heat-release split from a :class:`PreChamberResult`.

    The jet interval is the positive-$\\Delta p$ interval that contains the largest
    $\\Delta p$ inside ``window_deg``.
    """
    theta = result.crank_angle_deg
    dp = result.pressure_difference_pa
    x = result.transfer_rate_kg_per_deg
    lo, hi = max(window_deg[0], float(theta[0])), min(window_deg[1], float(theta[-1]))
    win = (theta >= lo) & (theta <= hi)
    if int(np.count_nonzero(win)) < 5:
        raise ValueError("window_deg does not overlap the trace.")
    idx = np.flatnonzero(win)
    i_max = int(idx[np.argmax(dp[win])])
    i_min = int(idx[np.argmin(dp[win])])
    nan = math.nan

    jet_start = jet_end = jet_mass = jet_fraction = choked_fraction = nan
    if dp[i_max] > 0.0:
        i0 = i_max
        while i0 > 0 and dp[i0 - 1] > 0.0:
            i0 -= 1
        i1 = i_max
        while i1 < theta.size - 1 and dp[i1 + 1] > 0.0:
            i1 += 1
        jet_start = _crossing(theta, dp, i0 - 1) if i0 > 0 else float(theta[0])
        jet_end = _crossing(theta, dp, i1) if i1 < theta.size - 1 else float(theta[-1])
        seg = slice(i0, i1 + 1)
        jet_mass = float(trapezoid(np.maximum(x[seg], 0.0), theta[seg]))
        jet_fraction = jet_mass / float(result.pre_chamber_mass_kg[i0])
        choked_fraction = float(np.mean(result.choked[seg]))

    def cumulative(dq: FloatArray) -> float:
        return float(trapezoid(dq[win], theta[win]))

    q_pc = cumulative(result.dq_pc_gross_j_per_deg)
    q_mc = cumulative(result.dq_mc_gross_j_per_deg)

    def burn(dq: FloatArray) -> tuple[float, float, float]:
        try:
            b = extract_burn_fractions(theta, dq, start_deg=lo, end_deg=hi)
        except (PhysicalBoundsError, ValueError):
            return nan, nan, nan
        return b.ca10_deg, b.ca50_deg, b.ca90_deg

    _, pc50, _ = burn(result.dq_pc_gross_j_per_deg)
    mc10, mc50, mc90 = burn(result.dq_mc_gross_j_per_deg)
    total = q_pc + q_mc
    return PreChamberMetrics(
        dp_max_pa=float(dp[i_max]),
        theta_dp_max_deg=float(theta[i_max]),
        dp_min_pa=float(dp[i_min]),
        theta_dp_min_deg=float(theta[i_min]),
        jet_start_deg=jet_start,
        jet_end_deg=jet_end,
        jet_mass_kg=jet_mass,
        jet_mass_fraction=jet_fraction,
        peak_jet_rate_kg_per_s=float(np.max(np.maximum(x[win], 0.0))) * 6.0 * rpm,
        choked_fraction=choked_fraction,
        q_pc_gross_j=q_pc,
        q_mc_gross_j=q_mc,
        pc_energy_share=q_pc / total if total > 0.0 else nan,
        pc_ca50_deg=pc50,
        mc_ca10_deg=mc10,
        mc_ca50_deg=mc50,
        mc_ca90_deg=mc90,
    )
