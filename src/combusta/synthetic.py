"""Manufactured engine cycles for exercising an analysis pipeline before real data exist.

:func:`manufactured_cycle` integrates the closed-cycle single-zone pressure equation

$$\\frac{dP}{d\\theta} = \\frac{(\\gamma - 1)\\left(\\frac{dQ_g}{d\\theta} -
  \\frac{dQ_w}{d\\theta}\\right) - \\gamma P \\frac{dV}{d\\theta}}{V}$$

forward in crank angle, with $\\gamma(T, x_b)$ from the unburned and complete-combustion
mixtures of the fuel (:func:`combusta.thermo.fuels.blend_burned_unburned`), a Wiebe
burn law $x_b(\\theta)$, and an optional prescribed wall heat loss. Because the heat
release that produced the pressure trace is known exactly, the inverse analysis of
:mod:`combusta.compare` must recover it; this is code verification by the method of
manufactured solutions (Roache, 1998), not a model of a real engine.

:func:`write_case_directory` writes such a case as ``case.toml`` plus a trace CSV so that
the file readers, case loader, and report can be tested end to end. Every case produced
here has ``source = "synthetic"``.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
from typing import Literal

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

from combusta.case import CycleTrace, EngineCase, OperatingPoint
from combusta.exceptions import PhysicalBoundsError
from combusta.prechamber import PreChamberGeometry, transfer_rate
from combusta.thermo.fuels import Fuel, blend_burned_unburned, make_product_mixture, make_reactant_mixture
from combusta.thermo.kinematics import SliderCrankGeometry

FloatArray = NDArray[np.float64]

WIEBE_EFFICIENCY_A: float = 6.908
"""Wiebe efficiency parameter $a = -\\ln(0.001)$: 99.9 % burned at the end of combustion."""


@dataclass(frozen=True, slots=True)
class ManufacturedCycle:
    """A manufactured cycle and the exact heat release that produced it.

    Attributes
    ----------
    case : EngineCase
        Synthetic case (``source = "synthetic"``) with pressure, volume, temperature,
        and wall heat-loss traces.
    gross_heat_release_j_per_deg : NDArray[np.float64]
        Prescribed gross heat-release rate $dQ_g/d\\theta$ [J/deg].
    mass_fraction_burned : NDArray[np.float64]
        Prescribed $x_b(\\theta)$ [-], normalized to reach 1 at the end of combustion.
    total_heat_release_j : float
        $\\int dQ_g$ [J] = fuel mass times LHV times ``combustion_efficiency``.
    """

    case: EngineCase
    gross_heat_release_j_per_deg: FloatArray
    mass_fraction_burned: FloatArray
    total_heat_release_j: float


def _wiebe(theta: FloatArray, start: float, duration: float, m: float) -> FloatArray:
    u = np.clip((theta - start) / duration, 0.0, 1.0)
    return (1.0 - np.exp(-WIEBE_EFFICIENCY_A * u ** (m + 1.0))) / (1.0 - math.exp(-WIEBE_EFFICIENCY_A))


def manufactured_cycle(
    name: str,
    geometry: SliderCrankGeometry,
    fuel: Fuel,
    *,
    lambda_air: float,
    rpm: float,
    p_ivc_pa: float = 1.0e5,
    t_ivc_k: float = 330.0,
    combustion_start_deg: float = -15.0,
    combustion_duration_deg: float = 45.0,
    wiebe_form_factor: float = 2.0,
    combustion_efficiency: float = 1.0,
    wall_heat_loss_j: float = 0.0,
    theta_range_deg: tuple[float, float] = (-180.0, 180.0),
    resolution_deg: float = 0.5,
    ignition_timing_deg: float | None = None,
) -> ManufacturedCycle:
    """Integrate a single-zone cycle with a prescribed Wiebe heat release.

    Parameters
    ----------
    name : str
        Case name.
    geometry : SliderCrankGeometry
        Cylinder geometry.
    fuel : Fuel
        Fuel; the charge is fuel and air at ``lambda_air`` (no residual gas).
    lambda_air : float
        Air-excess ratio $\\lambda$ [-].
    rpm : float
        Engine speed [rev/min] (only stored; the closed cycle does not depend on it).
    p_ivc_pa, t_ivc_k : float
        Pressure [Pa] and temperature [K] at the first crank angle.
    combustion_start_deg, combustion_duration_deg : float
        Wiebe start [deg ATDC] and duration [deg].
    wiebe_form_factor : float, default=2.0
        Wiebe form factor $m$.
    combustion_efficiency : float, default=1.0
        Fraction of the fuel energy released as heat.
    wall_heat_loss_j : float, default=0.0
        Total wall heat loss [J], released as a Gaussian centred 10 deg after the
        combustion start with a 25 deg standard deviation (a prescribed shape, not a
        heat-transfer model).
    theta_range_deg : tuple[float, float], default=(-180, 180)
        Closed-cycle crank-angle range [deg ATDC].
    resolution_deg : float, default=0.5
        Output sample spacing [deg].
    ignition_timing_deg : float | None, optional
        Spark timing [deg ATDC] stored in the operating point. It does not drive the
        burn law, which is set by ``combustion_start_deg`` and ``combustion_duration_deg``;
        a synthetic spark sweep chooses both consistently (e.g. a fixed delay).
    """
    if not 0.0 < combustion_efficiency <= 1.0:
        raise PhysicalBoundsError("combustion_efficiency must lie in (0, 1].")
    if combustion_duration_deg <= 0.0 or resolution_deg <= 0.0 or wall_heat_loss_j < 0.0:
        raise PhysicalBoundsError("Duration and resolution must be > 0 and wall heat loss >= 0.")
    theta = np.arange(theta_range_deg[0], theta_range_deg[1] + 0.5 * resolution_deg, resolution_deg)
    kin = geometry.evaluate(theta)
    unburned = make_reactant_mixture(fuel, lambda_air)
    burned = make_product_mixture(fuel, lambda_air)
    mass = p_ivc_pa * float(kin.volume_m3[0]) / (unburned.specific_gas_constant * t_ivc_k)
    fuel_mass = mass / (1.0 + lambda_air * fuel.stoichiometric_afr)
    q_total = combustion_efficiency * fuel_mass * fuel.lhv_j_per_kg

    a, m = WIEBE_EFFICIENCY_A, wiebe_form_factor
    start, dur = combustion_start_deg, combustion_duration_deg
    norm = 1.0 - math.exp(-a)
    w_mu, w_sigma = start + 10.0, 25.0
    w_scale = wall_heat_loss_j / (w_sigma * math.sqrt(2.0 * math.pi))

    def x_b(th: float) -> float:
        u = min(max((th - start) / dur, 0.0), 1.0)
        return (1.0 - math.exp(-a * u ** (m + 1.0))) / norm

    def dq_g(th: float) -> float:
        u = (th - start) / dur
        if u <= 0.0 or u >= 1.0:
            return 0.0
        return q_total * a * (m + 1.0) / dur * u**m * math.exp(-a * u ** (m + 1.0)) / norm

    def dq_w(th: float) -> float:
        return w_scale * math.exp(-0.5 * ((th - w_mu) / w_sigma) ** 2)

    r_u, r_b = unburned.specific_gas_constant, burned.specific_gas_constant

    def rhs(th: float, y: FloatArray) -> list[float]:
        k = geometry.evaluate(np.array([th]))
        v, dv = float(k.volume_m3[0]), float(k.dv_dtheta_m3_per_deg[0])
        xb = x_b(th)
        t = float(y[0]) * v / (mass * ((1.0 - xb) * r_u + xb * r_b))
        g = float(blend_burned_unburned(np.array([t]), np.array([xb]), unburned, burned).gamma[0])
        return [((g - 1.0) * (dq_g(th) - dq_w(th)) - g * float(y[0]) * dv) / v]

    sol = solve_ivp(
        rhs, (float(theta[0]), float(theta[-1])), [p_ivc_pa], t_eval=theta,
        rtol=1e-10, atol=1e-6, max_step=0.1,
    )
    if not sol.success:  # pragma: no cover - integration of a smooth problem
        raise RuntimeError(f"Forward integration failed: {sol.message}")
    pressure = np.asarray(sol.y[0], dtype=np.float64)
    xb_arr = _wiebe(theta, start, dur, m)
    r_arr = (1.0 - xb_arr) * r_u + xb_arr * r_b
    temperature = pressure * kin.volume_m3 / (mass * r_arr)
    dqg = np.array([dq_g(float(th)) for th in theta])
    dqw = w_scale * np.exp(-0.5 * ((theta - w_mu) / w_sigma) ** 2)

    case = EngineCase(
        name=name,
        geometry=geometry,
        operating_point=OperatingPoint(
            rpm=rpm,
            fuel=fuel,
            lambda_air=lambda_air,
            fuel_mass_per_cycle_kg=fuel_mass,
            trapped_mass_kg=mass,
            intake_pressure_pa=p_ivc_pa,
            intake_temperature_k=t_ivc_k,
            ignition_timing_deg=ignition_timing_deg,
        ),
        trace=CycleTrace(
            crank_angle_deg=theta,
            pressure_pa=pressure,
            volume_m3=kin.volume_m3,
            temperature_k=temperature,
            wall_heat_loss_j_per_deg=dqw,
        ),
        source="synthetic",
        metadata={
            "generator": "combusta.synthetic.manufactured_cycle",
            "wiebe": f"start={start:g} deg, duration={dur:g} deg, m={m:g}, a={a:g}",
        },
    )
    return ManufacturedCycle(case=case, gross_heat_release_j_per_deg=dqg,
                             mass_fraction_burned=xb_arr, total_heat_release_j=q_total)


def _toml_value(value: object) -> str:
    if isinstance(value, str):
        return '"' + value.replace("\\", "\\\\").replace('"', '\\"') + '"'
    if isinstance(value, float):
        return repr(value)
    return str(value)


def write_case_directory(
    case: EngineCase,
    directory: str | Path,
    *,
    trace_format: Literal["forte_thermo", "csv"] = "forte_thermo",
    fuel_key: str | None = None,
    wall_heat_loss_sign: float = 1.0,
    uncertainty: dict[str, float] | None = None,
    configuration: str | None = None,
) -> Path:
    """Write ``case.toml`` and a trace CSV for a case; return the case-file path.

    ``trace_format="forte_thermo"`` writes a ``thermo.csv`` with inline units (MPa,
    cm^3, K, J/deg); ``"csv"`` writes ``trace.csv`` with a unit row (bar, no volume or
    temperature), as a test-bench export would. The wall heat-loss column is written
    multiplied by ``wall_heat_loss_sign`` and the case file declares the same sign, so
    loading restores the Combusta convention.

    Parameters
    ----------
    case : EngineCase
        Case to write. Its fuel must be in :data:`combusta.casefile.FUEL_ALIASES`
        unless ``fuel_key`` names one.
    directory : str | Path
        Target directory (created if needed).
    fuel_key : str | None
        Fuel name to write; defaults to ``case.operating_point.fuel.name``.
    configuration : str | None
        Configuration name for grouping sweep cases.
    """
    if wall_heat_loss_sign not in (1.0, -1.0):
        raise ValueError("wall_heat_loss_sign must be +1.0 or -1.0.")
    out = Path(directory)
    out.mkdir(parents=True, exist_ok=True)
    tr = case.trace
    if trace_format == "forte_thermo":
        if tr.volume_m3 is None:
            raise ValueError("forte_thermo format needs a volume trace.")
        cols = ["Crank Angle [deg]", "Pressure [MPa]", "Volume [cm^3]"]
        data = [tr.crank_angle_deg, tr.pressure_pa * 1e-6, tr.volume_m3 * 1e6]
        if tr.temperature_k is not None:
            cols.append("Temperature [K]")
            data.append(tr.temperature_k)
        if tr.wall_heat_loss_j_per_deg is not None:
            cols.append("Wall Heat Loss [J/deg]")
            data.append(wall_heat_loss_sign * tr.wall_heat_loss_j_per_deg)
        trace_name = "thermo.csv"
        lines = [",".join(cols)]
    elif trace_format == "csv":
        cols = ["Crank Angle", "Cylinder Pressure"]
        units = ["deg", "bar"]
        data = [tr.crank_angle_deg, tr.pressure_pa * 1e-5]
        if tr.wall_heat_loss_j_per_deg is not None:
            cols.append("Wall Heat Loss")
            units.append("J/deg")
            data.append(wall_heat_loss_sign * tr.wall_heat_loss_j_per_deg)
        trace_name = "trace.csv"
        lines = [",".join(cols), ",".join(units)]
    else:
        raise ValueError(f"Unknown trace_format {trace_format!r}.")
    matrix = np.column_stack(data)
    lines += [",".join(f"{v:.12g}" for v in row) for row in matrix]
    (out / trace_name).write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")

    geo, op = case.geometry, case.operating_point
    toml = [
        f"name = {_toml_value(case.name)}",
        f"source = {_toml_value(case.source)}",
        *([f"configuration = {_toml_value(configuration)}"] if configuration else []),
        "",
        "[geometry]",
        f"bore_m = {geo.bore_m!r}",
        f"stroke_m = {geo.stroke_m!r}",
        f"conrod_m = {geo.conrod_m!r}",
        f"compression_ratio = {geo.compression_ratio!r}",
        "",
        "[operating_point]",
        f"rpm = {float(op.rpm)!r}",
        f"fuel = {_toml_value(fuel_key or op.fuel.name)}",
    ]
    for key in ("lambda_air", "fuel_mass_per_cycle_kg", "trapped_mass_kg", "ignition_timing_deg",
                "intake_pressure_pa", "intake_temperature_k", "wall_temperature_k"):
        value = getattr(op, key)
        if value is not None:
            toml.append(f"{key} = {float(value)!r}")
    toml += [
        "",
        "[trace]",
        f"file = {_toml_value(trace_name)}",
        f"format = {_toml_value(trace_format)}",
        f"wall_heat_loss_sign = {wall_heat_loss_sign!r}",
    ]
    if uncertainty:
        toml += ["", "[uncertainty]"] + [f"{k} = {float(v)!r}" for k, v in uncertainty.items()]
    if case.metadata:
        toml += ["", "[metadata]"] + [f"{_toml_value(k)} = {_toml_value(v)}" for k, v in case.metadata.items()]
    path = out / "case.toml"
    path.write_text("\n".join(toml) + "\n", encoding="utf-8", newline="\n")
    return path


@dataclass(frozen=True, slots=True)
class ManufacturedPreChamberCycle:
    """A manufactured two-volume pre-chamber cycle with its exact heat release.

    All arrays share ``crank_angle_deg``. Heat-release rates are those prescribed in
    the forward integration [J/deg]; ``transfer_rate_kg_per_deg`` is the orifice mass
    transfer from pre-chamber to main chamber [kg/deg].
    """

    crank_angle_deg: FloatArray
    pre_chamber_pressure_pa: FloatArray
    main_chamber_pressure_pa: FloatArray
    main_chamber_volume_m3: FloatArray
    main_chamber_dv_dtheta_m3_per_deg: FloatArray
    pre_chamber_mass_kg: FloatArray
    main_chamber_mass_kg: FloatArray
    pre_chamber_temperature_k: FloatArray
    main_chamber_temperature_k: FloatArray
    transfer_rate_kg_per_deg: FloatArray
    dq_pc_j_per_deg: FloatArray
    dq_mc_j_per_deg: FloatArray
    q_pc_total_j: float
    q_mc_total_j: float
    trapped_mass_kg: float
    rpm: float


def manufactured_pre_chamber_cycle(
    geometry: SliderCrankGeometry,
    pre_chamber: PreChamberGeometry,
    fuel: Fuel,
    *,
    lambda_air: float,
    rpm: float,
    p_ivc_pa: float = 1.0e5,
    t_ivc_k: float = 330.0,
    pc_energy_fraction: float | None = None,
    pc_start_deg: float = -18.0,
    pc_duration_deg: float = 15.0,
    mc_start_deg: float = -8.0,
    mc_duration_deg: float = 30.0,
    theta_range_deg: tuple[float, float] = (-180.0, 180.0),
    resolution_deg: float = 0.25,
) -> ManufacturedPreChamberCycle:
    """Integrate a two-volume cycle with prescribed Wiebe heat release in each volume.

    The geometry describes the main chamber; the pre-chamber volume is added to it.
    Both volumes use the unburned-charge gas model of ``fuel`` at ``lambda_air`` (frozen
    composition, the same assumption as :func:`combusta.prechamber.analyze_pre_chamber`),
    and exchange mass through the orifice model of :mod:`combusta.prechamber`. A fraction
    ``pc_energy_fraction`` of the fuel energy is released in the pre-chamber; by default
    the pre-chamber share of the charge at top dead centre for a uniform temperature,
    V_pc / (V_pc + V_c), with V_c the main-chamber clearance volume. The result is a
    manufactured solution for code verification, not a model of a real engine.
    """
    if pc_energy_fraction is None:
        v_c = geometry.displaced_volume_m3 / (geometry.compression_ratio - 1.0)
        pc_energy_fraction = pre_chamber.volume_m3 / (pre_chamber.volume_m3 + v_c)
    if not 0.0 < pc_energy_fraction < 1.0:
        raise PhysicalBoundsError("pc_energy_fraction must lie in (0, 1).")
    theta = np.arange(theta_range_deg[0], theta_range_deg[1] + 0.5 * resolution_deg, resolution_deg)
    kin = geometry.evaluate(theta)
    mix = make_reactant_mixture(fuel, lambda_air)
    r = mix.specific_gas_constant
    v_pc = pre_chamber.volume_m3
    m_total = p_ivc_pa * (float(kin.volume_m3[0]) + v_pc) / (r * t_ivc_k)
    m_pc0 = p_ivc_pa * v_pc / (r * t_ivc_k)
    fuel_mass = m_total / (1.0 + lambda_air * fuel.stoichiometric_afr)
    q_total = fuel_mass * fuel.lhv_j_per_kg
    q_pc, q_mc = pc_energy_fraction * q_total, (1.0 - pc_energy_fraction) * q_total
    a = WIEBE_EFFICIENCY_A
    norm = 1.0 - math.exp(-a)

    def wiebe_rate(th: float, start: float, dur: float, total: float) -> float:
        u = (th - start) / dur
        if u <= 0.0 or u >= 1.0:
            return 0.0
        return total * a * 3.0 / dur * u**2 * math.exp(-a * u**3) / norm

    def rhs(th: float, y: FloatArray) -> list[float]:
        p_pc, p_mc, m_pc = float(y[0]), float(y[1]), float(y[2])
        k = geometry.evaluate(np.array([th]))
        v_mc, dv = float(k.volume_m3[0]), float(k.dv_dtheta_m3_per_deg[0])
        m_mc = m_total - m_pc
        t_pc, t_mc = p_pc * v_pc / (m_pc * r), p_mc * v_mc / (m_mc * r)
        x = float(transfer_rate(p_pc, t_pc, p_mc, t_mc, mixture_pc=mix, mixture_mc=mix,
                                geometry=pre_chamber, rpm=rpm))
        h_flow = float(mix.enthalpy(t_pc if x >= 0.0 else t_mc))
        g_pc, g_mc = float(mix.gamma(t_pc)), float(mix.gamma(t_mc))
        b_pc = float(mix.internal_energy(t_pc) - mix.cv(t_pc) * t_pc) - h_flow
        b_mc = float(mix.internal_energy(t_mc) - mix.cv(t_mc) * t_mc) - h_flow
        dq_pc = wiebe_rate(th, pc_start_deg, pc_duration_deg, q_pc)
        dq_mc = wiebe_rate(th, mc_start_deg, mc_duration_deg, q_mc)
        dp_pc = (g_pc - 1.0) / v_pc * (dq_pc - b_pc * (-x))
        dp_mc = (g_mc - 1.0) / v_mc * (dq_mc - g_mc / (g_mc - 1.0) * p_mc * dv - b_mc * x)
        return [dp_pc, dp_mc, -x]

    sol = solve_ivp(rhs, (float(theta[0]), float(theta[-1])), [p_ivc_pa, p_ivc_pa, m_pc0],
                    t_eval=theta, method="LSODA", rtol=1e-10, atol=[1e-4, 1e-4, 1e-15], max_step=0.05)
    if not sol.success:  # pragma: no cover
        raise RuntimeError(f"Forward integration failed: {sol.message}")
    p_pc, p_mc, m_pc = (np.asarray(sol.y[i], dtype=np.float64) for i in range(3))
    m_mc = m_total - m_pc
    t_pc = p_pc * v_pc / (m_pc * r)
    t_mc = p_mc * kin.volume_m3 / (m_mc * r)
    x = transfer_rate(p_pc, t_pc, p_mc, t_mc, mixture_pc=mix, mixture_mc=mix, geometry=pre_chamber, rpm=rpm)
    dq_pc = np.array([wiebe_rate(float(t), pc_start_deg, pc_duration_deg, q_pc) for t in theta])
    dq_mc = np.array([wiebe_rate(float(t), mc_start_deg, mc_duration_deg, q_mc) for t in theta])
    return ManufacturedPreChamberCycle(
        crank_angle_deg=theta,
        pre_chamber_pressure_pa=p_pc,
        main_chamber_pressure_pa=p_mc,
        main_chamber_volume_m3=kin.volume_m3,
        main_chamber_dv_dtheta_m3_per_deg=kin.dv_dtheta_m3_per_deg,
        pre_chamber_mass_kg=m_pc,
        main_chamber_mass_kg=m_mc,
        pre_chamber_temperature_k=t_pc,
        main_chamber_temperature_k=t_mc,
        transfer_rate_kg_per_deg=x,
        dq_pc_j_per_deg=dq_pc,
        dq_mc_j_per_deg=dq_mc,
        q_pc_total_j=q_pc,
        q_mc_total_j=q_mc,
        trapped_mass_kg=m_total,
        rpm=rpm,
    )


def write_pre_chamber_case_directory(
    cycle: ManufacturedPreChamberCycle,
    geometry: SliderCrankGeometry,
    pre_chamber: PreChamberGeometry,
    directory: str | Path,
    *,
    name: str,
    fuel_key: str,
    lambda_air: float,
    include_temperatures: bool = True,
    configuration: str | None = None,
) -> Path:
    """Write a manufactured pre-chamber cycle as ``case.toml``, ``thermo.csv``, and ``prechamber.csv``.

    ``thermo.csv`` holds the main chamber (crank angle, pressure in MPa, volume in cm^3,
    and optionally temperature), ``prechamber.csv`` the pre-chamber (crank angle,
    pressure in bar, optionally temperature). Without temperatures the analysis uses
    the orifice model with the trapped mass written to the case file.
    """
    out = Path(directory)
    out.mkdir(parents=True, exist_ok=True)
    th = cycle.crank_angle_deg
    main_cols = ["Crank Angle [deg]", "Pressure [MPa]", "Volume [cm^3]"]
    main_data = [th, cycle.main_chamber_pressure_pa * 1e-6, cycle.main_chamber_volume_m3 * 1e6]
    pc_cols = ["Crank Angle [deg]", "PC Pressure [bar]"]
    pc_data = [th, cycle.pre_chamber_pressure_pa * 1e-5]
    if include_temperatures:
        main_cols.append("Temperature [K]")
        main_data.append(cycle.main_chamber_temperature_k)
        pc_cols.append("PC Temperature [K]")
        pc_data.append(cycle.pre_chamber_temperature_k)
    for fname, cols, data in (("thermo.csv", main_cols, main_data), ("prechamber.csv", pc_cols, pc_data)):
        rows = [",".join(cols)] + [",".join(f"{v:.12g}" for v in row) for row in np.column_stack(data)]
        (out / fname).write_text("\n".join(rows) + "\n", encoding="utf-8", newline="\n")

    toml = [
        f"name = {_toml_value(name)}",
        'source = "synthetic"',
        *([f"configuration = {_toml_value(configuration)}"] if configuration else []),
        "",
        "[geometry]",
        f"bore_m = {geometry.bore_m!r}",
        f"stroke_m = {geometry.stroke_m!r}",
        f"conrod_m = {geometry.conrod_m!r}",
        f"compression_ratio = {geometry.compression_ratio!r}",
        "",
        "[operating_point]",
        f"rpm = {float(cycle.rpm)!r}",
        f"fuel = {_toml_value(fuel_key)}",
        f"lambda_air = {float(lambda_air)!r}",
        f"trapped_mass_kg = {cycle.trapped_mass_kg!r}",
        "",
        "[trace]",
        'file = "thermo.csv"',
        'format = "forte_thermo"',
        "",
        "[pre_chamber]",
        f"volume_m3 = {pre_chamber.volume_m3!r}",
        f"orifice_count = {pre_chamber.orifice_count}",
        f"orifice_diameter_m = {pre_chamber.orifice_diameter_m!r}",
        f"discharge_coefficient = {pre_chamber.discharge_coefficient!r}",
        'file = "prechamber.csv"',
        "",
        "[pre_chamber.columns]",
        'pressure = "PC Pressure"',
        *(['temperature = "PC Temperature"'] if include_temperatures else []),
        "",
        "[metadata]",
        '"generator" = "combusta.synthetic.manufactured_pre_chamber_cycle"',
    ]
    path = out / "case.toml"
    path.write_text("\n".join(toml) + "\n", encoding="utf-8", newline="\n")
    return path
