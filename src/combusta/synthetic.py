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
    (out / trace_name).write_text("\n".join(lines) + "\n", encoding="utf-8")

    geo, op = case.geometry, case.operating_point
    toml = [
        f"name = {_toml_value(case.name)}",
        f"source = {_toml_value(case.source)}",
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
    path.write_text("\n".join(toml) + "\n", encoding="utf-8")
    return path
