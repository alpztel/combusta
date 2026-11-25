"""Aero-diesel altitude sweep: an end-to-end, fully synthetic demonstration.

Every number in this script comes from textbook models; no engine, CFD, or test-bench
data are used. The script shows how the Combusta building blocks chain together and
produces the figures in ``assets/``. Its outputs are illustrations, not engineering
results.

Model chain, evaluated at 0 ft, 15,000 ft and 30,000 ft (ISO 2533 atmosphere):

1. Intake: a two-stage turbocharger holds a target manifold pressure until its maximum
   overall pressure ratio is reached (critical altitude of roughly 21,000 ft), followed
   by charge-air cooling toward ambient temperature.
2. Motored compression: closed-cycle single-zone First-Law ODE with NASA dry-air
   properties and Woschni (1967) wall heat transfer.
3. Ignition: Livengood-Wu (1955) integral of the Assanis et al. (2003) ignition-delay
   correlation, evaluated on the motored state from start of injection; its root gives
   the start of combustion (SOC).
4. Fired cycle: Wiebe heat release starting at the predicted SOC, integrated with the
   same ODE; synthetic 3.5 kPa transducer noise is added and removed with a zero-phase
   Savitzky-Golay filter.
5. Inverse analysis: gross heat release, CA10/50/90 and IMEP; spray tip penetration at
   the in-cylinder density at injection; Second-Law exergy budget with an altitude
   dead state, including its independent closure check.
6. GCI: the three-grid procedure applied to a *manufactured* power-law triplet, to
   illustrate the method and the figure. It is not a mesh study.
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Final, TypeAlias

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

# Allow direct execution via ``python examples/altitude_sweep_demo.py``
_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT / "src"))

from combusta import (  # noqa: E402
    DRY_AIR,
    FEET_TO_M,
    JET_A1_EXERGY_FACTOR,
    JET_A1_SURROGATE,
    PREDEFINED_FUELS,
    AltitudeSprayPlotCase,
    BurnFractionMetrics,
    ExergyCycleBudget,
    DeadStateConditions,
    HeatReleaseResult,
    IMEPResult,
    LivengoodWuResult,
    SliderCrankGeometry,
    ThermoOverlayTrace,
    WoschniCorrelation,
    compute_exergy_cycle_budget,
    compute_exergy_rates,
    compute_gci,
    compute_grid_size,
    compute_heat_release,
    compute_imep,
    compute_isa_conditions,
    evaluate_altitude_spray_risk,
    extract_burn_fractions,
    hiroyasu_arai_penetration,
    livengood_wu_integral,
    naber_siebers_penetration,
    plot_altitude_spray_comparison,
    plot_gci_convergence,
    plot_pv_log_diagram,
    plot_thermo_diagnostics,
    savgol_zero_phase,
)

FloatArray: TypeAlias = NDArray[np.float64]

# --- Engine and operating assumptions (illustrative, not data of any specific engine) ---
GEOMETRY: Final = SliderCrankGeometry(bore_m=0.085, stroke_m=0.090, conrod_m=0.155, compression_ratio=16.5)
ENGINE_SPEED_RPM: Final[float] = 2400.0
WALL_TEMPERATURE_K: Final[float] = 460.0
THETA_SOI_DEG: Final[float] = -8.0
BURN_DURATION_DEG: Final[float] = 40.0
WIEBE_A: Final[float] = 6.908  # 99.9 % burned at the end of the duration
WIEBE_M: Final[float] = 2.0
NOMINAL_FUEL_KG: Final[float] = 2.45e-5
LAMBDA_MIN: Final[float] = 1.45  # smoke-limited minimum air-excess ratio

TARGET_MANIFOLD_PRESSURE_PA: Final[float] = 2.4 * 101_325.0
MAX_PRESSURE_RATIO: Final[float] = 5.5  # two-stage turbocharger
COMPRESSOR_EFFICIENCY: Final[float] = 0.72
CHARGE_AIR_COOLER_EFFECTIVENESS: Final[float] = 0.80
INTAKE_HEATING_K: Final[float] = 15.0  # wall and residual-gas heating up to IVC

INJECTION_DELTA_P_PA: Final[float] = 100.0e6
FUEL_DENSITY_KG_M3: Final[float] = 810.0
NOZZLE_DIAMETER_M: Final[float] = 90.0e-6
SPRAY_CONE_DEG: Final[float] = 19.0
BOWL_WALL_DISTANCE_M: Final[float] = 0.038

THETA: Final[FloatArray] = np.linspace(-180.0, 180.0, 3601, dtype=np.float64)


@dataclass(frozen=True, slots=True)
class AltitudeCase:
    """Results of one synthetic altitude case."""

    label: str
    altitude_ft: float
    p_ambient_pa: float
    p_ivc_pa: float
    t_ivc_k: float
    rho_soi_kg_m3: float
    fuel_mass_kg: float
    ignition: LivengoodWuResult
    pressure_pa: FloatArray
    volume_m3: FloatArray
    hrr: HeatReleaseResult
    burn: BurnFractionMetrics
    imep: IMEPResult
    exergy: ExergyCycleBudget


def intake_state(altitude_m: float) -> tuple[float, float, float]:
    """Return ambient pressure, IVC pressure and IVC temperature for a boost-controlled engine."""
    isa = compute_isa_conditions(altitude_m)
    p_amb = float(isa.pressure_pa[0])
    t_amb = float(isa.temperature_k[0])
    pressure_ratio = min(TARGET_MANIFOLD_PRESSURE_PA / p_amb, MAX_PRESSURE_RATIO)
    t_compressor = t_amb * (1.0 + (pressure_ratio ** (0.4 / 1.4) - 1.0) / COMPRESSOR_EFFICIENCY)
    t_manifold = t_compressor - CHARGE_AIR_COOLER_EFFECTIVENESS * (t_compressor - t_amb)
    return p_amb, pressure_ratio * p_amb, t_manifold + INTAKE_HEATING_K


def wiebe_rate(theta_deg: float, q_total_j: float, start_deg: float) -> float:
    """Wiebe heat-release rate dQ/dθ [J/deg]."""
    u = (theta_deg - start_deg) / BURN_DURATION_DEG
    if u <= 0.0 or u >= 1.0:
        return 0.0
    return (
        q_total_j
        * WIEBE_A
        * (WIEBE_M + 1.0)
        / BURN_DURATION_DEG
        * u**WIEBE_M
        * math.exp(-WIEBE_A * u ** (WIEBE_M + 1.0))
    )


def single_zone_cycle(p_ivc_pa: float, t_ivc_k: float, q_total_j: float, soc_deg: float) -> FloatArray:
    """Integrate the closed-cycle single-zone pressure equation with Woschni wall heat loss."""
    kin_ivc = GEOMETRY.evaluate(THETA[:1])
    v_ivc = float(kin_ivc.volume_m3[0])
    r_air = DRY_AIR.specific_gas_constant
    mass = p_ivc_pa * v_ivc / (r_air * t_ivc_k)
    s_p = 2.0 * GEOMETRY.stroke_m * ENGINE_SPEED_RPM / 60.0
    omega = 6.0 * ENGINE_SPEED_RPM

    def rhs(th: float, state: FloatArray) -> list[float]:
        p = max(float(state[0]), 1.0e3)
        k = GEOMETRY.evaluate(np.array([th]))
        v, dv, area = float(k.volume_m3[0]), float(k.dv_dtheta_m3_per_deg[0]), float(k.wall_area_m2[0])
        t = min(max(p * v / (mass * r_air), 210.0), 5800.0)
        gamma = float(DRY_AIR.gamma(np.array([t]))[0])
        p_motored = p_ivc_pa * (v_ivc / v) ** 1.35
        w = 2.28 * s_p + 3.24e-3 * (GEOMETRY.displaced_volume_m3 * t_ivc_k / (p_ivc_pa * v_ivc)) * max(p - p_motored, 0.0)
        h = 130.0 * GEOMETRY.bore_m**-0.2 * (p * 1.0e-5) ** 0.8 * t**-0.53 * w**0.8
        dq_wall = h * area * (t - WALL_TEMPERATURE_K) / omega
        dq_fuel = wiebe_rate(th, q_total_j, soc_deg)
        return [((gamma - 1.0) * (dq_fuel - dq_wall) - gamma * p * dv) / v]

    sol = solve_ivp(rhs, (THETA[0], THETA[-1]), [p_ivc_pa], t_eval=THETA, rtol=1.0e-8, atol=1.0e-4)
    if not sol.success:
        raise RuntimeError(sol.message)
    return np.asarray(sol.y[0], dtype=np.float64)


def run_altitude_case(altitude_ft: float, rng: np.random.Generator) -> AltitudeCase:
    """Synthesize and analyze one altitude case."""
    altitude_m = altitude_ft * FEET_TO_M
    p_amb, p_ivc, t_ivc = intake_state(altitude_m)
    kin = GEOMETRY.evaluate(THETA)
    r_air = DRY_AIR.specific_gas_constant
    mass = p_ivc * float(kin.volume_m3[0]) / (r_air * t_ivc)

    fuel_mass = min(NOMINAL_FUEL_KG, mass / (LAMBDA_MIN * JET_A1_SURROGATE.stoichiometric_afr))
    q_fuel = fuel_mass * JET_A1_SURROGATE.lhv_j_per_kg

    # Ignition predicted on the motored state (before combustion the charge is motored)
    p_motored = single_zone_cycle(p_ivc, t_ivc, 0.0, soc_deg=1.0e3)
    t_motored = p_motored * kin.volume_m3 / (mass * r_air)
    phi_global = fuel_mass * JET_A1_SURROGATE.stoichiometric_afr / mass
    ignition = livengood_wu_integral(
        THETA,
        p_motored,
        t_motored,
        ENGINE_SPEED_RPM,
        THETA_SOI_DEG,
        equivalence_ratio=phi_global,
        cetane_number=PREDEFINED_FUELS["JET_A1_JP8"].cetane_number,
    )
    if ignition.misfire_detected:
        raise RuntimeError(f"No auto-ignition predicted at {altitude_ft:.0f} ft.")

    # Fired cycle starting at the predicted SOC, then synthetic measurement noise + filtering
    p_fired = single_zone_cycle(p_ivc, t_ivc, q_fuel, soc_deg=ignition.theta_soc_deg)
    p_measured = p_fired + rng.normal(0.0, 3_500.0, size=p_fired.shape)
    p_filtered = savgol_zero_phase(p_measured, window_length=31, polyorder=3)
    t_filtered = p_filtered * kin.volume_m3 / (mass * r_air)

    v_ivc = float(kin.volume_m3[0])
    heat_transfer = WoschniCorrelation(
        bore_m=GEOMETRY.bore_m,
        stroke_m=GEOMETRY.stroke_m,
        rpm=ENGINE_SPEED_RPM,
        wall_temp_k=WALL_TEMPERATURE_K,
        displaced_volume_m3=GEOMETRY.displaced_volume_m3,
    ).evaluate(
        p_filtered,
        t_filtered,
        kin.wall_area_m2,
        motored_pressure_pa=p_ivc * (v_ivc / kin.volume_m3) ** 1.35,
        p_ref_pa=p_ivc,
        v_ref_m3=v_ivc,
        t_ref_k=t_ivc,
    )
    hrr = compute_heat_release(
        THETA,
        p_filtered,
        kin.volume_m3,
        temperature_k=t_filtered,
        dv_dtheta_m3_per_deg=kin.dv_dtheta_m3_per_deg,
        dq_ht_j_per_deg=heat_transfer.dq_ht_j_per_deg,
        mixture=DRY_AIR,
    )
    burn = extract_burn_fractions(THETA, hrr.dq_gross_j_per_deg, start_deg=-15.0, end_deg=70.0)
    imep = compute_imep(
        THETA,
        p_filtered,
        kin.volume_m3,
        displaced_volume_m3=GEOMETRY.displaced_volume_m3,
        dv_dtheta_m3_per_deg=kin.dv_dtheta_m3_per_deg,
    )

    rates = compute_exergy_rates(
        THETA,
        p_filtered,
        t_filtered,
        kin.volume_m3,
        kin.dv_dtheta_m3_per_deg,
        hrr.dq_gross_j_per_deg,
        heat_transfer.dq_ht_j_per_deg,
        DeadStateConditions.from_altitude(altitude_m),
        mass,
        wall_temperature_k=WALL_TEMPERATURE_K,
        phi_factor=JET_A1_EXERGY_FACTOR,
        gas_mixture=DRY_AIR,
    )
    exergy = compute_exergy_cycle_budget(
        rates,
        fuel_mass_kg=fuel_mass,
        lhv_j_per_kg=JET_A1_SURROGATE.lhv_j_per_kg,
        phi_factor=JET_A1_EXERGY_FACTOR,
    )

    soi_index = int(np.argmin(np.abs(THETA - THETA_SOI_DEG)))
    label = "Sea level (0 ft)" if altitude_ft == 0.0 else f"{altitude_ft:,.0f} ft"
    return AltitudeCase(
        label=label,
        altitude_ft=altitude_ft,
        p_ambient_pa=p_amb,
        p_ivc_pa=p_ivc,
        t_ivc_k=t_ivc,
        rho_soi_kg_m3=float(mass / kin.volume_m3[soi_index]),
        fuel_mass_kg=fuel_mass,
        ignition=ignition,
        pressure_pa=p_filtered,
        volume_m3=kin.volume_m3,
        hrr=hrr,
        burn=burn,
        imep=imep,
        exergy=exergy,
    )


def spray_case(case: AltitudeCase, rho_sea_level: float) -> tuple[AltitudeSprayPlotCase, float | None]:
    """Spray tip penetration curves and the time the tip reaches the bowl wall."""
    t_spray = np.linspace(1.0e-6, 0.85e-3, 300, dtype=np.float64)
    common = dict(
        delta_p_pa=INJECTION_DELTA_P_PA,
        rho_liquid_kg_m3=FUEL_DENSITY_KG_M3,
        nozzle_diameter_m=NOZZLE_DIAMETER_M,
    )
    ha = hiroyasu_arai_penetration(t_spray, rho_gas_kg_m3=case.rho_soi_kg_m3, **common)
    ns = naber_siebers_penetration(
        t_spray, rho_gas_kg_m3=case.rho_soi_kg_m3, spray_angle_deg=SPRAY_CONE_DEG, **common
    )
    risk = evaluate_altitude_spray_risk(
        t_spray,
        rho_gas_sea_level_kg_m3=rho_sea_level,
        rho_gas_altitude_kg_m3=case.rho_soi_kg_m3,
        wall_distance_m=BOWL_WALL_DISTANCE_M,
        model="hiroyasu_arai",
        **common,
    )
    plot_case = AltitudeSprayPlotCase(
        label=case.label,
        rho_gas_kg_m3=case.rho_soi_kg_m3,
        time_s=t_spray,
        hiroyasu_penetration_m=ha.penetration_m,
        siebers_penetration_m=ns.penetration_m,
        impingement_time_s=risk.impingement_time_s,
    )
    return plot_case, risk.impingement_time_s


def manufactured_gci_example(p_max_pa: float):  # noqa: ANN201 - returns GCIResult
    """Three-grid GCI on a manufactured second-order triplet (illustration only)."""
    grid = compute_grid_size(cell_counts=(1_920_000.0, 680_000.0, 235_000.0))
    order = 1.92
    coefficient = -0.0085 * p_max_pa / grid.h1**order
    phi = [p_max_pa + coefficient * h**order for h in (grid.h1, grid.h2, grid.h3)]
    return compute_gci(phi[0], phi[1], phi[2], grid=grid, formal_order=2.0)


def print_summary(cases: list[AltitudeCase], tip_times: list[float | None]) -> None:
    """Print the synthetic results table (Markdown)."""
    print("\nSYNTHETIC DEMONSTRATION - textbook models only, not engineering results\n")
    print(
        "| Altitude | p_amb [kPa] | p_IVC [kPa] | rho_SOI [kg/m3] | Ignition delay [deg / ms] "
        "| SOC [deg ATDC] | p_max [MPa] | IMEP_net [bar] | CA50 [deg ATDC] | Tip at 38 mm [ms] "
        "| eta_II [%] | Exergy closure error [%] |"
    )
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for case, tip in zip(cases, tip_times):
        ign = case.ignition
        tip_str = f"{tip * 1e3:.2f}" if tip is not None else "not reached"
        closure = 100.0 * case.exergy.closure_error_j / case.exergy.e_released_j
        print(
            f"| {case.label} | {case.p_ambient_pa * 1e-3:.1f} | {case.p_ivc_pa * 1e-3:.1f} "
            f"| {case.rho_soi_kg_m3:.1f} | {ign.ignition_delay_deg:.1f} / {ign.ignition_delay_s * 1e3:.2f} "
            f"| {ign.theta_soc_deg:+.2f} | {float(np.max(case.pressure_pa)) * 1e-6:.2f} "
            f"| {case.imep.imep_net_pa * 1e-5:.2f} | {case.burn.ca50_deg:.2f} | {tip_str} "
            f"| {case.exergy.eta_second_law * 100:.1f} | {closure:.2f} |"
        )
    print()


def main() -> None:
    """Run the three altitude cases and export the figures to ``assets/``."""
    rng = np.random.default_rng(42)
    cases = [run_altitude_case(alt, rng) for alt in (0.0, 15_000.0, 30_000.0)]
    sprays = [spray_case(c, rho_sea_level=cases[0].rho_soi_kg_m3) for c in cases]
    print_summary(cases, [tip for _, tip in sprays])

    gci = manufactured_gci_example(float(np.max(cases[0].pressure_pa)))
    print(
        f"Manufactured GCI example: p = {gci.apparent_order_p:.3f} (order used {gci.order_used_p:.2f}), "
        f"GCI_fine = {gci.gci_fine_21 * 100:.3f} %\n"
    )

    assets = _REPO_ROOT / "assets"
    assets.mkdir(parents=True, exist_ok=True)
    sl, mid, high = cases
    tag = " (synthetic)"
    plot_thermo_diagnostics(
        THETA,
        sl.pressure_pa,
        sl.hrr,
        sl.burn,
        primary_label=sl.label + tag,
        overlay_traces=[
            ThermoOverlayTrace(c.label + tag, THETA, c.pressure_pa, c.hrr, c.burn) for c in (mid, high)
        ],
        save_path=assets / "thermo_diagnostics.png",
    )
    plot_pv_log_diagram(
        sl.volume_m3,
        sl.pressure_pa,
        THETA,
        primary_label=sl.label + tag,
        overlay_loops=[(c.label + tag, c.volume_m3, c.pressure_pa) for c in (mid, high)],
        save_path=assets / "pv_log_diagram.png",
    )
    plot_altitude_spray_comparison(
        [s for s, _ in sprays],
        wall_distance_m=BOWL_WALL_DISTANCE_M,
        save_path=assets / "altitude_spray_comparison.png",
    )
    plot_gci_convergence(
        gci,
        metric_name=r"Manufactured example: $P_{\max}$",
        metric_unit="MPa",
        scale_factor=1.0e-6,
        save_path=assets / "gci_convergence.png",
    )
    for name in ("thermo_diagnostics", "pv_log_diagram", "altitude_spray_comparison", "gci_convergence"):
        print(f"  saved assets/{name}.png")


if __name__ == "__main__":
    main()
