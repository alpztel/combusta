"""Source-independent engine case model.

An :class:`EngineCase` bundles everything an analysis needs about one operating
point of one engine configuration: geometry, operating conditions, fuel, the
crank-angle-resolved trace, and scalar results such as end-of-cycle emissions.
Readers (ANSYS Forte, test-bench files, 1D models) only have to produce an
``EngineCase``; every analysis consumes it, so no analysis depends on a file format.

The ``source`` field records provenance. Cases marked ``"synthetic"`` come from
the analytical generators used for code verification and must never be reported
as engineering results.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Literal, TypeAlias

import numpy as np
from numpy.typing import ArrayLike, NDArray

from combusta.exceptions import PhysicalBoundsError
from combusta.io.forte_parser import ForteThermoData
from combusta.thermo.fuels import Fuel, make_reactant_mixture
from combusta.thermo.kinematics import KinematicsState, SliderCrankGeometry

FloatArray: TypeAlias = NDArray[np.float64]
CaseSource: TypeAlias = Literal["cfd", "test_bench", "model_1d", "synthetic"]


def _optional_positive(value: float | None, name: str) -> None:
    if value is not None and (not np.isfinite(value) or value <= 0.0):
        raise PhysicalBoundsError(f"{name} must be finite and > 0 when provided, got {value}.")


@dataclass(frozen=True, slots=True)
class OperatingPoint:
    """Operating conditions of one engine case (strict SI units).

    Parameters
    ----------
    rpm : float
        Engine speed $N$ [rev/min].
    fuel : Fuel
        Fuel definition.
    lambda_air : float | None
        Global air-excess ratio $\\lambda$ [-].
    fuel_mass_per_cycle_kg : float | None
        Fuel mass delivered per cycle and cylinder $m_f$ [kg].
    trapped_mass_kg : float | None
        In-cylinder trapped mass at intake valve closing [kg].
    ignition_timing_deg : float | None
        Spark timing or start of injection [deg ATDC].
    intake_pressure_pa, intake_temperature_k : float | None
        Intake manifold state [Pa], [K].
    wall_temperature_k : float | None
        Representative combustion-chamber wall temperature [K].
    """

    rpm: float
    fuel: Fuel
    lambda_air: float | None = None
    fuel_mass_per_cycle_kg: float | None = None
    trapped_mass_kg: float | None = None
    ignition_timing_deg: float | None = None
    intake_pressure_pa: float | None = None
    intake_temperature_k: float | None = None
    wall_temperature_k: float | None = None

    def __post_init__(self) -> None:
        if not np.isfinite(self.rpm) or self.rpm <= 0.0:
            raise PhysicalBoundsError(f"Engine speed rpm must be finite and > 0, got {self.rpm}.")
        _optional_positive(self.lambda_air, "lambda_air")
        _optional_positive(self.fuel_mass_per_cycle_kg, "fuel_mass_per_cycle_kg")
        _optional_positive(self.trapped_mass_kg, "trapped_mass_kg")
        _optional_positive(self.intake_pressure_pa, "intake_pressure_pa")
        _optional_positive(self.intake_temperature_k, "intake_temperature_k")
        _optional_positive(self.wall_temperature_k, "wall_temperature_k")
        if self.ignition_timing_deg is not None and not np.isfinite(self.ignition_timing_deg):
            raise PhysicalBoundsError("ignition_timing_deg must be finite when provided.")

    @property
    def fuel_energy_per_cycle_j(self) -> float:
        """Fuel chemical energy per cycle $m_f \\cdot \\mathrm{LHV}$ [J]."""
        if self.fuel_mass_per_cycle_kg is None:
            raise ValueError("fuel_mass_per_cycle_kg is required for the fuel energy.")
        return self.fuel_mass_per_cycle_kg * self.fuel.lhv_j_per_kg

    @property
    def equivalence_ratio(self) -> float:
        """Fuel-air equivalence ratio $\\phi = 1 / \\lambda$ [-]."""
        if self.lambda_air is None:
            raise ValueError("lambda_air is required for the equivalence ratio.")
        return 1.0 / self.lambda_air


@dataclass(frozen=True, slots=True)
class CycleTrace:
    """Crank-angle-resolved in-cylinder signals in strict SI units.

    Parameters
    ----------
    crank_angle_deg : ArrayLike
        Strictly increasing crank angle $\\theta$ [deg ATDC].
    pressure_pa : ArrayLike
        Cylinder pressure $P(\\theta)$ [Pa].
    volume_m3, temperature_k, wall_heat_loss_j_per_deg : ArrayLike | None
        Optional volume [m^3], mass-averaged temperature [K], and wall heat-loss rate
        [J/deg, positive = heat leaving the gas] on the same crank-angle grid.
    channels : Mapping[str, ArrayLike]
        Additional named signals on the same grid (e.g. ``"NO_ppm"``).
    """

    crank_angle_deg: FloatArray
    pressure_pa: FloatArray
    volume_m3: FloatArray | None = None
    temperature_k: FloatArray | None = None
    wall_heat_loss_j_per_deg: FloatArray | None = None
    channels: Mapping[str, FloatArray] = field(default_factory=dict)

    def __post_init__(self) -> None:
        theta = np.asarray(self.crank_angle_deg, dtype=np.float64)
        if theta.ndim != 1 or theta.size < 5:
            raise ValueError("crank_angle_deg must be a 1D array with at least 5 samples.")
        if not np.all(np.isfinite(theta)) or np.any(np.diff(theta) <= 0.0):
            raise ValueError("crank_angle_deg must be finite and strictly increasing.")
        object.__setattr__(self, "crank_angle_deg", theta)

        def checked(values: ArrayLike, name: str, *, positive: bool) -> FloatArray:
            arr = np.asarray(values, dtype=np.float64)
            if arr.shape != theta.shape:
                raise ValueError(f"{name} must match crank_angle_deg shape {theta.shape}.")
            if not np.all(np.isfinite(arr)):
                raise PhysicalBoundsError(f"{name} contains non-finite values.")
            if positive and np.any(arr <= 0.0):
                raise PhysicalBoundsError(f"{name} must be strictly positive.")
            return arr

        object.__setattr__(self, "pressure_pa", checked(self.pressure_pa, "pressure_pa", positive=True))
        if self.volume_m3 is not None:
            object.__setattr__(self, "volume_m3", checked(self.volume_m3, "volume_m3", positive=True))
        if self.temperature_k is not None:
            object.__setattr__(
                self, "temperature_k", checked(self.temperature_k, "temperature_k", positive=True)
            )
        if self.wall_heat_loss_j_per_deg is not None:
            object.__setattr__(
                self,
                "wall_heat_loss_j_per_deg",
                checked(self.wall_heat_loss_j_per_deg, "wall_heat_loss_j_per_deg", positive=False),
            )
        object.__setattr__(
            self,
            "channels",
            {key: checked(val, f"channel '{key}'", positive=False) for key, val in self.channels.items()},
        )


@dataclass(frozen=True, slots=True)
class EngineCase:
    """One engine configuration at one operating point, independent of data source.

    Parameters
    ----------
    name : str
        Case label used in comparisons and reports.
    geometry : SliderCrankGeometry
        Cylinder geometry.
    operating_point : OperatingPoint
        Speed, fuel, mixture, and boundary conditions.
    trace : CycleTrace
        Crank-angle-resolved signals.
    source : {"cfd", "test_bench", "model_1d", "synthetic"}
        Data provenance. ``"synthetic"`` marks verification fixtures only.
    scalars : Mapping[str, float]
        Scalar results not resolved in crank angle (e.g. ``"NOx_ppm_evo"``).
    metadata : Mapping[str, str]
        Free-form provenance (solver version, mesh, file path, ...).
    """

    name: str
    geometry: SliderCrankGeometry
    operating_point: OperatingPoint
    trace: CycleTrace
    source: CaseSource
    scalars: Mapping[str, float] = field(default_factory=dict)
    metadata: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.source not in ("cfd", "test_bench", "model_1d", "synthetic"):
            raise ValueError(f"Unknown case source '{self.source}'.")
        object.__setattr__(self, "scalars", {k: float(v) for k, v in self.scalars.items()})
        object.__setattr__(self, "metadata", dict(self.metadata))

    @property
    def is_synthetic(self) -> bool:
        """``True`` for analytical verification fixtures (never engineering results)."""
        return self.source == "synthetic"

    def kinematics(self) -> KinematicsState:
        """Slider-crank kinematics evaluated on the trace crank-angle grid."""
        return self.geometry.evaluate(self.trace.crank_angle_deg)

    def volume_m3(self) -> FloatArray:
        """Cylinder volume [m^3]: the traced volume if present, else from geometry."""
        if self.trace.volume_m3 is not None:
            return self.trace.volume_m3
        return self.kinematics().volume_m3

    def volume_consistency_error(self) -> float:
        """Maximum relative deviation between traced and geometric volume [-].

        A large value usually means a wrong compression ratio, connecting-rod length,
        or crank-angle offset between the data source and the declared geometry.
        """
        if self.trace.volume_m3 is None:
            raise ValueError("Trace has no volume; nothing to compare with the geometry.")
        v_geom = self.kinematics().volume_m3
        return float(np.max(np.abs(self.trace.volume_m3 - v_geom) / v_geom))

    def temperature_k(self) -> FloatArray:
        """Mass-averaged gas temperature [K].

        Uses the traced temperature when available; otherwise applies the ideal-gas law
        $T = P V / (m R_u)$ with the trapped mass and the unburned-charge gas constant.
        """
        if self.trace.temperature_k is not None:
            return self.trace.temperature_k
        op = self.operating_point
        if op.trapped_mass_kg is None or op.lambda_air is None:
            raise ValueError(
                "Temperature needs a traced temperature or both trapped_mass_kg and lambda_air."
            )
        r_spec = make_reactant_mixture(op.fuel, op.lambda_air).specific_gas_constant
        return np.asarray(
            self.trace.pressure_pa * self.volume_m3() / (op.trapped_mass_kg * r_spec),
            dtype=np.float64,
        )

    @classmethod
    def from_forte_thermo(
        cls,
        name: str,
        thermo: ForteThermoData,
        *,
        geometry: SliderCrankGeometry,
        operating_point: OperatingPoint,
        scalars: Mapping[str, float] | None = None,
        metadata: Mapping[str, str] | None = None,
        wall_heat_loss_sign: float = 1.0,
    ) -> EngineCase:
        """Build a CFD case from a parsed ANSYS Forte ``thermo.csv``.

        ``wall_heat_loss_sign`` converts the exported wall heat-transfer column to the
        Combusta convention (positive = heat leaving the gas). Pass ``-1.0`` if the
        export reports heat loss as negative; check it on the first real case.
        """
        if wall_heat_loss_sign not in (1.0, -1.0):
            raise ValueError("wall_heat_loss_sign must be +1.0 or -1.0.")
        meta = dict(metadata or {})
        if thermo.source_units:
            meta.update({f"unit:{k}": v for k, v in thermo.source_units.items()})
        trace = CycleTrace(
            crank_angle_deg=thermo.crank_angle_deg,
            pressure_pa=thermo.pressure_pa,
            volume_m3=thermo.volume_m3,
            temperature_k=thermo.temperature_k,
            wall_heat_loss_j_per_deg=(
                None
                if thermo.wall_heat_loss_j_per_deg is None
                else wall_heat_loss_sign * thermo.wall_heat_loss_j_per_deg
            ),
        )
        meta["wall_heat_loss_sign"] = f"{wall_heat_loss_sign:+.0f}"
        return cls(
            name=name,
            geometry=geometry,
            operating_point=operating_point,
            trace=trace,
            source="cfd",
            scalars=dict(scalars or {}),
            metadata=meta,
        )
