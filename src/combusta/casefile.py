"""Engine case definitions in TOML files (``case.toml``).

A case file describes one :class:`~combusta.case.EngineCase` without code: geometry,
operating point, the trace file and how to read it, scalar results, emissions, and the
numerical uncertainties to use in comparisons. Keys are the field names of the
corresponding Combusta dataclasses and carry SI units in their names (``bore_m``,
``intake_pressure_pa``); unknown keys are rejected so that a typing error cannot pass
silently. ``docs/case-file.md`` documents the format.

Example::

    name = "open_chamber"
    source = "cfd"

    [geometry]
    bore_m = 0.0920
    stroke_m = 0.0850
    conrod_m = 0.1450
    compression_ratio = 11.0

    [operating_point]
    rpm = 2000.0
    fuel = "hydrogen"
    lambda_air = 2.0
    fuel_mass_per_cycle_kg = 1.0e-5

    [trace]
    file = "thermo.csv"
    format = "forte_thermo"
    wall_heat_loss_sign = 1.0

    [trace.columns]
    pressure = { column = "Pressure", unit = "MPa" }

    [uncertainty]
    p_max = 5.0e4
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, fields
import hashlib
import math
from pathlib import Path
import tomllib
from typing import Any, Final, TypeAlias

import numpy as np
from numpy.typing import NDArray

from combusta.case import CycleTrace, EngineCase, OperatingPoint
from combusta.exceptions import CaseFileError, CombustaError
from combusta.io.forte_parser import THERMO_FIELDS, ColumnSpec, parse_thermo_csv
from combusta.io.trace_csv import TRACE_FIELDS, parse_trace_csv
from combusta.prechamber import PreChamberGeometry
from combusta.thermo.fuels import (
    DIESEL_SURROGATE,
    FUEL_LIBRARY,
    HYDROGEN,
    ISO_OCTANE,
    JET_A1_SURROGATE,
    METHANE,
    N_DODECANE,
    Fuel,
)
from combusta.thermo.kinematics import SliderCrankGeometry

FloatArray: TypeAlias = NDArray[np.float64]

CASE_FILE_NAME: Final[str] = "case.toml"
"""File name looked up when a directory is given instead of a case file."""

TRACE_FORMATS: Final[tuple[str, ...]] = ("forte_thermo", "csv")
"""Supported trace formats: ANSYS Forte ``thermo.csv`` and a generic crank-angle table."""

FUEL_ALIASES: Final[Mapping[str, Fuel]] = {
    **{name.lower(): fuel for name, fuel in FUEL_LIBRARY.items()},
    "h2": HYDROGEN,
    "ch4": METHANE,
    "isooctane": ISO_OCTANE,
    "c8h18": ISO_OCTANE,
    "c12h26": N_DODECANE,
    "diesel": DIESEL_SURROGATE,
    "jet-a1": JET_A1_SURROGATE,
    "jp-8": JET_A1_SURROGATE,
}
"""Fuel names accepted by ``operating_point.fuel`` (case-insensitive)."""

_TOP_LEVEL_KEYS: Final[frozenset[str]] = frozenset({
    "name", "source", "description", "configuration", "geometry", "operating_point", "trace",
    "scalars", "emissions_kg_per_cycle", "uncertainty", "metadata", "pre_chamber",
})
_TRACE_KEYS: Final[frozenset[str]] = frozenset({"file", "format", "wall_heat_loss_sign", "columns"})
_PRE_CHAMBER_KEYS: Final[frozenset[str]] = frozenset({
    "volume_m3", "orifice_count", "orifice_diameter_m", "discharge_coefficient", "file", "columns",
})


@dataclass(frozen=True, slots=True)
class PreChamberInput:
    """Pre-chamber geometry and signals from the ``[pre_chamber]`` table of a case file.

    Attributes
    ----------
    geometry : PreChamberGeometry
        Volume and orifices.
    pressure_pa : NDArray[np.float64]
        Pre-chamber pressure on the case's crank-angle grid [Pa].
    temperature_k : NDArray[np.float64] | None
        Pre-chamber temperature [K], if exported.
    wall_heat_loss_j_per_deg : NDArray[np.float64] | None
        Pre-chamber wall heat loss [J/deg] in the Combusta sign convention, if exported.
    file : Path
        File the signals were read from.
    sha256 : str
        SHA-256 digest of that file.
    """

    geometry: PreChamberGeometry
    pressure_pa: FloatArray
    temperature_k: FloatArray | None
    wall_heat_loss_j_per_deg: FloatArray | None
    file: Path
    sha256: str


@dataclass(frozen=True, slots=True)
class LoadedCase:
    """An :class:`EngineCase` together with the provenance of its case file.

    Attributes
    ----------
    case : EngineCase
        The case. Its metadata records the case file, trace file, and trace SHA-256.
    uncertainties : Mapping[str, float]
        Absolute numerical uncertainty $U$ [SI] per metric key, from ``[uncertainty]``.
    case_file : Path
        Resolved path of the case file.
    trace_file : Path
        Resolved path of the trace file.
    trace_sha256 : str
        SHA-256 digest of the trace file, so a report identifies its input exactly.
    configuration : str | None
        Engine configuration the case belongs to (``configuration`` key), used to group
        the cases of a spark-timing sweep.
    pre_chamber : PreChamberInput | None
        Pre-chamber data from ``[pre_chamber]``; the case trace then describes the main
        chamber only.
    """

    case: EngineCase
    uncertainties: Mapping[str, float]
    case_file: Path
    trace_file: Path
    trace_sha256: str
    configuration: str | None = None
    pre_chamber: PreChamberInput | None = None


def _table(data: Mapping[str, Any], key: str, where: str, *, required: bool) -> dict[str, Any]:
    value = data.get(key)
    if value is None:
        if required:
            raise CaseFileError(f"{where}: missing [{key}] table.")
        return {}
    if not isinstance(value, dict):
        raise CaseFileError(f"{where}: '{key}' must be a table.")
    return value


def _check_keys(table: Mapping[str, Any], allowed: frozenset[str] | set[str], where: str) -> None:
    unknown = sorted(set(table) - set(allowed))
    if unknown:
        raise CaseFileError(f"{where}: unknown key(s) {unknown}; allowed: {sorted(allowed)}.")


def _numbers(table: Mapping[str, Any], where: str) -> dict[str, float]:
    out: dict[str, float] = {}
    for key, value in table.items():
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise CaseFileError(f"{where}: '{key}' must be a number, got {value!r}.")
        if not math.isfinite(float(value)):
            raise CaseFileError(f"{where}: '{key}' must be finite.")
        out[key] = float(value)
    return out


def _build(cls: type, table: Mapping[str, Any], where: str, **extra: Any) -> Any:
    allowed = {f.name for f in fields(cls)} - set(extra)
    _check_keys(table, allowed, where)
    try:
        return cls(**table, **extra)
    except (TypeError, ValueError, CombustaError) as exc:
        raise CaseFileError(f"{where}: {exc}") from exc


def _fuel(value: Any, where: str) -> Fuel:
    if isinstance(value, str):
        fuel = FUEL_ALIASES.get(value.strip().lower())
        if fuel is None:
            raise CaseFileError(
                f"{where}: unknown fuel '{value}'. Use one of {sorted(FUEL_ALIASES)} "
                "or define it as a table with the fields of combusta.thermo.fuels.Fuel."
            )
        return fuel
    if isinstance(value, dict):
        return _build(Fuel, value, f"{where}.fuel")
    raise CaseFileError(f"{where}: 'fuel' must be a fuel name or a table.")


def _column_map(
    table: Mapping[str, Any], allowed: Mapping[str, tuple[str, ...]], where: str
) -> dict[str, ColumnSpec]:
    _check_keys(table, set(allowed), where)
    out: dict[str, ColumnSpec] = {}
    for key, value in table.items():
        if isinstance(value, str):
            out[key] = ColumnSpec(column=value)
        elif isinstance(value, dict):
            out[key] = _build(ColumnSpec, value, f"{where}.{key}")
        else:
            raise CaseFileError(f"{where}.{key}: expected a header name or {{ column, unit }}.")
    return out


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def resolve_case_path(path: str | Path) -> Path:
    """Return the case file for ``path``: the file itself or ``<dir>/case.toml``."""
    p = Path(path)
    if p.is_dir():
        p = p / CASE_FILE_NAME
    if not p.is_file():
        raise CaseFileError(f"Case file not found: {p}")
    return p.resolve()


def load_case(path: str | Path) -> LoadedCase:
    """Load an engine case from a ``case.toml`` file or a directory containing one.

    The trace file path is resolved relative to the case file. Trace format
    ``"forte_thermo"`` uses :func:`~combusta.io.forte_parser.parse_thermo_csv` (volume
    required), format ``"csv"`` uses :func:`~combusta.io.trace_csv.parse_trace_csv`
    (volume optional; taken from the geometry when absent).

    Raises
    ------
    CaseFileError
        Missing or unknown keys, wrong types, unknown fuel or trace format, or values
        rejected by the Combusta data classes.
    ForteParseError, UnitDetectionError, PhysicalBoundsError
        Problems in the trace file itself.
    """
    case_file = resolve_case_path(path)
    where = str(case_file)
    try:
        data = tomllib.loads(case_file.read_text(encoding="utf-8"))
    except tomllib.TOMLDecodeError as exc:
        raise CaseFileError(f"{where}: invalid TOML: {exc}") from exc
    _check_keys(data, _TOP_LEVEL_KEYS, where)

    name = data.get("name")
    if not isinstance(name, str) or not name.strip():
        raise CaseFileError(f"{where}: 'name' must be a non-empty string.")
    configuration = data.get("configuration")
    if configuration is not None and (not isinstance(configuration, str) or not configuration.strip()):
        raise CaseFileError(f"{where}: 'configuration' must be a non-empty string.")
    source = data.get("source")
    if source not in ("cfd", "test_bench", "model_1d", "synthetic"):
        raise CaseFileError(
            f"{where}: 'source' must be one of cfd, test_bench, model_1d, synthetic; got {source!r}."
        )

    geo_where = f"{where} [geometry]"
    geometry = _build(
        SliderCrankGeometry, _numbers(_table(data, "geometry", where, required=True), geo_where), geo_where
    )

    op_table = dict(_table(data, "operating_point", where, required=True))
    if "fuel" not in op_table:
        raise CaseFileError(f"{where} [operating_point]: 'fuel' is required.")
    fuel = _fuel(op_table.pop("fuel"), f"{where} [operating_point]")
    operating_point = _build(
        OperatingPoint, _numbers(op_table, f"{where} [operating_point]"),
        f"{where} [operating_point]", fuel=fuel,
    )

    trace_table = _table(data, "trace", where, required=True)
    _check_keys(trace_table, _TRACE_KEYS, f"{where} [trace]")
    fmt = trace_table.get("format", "forte_thermo")
    if fmt not in TRACE_FORMATS:
        raise CaseFileError(f"{where} [trace]: 'format' must be one of {TRACE_FORMATS}, got {fmt!r}.")
    rel = trace_table.get("file")
    if not isinstance(rel, str) or not rel:
        raise CaseFileError(f"{where} [trace]: 'file' is required.")
    trace_file = (case_file.parent / rel).resolve()
    if not trace_file.is_file():
        raise CaseFileError(f"{where} [trace]: trace file not found: {trace_file}")
    sign = trace_table.get("wall_heat_loss_sign", 1.0)
    if isinstance(sign, bool) or sign not in (1, -1):
        raise CaseFileError(f"{where} [trace]: 'wall_heat_loss_sign' must be +1 or -1.")
    sign = float(sign)
    columns_table = _table(trace_table, "columns", f"{where} [trace]", required=False)

    if fmt == "forte_thermo":
        cmap = _column_map(columns_table, THERMO_FIELDS, f"{where} [trace.columns]")
        parsed = parse_thermo_csv(trace_file, columns=cmap)
        volume, temperature, wall = parsed.volume_m3, parsed.temperature_k, parsed.wall_heat_loss_j_per_deg
        theta, pressure, units = parsed.crank_angle_deg, parsed.pressure_pa, parsed.source_units or {}
    else:
        cmap = _column_map(columns_table, TRACE_FIELDS, f"{where} [trace.columns]")
        table = parse_trace_csv(trace_file, columns=cmap)
        volume, temperature, wall = table.volume_m3, table.temperature_k, table.wall_heat_loss_j_per_deg
        theta, pressure, units = table.crank_angle_deg, table.pressure_pa, table.source_units

    try:
        trace = CycleTrace(
            crank_angle_deg=theta,
            pressure_pa=pressure,
            volume_m3=volume,
            temperature_k=temperature,
            wall_heat_loss_j_per_deg=None if wall is None else sign * wall,
        )
    except (ValueError, CombustaError) as exc:
        raise CaseFileError(f"{trace_file}: {exc}") from exc

    metadata_table = _table(data, "metadata", where, required=False)
    metadata = {k: str(v) for k, v in metadata_table.items()}
    if isinstance(data.get("description"), str):
        metadata["description"] = data["description"]
    if configuration is not None:
        metadata["configuration"] = configuration.strip()
    digest = _sha256(trace_file)
    metadata.update({f"unit:{k}": v for k, v in units.items()})
    metadata.update({
        "case_file": str(case_file),
        "trace_file": str(trace_file),
        "trace_format": fmt,
        "trace_sha256": digest,
        "wall_heat_loss_sign": f"{sign:+.0f}",
    })

    uncertainties = _numbers(_table(data, "uncertainty", where, required=False), f"{where} [uncertainty]")
    if any(u < 0.0 for u in uncertainties.values()):
        raise CaseFileError(f"{where} [uncertainty]: uncertainties must be >= 0.")

    try:
        case = EngineCase(
            name=name.strip(),
            geometry=geometry,
            operating_point=operating_point,
            trace=trace,
            source=source,
            scalars=_numbers(_table(data, "scalars", where, required=False), f"{where} [scalars]"),
            emissions_kg_per_cycle=_numbers(
                _table(data, "emissions_kg_per_cycle", where, required=False),
                f"{where} [emissions_kg_per_cycle]",
            ),
            metadata=metadata,
        )
    except (ValueError, CombustaError) as exc:
        raise CaseFileError(f"{where}: {exc}") from exc

    pre_chamber = _pre_chamber(data, where, case_file, trace_file, theta, sign)
    if pre_chamber is not None:
        metadata_extra = {"pre_chamber_file": str(pre_chamber.file), "pre_chamber_sha256": pre_chamber.sha256}
        case = EngineCase(
            name=case.name, geometry=case.geometry, operating_point=case.operating_point, trace=case.trace,
            source=case.source, scalars=case.scalars, emissions_kg_per_cycle=case.emissions_kg_per_cycle,
            metadata={**case.metadata, **metadata_extra},
        )

    return LoadedCase(
        case=case,
        uncertainties=uncertainties,
        case_file=case_file,
        trace_file=trace_file,
        trace_sha256=digest,
        configuration=None if configuration is None else configuration.strip(),
        pre_chamber=pre_chamber,
    )


def _pre_chamber(
    data: Mapping[str, Any],
    where: str,
    case_file: Path,
    trace_file: Path,
    theta: FloatArray,
    sign: float,
) -> PreChamberInput | None:
    """Read the optional ``[pre_chamber]`` table and its signals onto the case grid."""
    table = _table(data, "pre_chamber", where, required=False)
    if not table:
        return None
    pc_where = f"{where} [pre_chamber]"
    _check_keys(table, _PRE_CHAMBER_KEYS, pc_where)
    geo_values = {k: table[k] for k in ("volume_m3", "orifice_count", "orifice_diameter_m", "discharge_coefficient")
                  if k in table}
    if isinstance(geo_values.get("orifice_count"), float) and float(geo_values["orifice_count"]).is_integer():
        geo_values["orifice_count"] = int(geo_values["orifice_count"])
    geometry = _build(PreChamberGeometry, geo_values, pc_where)
    rel = table.get("file")
    pc_file = trace_file if rel is None else (case_file.parent / str(rel)).resolve()
    if not pc_file.is_file():
        raise CaseFileError(f"{pc_where}: file not found: {pc_file}")
    columns = _column_map(_table(table, "columns", pc_where, required=False), TRACE_FIELDS, f"{pc_where}.columns")
    if pc_file == trace_file and "pressure" not in columns:
        raise CaseFileError(
            f"{pc_where}: the pre-chamber signals come from the trace file, so "
            "[pre_chamber.columns] must name the pre-chamber pressure column explicitly."
        )
    signals = parse_trace_csv(pc_file, columns=columns)
    th = signals.crank_angle_deg
    if th[0] > theta[0] + 1e-9 or th[-1] < theta[-1] - 1e-9:
        raise CaseFileError(
            f"{pc_where}: pre-chamber signals cover [{th[0]:g}, {th[-1]:g}] deg but the case trace "
            f"covers [{theta[0]:g}, {theta[-1]:g}] deg."
        )

    def on_grid(values: FloatArray | None) -> FloatArray | None:
        if values is None:
            return None
        if th.shape == theta.shape and np.allclose(th, theta, rtol=0.0, atol=1e-9):
            return values
        return np.asarray(np.interp(theta, th, values), dtype=np.float64)

    pressure = on_grid(signals.pressure_pa)
    assert pressure is not None
    wall = on_grid(signals.wall_heat_loss_j_per_deg)
    return PreChamberInput(
        geometry=geometry,
        pressure_pa=pressure,
        temperature_k=on_grid(signals.temperature_k),
        wall_heat_loss_j_per_deg=None if wall is None else sign * wall,
        file=pc_file,
        sha256=_sha256(pc_file),
    )
