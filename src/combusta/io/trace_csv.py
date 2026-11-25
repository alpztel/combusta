"""Generic crank-angle trace CSV reader for test-bench and 1D-model data.

Test-bench indicating systems and 1D engine models export crank-angle tables whose
only common content is crank angle and cylinder pressure; volume usually comes from
the geometry. :func:`parse_trace_csv` reads such a table with the same header, unit,
and column-mapping rules as the ANSYS Forte readers (:mod:`combusta.io.forte_parser`),
but requires only crank angle and pressure.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import IO, Final

import numpy as np

from combusta.constants import ANGLE_TO_DEG, HEAT_RATE_TO_J_PER_DEG, PRESSURE_TO_PA, VOLUME_TO_M3
from combusta.exceptions import PhysicalBoundsError
from combusta.io.forte_parser import (
    THERMO_FIELDS,
    ColumnMap,
    FloatArray,
    _check_column_map,
    _read_raw_columns,
    _resolve_scaled,
    _resolve_temperature,
)

TRACE_FIELDS: Final[Mapping[str, tuple[str, ...]]] = {
    key: THERMO_FIELDS[key]
    for key in ("crank_angle", "pressure", "volume", "temperature", "wall_heat_loss")
}
"""Canonical fields of a generic trace table and the header aliases searched for each."""


@dataclass(frozen=True, slots=True)
class TraceTableData:
    """SI-normalized crank-angle trace from a generic CSV table.

    Attributes
    ----------
    crank_angle_deg : NDArray[np.float64]
        Crank angle $\\theta$ [deg].
    pressure_pa : NDArray[np.float64]
        Absolute cylinder pressure $P$ [Pa].
    volume_m3, temperature_k, wall_heat_loss_j_per_deg : NDArray[np.float64] | None
        Optional cylinder volume [m^3], gas temperature [K], and wall heat-loss rate
        [J/deg] as exported (sign convention of the source).
    source_units : dict[str, str]
        Unit tokens found in the file, per canonical field.
    """

    crank_angle_deg: FloatArray
    pressure_pa: FloatArray
    volume_m3: FloatArray | None
    temperature_k: FloatArray | None
    wall_heat_loss_j_per_deg: FloatArray | None
    source_units: dict[str, str]


def parse_trace_csv(
    source: Path | str | IO[str],
    *,
    columns: ColumnMap | None = None,
) -> TraceTableData:
    """Read a crank-angle trace table and convert it to SI units.

    Parameters
    ----------
    source : Path | str | IO[str]
        Filesystem path or text stream. One header row with inline units
        (``"Pressure [bar]"``) or a header row followed by a unit row.
    columns : ColumnMap | None, optional
        Explicit header (and unit) per canonical field of :data:`TRACE_FIELDS`.

    Returns
    -------
    TraceTableData
        SI trace; crank angle and pressure are always present.

    Raises
    ------
    ForteParseError
        Missing crank-angle or pressure column, or non-numeric data.
    UnitDetectionError
        Missing or unknown unit of a used column.
    PhysicalBoundsError
        Non-positive pressure, volume, or temperature.
    """
    spec = _check_column_map(columns, TRACE_FIELDS)
    data = _read_raw_columns(source)

    ca_deg, ca_unit = _resolve_scaled(
        data, TRACE_FIELDS["crank_angle"], ANGLE_TO_DEG, "crank_angle",
        override=spec.get("crank_angle"),
    )
    p_pa, p_unit = _resolve_scaled(
        data, TRACE_FIELDS["pressure"], PRESSURE_TO_PA, "pressure", override=spec.get("pressure")
    )
    v_m3, v_unit = _resolve_scaled(
        data, TRACE_FIELDS["volume"], VOLUME_TO_M3, "volume", required=False,
        override=spec.get("volume"),
    )
    t_k, t_unit = _resolve_temperature(data, required=False, override=spec.get("temperature"))
    q_ht, q_unit = _resolve_scaled(
        data, TRACE_FIELDS["wall_heat_loss"], HEAT_RATE_TO_J_PER_DEG, "wall_heat_loss",
        required=False, override=spec.get("wall_heat_loss"),
    )
    assert ca_deg is not None and p_pa is not None and ca_unit is not None and p_unit is not None

    if np.any(p_pa <= 0.0):
        raise PhysicalBoundsError("Converted cylinder pressure contains values <= 0 Pa.")
    if v_m3 is not None and np.any(v_m3 <= 0.0):
        raise PhysicalBoundsError("Converted cylinder volume contains values <= 0 m^3.")
    if t_k is not None and np.any(t_k <= 0.0):
        raise PhysicalBoundsError("Converted gas temperature contains values <= 0 K.")

    units = {"crank_angle": ca_unit, "pressure": p_unit}
    for key, unit in (("volume", v_unit), ("temperature", t_unit), ("wall_heat_loss", q_unit)):
        if unit is not None:
            units[key] = unit
    return TraceTableData(
        crank_angle_deg=ca_deg,
        pressure_pa=p_pa,
        volume_m3=v_m3,
        temperature_k=t_k,
        wall_heat_loss_j_per_deg=q_ht,
        source_units=units,
    )
