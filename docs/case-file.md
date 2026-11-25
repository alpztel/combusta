# Case file reference (`case.toml`)

A case file describes one engine configuration at one operating point: geometry, operating conditions, the crank-angle trace and how to read it, and optional scalar results, emissions, and numerical uncertainties. `combusta.load_case` turns it into an `EngineCase`; `combusta compare` turns several of them into a report.

Keys are the field names of the Combusta data classes and carry their SI unit in the name (`bore_m`, `intake_pressure_pa`). Unknown keys are rejected, so a typing error like `bore_mm` stops the analysis instead of being silently ignored.

## Minimal example

```toml
name = "open_chamber"
source = "cfd"                    # cfd | test_bench | model_1d | synthetic

[geometry]
bore_m = 0.086
stroke_m = 0.086
conrod_m = 0.145
compression_ratio = 12.0

[operating_point]
rpm = 2000.0
fuel = "hydrogen"
lambda_air = 2.0
fuel_mass_per_cycle_kg = 6.93e-6

[trace]
file = "thermo.csv"               # relative to this file
format = "forte_thermo"           # forte_thermo | csv
```

## Sections

### Top level

| Key | Type | Required | Meaning |
| :--- | :--- | :---: | :--- |
| `name` | string | yes | Case label in tables and figures; must be unique within one comparison. |
| `source` | string | yes | Provenance: `cfd`, `test_bench`, `model_1d`, or `synthetic`. Reports with a `synthetic` case carry a warning banner. |
| `description` | string | no | Stored in the case metadata. |

### `[geometry]`

Fields of `combusta.thermo.kinematics.SliderCrankGeometry`.

| Key | Unit | Required |
| :--- | :--- | :---: |
| `bore_m`, `stroke_m`, `conrod_m` | m | yes |
| `compression_ratio` | – | yes |
| `pin_offset_m` | m | no (0) |
| `crown_area_multiplier` | – | no (1) |

### `[operating_point]`

Fields of `combusta.case.OperatingPoint`.

| Key | Unit | Notes |
| :--- | :--- | :--- |
| `rpm` | rev/min | required |
| `fuel` | – | required; a name (below) or a table |
| `lambda_air` | – | needed for composition-dependent gas properties |
| `fuel_mass_per_cycle_kg` | kg | needed for efficiency, ISFC, and emission indices |
| `trapped_mass_kg` | kg | needed when the trace has no temperature |
| `ignition_timing_deg` | deg ATDC | spark timing or start of injection |
| `intake_pressure_pa`, `intake_temperature_k`, `wall_temperature_k` | Pa, K | optional |

Fuel names (case-insensitive): `hydrogen`/`H2`, `methane`/`CH4`, `iso-octane`/`isooctane`/`C8H18`, `n-dodecane`/`C12H26`, `diesel`, `jet-a1`/`jp-8`, or the full names in `combusta.FUEL_LIBRARY`. Any other C<sub>x</sub>H<sub>y</sub>O<sub>z</sub> fuel is defined as a table with the fields of `combusta.thermo.fuels.Fuel`:

```toml
[operating_point.fuel]
name = "ethanol"
carbon_atoms = 2.0
hydrogen_atoms = 6.0
oxygen_atoms = 1.0
lhv_j_per_kg = 26.8e6             # take the LHV from your fuel specification
```

### `[trace]`

| Key | Default | Meaning |
| :--- | :--- | :--- |
| `file` | – | Trace file, relative to the case file. |
| `format` | `forte_thermo` | `forte_thermo`: ANSYS Forte cylinder-average table; crank angle, pressure, and volume are required. `csv`: generic crank-angle table (test bench, 1D model); only crank angle and pressure are required, volume comes from the geometry when absent. |
| `wall_heat_loss_sign` | `1.0` | `+1` if the file reports heat leaving the gas as positive, `-1` if as negative. Combusta uses positive = heat leaving the gas. |

Both formats accept CSV files with units inline (`Pressure [MPa]`, `Volume (cm^3)`) or in a second header row. Headers are matched case-insensitively, with runs of non-alphanumeric characters treated as one underscore. A header that appears twice gets a numeric suffix in file order (`pressure`, `pressure_2`), so no column is dropped.

### `[trace.columns]`

Explicit header (and unit) per canonical field, for files whose headers are not recognized or carry no unit. Unmapped fields are still searched by their aliases.

```toml
[trace.columns]
crank_angle = { column = "Crank Angle", unit = "deg" }
pressure = { column = "Mean Pressure", unit = "MPa" }
volume = "Total Volume"           # header name only; the unit comes from the file
```

Canonical fields: `crank_angle`, `pressure`, `volume`, `temperature`, `wall_heat_loss`; `time` and `trapped_mass` are also accepted for `forte_thermo`. Run `combusta inspect <file>` to see the headers, units, and current matches of a file. Accepted unit tokens are listed in `combusta.constants` (`PRESSURE_TO_PA`, `VOLUME_TO_M3`, `TEMP_TO_K`, `ANGLE_TO_DEG`, `HEAT_RATE_TO_J_PER_DEG`, ...).

### `[scalars]`

Scalar results compared as given, for example a value read from a solver summary:

```toml
[scalars]
nox_ppm_evo = 118.0
```

### `[emissions_kg_per_cycle]`

Exhaust species mass per cycle and cylinder [kg]. NO<sub>x</sub> is given as NO<sub>2</sub> equivalent. They become specific emissions (`is_<species>`, g/kWh), emission indices (`ei_<species>`, g/kg fuel), and the combustion efficiency. Convert ppm values with `combusta.species_mass_from_mole_fraction`.

```toml
[emissions_kg_per_cycle]
NOx = 2.1e-8
H2 = 4.0e-8
```

### `[uncertainty]`

Absolute numerical uncertainty $U$ in SI units per metric key (`p_max` in Pa, `ca50` in deg, `imep_gross` in Pa, ...). `compare_cases` flags a difference as significant when it exceeds $\sqrt{U_a^2 + U_b^2}$; the test is made only when both cases give an uncertainty for that metric. Typical sources are a three-grid GCI study (`uncertainty_from_gci`) or a time-step study. Keys that match no reported metric are listed as a warning in the report.

```toml
[uncertainty]
p_max = 3.0e4
ca50 = 0.3
```

### `[metadata]`

Free-form strings stored with the case (solver version, mesh name, run date). The loader adds `case_file`, `trace_file`, `trace_format`, `trace_sha256`, the source unit of every column, and `wall_heat_loss_sign`.

## Command line

```bash
combusta inspect open/thermo.csv
combusta compare open prechamber --out report --baseline open_chamber
```

`compare` accepts case files or directories containing `case.toml` and writes `report.md` (provenance with SHA-256 of every trace, data checks, comparison table, analysis settings and notes), `metrics.csv`, and `pressure_heat_release.png`. Options: `--combustion-window START END`, `--gross-window START END`, `--revolutions 1|2`, `--no-figure`.
