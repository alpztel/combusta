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
| `configuration` | string | no | Engine configuration the case belongs to (for example `open_chamber`). Cases with the same value form one spark-timing sweep for `combusta matched`; required there. |

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

### `[measurement]`

Standard uncertainties of the cylinder-pressure measurement, used by `combusta measurement` (command-line options override them). All keys are optional; each error is treated as normally distributed and independent of the others.

| Key | Unit | Meaning |
| :--- | :--- | :--- |
| `tdc_offset_deg` | deg | Crank-angle phasing (TDC) error |
| `pegging_offset_pa` | Pa | Absolute pressure referencing error |
| `gain_relative` | – | Relative transducer sensitivity error (0.01 = 1 %) |
| `noise_pa` | Pa | Standard deviation of uncorrelated sample noise |

```toml
[measurement]
tdc_offset_deg = 0.2
pegging_offset_pa = 5000.0
gain_relative = 0.01
noise_pa = 2000.0
```

```bash
combusta measurement bench_case --samples 500 --savgol 21 --out meas
```

The report gives the contribution of each systematic source (linear propagation) and the Monte Carlo standard deviation of all sources including noise. Use `--savgol WINDOW` when the real data are filtered, so that the noise is filtered the same way. The result is the measurement part of the experimental uncertainty; combine it with the cycle-to-cycle (Type A) uncertainty of the mean (`experimental_uncertainty`) before using it as $u_D$ in a validation.

### `[pre_chamber]`

Declares a pre-chamber engine for `combusta prechamber`. The case trace then describes the **main chamber only** (its pressure, volume, and optionally temperature and wall heat loss); the pre-chamber volume must not be part of it.

| Key | Unit | Required | Meaning |
| :--- | :--- | :---: | :--- |
| `volume_m3` | m³ | yes | Pre-chamber volume including the orifices |
| `orifice_count` | – | yes | Number of orifices |
| `orifice_diameter_m` | m | yes | Orifice diameter |
| `discharge_coefficient` | – | yes | $C_D$ in (0, 1]; no default, take it from measurement, CFD, or the fitted value in the report |
| `file` | – | no | File with the pre-chamber signals; defaults to the trace file |

`[pre_chamber.columns]` maps the pre-chamber signals with the same syntax as `[trace.columns]` (fields `crank_angle`, `pressure`, `temperature`, `wall_heat_loss`). If the signals are in the trace file itself, `pressure` must be mapped explicitly, otherwise the main-chamber pressure would be found by its alias. Signals on a different crank-angle grid are interpolated linearly onto the trace grid; they must cover it.

Zone masses come from the pre-chamber and main-chamber temperatures when both are exported; otherwise the orifice model is integrated with `operating_point.trapped_mass_kg`, which is then the total mass of both volumes. `lambda_air` is required (gas model of both volumes).

```toml
[pre_chamber]
volume_m3 = 1.0e-6
orifice_count = 6
orifice_diameter_m = 1.2e-3
discharge_coefficient = 0.7
file = "prechamber.csv"

[pre_chamber.columns]
pressure = "PC Pressure"
temperature = "PC Temperature"
```

### `[metadata]`

Free-form strings stored with the case (solver version, mesh name, run date). The loader adds `case_file`, `trace_file`, `trace_format`, `trace_sha256`, the source unit of every column, and `wall_heat_loss_sign`.

## Command line

```bash
combusta prechamber prechamber_case --out pc_report   # case with a [pre_chamber] table
combusta inspect open/thermo.csv
combusta compare open prechamber --out report --baseline open_chamber
```

`compare` accepts case files or directories containing `case.toml` and writes `report.md` (provenance with SHA-256 of every trace, data checks, comparison table, analysis settings and notes), `metrics.csv`, and `pressure_heat_release.png`. Options: `--combustion-window START END`, `--gross-window START END`, `--revolutions 1|2`, `--no-figure`.

### Matched comparison of spark-timing sweeps

```bash
combusta matched open/*/ prechamber/*/ --ca50 8 --out matched      # at a common CA50
combusta matched open/*/ prechamber/*/ --mbt --out matched_mbt     # each at its MBT timing
```

`matched` groups the cases by `configuration`; every configuration needs at least two cases, and an MBT comparison needs `ignition_timing_deg` in every case. All metrics are interpolated within each sweep with a monotone cubic (PCHIP); targets outside a sweep are rejected. The MBT timing is the vertex of a quadratic least-squares fit of gross IMEP over the `--fit-points` timings (default 5) around the sweep maximum; it is reported with its standard uncertainty when more than three points are fitted. Numerical uncertainties from `[uncertainty]` are interpolated linearly between the two bracketing cases. Output: `matched_report.md`, `matched_metrics.csv` (including `interpolation_spread_si`, the difference between cubic and linear interpolation, as an indicator of sweep resolution), and `sweep_overview.png`. The same analysis options as for `compare` apply, plus `--baseline`.
