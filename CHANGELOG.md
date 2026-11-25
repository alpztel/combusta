# Changelog

All notable changes to **`combusta`** are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Documentation
- README: the GCI sign function uses `\mathrm{sgn}` (GitHub math does not allow `\operatorname`, so the equation block did not render); scope, Turkish summary, quick start (measurement uncertainty), package table, limitations, development notes, and citation version brought up to date.
- CONTRIBUTING: module map, SI boundary rules for case files and reports, test-speed guidance; SECURITY: scope covers case files, generic traces, and written reports, and states that case files are trusted input; CODE_OF_CONDUCT: private reporting path clarified.
- Issue and pull-request templates updated; a contact link to the case-file reference added.

## [0.4.0] - 2026-11-26

From exported data to a comparison report, and the analyses needed to compare an open-chamber and a pre-chamber engine: case files and a command line, comparison at matched combustion phasing or MBT timing, a two-volume pre-chamber analysis, and the uncertainty of pressure-derived metrics from measurement errors. All new models are verified against closed-form results and manufactured solutions; the library has **still not been validated against engine or CFD measurements**.

### Added (measurement uncertainty)
- `combusta.measurement`: `PressureMeasurementUncertainty` (TDC phasing, pegging offset, gain, noise), `perturb_pressure`, and `pressure_measurement_uncertainty` with a linear uncertainty budget (central differences, JCGM 100:2008) and a Monte Carlo propagation (JCGM 101:2008) reported with the standard error of its standard deviation. A traced temperature is replaced by the ideal-gas reconstruction when the trapped mass is known.
- `[measurement]` table in case files, `combusta measurement` (with `--savgol` to reproduce the filtering of real data), `write_measurement_report`.
- `tests/test_15_measurement_uncertainty.py`.
- Documentation: README Section 8 (measurement uncertainty), Brunt & Emtage (1996), Randolph (1990), and JCGM 101:2008 references; case-file reference.

### Added (pre-chamber)
- `combusta.prechamber`: isentropic orifice flow with choking (`orifice_mass_flow`, `critical_pressure_ratio`, `transfer_rate`), open-system apparent heat release with mass exchange (`zone_apparent_heat_release`), two-volume analysis `analyze_pre_chamber` (zone masses from exported masses, exported temperatures, or the orifice model integrated from the trapped mass), `pre_chamber_metrics` (pressure difference, jet start/end, jet mass, peak jet flow, choked share, heat release of both volumes and pre-chamber share, burn fractions), and `estimate_discharge_coefficient`.
- `[pre_chamber]` table in case files (`PreChamberInput`), `combusta prechamber`, `write_pre_chamber_report`, `plot_pre_chamber`.
- `manufactured_pre_chamber_cycle` and `write_pre_chamber_case_directory` for verification.
- `tests/test_14_prechamber.py`.
- Documentation: README Section 10 (later sections renumbered), Toulson et al. (2010) reference, case-file and checklist updates.

### Added (matched comparison)
- `combusta.sweep`: `build_sweep`, `metrics_at` (PCHIP interpolation to a target CA50 or spark timing, no extrapolation), `find_mbt` (quadratic least-squares vertex of gross IMEP with delta-method uncertainty), and `compare_matched` (configurations at matched CA50 or at their own MBT timing, with the significance test of `compare_cases`).
- `configuration` key in `case.toml`; `combusta matched CASE ... (--ca50 DEG | --mbt)`; `write_matched_report` and `plot_sweep_overview`.
- `ignition_timing` metric in cycle metrics; `manufactured_cycle(ignition_timing_deg=...)`.
- `tests/test_13_matched_comparison.py`.
- Documentation: matched comparison in the README (Section 9, quick start, verification table), `configuration` and `combusta matched` in `docs/case-file.md`, sweep planning in `docs/forte-output-checklist.md`, and `examples/matched_sweep_demo.py`.

### Added (data intake and command line)
- **Case files** (`combusta.casefile`): `load_case` builds an `EngineCase` from a `case.toml` (geometry, operating point, fuel by name or as a table, trace file and format, column map, wall heat-loss sign, scalars, emissions, uncertainties, metadata). Unknown keys are rejected. The trace file's SHA-256 is recorded. Format reference in `docs/case-file.md`.
- **Column mapping** for all readers: `ColumnSpec(column, unit)` and a `columns=` argument on `parse_thermo_csv`, `parse_spray_csv`, and `parse_dynamic_csv`; `inspect_columns` lists headers, units, and matched fields. Error messages name the available columns. Repeated headers get a numeric suffix instead of overwriting each other.
- `parse_trace_csv` (`combusta.io.trace_csv`): generic crank-angle table reader for test-bench and 1D-model data; only crank angle and pressure are required.
- **Reports and command line**: `write_comparison_report` (`combusta.report`) and the `combusta` command with `inspect` and `compare`; the report contains provenance, data checks (volume consistency, missing inputs, unused uncertainty keys), the comparison table, analysis settings and notes, a metrics CSV, and a pressure / heat-release overlay (`plot_case_overlay`).
- `combusta.synthetic`: `manufactured_cycle` (forward-integrated single-zone cycle with Wiebe heat release, composition-dependent gas properties, and optional wall heat loss) and `write_case_directory`, for testing an analysis chain before real data exist.
- `docs/forte-output-checklist.md`: outputs a CFD run must write for case comparison, grid studies, and pre-chamber analysis.
- `examples/make_example_cases.py` and two manufactured example cases in `examples/cases/`.
- `tests/test_12_data_intake.py`.

### Changed
- Manufactured cycles (`manufactured_cycle`, `manufactured_pre_chamber_cycle`) evaluate their right-hand side with plain-float NASA polynomials and slider-crank kinematics instead of the array APIs; the same equations, about 10 and 35 times faster. The test suite uses the same helpers in its own forward models and runs pyrefly once per session instead of once per package; the local suite time drops from about 72 s to 41 s.
- CI no longer runs the altitude demo as a separate step (the test suite already runs it) and prints the 15 slowest tests.
- CI runs on pull requests and on pushes to `main` only (previously on every push as well, which ran the suite twice per pull request), cancels superseded pull-request runs, and can be started manually.

### Fixed
- `write_case_directory` writes LF line endings on every platform, so the SHA-256 of a written trace no longer depends on the operating system (#17).
- Without a temperature trace, `compute_cycle_metrics` reconstructed the temperature with the unburned gas constant. It now passes the trapped mass to `compute_heat_release(trapped_mass_kg=...)`, which evaluates $T = PV/(m R(x_b))$ inside the $x_b$ iteration. For hydrogen at λ = 2 the old behaviour overstated the gross heat release of a manufactured case by 0.4 %.

## [0.3.0] - 2026-11-25

First release under the name **Combusta**. Versions 0.1.0 and 0.2.0 were published as `ice-cfd-validation-assistant`. This release has been verified against closed-form solutions and published reference data, but **not yet validated against engine or CFD measurements**.

### Repository
- Security policy with private vulnerability reporting, a code of conduct (Contributor Covenant 2.1), rewritten contribution guidelines, issue forms (bug report, physics/reference issue, feature request), and a pull request template.
- Dependabot for Python dependencies and GitHub Actions; CodeQL code scanning.
- Release workflow that builds the sdist and wheel, checks them with `twine`, smoke-tests the wheel outside the source tree, attaches them to the GitHub release, and generates build provenance attestations.
- Package classifier changed from "Beta" to "Alpha" to reflect the validation status.

### Documentation
- README rewritten: scope, quick start, package overview, a complete mathematical formulation, a verification table that names the reference of every check, assumptions and limitations, and a verified reference list. Third-party license notices moved from `LICENSE` (now the plain MIT text) into the README and corrected.
- Corrected citations: Gatowski et al. (1984) is SAE 841359; Caton (2000) is SAE 2000-01-1081; Rakopoulos & Giakoumis (2006) is a review in *Progress in Energy and Combustion Science*; the 2008 TDC paper by Pipitone et al. is SAE 2008-36-0059.
- `CITATION.cff` rewritten and validated against the CFF 1.2.0 schema.

### Changed (example)
- `examples/altitude_sweep_demo.py` rewritten: start of combustion now comes from the Livengood-Wu prediction on the motored state (it previously started before the predicted ignition), the turbocharger holds manifold pressure up to a critical altitude, the GCI figure is labeled as a manufactured example, and synthetic spray samples are no longer labeled as Forte output.

### Fixed
- **Hardenberg–Hase ignition delay** now uses the published form, in which the temperature and pressure terms multiply inside the exponent (Hardenberg & Hase, 1979; Heywood, 2018, Eq. 10.37). The former default added them, which lengthened delays by a factor of about 2.3 at typical conditions (0.46 ms vs. 0.20 ms at 900 K, 60 bar, CN 50). `pressure_coupling="additive"` remains for reproducing old results and emits a `DeprecationWarning`.
- The 1.24 MPa Hardenberg–Hase limit is described as the validity limit of the correlation (its pressure term is singular), no longer as a physical auto-ignition or misfire threshold.
- **Exergy budget**: exergy destruction is now computed from entropy generation (Gouy–Stodola) instead of as the balance residual, so `closure_error_j` is a real consistency check of the inputs (it was zero by construction). New `balance_residual_rate_j_per_deg`, `e_released_j`, `e_unreleased_j`, and `unreleased_share_pct`. Thermomechanical exergy is no longer clipped at zero.
- **GCI**: divergent triplets emit a `UserWarning`; the new `formal_order` argument caps an apparent order above the formal order of the scheme (Roache, 1998), reported via `order_used_p` and `order_limited`.
- Naber & Siebers (1996) penetration was mislabeled "Siebers (1999)"; Siebers (1999) is the liquid-length scaling, which is not implemented. The altitude spray risk docstring now states that tip penetration is not a liquid wall-wetting criterion.
- `compute_knock_metrics` documents that MAPO needs a local pressure signal; spatially averaged CFD pressure is not suitable.

### Changed
- Project renamed from `ice-cfd-validation-assistant` / `ice_cfd_validator` to **Combusta** (`combusta`), reflecting the move to a fuel-agnostic engine analysis, case comparison, and V&V toolkit. The GitHub repository is now `alpztel/combusta`; old URLs redirect.
- Base exception renamed to `CombustaError`; `ValidationAssistantError` remains as an alias.
- Removed the AST "no explicit loops" test audits; the pyrefly type check is kept.

### Added
- `combusta.thermo.fuels`: `Fuel` (any C<sub>x</sub>H<sub>y</sub>O<sub>z</sub>, including pure H<sub>2</sub>), stoichiometric air-fuel ratio, reactant/product species and NASA mixtures for any air-excess ratio λ (rich mixtures for carbon-free fuels), element bookkeeping, and the fuel library `HYDROGEN`, `METHANE`, `ISO_OCTANE`, `N_DODECANE`, `DIESEL_SURROGATE`, `JET_A1_SURROGATE`.
- GRI-Mech 3.0 NASA polynomials for H<sub>2</sub> and CH<sub>4</sub>.
- `blend_burned_unburned`: mass-weighted unburned/burned gas properties at mass fraction burned x<sub>b</sub>.
- `compute_heat_release(..., burned_mixture=...)`: x<sub>b</sub>-dependent γ resolved by fixed-point iteration (Gatowski et al., 1984; Brunt et al., 1998). Results expose `mass_fraction_burned` and `property_iterations`.
- `combusta.case`: source-independent `EngineCase`, `OperatingPoint`, and `CycleTrace`, with `EngineCase.from_forte_thermo`, provenance (`source`, including `"synthetic"` for verification fixtures), ideal-gas temperature reconstruction, and a geometry/volume consistency check that detects crank-angle offsets.
- `combusta.compare`: `compute_cycle_metrics` (peak pressure and phasing, max. pressure-rise rate, gross/net IMEP, gross indicated power, efficiency and ISFC, CA10/50/90, heat release) and `compare_cases` with baseline differences, significance against the combined numerical uncertainty U<sub>Δ</sub> = √(U<sub>a</sub>² + U<sub>b</sub>²), Markdown/DataFrame output, and a banner for synthetic data. `uncertainty_from_gci` converts a GCI result into an absolute uncertainty and rejects divergent triplets.
- `combusta.emissions`: exhaust amount and water fraction of complete combustion, ppm (wet/dry) to mass per cycle with NOx as NO<sub>2</sub> equivalent, specific emissions (g/kWh), emission indices (g/kg fuel), and combustion efficiency from exhaust chemical energy (Heywood, 2018, Eq. 4.69).
- `EngineCase.emissions_kg_per_cycle`: emission masses are turned into `is_<species>`, `ei_<species>`, and `combustion_efficiency` metrics and appear in `compare_cases` tables.
- `combusta.vv.validation` (ASME V&V 20-2009): `validate` (comparison error E = S − D, validation uncertainty u<sub>val</sub>, model-error interval E ± k·u<sub>val</sub>, scalar or crank-angle-resolved), `numerical_uncertainty_from_gci`, `sensitivity_coefficient`, `input_uncertainty` (first-order propagation), and `experimental_uncertainty` (GUM Type A + Type B).
- `tests/test_10_emissions_and_validation.py`: exhaust mass conservation, NIST molar masses, wet/dry conversion, Heywood combustion efficiency, closed-form GUM statistics, and exactness of central differences for quadratic responses.
- CI fix: `cumulative_trapezoid(..., initial=0)` for scipy-stubs 1.18.
- Combustion plausibility check (net vs. absolute apparent heat release) so motored or misfiring traces do not yield spurious CA50 values.
- `EngineCase.from_forte_thermo(..., wall_heat_loss_sign=...)` to align the exported wall heat-loss sign with the Combusta convention (positive = heat leaving the gas).
- `tests/test_09_case_comparison.py`: Otto-cycle efficiency limit, analytical Wiebe CA10/50/90, zero work of a reversible motored cycle, and the GCI bound on a power-law triplet.
- `tests/test_08_fuels_and_case.py`: checks against NIST-JANAF/ATcT data and Heywood (2018) stoichiometry, element conservation, and a manufactured-solution verification of the x<sub>b</sub>-dependent heat release.

## [0.2.0] - 2026-11-08

### Added
- **`combusta.ignition` (Step 6 — High-Altitude Auto-Ignition Kinetics & CFD Diagnostics)**:
  - `correlations.py`: Arrhenius-type empirical auto-ignition models (`arrhenius_state_delay`, `compute_apparent_activation_energy`, `scale_activation_energy`) with verified Wolfer (1938) and Assanis et al. (2003) parameters for EN590 Diesel and Jet-A1 / JP-8 / F-34.
  - Hardenberg & Hase (1979, SAE 790493) global ignition delay correlation (`hardenberg_hase_delay`, `HardenbergHaseResult`) supporting both classic additive and multiplicative pressure coupling, apparent activation energy scaling from Cetane Number ($E_a = 618840 / (\text{CN} + 25)$), and strict $P_{\text{TDC}} \le 1.24\ \text{MPa}$ (`12.4 bar`) singularity guards.
  - `livengood_wu.py`: Vectorized Livengood-Wu (1955) auto-ignition integral solver (`solve_livengood_wu`, `livengood_wu_integral`, `LivengoodWuResult`) on dynamic in-cylinder pressure and temperature histories. Features Hermite cubic sub-grid interpolation for off-node SOI injection angles, exact machine-precision root solving ($I(\theta_{\text{SOC}}) = 1.0$), and high-altitude misfire detection ($I(\theta) < 1.0$).
  - `cfd_detector.py`: Multi-criteria chronological ignition detector (`extract_cfd_ignition_metrics`, `CFDIgnitionMetrics`) resolving evaporative cooling trough energy ($Q_{\text{evap}}$ [J]), zero-crossing $\theta_{\text{SOC}}^{\text{AHRR0}}$, $\text{CA01}$, $\text{CA02}$, $\text{CA05}$, and maximum pressure curvature $\theta_{\max(d^2P/d\theta^2)}$.
  - Two-stage ignition cool-flame separator (`detect_two_stage_ignition`, `TwoStageIgnitionResult`) decomposing Low-Temperature Heat Release (LTHR), Negative Temperature Coefficient (NTC) dwell, and High-Temperature Heat Release (HTHR) peaks, with defensive single-stage false-positive prevention.
- **`combusta.exergy` (Step 7 — Second-Law Thermodynamics & Exergy Budget)**:
  - `dead_state.py`: Altitude-dependent thermodynamic dead-state generator (`DeadStateConditions.from_altitude`) via ISO 2533 extending from sea level through the $11\ \text{km}$ tropopause to the $15\ \text{km}$ isothermal stratosphere.
  - Kotas (1995) fuel chemical exergy calculator (`compute_fuel_chemical_exergy`, `FuelChemicalExergy`) supporting mass ratio $h/c$, atomic ratio $(H/C)$ with NIST IUPAC atomic weights ($M_H = 1.00794, M_C = 12.0107$), and predefined fuel factors (`DIESEL_EN590_EXERGY_FACTOR`, `JET_A1_EXERGY_FACTOR`, `HYDROGEN_EXERGY_FACTOR`).
  - `second_law.py`: Vectorized thermomechanical exergy engine (`compute_thermomechanical_exergy`) strictly positive-definite over $50 \times 50$ $(T, P)$ Gibbs convexity grids, with NASA RP-1311 14-coefficient polynomial entropy continuity ($|\Delta e_{\text{tm}}| < 0.01\ \text{J/kg}$ across $T = 1000\ \text{K}$) and exact Maxwell partial derivative identities.
  - Closed-cycle Second-Law exergy rate and cumulative budget solver (`compute_exergy_rates`, `compute_exergy_cycle_budget`, `ExergyRatesResult`, `ExergyCycleBudget`) proving Gouy-Stodola theorem ($E_{\text{dest}} = T_0 S_{\text{gen}}$) closure within $5 \times 10^{-7}\ \text{J}$, 100.0000% closed-cycle budget closure with $< 10^{-9}\ \text{J}$ error, and high-altitude exhaust blowdown exergy recovery potential ($+24.87\%$ gain at $30{,}000\ \text{ft}$).
- **Centralized Test Architecture & Pytest Markers (`tests/`)**:
  - Consolidated 11 fragmented test scripts into 7 modular, standard-compliant test suites (`test_01` through `test_07`) sharing central fixtures (`conftest.py`: `standard_engine_geom`, `generate_synthetic_rk45_cycle`, `assert_package_vectorized_and_typed`).
  - Configured custom pytest markers `unit` and `strict_audit` in `pyproject.toml` enabling fast unit runs (`pytest -m unit` in 4.35s) and rigorous analytical torture audits (`pytest -m strict_audit` in 19.62s).
  - All 65 tests passing with zero explicit element-wise loops in numerical packages and 0 errors in `pyrefly check`.

## [0.1.0] - 2026-11-07

### Added
- **`combusta.constants`**:
  - Exact CODATA 2018 / CIPM-2007 / NIST SP 811 fundamental constants (`R_UNIVERSAL`, `MOLAR_MASS_DRY_AIR`, `R_SPECIFIC_DRY_AIR`, `STD_GRAVITY_M_S2`, `STD_ATM_PA`, `ZERO_CELSIUS_K`, `FEET_TO_M`).
  - Eight strict-SI conversion dictionaries (`PRESSURE_TO_PA`, `VOLUME_TO_M3`, `LENGTH_TO_M`, `MASS_TO_KG`, `TIME_TO_S`, `ANGLE_TO_DEG`, `TEMP_TO_K`, `HEAT_RATE_TO_J_PER_DEG`) supporting metric, CGS (`dyne/cm^2`, `erg/CA`), and affine Celsius (`degC -> K`) tokens.
  - Vectorized ISO 2533:1975 International Standard Atmosphere model (`compute_isa_conditions`, `ISAConditions`) covering both the linear-lapse troposphere (`0–11,000 m`) and isothermal lower stratosphere (`11,000–20,000 m`).
- **`combusta.io`**:
  - `forte_parser.py`: Automated 1-row and 2-row header parsers (`parse_thermo_csv`, `parse_spray_csv`, `parse_dynamic_csv`) returning immutable `ForteThermoData`, `ForteSprayData`, and `ForteDynamicData` containers with strict positivity bounds and `source_units` audit provenance.
  - `ecn_loader.py`: Sandia Engine Combustion Network (ECN Spray A / Spray D) experimental loader (`load_ecn_spray_data`, `ECNSprayDataset`) supporting comma- and whitespace-delimited tables, `#`/`%` comment headers, pre-SOI timestamp trimming, and experimental standard-deviation (`u_D`) extraction.
- **`combusta.signal`**:
  - `savgol_zero_phase` and `butterworth_zero_phase` (`sosfiltfilt`) zero-phase pressure smoothing and acoustic bandpass filters preventing artificial `CA50` crank-angle phase lag.
  - Loop-free `compute_derivative` combining 5-point $\mathcal{O}(\Delta\theta^4)$ central/Fornberg boundary stencils on uniform grids with analytical `CubicSpline` differentiation on non-uniform adaptive CFD grids.
- **`combusta.thermo`**:
  - `kinematics.py`: Centered and wrist-pin-offset slider-crank solver (`SliderCrankGeometry`, `KinematicsState`, `compute_kinematics`, `compute_slider_crank_kinematics`) with exact analytical $dV/d\theta$ (`rad -> deg` Jacobian) and sub-degree motored TDC loss-angle auditor (`verify_tdc_phasing`, `TDCPhasingResult`).
  - `nasa_polynomials.py`: 14-coefficient (two-interval `200–1000 K` & `1000–6000 K`) NASA RP-1311 / GRI-Mech 3.0 ideal-gas property engine (`NASAGasMixture`, `DRY_AIR`, `DIESEL_BURNED_PRODUCTS`, `STOICHIOMETRIC_PRODUCTS`, `mix_nasa_species`, `make_combustion_mixture`, `compute_cp`, `compute_cv`, `compute_gamma`).
  - `heat_transfer.py`: Vectorized Woschni (1967, SAE 670931) and Hohenberg (1979, SAE 790825) in-cylinder wall heat-transfer models (`WoschniCorrelation`, `HohenbergCorrelation`, `compute_woschni_heat_transfer`, `compute_hohenberg_heat_transfer`).
  - `heat_release.py`: Heywood (2018) single-zone First-Law net and gross AHRR (`compute_heat_release`), sub-degree cubic-spline `CA10`/`CA50`/`CA90` MFB solver (`extract_burn_fractions`), Gross/Net `IMEP` & `PMEP` (`compute_imep`), and `MAPO` + Eng (2002) acoustic Ringing Intensity (`compute_knock_metrics`).
- **`combusta.spray`**:
  - Vectorized continuous ($C^0$) Hiroyasu & Arai (1990, SAE 900475) two-regime spray penetration (`hiroyasu_arai_penetration`).
  - Vectorized Newton-exact Naber & Siebers (1996, SAE 960034) / Siebers (1999) momentum scaling (`naber_siebers_penetration`).
  - High-altitude gas-density over-penetration and sub-step piston-bowl wall-impingement auditor (`evaluate_altitude_spray_risk`) plus asynchronous CFD error norms (`compute_spray_validation_metrics`: `RMSE`, `MAPE`, $L_2$, $L_\infty$).
- **`combusta.vv`**:
  - Roache (1998), Celik et al. (2008), and ASME V&V 20-2009 three-mesh verification module (`compute_grid_size`, `compute_gci`, `GridSizeResult`, `GCIResult`) featuring 3D cube-root grid sizing, non-uniform $q(p)$ fixed-point / Brent root-finding, and `NaN`-proof Adaptive Mesh Refinement (AMR) oscillatory convergence handling.
- **`combusta.plotting`**:
  - Headless Agg 300-DPI SAE/IEEE STIX-mathtext diagnostic figure generators (`plot_thermo_diagnostics`, `plot_pv_log_diagram`, `plot_altitude_spray_comparison`, `plot_gci_convergence`).
- **Verification, Showcase & CI/CD**:
  - End-to-end non-proprietary Aero-Diesel altitude sweep (`examples/altitude_sweep_demo.py`: `0 ft`, `15,000 ft`, `30,000 ft`) and exported 300-DPI figures in `assets/`.
  - 39 automated unit, analytical roundtrip, AST vectorization, memory-leak, and `pyrefly` static-typing audit tests under `tests/`.
  - Multi-version Ubuntu GitHub Actions workflow (`.github/workflows/ci_validation.yml`) for Python 3.11, 3.12, and 3.13.

[Unreleased]: https://github.com/alpztel/combusta/compare/v0.4.0...HEAD
[0.4.0]: https://github.com/alpztel/combusta/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/alpztel/combusta/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/alpztel/combusta/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/alpztel/combusta/releases/tag/v0.1.0
