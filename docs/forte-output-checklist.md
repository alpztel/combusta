# CFD output checklist (ANSYS Forte and similar solvers)

An analysis can only use what the simulation wrote out, and missing output usually means repeating the run. This checklist lists the quantities Combusta needs for a case comparison, a grid study, and a pre-chamber analysis. Decide on them **before** the simulations are started.

> [!NOTE]
> This document lists *what* is needed, not *where to click*. Output variable names and the way regions or monitors are defined differ between solver versions and setups; check them in the documentation of your solver version. Combusta's readers do not depend on particular header names: `combusta inspect` shows what a file contains, and `[trace.columns]` in the case file maps any header (see [case-file.md](case-file.md)).

## 1. Every case

Cylinder-averaged quantities over the closed part of the cycle (at least intake valve closing to exhaust valve opening), on one crank-angle grid:

| Quantity | Why | Required |
| :--- | :--- | :---: |
| Crank angle | Common abscissa; record the reference (deg ATDC firing) | yes |
| Pressure (absolute) | Heat release, IMEP, peak pressure, pressure-rise rate | yes |
| Volume | Work; checked against the declared geometry (a crank-angle offset shows up here) | strongly recommended |
| Mass-averaged temperature | Gas properties; otherwise reconstructed from trapped mass | recommended |
| Total mass in the cylinder | Trapped mass; needed if temperature is missing; reveals crevice/blow-by models | recommended |
| Wall heat-transfer rate | Gross heat release; **note the sign convention** | yes |
| Heat-release rate computed by the solver | Independent cross-check of Combusta's apparent heat release | recommended |
| Species masses or mole fractions at exhaust valve opening: NO, NO<sub>2</sub>, CO, CO<sub>2</sub>, unburned fuel (H<sub>2</sub> for hydrogen), H<sub>2</sub>O, O<sub>2</sub> | Specific emissions, emission indices, combustion efficiency | yes for emissions |

Case data that are not in the trace but belong in `case.toml`: bore, stroke, connecting-rod length, compression ratio, speed, fuel and λ, fuel mass per cycle, spark timing, intake pressure and temperature, wall temperatures.

## 2. Grid and time-step studies

A Grid Convergence Index needs the **same outputs** from at least three systematically refined meshes (refinement ratio above 1.3 per Celik et al., 2008), plus:

- the cell count history (or mean cell size) of each run,
- the base grid size and any adaptive-refinement and embedding settings.

Adaptive mesh refinement does not give a systematic refinement by itself; keep the refinement criteria proportional between the three runs and record them. The resulting uncertainty goes into `[uncertainty]` of each case.

## 3. Pre-chamber (turbulent jet ignition) cases

The pre-chamber analysis planned for Combusta 0.4 treats the pre-chamber and the main chamber as two control volumes connected by the orifices. It needs, **separately for each volume**:

| Quantity | Why |
| :--- | :--- |
| Volume-averaged pressure | Pressure difference across the orifices, jet timing |
| Mass-averaged temperature | Gas properties, orifice flow |
| Mass in the volume | Mass exchange through the orifices |
| Volume (constant for the pre-chamber) | First Law per volume |
| Heat-release rate and wall heat transfer, if the solver can report them per region | Cross-check and gross heat release per volume |
| Mass flow rate through the orifices, if available | Direct check of the orifice flow model |

If the pre-chamber is not defined as a separate region or monitor volume in the setup, its pressure cannot be recovered afterwards from cylinder averages. Record the orifice number, diameter, and length in the case metadata.

## 4. Simulation plan for a fair comparison

Comparing two combustion systems at the same spark timing compares them at different combustion phasing, because their burn rates differ. Plan each configuration as a **spark-timing sweep** (several timings bracketing the timing of maximum IMEP) at the same fuel mass and λ, so that the configurations can be compared at matched CA50 or at their own maximum-brake-torque timing. Comparison at matched conditions is planned for Combusta 0.4.

## 5. Export and first use

1. Export CSV tables with units in the headers (or a unit row), one file per case and quantity group.
2. Record the solver version, mesh, and run date in `[metadata]`.
3. Never commit licensed files (solver projects, restart files, tutorial geometry) to a public repository.
4. On the first exported file:
   - run `combusta inspect thermo.csv` and map unrecognized headers in `[trace.columns]`,
   - check the wall heat-loss sign: with the correct sign the cumulative wall heat loss is positive and the gross heat release is close to fuel energy times combustion efficiency; with the wrong sign the gross heat release is low by about twice the wall heat loss,
   - check the volume warning of the report (traced vs. declared geometry, threshold 1 %),
   - compare Combusta's gross heat release with the solver's own heat-release output.
