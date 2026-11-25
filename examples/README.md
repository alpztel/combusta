# Examples

## `altitude_sweep_demo.py` — synthetic aero-diesel altitude sweep

```bash
python examples/altitude_sweep_demo.py
```

A fully synthetic, end-to-end demonstration of the Combusta workflow for a single-cylinder aero-diesel at 0 ft, 15,000 ft, and 30,000 ft. All inputs are textbook models and illustrative assumptions; the outputs are not engineering results.

The script chains:

1. ISO 2533 ambient conditions and a boost-controlled two-stage turbocharger (critical altitude of about 21,000 ft) with charge-air cooling.
2. A motored closed-cycle single-zone simulation with NASA dry-air properties and Woschni heat transfer.
3. The Livengood-Wu integral of the Assanis et al. (2003) ignition delay on the motored state, which sets the start of combustion.
4. A fired cycle with Wiebe heat release from the predicted start of combustion, synthetic transducer noise, and zero-phase filtering.
5. Inverse heat-release analysis, burn fractions, IMEP, spray tip penetration, and the Second-Law exergy budget with its independent closure check.
6. The Grid Convergence Index procedure applied to a manufactured three-grid triplet (an illustration of the method, not a mesh study).

It prints a Markdown summary table and writes four figures to `../assets/`. The README section [Example: synthetic altitude sweep](../README.md#example-synthetic-altitude-sweep) discusses the results.

## `make_example_cases.py` and `cases/` — case files and the command line

```bash
python examples/make_example_cases.py
combusta compare examples/cases/synthetic_slow_burn examples/cases/synthetic_fast_burn --out examples/output
```

Two manufactured hydrogen cases (λ = 2, 2000 rpm) whose pressure traces were integrated forward from a prescribed Wiebe heat release (`combusta.synthetic.manufactured_cycle`). They differ only in burn start and duration and do not represent any particular combustion system. `synthetic_slow_burn` uses the ANSYS Forte cylinder-average layout (MPa, cm³, wall heat loss reported as negative); `synthetic_fast_burn` is a test-bench style table (bar, crank angle and pressure only). Both recover the prescribed heat release of 831.83 J, which makes them a check of unit conversion, sign handling, and temperature reconstruction. The uncertainties in their case files are illustrative, not the result of a grid study. The `cases/` directory is committed; regenerate it with the script. The report goes to `examples/output/`, which is not tracked.
