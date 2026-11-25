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
