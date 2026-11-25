"""Write two manufactured hydrogen cases for trying the case-file workflow.

Usage::

    python examples/make_example_cases.py
    combusta compare examples/cases/synthetic_slow_burn examples/cases/synthetic_fast_burn --out examples/output

Both cases are manufactured solutions (``source = "synthetic"``): the pressure traces
are integrated forward from a prescribed Wiebe heat release, so the analysis must return
that heat release. They differ only in burn rate and phasing, not in any modeled
combustion system, and are not engineering results.

``synthetic_slow_burn`` is written in the ANSYS Forte cylinder-average layout (MPa,
cm^3, wall heat loss reported as negative), ``synthetic_fast_burn`` as a test-bench
table (bar, crank angle and pressure only), so the example also exercises unit
conversion, sign handling, and volume and temperature reconstruction.
"""

from __future__ import annotations

from pathlib import Path

from combusta import HYDROGEN, SliderCrankGeometry, manufactured_cycle, write_case_directory

OUT = Path(__file__).resolve().parent / "cases"
GEOMETRY = SliderCrankGeometry(bore_m=0.086, stroke_m=0.086, conrod_m=0.145, compression_ratio=12.0)
# Illustrative numerical uncertainties (as if from a grid study); not derived from any study.
UNCERTAINTY = {"p_max": 3.0e4, "ca50": 0.3, "imep_gross": 2.0e3}


def main() -> None:
    slow = manufactured_cycle(
        "synthetic_slow_burn", GEOMETRY, HYDROGEN, lambda_air=2.0, rpm=2000.0,
        combustion_start_deg=-15.0, combustion_duration_deg=45.0, wall_heat_loss_j=60.0,
    )
    fast = manufactured_cycle(
        "synthetic_fast_burn", GEOMETRY, HYDROGEN, lambda_air=2.0, rpm=2000.0,
        combustion_start_deg=-8.0, combustion_duration_deg=28.0, wall_heat_loss_j=60.0,
    )
    for cycle, fmt, sign in ((slow, "forte_thermo", -1.0), (fast, "csv", 1.0)):
        path = write_case_directory(
            cycle.case, OUT / cycle.case.name, trace_format=fmt, fuel_key="hydrogen",
            wall_heat_loss_sign=sign, uncertainty=UNCERTAINTY,
        )
        print(f"wrote {path} (exact heat release {cycle.total_heat_release_j:.2f} J)")


if __name__ == "__main__":
    main()
