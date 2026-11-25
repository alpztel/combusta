"""Matched-CA50 and MBT comparison of two manufactured spark-timing sweeps.

Usage::

    python examples/matched_sweep_demo.py

Writes two synthetic sweeps of case files to ``examples/output/sweeps/`` and the matched
reports to ``examples/output/matched_ca50/`` and ``examples/output/matched_mbt/``
(``examples/output/`` is not tracked).

Both configurations use the same fuel, λ, speed, and fuel mass. ``slow_burn`` burns over
45 deg with a 6 deg ignition delay, ``fast_burn`` over 28 deg with a 4 deg delay. They are
manufactured solutions with a prescribed wall heat loss whose shape does not depend on
phasing; their MBT phasing (CA50 near TDC) is therefore not representative of real
engines, where heat transfer moves MBT CA50 later. The example demonstrates the method,
not a combustion-system result.
"""

from __future__ import annotations

from pathlib import Path
import sys

from combusta import HYDROGEN, SliderCrankGeometry, manufactured_cycle, write_case_directory, write_matched_report

ROOT = Path(__file__).resolve().parent / "output"
GEOMETRY = SliderCrankGeometry(bore_m=0.086, stroke_m=0.086, conrod_m=0.145, compression_ratio=12.0)
# Illustrative numerical uncertainties; not derived from a grid study.
UNCERTAINTY = {"imep_gross": 2.0e3, "eta_indicated_gross": 2.0e-3, "p_max": 3.0e4}
SWEEPS = {
    # name: (spark timings [deg ATDC], ignition delay [deg], burn duration [deg])
    "slow_burn": ((-36.0, -30.0, -24.0, -18.0, -12.0, -6.0), 6.0, 45.0),
    "fast_burn": ((-24.0, -18.0, -12.0, -6.0, 0.0, 6.0), 4.0, 28.0),
}


def main() -> None:
    reconfigure = getattr(sys.stdout, "reconfigure", None)
    if reconfigure is not None:  # tables contain Δ; legacy Windows consoles use cp125x
        reconfigure(errors="replace")
    case_dirs = []
    for name, (timings, delay, duration) in SWEEPS.items():
        for spark in timings:
            cycle = manufactured_cycle(
                f"{name}_{spark:+.0f}", GEOMETRY, HYDROGEN, lambda_air=2.0, rpm=2000.0,
                combustion_start_deg=spark + delay, combustion_duration_deg=duration,
                wall_heat_loss_j=60.0, ignition_timing_deg=spark,
            )
            path = write_case_directory(
                cycle.case, ROOT / "sweeps" / name / f"{spark:+.0f}", fuel_key="hydrogen",
                configuration=name, uncertainty=UNCERTAINTY,
            )
            case_dirs.append(path.parent)
        print(f"wrote sweep '{name}' ({len(timings)} cases)")

    ca50 = write_matched_report(case_dirs, ROOT / "matched_ca50", ca50_deg=8.0)
    mbt = write_matched_report(case_dirs, ROOT / "matched_mbt", at_mbt=True)
    print(ca50.comparison.to_markdown())
    print()
    print(mbt.comparison.to_markdown())
    print(f"\nReports: {ca50.report_path}, {mbt.report_path}")


if __name__ == "__main__":
    main()
