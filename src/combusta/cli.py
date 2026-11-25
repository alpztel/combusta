"""Command-line interface: ``combusta inspect``, ``combusta compare``, ``combusta matched``.

``combusta inspect FILE`` lists the columns of a CSV export, their units, and the
canonical fields they are matched to; use it on the first file of a new format and
put the result into ``[trace.columns]`` of the case file.

``combusta compare CASE CASE [...]`` loads case files, compares them, and writes a
Markdown report, a metrics CSV, and an overlay figure (see :mod:`combusta.report`).

``combusta matched CASE ... (--ca50 DEG | --mbt)`` groups case files by their
``configuration`` key into spark-timing sweeps and compares the configurations at
matched CA50 or at their own MBT timing (see :mod:`combusta.sweep`).
"""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from pathlib import Path
import sys

from combusta.exceptions import CombustaError


def _inspect(args: argparse.Namespace) -> int:
    from combusta.io.forte_parser import DYNAMIC_FIELDS, SPRAY_FIELDS, THERMO_FIELDS, inspect_columns
    from combusta.io.trace_csv import TRACE_FIELDS

    tables = {"thermo": THERMO_FIELDS, "spray": SPRAY_FIELDS, "dynamic": DYNAMIC_FIELDS, "trace": TRACE_FIELDS}
    fields = None if args.kind == "all" else tables[args.kind]
    infos = inspect_columns(args.file, fields)
    width = max(len(i.header) for i in infos)
    print(f"{len(infos)} columns in {args.file}\n")
    print(f"{'header':<{width}}  {'name':<24} {'unit':<10} matched field(s)")
    for info in infos:
        matched = ", ".join(info.matched_fields) or "-"
        print(f"{info.header:<{width}}  {info.name:<24} {info.unit or '-':<10} {matched}")
    print("\nUnmatched or unit-less columns can be mapped in case.toml, e.g.\n"
          "[trace.columns]\npressure = { column = \"<header>\", unit = \"MPa\" }")
    return 0


def _compare(args: argparse.Namespace) -> int:
    from combusta.report import write_comparison_report

    baseline: int | str = args.baseline if args.baseline is not None else 0
    report = write_comparison_report(
        args.cases,
        args.out,
        baseline=baseline,
        figure=not args.no_figure,
        combustion_window_deg=(float(args.combustion_window[0]), float(args.combustion_window[1])),
        gross_window_deg=(float(args.gross_window[0]), float(args.gross_window[1])),
        revolutions_per_cycle=int(args.revolutions),
    )
    print(report.comparison.to_markdown())
    if report.warnings:
        print("\nData checks:")
        for w in report.warnings:
            print(f"  - {w}")
    print(f"\nWrote {report.report_path}, {report.metrics_path}"
          + ("" if args.no_figure else f", {report.figure_path}"))
    return 0


def _matched(args: argparse.Namespace) -> int:
    from combusta.report import write_matched_report

    baseline: int | str = args.baseline if args.baseline is not None else 0
    report = write_matched_report(
        args.cases,
        args.out,
        ca50_deg=args.ca50,
        at_mbt=args.mbt,
        baseline=baseline,
        fit_points=args.fit_points,
        figure=not args.no_figure,
        combustion_window_deg=(float(args.combustion_window[0]), float(args.combustion_window[1])),
        gross_window_deg=(float(args.gross_window[0]), float(args.gross_window[1])),
        revolutions_per_cycle=int(args.revolutions),
    )
    print(report.comparison.to_markdown())
    if report.warnings:
        print("\nData checks:")
        for w in report.warnings:
            print(f"  - {w}")
    print(f"\nWrote {report.report_path}, {report.metrics_path}"
          + ("" if args.no_figure else f", {report.figure_path}"))
    return 0


def _add_analysis_options(p: argparse.ArgumentParser) -> None:
    p.add_argument("--baseline", help="Baseline name (default: the first).")
    p.add_argument("--combustion-window", nargs=2, type=float, default=(-40.0, 90.0),
                   metavar=("START", "END"), help="Burn-fraction window [deg ATDC].")
    p.add_argument("--gross-window", nargs=2, type=float, default=(-180.0, 180.0),
                   metavar=("START", "END"), help="Gross indicated work window [deg ATDC].")
    p.add_argument("--revolutions", type=int, choices=(1, 2), default=2,
                   help="Crank revolutions per cycle (2 = four-stroke).")
    p.add_argument("--no-figure", action="store_true", help="Do not write the figure.")


def build_parser() -> argparse.ArgumentParser:
    """Return the argument parser of the ``combusta`` command."""
    from combusta import __version__

    parser = argparse.ArgumentParser(prog="combusta", description="Combusta engine case analysis.")
    parser.add_argument("--version", action="version", version=f"combusta {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    p_inspect = sub.add_parser("inspect", help="List the columns of a CSV export and their matches.")
    p_inspect.add_argument("file", type=Path)
    p_inspect.add_argument("--kind", choices=("all", "thermo", "spray", "dynamic", "trace"), default="all",
                           help="Field table to match against (default: all).")
    p_inspect.set_defaults(func=_inspect)

    p_compare = sub.add_parser("compare", help="Compare engine cases defined in case.toml files.")
    p_compare.add_argument("cases", nargs="+", type=Path, help="case.toml files or directories containing one.")
    p_compare.add_argument("--out", type=Path, default=Path("combusta_report"), help="Output directory.")
    _add_analysis_options(p_compare)
    p_compare.set_defaults(func=_compare)

    p_matched = sub.add_parser(
        "matched",
        help="Compare configurations (spark-timing sweeps) at matched CA50 or at MBT timing.",
    )
    p_matched.add_argument("cases", nargs="+", type=Path,
                           help="case.toml files or directories; grouped by their 'configuration' key.")
    mode = p_matched.add_mutually_exclusive_group(required=True)
    mode.add_argument("--ca50", type=float, metavar="DEG", help="Compare at this CA50 [deg ATDC].")
    mode.add_argument("--mbt", action="store_true", help="Compare each configuration at its MBT timing.")
    p_matched.add_argument("--fit-points", type=int, default=5,
                           help="Sweep points around the IMEP maximum used for the MBT fit (default 5).")
    p_matched.add_argument("--out", type=Path, default=Path("combusta_matched"), help="Output directory.")
    _add_analysis_options(p_matched)
    p_matched.set_defaults(func=_matched)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Entry point of the ``combusta`` command; returns the process exit code."""
    for stream in (sys.stdout, sys.stderr):
        # Tables contain non-ASCII symbols (Δ, λ); legacy Windows consoles default to cp1252.
        reconfigure = getattr(stream, "reconfigure", None)
        if reconfigure is not None:
            reconfigure(errors="replace")
    parser = build_parser()
    args = parser.parse_args(argv)
    if getattr(args, "command", None) == "compare" and len(args.cases) < 2:
        parser.error("compare needs at least two cases.")
    if getattr(args, "command", None) == "matched" and len(args.cases) < 4:
        parser.error("matched needs at least two configurations with two cases each.")
    try:
        return int(args.func(args))
    except (CombustaError, ValueError, KeyError, OSError) as exc:
        print(f"combusta: error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
