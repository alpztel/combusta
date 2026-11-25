# Accessibility

Combusta is a Python library with a command-line interface. It has no graphical user interface; its outputs are text, Markdown reports, CSV tables, and PNG figures.

## What is supported

- **Text first.** All scalar results are available as text: the command line prints the result tables, and each report is accompanied by a CSV file with the values in SI and display units. The crank-angle curves shown in the figures come from the input traces (CSV files) and are available as arrays from the Python API (for example `CycleMetrics.heat_release` and `PreChamberResult`).
- **Plain formats.** Reports are Markdown with standard headings and tables; figures are embedded with a short descriptive alternative text.
- **Colour.** Figures use the colours #004488, #DDAA33, and #BB5566 for the first three series, which form Paul Tol's high-contrast scheme designed to remain distinguishable for readers with colour-vision deficiencies. Every figure has a legend, and markers (circles in sweep plots) or line styles (dotted CA50 and MBT markers, dashed fits) carry additional meaning.
- **Documentation.** The README and the files in `docs/` are plain Markdown; equations are written in LaTeX source, which remains readable as text where it is not rendered.

## Known limitations

- Figures with more than three series rely on colour and legend order to distinguish the fourth and later series.
- Alternative texts of figures in reports are short labels, not descriptions of the curves.
- The figure style (font sizes, line widths) is fixed for printed publications and is not configurable yet.

## Reporting a barrier

If something in Combusta's output or documentation is hard or impossible for you to use, please open an [issue](https://github.com/alpztel/combusta/issues/new/choose) describing the output, what you tried, and the assistive technology you use. Accessibility problems are treated as bugs.
