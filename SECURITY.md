# Security Policy

## Supported versions

Combusta is maintained by a single developer. Security fixes are made only on the latest minor release line.

| Version | Supported |
| :--- | :---: |
| 0.4.x | Yes |
| 0.3.x | No |
| < 0.3 (released as `ice-cfd-validation-assistant`) | No |

## Reporting a vulnerability

Please **do not** report security problems in public issues, discussions, or pull requests.

Report them privately through GitHub's [private vulnerability reporting](https://github.com/alpztel/combusta/security/advisories/new) ("Security" tab → "Report a vulnerability"). A useful report contains:

- the affected version or commit,
- a minimal reproduction built from **synthetic** data (see [Data hygiene](#data-hygiene)),
- the observed and the expected behaviour, and
- your assessment of the impact.

This is a volunteer-maintained research project. Reports are handled on a best-effort basis: the maintainer aims to acknowledge a report within 14 days and to agree on a disclosure timeline with the reporter once the issue is confirmed. Confirmed vulnerabilities are fixed in a patch release and disclosed through a GitHub Security Advisory that credits the reporter, unless the reporter prefers to remain anonymous.

## Scope

Combusta is a numerical library with a command-line front end. The library and the `combusta` command open no network connections, run no subprocesses, and do not deserialize pickled or otherwise executable data. The attack surface consists of the files Combusta reads and writes:

- case files (`case.toml`, parsed with the standard-library `tomllib` in `combusta.casefile`),
- ANSYS Forte CSV exports and generic crank-angle tables (`combusta.io.forte_parser`, `combusta.io.trace_csv`),
- Engine Combustion Network text tables (`combusta.io.ecn_loader`),
- reports, CSV files, and figures written by `combusta compare`, `matched`, `prechamber`, and `measurement` into the output directory given on the command line (existing files of the same name are overwritten).

Case files are treated as trusted input: a case file may point to any file the user can read, relative to its own location. Do not run Combusta on case files from untrusted sources without reading them first.

In scope are, for example, crashes, excessive memory or CPU use, or path handling problems triggered by a crafted input file, and vulnerabilities in how Combusta uses its dependencies.

Out of scope:

- vulnerabilities in NumPy, SciPy, pandas, or Matplotlib themselves; report those upstream,
- **numerically wrong results**. A wrong heat release, uncertainty, or correlation value is a correctness bug, not a security issue. Please open a regular [bug report](https://github.com/alpztel/combusta/issues/new/choose) and cite the reference you compare against.

Dependency vulnerabilities are tracked with Dependabot, and the code is scanned with GitHub CodeQL.

## Data hygiene

Never attach proprietary or confidential material to an issue, pull request, or advisory: engine geometry or meshes, OEM calibration data, test-cell logs, ANSYS Forte project or restart files, or licensed tutorial files. Reproduce problems with synthetic traces (the test suite contains generators for them) or with openly licensed data such as the [Engine Combustion Network](https://ecn.sandia.gov/) archives.
