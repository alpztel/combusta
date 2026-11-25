# Contributing to Combusta

Contributions are welcome: bug reports, verification cases, corrections to the physics or the references, documentation, and code. Combusta produces numbers that may end up in engineering reports and theses, so every change has to be traceable to a published source or a closed-form result.

Please read the [Code of Conduct](CODE_OF_CONDUCT.md) before taking part. Security problems are **not** reported through issues; see [SECURITY.md](SECURITY.md).

## Ways to contribute

- **Bug reports.** Use the [bug report form](https://github.com/alpztel/combusta/issues/new/choose). Include the version, a minimal reproduction on synthetic data, and the expected result.
- **Physics and reference corrections.** If a formula, constant, or citation is wrong, open a *Physics / reference issue* and cite the source (equation, table, or page) that shows the correct value.
- **Feature requests.** Describe the analysis you need and the literature it is based on.
- **Pull requests.** For anything larger than a small fix, open an issue first so the approach can be agreed before you write code.

## Ground rules

### No proprietary data

Never commit, attach, or paste confidential material: engine geometry or meshes, OEM calibration maps, test-cell logs, ANSYS Forte project or restart files, or licensed tutorial files. Tests, examples, and benchmarks use only

1. closed-form or manufactured solutions (polytropic and Otto cycles, forward-integrated Wiebe cycles, power-law grid triplets), or
2. openly licensed data, such as the Sandia [Engine Combustion Network](https://ecn.sandia.gov/) archives, with their source and terms stated.

### SI units inside, conversions at the boundary

Every function in the numerical packages (`thermo`, `signal`, `spray`, `vv`, `ignition`, `exergy`, `emissions`, `compare`, `case`) takes and returns SI quantities. Units appear in parameter names (`pressure_pa`, `crank_angle_deg`, `dv_dtheta_m3_per_deg`) and in brackets in docstrings (`[Pa]`, `[J/deg]`). Unit conversion happens only in `combusta.io` (ingestion) and `combusta.plotting` (display), using the tables in `combusta.constants`.

### Every model has a source and a test

A new correlation, model, or procedure must

1. cite its primary source (SAE paper number, journal article, standard and clause, or textbook equation) in the module docstring and in the README reference list,
2. reject non-physical input with `PhysicalBoundsError` from `combusta.exceptions` (for example $P \le 0$, $T \le 0$, a compression ratio $\le 1$), and state the validity range of the correlation, and
3. come with a verification test against a closed-form solution, a published table, or a theoretical order of accuracy. Name the reference in the test docstring and add a row to the README verification table.

Heuristics that are not taken from the literature (for example a plausibility threshold) must be labeled as such in the docstring and the README.

### Typing and style

- Python 3.11 or newer; complete type annotations on public functions, including return types.
- Results are returned as `@dataclass(frozen=True, slots=True)` containers.
- NumPy-style docstrings.
- `pyrefly check` must report no errors.

## Development setup

```bash
git clone https://github.com/alpztel/combusta.git
cd combusta
python -m venv .venv
# Windows PowerShell: .\.venv\Scripts\Activate.ps1
# Linux / macOS:      source .venv/bin/activate
python -m pip install --upgrade pip
pip install -e ".[dev]"
```

## Before you open a pull request

Run these from the repository root; CI runs the same steps on Python 3.11, 3.12, and 3.13:

```bash
pyrefly check                        # static type check, 0 errors
pytest                               # full test suite
python examples/altitude_sweep_demo.py   # end-to-end smoke test
```

Then:

- add an entry under `[Unreleased]` in [CHANGELOG.md](CHANGELOG.md) (Keep a Changelog sections: Added, Changed, Fixed, Removed, Documentation),
- update the README if you changed a public API, an equation, or a reference,
- write commit messages in the Conventional Commits style (`feat(thermo): ...`, `fix(vv): ...`, `docs: ...`).

## Review and merging

`main` is protected. Pull requests are merged after the three CI jobs pass and the maintainer has reviewed the change. Review concentrates on physical correctness, traceability to the literature, and test coverage, in that order.

## Releases

Combusta follows [Semantic Versioning](https://semver.org/). While the version is 0.x, minor releases may change the public API; such changes are listed in the changelog. A release consists of a version bump in `pyproject.toml`, `src/combusta/__init__.py`, and `CITATION.cff`, a dated changelog section, an annotated `vX.Y.Z` tag, and a GitHub release whose wheel and source distribution are built by the release workflow.

## License

By contributing, you agree that your contributions are licensed under the [MIT License](LICENSE) of this project.
