"""Verification (GCI, Richardson extrapolation) and validation (ASME V&V 20-2009)."""

from __future__ import annotations

from combusta.vv.gci import (
    ASME_MIN_REFINEMENT_RATIO,
    DEFAULT_THREE_GRID_FS,
    ConvergenceRegime,
    GCIResult,
    GridSizeResult,
    compute_gci,
    compute_grid_size,
)
from combusta.vv.validation import (
    DEFAULT_GCI_COVERAGE_FACTOR,
    DEFAULT_VALIDATION_COVERAGE_FACTOR,
    ValidationResult,
    experimental_uncertainty,
    input_uncertainty,
    numerical_uncertainty_from_gci,
    sensitivity_coefficient,
    validate,
)

__all__ = [
    "DEFAULT_GCI_COVERAGE_FACTOR",
    "DEFAULT_VALIDATION_COVERAGE_FACTOR",
    "ValidationResult",
    "experimental_uncertainty",
    "input_uncertainty",
    "numerical_uncertainty_from_gci",
    "sensitivity_coefficient",
    "validate",
    "ASME_MIN_REFINEMENT_RATIO",
    "DEFAULT_THREE_GRID_FS",
    "ConvergenceRegime",
    "GCIResult",
    "GridSizeResult",
    "compute_gci",
    "compute_grid_size",
]
