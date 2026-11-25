"""Domain-specific exceptions for Combusta."""

from __future__ import annotations


class CombustaError(Exception):
    """Base exception for all ``combusta`` errors."""


# Backward-compatible alias (pre-0.3 name)
ValidationAssistantError = CombustaError


class ForteParseError(CombustaError):
    """Raised when an ANSYS Forte CSV file is malformed or lacks required columns."""


class UnitDetectionError(ForteParseError):
    """Raised when a physical header unit is missing, unknown, or dimensionally invalid."""


class PhysicalBoundsError(CombustaError):
    """Raised when SI-converted arrays violate fundamental thermodynamic/kinematic bounds."""


class TDCPhaseError(CombustaError):
    """Raised when motored peak pressure violates thermodynamic TDC phasing tolerance."""


class CaseFileError(CombustaError, ValueError):
    """Raised when a ``case.toml`` file is missing, malformed, or inconsistent."""
