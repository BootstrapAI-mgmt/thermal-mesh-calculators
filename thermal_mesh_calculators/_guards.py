"""
Input range guards shared by the calculators.

A physically impossible input (a negative conductivity, an emissivity
above one, a temperature at or below absolute zero) fails where it enters,
with a ValueError naming the parameter, the range it must lie in and the
value given, instead of coming back as a negative or infinite element size.
A value that is not a number at all raises TypeError.
"""

import math


def _number(name: str, value) -> float:
    # A string such as "45" converts, but the calculators' arithmetic
    # would then fail far from here; a bool is a number by accident.
    if isinstance(value, (bool, str, bytes)):
        raise TypeError(f"{name} must be a number, got {value!r}")
    try:
        return float(value)
    except (TypeError, ValueError):
        raise TypeError(f"{name} must be a number, got {value!r}") from None


def require_positive(name: str, value, unit: str = "") -> None:
    """Raise ValueError unless value is a finite number > 0."""
    number = _number(name, value)
    if not (math.isfinite(number) and number > 0.0):
        suffix = " " + unit if unit else ""
        raise ValueError(f"{name} must be > 0{suffix}, got {value!r}")


def require_non_negative(name: str, value, unit: str = "") -> None:
    """Raise ValueError unless value is a finite number >= 0."""
    number = _number(name, value)
    if not (math.isfinite(number) and number >= 0.0):
        suffix = " " + unit if unit else ""
        raise ValueError(f"{name} must be >= 0{suffix}, got {value!r}")


def require_fraction(name: str, value) -> None:
    """Raise ValueError unless 0 <= value <= 1 (emissivity, view factor)."""
    number = _number(name, value)
    if not (0.0 <= number <= 1.0):
        raise ValueError(f"{name} must lie in [0, 1], got {value!r}")


def require_temperatures(**temperatures) -> None:
    """Every temperature (K) must be > 0."""
    for name, value in temperatures.items():
        require_positive(name, value, "K")
