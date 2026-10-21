"""
The results that were there before the internal-flow correlations, held where
they were.

golden_outputs.py was captured from the package before estimate_h() gained its
duct argument: the eleven worked automotive examples and the batch example,
byte for byte as they print, and estimate_h() over a grid of external-flow
inputs with its floats written to 12 significant digits.  These tests passed on
the package as it was and pass on it now.
"""

import contextlib
import importlib
import io
from pathlib import Path

import pytest

from tests import golden_outputs
from thermal_mesh_calculators.h_estimator import estimate_h

REPO_ROOT = Path(__file__).resolve().parent.parent


def _example_stdout(name):
    """What examples/<name>.py prints.  The sdist ships no examples/, so
    there the test skips; in a checkout a missing example fails it."""
    try:
        module = importlib.import_module(f"examples.{name}")
    except ImportError:
        if (REPO_ROOT / ".git").exists():
            raise
        pytest.skip("examples/ is not in this tree (the sdist ships none)")
    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer):
        module.main()
    return buffer.getvalue()


# The grid the capture ran (golden_outputs.ESTIMATE_H_GRID_TEXT), unchanged.
ESTIMATE_H_GRID = [
    (velocity, t_surf, t_fluid, char_length, orientation)
    for velocity in (0.0, 0.5, 4.47, 15.0)
    for (t_surf, t_fluid) in ((400.0, 300.0), (900.0, 330.0))
    for char_length in (0.05, 0.4)
    for orientation in ("vertical", "horizontal_up", "horizontal_down")
]


def _canonical(value):
    """A value as text: floats to 12 significant digits, so a last-place
    difference between platforms' maths libraries cannot show."""
    if isinstance(value, bool) or value is None or isinstance(value, str):
        return repr(value)
    if isinstance(value, float):
        return format(value, ".12g")
    if isinstance(value, int):
        return repr(value)
    if isinstance(value, dict):
        return "{" + ", ".join(f"{k}: {_canonical(value[k])}"
                               for k in sorted(value)) + "}"
    raise TypeError(type(value).__name__)


def _estimate_h_grid_text():
    lines = []
    for velocity, t_surf, t_fluid, char_length, orientation in ESTIMATE_H_GRID:
        result = estimate_h(velocity, t_surf, t_fluid, char_length, orientation)
        lines.append(f"{velocity} {t_surf} {t_fluid} {char_length} {orientation}: "
                     + _canonical(result))
    return "\n".join(lines) + "\n"


def test_the_eleven_worked_examples_print_what_they_printed():
    assert _example_stdout("automotive_examples") == \
        golden_outputs.AUTOMOTIVE_EXAMPLES_STDOUT


def test_the_batch_example_prints_what_it_printed():
    # process_batch() reaches estimate_h() through the convection zones.
    assert _example_stdout("batch_example") == \
        golden_outputs.BATCH_EXAMPLE_STDOUT


def test_estimate_h_without_a_duct_returns_what_it_returned():
    assert len(ESTIMATE_H_GRID) == 48
    assert _estimate_h_grid_text() == golden_outputs.ESTIMATE_H_GRID_TEXT
