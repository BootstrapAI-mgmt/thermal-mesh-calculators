"""
Every result the suite produces is checked against the result schema.

process_part(), process_part_from_props() and process_batch() are wrapped
here, before any test module imports them, so every result they return in
any test is passed to schema.validate_result(): a result that does not
conform fails the test that produced it, with the problems named.  The
wrappers count what they check; test_result_schema.py reads the count and
plants a result the hook must refuse.
"""

import functools

import pytest

import thermal_mesh_calculators
import thermal_mesh_calculators.batch as batch
from thermal_mesh_calculators.schema import validate_result

CHECKED = {"results": 0}


def schema_checked(func, returns_list=False):
    """func, with every result it returns checked by validate_result()."""
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        returned = func(*args, **kwargs)
        for result in (returned if returns_list else [returned]):
            problems = validate_result(result)
            assert not problems, (
                f"{func.__name__}() returned a result that does not conform "
                f"to the result schema: {problems}")
            CHECKED["results"] += 1
        return returned
    wrapper.schema_checked = True
    return wrapper


for _name, _returns_list in (("process_part", False),
                             ("process_part_from_props", False),
                             ("process_batch", True)):
    _checked = schema_checked(getattr(batch, _name), _returns_list)
    setattr(batch, _name, _checked)
    setattr(thermal_mesh_calculators, _name, _checked)


@pytest.fixture
def result_checks():
    """The wrapper factory and the running count of checked results."""
    return schema_checked, CHECKED
