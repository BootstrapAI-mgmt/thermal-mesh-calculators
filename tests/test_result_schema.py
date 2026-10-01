"""
Every result validates against the versioned result schema
==========================================================

conftest.py passes every result the suite produces to validate_result();
these tests show that check can fail: a result missing schema_version, an
undeclared key, an unregistered warning code, a governing size that is not
the smallest candidate are each reported, and the hook fails a test that
produces one.  The registry of warning codes is held to the codes the
source emits, the documentation's key and warning tables to the schemas,
and the integrator guide's examples run as doctests.  The package imports
nothing outside the standard library.
"""

import ast
import doctest
import re
from pathlib import Path

import pytest

import thermal_mesh_calculators
from thermal_mesh_calculators.batch import process_batch, process_part
from thermal_mesh_calculators.schema import (
    PART_SCHEMA,
    PROBLEM_CODES,
    PROJECT_SCHEMA,
    RESULT_SCHEMA_VERSION,
    WARNING_CODES,
    validate_result,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
PACKAGE = REPO_ROOT / "thermal_mesh_calculators"

PROJECT = {"t_fluid_K": 353.15, "t_surr_K": 353.15, "t_exh_K": 1073.15}
PART = {"part_id": "RS-1", "material": "steel_mild",
        "component_class": "structural", "convection_zone": "engine_beside",
        "thickness_mm": 3.0, "t_surf_K": 473.15, "surface": "painted"}


def _codes(problems):
    return [(p["key"], p["code"]) for p in problems]


@pytest.fixture
def result():
    return process_part(PART, PROJECT)


# --- The validator, on planted results ---

class TestValidateResult:

    def test_a_result_conforms(self, result):
        assert result["schema_version"] == RESULT_SCHEMA_VERSION == 1
        assert validate_result(result) == []

    def test_a_result_missing_its_schema_version_fails(self, result):
        del result["schema_version"]
        assert _codes(validate_result(result)) == [
            ("schema_version", "MISSING_KEY")]

    def test_another_schema_version_fails(self, result):
        result["schema_version"] = 2
        assert _codes(validate_result(result)) == [
            ("schema_version", "UNSUPPORTED_VERSION")]

    def test_an_undeclared_key_fails(self, result):
        result["dx_mm"] = 1.0
        assert _codes(validate_result(result)) == [("dx_mm", "UNKNOWN_KEY")]

    def test_a_missing_required_key_fails(self, result):
        del result["warnings"]
        assert _codes(validate_result(result)) == [("warnings", "MISSING_KEY")]

    def test_a_wrong_type_fails(self, result):
        result["part_id"] = 3
        assert _codes(validate_result(result)) == [("part_id", "WRONG_TYPE")]

    def test_a_promised_sub_key_fails(self, result):
        del result["biot"]["mesh_type"]
        assert _codes(validate_result(result)) == [
            ("biot.mesh_type", "MISSING_KEY")]

    def test_an_unregistered_warning_code_fails(self, result):
        result["warnings"].append({"code": "MYSTERY", "severity": "caution",
                                   "message": "?"})
        problems = validate_result(result)
        assert len(problems) == 1 and problems[0]["code"] == "UNKNOWN_NAME"
        assert "unregistered warning code 'MYSTERY'" in problems[0]["reason"]

    def test_a_warning_with_the_wrong_severity_fails(self, result):
        result["warnings"].append({"code": "NO_FINITE_SIZE",
                                   "severity": "info", "message": "?"})
        problems = validate_result(result)
        assert len(problems) == 1 and "severity 'info'" in problems[0]["reason"]

    def test_a_warning_with_an_undeclared_field_fails(self, result):
        result["warnings"].append({"code": "NO_FINITE_SIZE", "severity":
                                   "warning", "message": "?", "hint": "?"})
        assert [p["code"] for p in validate_result(result)] == ["UNKNOWN_KEY"]

    def test_a_governing_size_that_is_not_the_smallest_fails(self, result):
        result["governing_dx_mm"] = result["all_constraints"][-1]["dx_mm"] + 1
        assert _codes(validate_result(result)) == [
            ("governing_dx_mm", "INCONSISTENT")]

    @pytest.mark.parametrize("dx", [0.0, -1.0, float("nan")])
    def test_a_governing_size_must_be_positive(self, result, dx):
        result["governing_dx_mm"] = dx
        result["all_constraints"] = [{"dx_mm": dx, "source": "conduction"}]
        result["governing_constraint"] = "conduction"
        assert ("governing_dx_mm", "OUT_OF_RANGE") in _codes(validate_result(result))

    def test_an_unbounded_part_may_be_infinite(self):
        shield = {"part_id": "SH", "material": "steel_mild",
                  "component_class": "shield",
                  "convection_zone": "engine_beside", "thickness_mm": 1.0,
                  "t_exh_K": 353.15}
        r = process_part(shield, PROJECT)
        assert r["governing_dx_mm"] == float("inf")
        assert validate_result(r) == []

    def test_not_a_dict(self):
        assert _codes(validate_result([])) == [(None, "NOT_AN_OBJECT")]


class TestErrorResults:

    def test_an_error_result_conforms(self):
        [r] = process_batch([dict(PART, thickness_mm=-1.0)], PROJECT)
        assert r["error"] and r["governing_dx_mm"] is None
        assert r["schema_version"] == RESULT_SCHEMA_VERSION
        assert _codes(r["problems"]) == [("thickness_mm", "OUT_OF_RANGE")]
        assert validate_result(r) == []

    def test_a_calculator_failure_is_reported(self, monkeypatch):
        import thermal_mesh_calculators.batch as batch_mod

        def broken(*args, **kwargs):
            raise ArithmeticError("planted")
        monkeypatch.setattr(batch_mod.BoundaryDrivenConductionCalculator,
                            "max_mesh_size", broken)
        [r] = process_batch([PART], PROJECT)
        assert r["error"] == "planted"
        assert r["problems"] == [{"key": None, "code": "CALCULATION_FAILED",
                                  "reason": "ArithmeticError: planted"}]

    def test_an_error_result_with_an_unregistered_code_fails(self):
        [r] = process_batch([dict(PART, thickness_mm=-1.0)], PROJECT)
        r["problems"][0]["code"] = "MYSTERY"
        assert _codes(validate_result(r)) == [("problems[0].code", "UNKNOWN_NAME")]

    def test_an_error_result_missing_its_problems_fails(self):
        [r] = process_batch([dict(PART, thickness_mm=-1.0)], PROJECT)
        del r["problems"]
        assert _codes(validate_result(r)) == [("problems", "MISSING_KEY")]


# --- The hook in conftest.py ---

def test_the_hook_wraps_the_three_entry_points():
    import thermal_mesh_calculators.batch as batch_mod
    for name in ("process_part", "process_part_from_props", "process_batch"):
        assert getattr(getattr(batch_mod, name), "schema_checked", False), name
        assert getattr(thermal_mesh_calculators, name) is getattr(batch_mod, name)


def test_the_hook_counts_and_fails_a_planted_result(result_checks):
    schema_checked, checked = result_checks
    before = checked["results"]
    process_part(PART, PROJECT)
    assert checked["results"] == before + 1

    def missing_version(part, project):
        r = process_part(part, project)
        del r["schema_version"]
        return r

    with pytest.raises(AssertionError, match="'schema_version'.*MISSING_KEY"):
        schema_checked(missing_version)(PART, PROJECT)


# --- The registry of warning codes ---

def test_every_registered_code_is_emitted_and_every_emitted_code_registered():
    source = (PACKAGE / "batch.py").read_text(encoding="utf-8")
    emitted = set(re.findall(r'"code": "([A-Z_]+)"', source))
    assert emitted - set(WARNING_CODES) - set(PROBLEM_CODES) == set()
    assert set(WARNING_CODES) - emitted == set()


# --- The package needs only the standard library ---

STDLIB_IMPORTS = {"csv", "io", "json", "math", "os", "re", "typing"}


def test_the_package_imports_only_the_standard_library():
    imported = set()
    for module in sorted(PACKAGE.glob("*.py")):
        for node in ast.walk(ast.parse(module.read_text(encoding="utf-8"))):
            if isinstance(node, ast.Import):
                imported |= {a.name.split(".")[0] for a in node.names}
            elif isinstance(node, ast.ImportFrom) and node.level == 0:
                imported.add(node.module.split(".")[0])
    imported.discard("thermal_mesh_calculators")
    assert imported <= STDLIB_IMPORTS, imported - STDLIB_IMPORTS


# --- The documentation states what the schemas say ---

def _doc(relpath):
    path = REPO_ROOT / relpath
    if not path.is_file():
        if (REPO_ROOT / ".git").exists():
            pytest.fail(f"{relpath} is missing from this checkout")
        pytest.skip(f"{relpath} is not in this tree (the sdist ships no docs/)")
    return path.read_text(encoding="utf-8")


def _table(text, name):
    """The rows between <!-- name:begin --> and <!-- name:end -->, each as
    its list of cells."""
    match = re.search(rf"<!-- {name}:begin -->(.*?)<!-- {name}:end -->",
                      text, re.S)
    assert match, f"no {name} table"
    rows = []
    for line in match.group(1).splitlines():
        if line.startswith("| `"):
            rows.append([c.strip() for c in line.strip().strip("|").split("|")])
    return rows


@pytest.mark.parametrize("table, schema", [
    ("part-keys", PART_SCHEMA), ("project-keys", PROJECT_SCHEMA),
])
def test_the_integrator_guide_lists_every_key(table, schema):
    rows = _table(_doc("docs/integration.md"), table)
    keys = [row[0].strip("`") for row in rows]
    assert sorted(keys) == sorted(schema)
    assert len(keys) == len(set(keys))


def test_the_quick_reference_lists_every_warning_code():
    rows = _table(_doc("docs/quick_reference.md"), "warning-codes")
    listed = {row[0].strip("`"): set(row[1].split(" or ")) for row in rows}
    assert listed == {code: set(entry["severities"])
                      for code, entry in WARNING_CODES.items()}


def test_the_integrator_guide_examples_run():
    _doc("docs/integration.md")
    outcome = doctest.testfile(
        str(REPO_ROOT / "docs" / "integration.md"), module_relative=False,
        optionflags=doctest.ELLIPSIS | doctest.NORMALIZE_WHITESPACE,
        verbose=False, report=True)
    assert outcome.attempted > 0
    assert outcome.failed == 0
