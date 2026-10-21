"""
The sizing hook for an external mesher: size_for_fea() reproduces
BoundaryDrivenConductionCalculator.size_wall() on a plate, its result is
plain JSON and conforms to its schema, and its request is checked the way the
package checks a part.
"""

import json
import math
import re
from pathlib import Path

import pytest

from thermal_mesh_calculators import (
    FEA_SIZING_REQUEST_SCHEMA,
    FEA_SIZING_RESULT_SCHEMA,
    FEA_SIZING_SCHEMA_VERSION,
    MATERIALS,
    BoundaryDrivenConductionCalculator,
    PartInputError,
    size_for_fea,
    validate_fea_sizing,
    validate_fea_sizing_request,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
REL = 1e-12

# --- The plate grid -----------------------------------------------------------

GRID_MATERIALS = ("steel_mild", "aluminium_6061", "cast_iron",
                  "plastic_pa66_gf30", "insulation_ceramic_blanket",
                  "rubber_epdm")
GRID_THICKNESS_MM = (0.5, 2.0, 6.0, 25.0)
GRID_H = (5.0, 25.0, 150.0, 800.0)
# (t_surf_K, t_fluid_K, t_surr_K): a hot part losing heat to cooler air and a
# cooler enclosure; a part heated by hotter air and a hotter enclosure.
GRID_TEMPERATURES = ((473.15, 353.15, 323.15), (373.15, 573.15, 600.0))
GRID_EPSILON = (0.0, 0.9)
GRID_MAX_DT = (2.0, 10.0)

GRID = [
    (m, t, h, temps, eps, dt)
    for m in GRID_MATERIALS
    for t in GRID_THICKNESS_MM
    for h in GRID_H
    for temps in GRID_TEMPERATURES
    for eps in GRID_EPSILON
    for dt in GRID_MAX_DT
]

# size_wall()'s keys and the hook's names for them.  Every numeric quantity
# the two share is compared to 1e-12 relative; the rest exactly.
SHARED = {
    "n_cells": "n_cells",
    "biot": "biot",
    "regime": "regime",
    "dx_used_mm": "char_size_mm",
    "dx_raw_mm": "dx_raw_mm",
    "dx_exceeds_part": "dx_exceeds_part",
    "dt_wall_K": "dt_wall_K",
    "h_effective": "h_effective",
    "q_conv": "q_conv",
    "q_rad": "q_rad",
    "q_total": "q_total",
    "rad_fraction": "rad_fraction",
    "recommendation": "recommendation",
    "warnings": "warnings",
}


def _request(material, thickness_mm, h, temps, epsilon, max_dt, **extra):
    t_surf, t_fluid, t_surr = temps
    req = {"material": material, "thickness_mm": thickness_mm, "h": h,
           "t_surf_K": t_surf, "t_fluid_K": t_fluid, "t_surr_K": t_surr,
           "epsilon": epsilon, "max_dt": max_dt}
    req.update(extra)
    return req


def _size_wall(k, thickness_mm, h, temps, epsilon, max_dt, max_cells=100):
    t_surf, t_fluid, t_surr = temps
    return BoundaryDrivenConductionCalculator.size_wall(
        k=k, h=h, t_surf=t_surf, t_fluid=t_fluid, epsilon=epsilon,
        t_surr=t_surr, max_dt=max_dt, thickness_m=thickness_mm / 1000.0,
        max_cells=max_cells)


def _mismatches(direct, hook):
    """The shared quantities on which the hook differs from size_wall()."""
    bad = []
    for wall_key, hook_key in SHARED.items():
        a, b = direct[wall_key], hook[hook_key]
        if isinstance(a, float) and not math.isfinite(a):
            same = b is None          # JSON has no inf: the hook says None
        elif isinstance(a, float):
            same = (isinstance(b, float)
                    and math.isclose(a, b, rel_tol=REL, abs_tol=0.0))
        else:
            same = (a == b and type(a) is type(b))
        if not same:
            bad.append(f"{hook_key}: size_wall {a!r}, hook {b!r}")
    return bad


class TestReproducesSizeWallOnAPlate:
    """The gate: the hook is size_wall() on every point of the grid."""

    def test_the_grid_is_what_it_says(self):
        assert len(GRID) == 6 * 4 * 4 * 2 * 2 * 2 == 768

    def test_every_grid_point_by_material_name(self):
        failures = []
        for m, t, h, temps, eps, dt in GRID:
            direct = _size_wall(MATERIALS[m]["k"], t, h, temps, eps, dt)
            hook = size_for_fea(_request(m, t, h, temps, eps, dt))
            bad = _mismatches(direct, hook)
            if bad:
                failures.append(f"{m} {t} mm h={h} {temps} eps={eps} "
                                f"max_dt={dt}: {bad}")
        assert not failures, f"{len(failures)} of {len(GRID)}: {failures[:3]}"

    def test_every_grid_point_by_explicit_properties(self):
        failures = []
        for m, t, h, temps, eps, dt in GRID:
            props = {key: MATERIALS[m][key] for key in ("k", "rho", "cp")}
            req = _request(None, t, h, temps, eps, dt, **props)
            del req["material"]
            hook = size_for_fea(req)
            direct = _size_wall(props["k"], t, h, temps, eps, dt)
            bad = _mismatches(direct, hook)
            if bad:
                failures.append(f"{m} by properties, {t} mm h={h}: {bad}")
        assert not failures, f"{len(failures)} of {len(GRID)}: {failures[:3]}"

    def test_the_cell_size_times_the_count_is_the_wall(self):
        for m, t, h, temps, eps, dt in GRID[::17]:
            hook = size_for_fea(_request(m, t, h, temps, eps, dt))
            assert math.isclose(hook["char_size_mm"] * hook["n_cells"], t,
                                rel_tol=REL)

    def test_max_cells_reaches_size_wall(self):
        temps = (1073.15, 353.15, 353.15)
        direct = _size_wall(0.2, 25.0, 150.0, temps, 0.9, 0.5, max_cells=40)
        req = _request(None, 25.0, 150.0, temps, 0.9, 0.5, k=0.2,
                       rho=1000.0, cp=1000.0, max_cells=40)
        del req["material"]
        hook = size_for_fea(req)
        assert not _mismatches(direct, hook)
        assert hook["n_cells"] == 40
        assert hook["binding_constraint"] == "max_cells"


# --- The binding constraint and the rationale --------------------------------

class TestBindingConstraint:

    def test_conduction_binds_a_resolved_wall(self):
        # PA66-GF30 (k = 0.25 W/m K in MATERIALS), 6 mm, h = 150,
        # 473.15 -> 353.15 K, eps 0: q = 150 * 120 = 18000 W/m^2;
        # raw dx = 0.25 * 10 / 18000 m = 0.13889 mm; 6 / 0.13889 = 43.2,
        # so 44 cells.
        r = size_for_fea(_request("plastic_pa66_gf30", 6.0, 150.0,
                                  (473.15, 353.15, 353.15), 0.0, 10.0))
        assert r["n_cells"] == 44
        assert r["binding_constraint"] == "conduction"
        assert "44" in r["rationale"] and "max_dt" in r["rationale"]

    def test_the_wall_binds_a_thermally_thin_plate(self):
        r = size_for_fea(_request("aluminium_6061", 2.0, 25.0,
                                  (473.15, 353.15, 353.15), 0.0, 10.0))
        assert r["n_cells"] == 1 and r["dx_exceeds_part"] is True
        assert r["binding_constraint"] == "wall_thickness"
        assert r["regime"] == "thermally_thin"

    def test_near_equilibrium_binds_when_the_fluxes_cancel(self):
        # Ts between a hotter enclosure and cooler air: convection out and
        # radiation in cancel at h (Ts - Tf) = eps sigma (Tsurr^4 - Ts^4).
        sigma = 5.670374419e-8
        ts, tf, tsurr, eps = 400.0, 350.0, 500.0, 0.9
        h = eps * sigma * (tsurr ** 4 - ts ** 4) / (ts - tf)
        r = size_for_fea(_request("steel_mild", 3.0, h, (ts, tf, tsurr),
                                  eps, 10.0))
        assert r["regime"] == "near_equilibrium"
        assert r["binding_constraint"] == "near_equilibrium"
        assert r["dx_raw_mm"] is None and r["n_cells"] == 1

    def test_an_undefined_biot_is_none(self):
        # Ts == Tf: radiation alone drives the flux; Biot is undefined.
        r = size_for_fea(_request("steel_mild", 3.0, 25.0,
                                  (400.0, 400.0, 300.0), 0.9, 10.0))
        assert r["biot"] is None and r["h_effective"] is None
        assert validate_fea_sizing(r) == []


# --- The result is plain JSON and conforms to its schema ----------------------

class TestResultSchema:

    @pytest.fixture
    def result(self):
        return size_for_fea(_request("steel_mild", 3.0, 25.0,
                                     (473.15, 353.15, 323.15), 0.9, 10.0))

    def test_every_grid_result_conforms_and_is_json(self):
        for m, t, h, temps, eps, dt in GRID[::5]:
            r = size_for_fea(_request(m, t, h, temps, eps, dt))
            assert validate_fea_sizing(r) == [], (m, t, h)
            assert json.loads(json.dumps(r, allow_nan=False)) == r

    def test_the_result_carries_its_version_and_material(self, result):
        assert result["schema_version"] == FEA_SIZING_SCHEMA_VERSION == 1
        assert result["material"]["name"] == "steel_mild"
        k, rho, cp = (MATERIALS["steel_mild"][key] for key in ("k", "rho", "cp"))
        assert result["material"]["k"] == k
        assert result["material"]["alpha_m2_s"] == pytest.approx(
            k / (rho * cp), rel=1e-15)

    def test_a_result_missing_its_schema_version_fails(self, result):
        del result["schema_version"]
        assert [p["code"] for p in validate_fea_sizing(result)] == \
            ["MISSING_KEY"]

    def test_another_schema_version_fails(self, result):
        result["schema_version"] = 2
        assert [p["code"] for p in validate_fea_sizing(result)] == \
            ["UNSUPPORTED_VERSION"]

    def test_an_undeclared_key_fails(self, result):
        result["mesh"] = "fine"
        assert [p["code"] for p in validate_fea_sizing(result)] == \
            ["UNKNOWN_KEY"]

    def test_an_infinite_value_fails(self, result):
        result["dx_raw_mm"] = float("inf")
        assert [p["code"] for p in validate_fea_sizing(result)] == \
            ["OUT_OF_RANGE"]

    def test_a_cell_size_that_does_not_fill_the_wall_fails(self, result):
        result["char_size_mm"] *= 1.5
        assert [p["code"] for p in validate_fea_sizing(result)] == \
            ["INCONSISTENT"]

    def test_an_unknown_binding_constraint_fails(self, result):
        result["binding_constraint"] = "guess"
        assert [p["code"] for p in validate_fea_sizing(result)] == \
            ["UNKNOWN_NAME"]

    def test_a_warning_that_is_not_a_string_fails(self, result):
        result["warnings"] = [{"code": "X"}]
        assert [p["code"] for p in validate_fea_sizing(result)] == \
            ["WRONG_TYPE"]

    def test_not_a_dict(self):
        assert [p["code"] for p in validate_fea_sizing([])] == \
            ["NOT_AN_OBJECT"]


# --- The request is checked before anything is computed -----------------------

GOOD = _request("steel_mild", 3.0, 25.0, (473.15, 353.15, 323.15), 0.9, 10.0)


class TestRequest:

    def test_a_good_request_has_no_problems(self):
        assert validate_fea_sizing_request(GOOD) == []

    def test_an_alias_names_the_material(self):
        r = size_for_fea(dict(GOOD, material="Mild steel"))
        assert r["material"]["name"] == "steel_mild"

    def test_caller_keys_are_carried_unread(self):
        assert validate_fea_sizing_request(dict(GOOD, x_mesher_id=7)) == []

    @pytest.mark.parametrize("change, code, words", [
        ({"thickness": 3.0}, "UNKNOWN_KEY", "thickness"),
        ({"thickness_mm": -3.0}, "OUT_OF_RANGE", "thickness_mm"),
        ({"epsilon": 1.2}, "OUT_OF_RANGE", "epsilon"),
        ({"h": "25"}, "WRONG_TYPE", "h"),
        ({"max_cells": 0}, "OUT_OF_RANGE", "max_cells"),
        ({"material": "unobtainium"}, "UNKNOWN_NAME", "steel_mild"),
        ({"schema_version": 2}, "UNSUPPORTED_VERSION", "schema_version"),
        ({"k": 45.0}, "INCONSISTENT", "both"),
    ])
    def test_a_bad_request_raises_with_the_problem(self, change, code, words):
        req = dict(GOOD, **change)
        problems = validate_fea_sizing_request(req)
        assert [p["code"] for p in problems] == [code]
        assert words in problems[0]["reason"]
        with pytest.raises(PartInputError) as raised:
            size_for_fea(req)
        assert raised.value.problems == problems

    def test_a_missing_key_is_named(self):
        req = dict(GOOD)
        del req["t_surr_K"]
        assert [(p["key"], p["code"]) for p in
                validate_fea_sizing_request(req)] == [("t_surr_K",
                                                       "MISSING_KEY")]

    def test_no_material_and_no_properties_is_named(self):
        req = dict(GOOD)
        del req["material"]
        problems = validate_fea_sizing_request(req)
        assert [p["code"] for p in problems] == ["MISSING_KEY"]
        assert "'material'" in problems[0]["reason"]

    def test_incomplete_properties_are_named(self):
        req = dict(GOOD, k=45.0, rho=7850.0)
        del req["material"]
        problems = validate_fea_sizing_request(req)
        assert [(p["key"], p["code"]) for p in problems] == [("cp",
                                                             "MISSING_KEY")]

    def test_every_problem_is_reported_at_once(self):
        req = dict(GOOD, thickness_mm=0.0, epsilon=-1.0, colour="red")
        assert len(validate_fea_sizing_request(req)) == 3

    def test_not_a_dict(self):
        with pytest.raises(PartInputError, match="dict"):
            size_for_fea([("material", "steel_mild")])


# --- The integrator guide states what the schemas say -------------------------

def _guide_table(name):
    path = REPO_ROOT / "docs" / "integration.md"
    if not path.is_file():
        if (REPO_ROOT / ".git").exists():
            pytest.fail("docs/integration.md is missing from this checkout")
        pytest.skip("docs/ is not in this tree (the sdist ships none)")
    text = path.read_text(encoding="utf-8")
    match = re.search(rf"<!-- {name}:begin -->(.*?)<!-- {name}:end -->",
                      text, re.S)
    assert match, f"no {name} table"
    return [line.split("|")[1].strip().strip("`")
            for line in match.group(1).splitlines() if line.startswith("| `")]


@pytest.mark.parametrize("table, schema", [
    ("fea-request-keys", FEA_SIZING_REQUEST_SCHEMA),
    ("fea-result-keys", FEA_SIZING_RESULT_SCHEMA),
])
def test_the_integrator_guide_lists_every_hook_key(table, schema):
    keys = _guide_table(table)
    assert sorted(keys) == sorted(schema)
    assert len(keys) == len(set(keys))
