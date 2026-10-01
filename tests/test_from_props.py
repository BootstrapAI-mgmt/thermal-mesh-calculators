"""
Sizing from explicit material properties, fluid regions, service limits
========================================================================

process_part_from_props(part, project, **MATERIALS[m]) must return what
process_part(part, project) returns for a part whose material is m, for
every material in the table and every sizing path: the properties path and
the table path are one computation.  The parts below carry every input
their path reads, and half the projects carry a time step, so the density
and the specific heat reach the transient constraint too.

A fluid region is sized by its boundary layer alone, checked against the
boundary-layer calculator called directly.  A material record that gives
a service temperature limit turns a temperature above it into a warning.
"""

import inspect
import math

import pytest

import thermal_mesh_calculators.batch as batch_mod
from thermal_mesh_calculators.batch import (
    MATERIALS,
    PartInputError,
    process_batch,
    process_part,
    process_part_from_props,
    summary_table,
)
from thermal_mesh_calculators.boundary_layer import BoundaryLayerCalculator
from thermal_mesh_calculators.schema import (
    MATERIAL_RECORD_SCHEMA,
    validate_material_record,
)
from thermal_mesh_calculators.zones import get_zone


PROJECTS = (
    {"t_fluid_K": 353.15, "t_surr_K": 353.15, "t_exh_K": 1073.15},
    {"t_fluid_K": 353.15, "t_surr_K": 353.15, "t_exh_K": 1073.15,
     "dt": 0.5, "fo_max": 0.5, "tau_bc": 20.0},
)

# One part per sizing path.  The material is filled in per case.
TEMPLATES = {
    "structural": {
        "part_id": "PROPS-STR", "component_class": "structural",
        "convection_zone": "engine_beside", "thickness_mm": 3.0,
        "t_surf_K": 473.15, "surface": "painted", "radius_mm": 25.0,
        "bl_regime": "external_forced", "velocity_ms": 4.0,
    },
    # No surface: the fallback emissivity is chosen by the material's name.
    "exhaust_defaulted_surface": {
        "part_id": "PROPS-EXH", "component_class": "exhaust",
        "convection_zone": "exhaust_internal", "thickness_mm": 6.0,
        "t_surf_K": 1173.15, "t_fluid_K": 1073.15,
    },
    "exhaust_adjacent": {
        "part_id": "PROPS-ADJ", "component_class": "exhaust_adjacent",
        "convection_zone": "exhaust_above", "thickness_mm": 2.0,
        "t_surf_K": 520.0, "epsilon": 0.6,
    },
    "shield": {
        "part_id": "PROPS-SH", "component_class": "shield",
        "convection_zone_in": "exhaust_beside",
        "convection_zone_out": "engine_beside", "thickness_mm": 0.8,
        "surface_in": "aluminised", "surface_out": "aluminised",
    },
    "multilayer_shield": {
        "part_id": "PROPS-ML", "component_class": "multilayer_shield",
        "convection_zone_in": "exhaust_beside",
        "convection_zone_out": "shield_gap_confined", "thickness_mm": 0.6,
        "surface_in": "aluminised", "surface_out": "painted",
        "surface_g1": "aluminised", "surface_g2": "aluminised",
        "h_gap": 12.0, "f12": 0.8,
    },
}


def _same(a, b):
    """Exact equality, recursively, with NaN equal to NaN."""
    if isinstance(a, float) and isinstance(b, float):
        return a == b or (math.isnan(a) and math.isnan(b))
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(_same(a[k], b[k]) for k in a)
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(_same(x, y) for x, y in zip(a, b))
    return type(a) is type(b) and a == b


# --- The gate: the properties path equals the table path ---

@pytest.mark.parametrize("name", sorted(MATERIALS))
def test_properties_path_equals_table_path(name):
    diverged = []
    for path, template in TEMPLATES.items():
        for i, project in enumerate(PROJECTS):
            part = dict(template, material=name)
            by_name = process_part(part, project)
            by_props = process_part_from_props(part, project,
                                               **MATERIALS[name])
            if not _same(by_props, by_name):
                diverged.append(f"{path} / project {i}")
    assert not diverged, f"{name}: the paths differ for {diverged}"


def test_all_43_materials_are_covered():
    assert len(MATERIALS) == 43


def test_a_part_without_a_label_echoes_none_and_sizes_alike():
    # Explicit surfaces and a metal: no rule reads the name, so only the
    # echoed material differs.
    part = dict(TEMPLATES["structural"], material="steel_mild")
    unlabelled = dict(part)
    del unlabelled["material"]
    project = PROJECTS[1]
    by_props = process_part_from_props(unlabelled, project,
                                       **MATERIALS["steel_mild"])
    assert by_props["material"] is None
    assert _same(dict(by_props, material="steel_mild"),
                 process_part(part, project))


def test_a_label_need_not_name_a_table_material():
    part = dict(TEMPLATES["structural"], material="our-own-grade-7")
    r = process_part_from_props(part, PROJECTS[0], k=40.0, rho=7800.0,
                                cp=470.0)
    assert r["material"] == "our-own-grade-7"
    assert 0 < r["governing_dx_mm"] < float("inf")


def test_the_properties_are_used():
    part = dict(TEMPLATES["structural"], material="steel_mild")
    steel = process_part_from_props(part, PROJECTS[0],
                                    **MATERIALS["steel_mild"])
    half_k = dict(MATERIALS["steel_mild"], k=MATERIALS["steel_mild"]["k"] / 2)
    softer = process_part_from_props(part, PROJECTS[0], **half_k)
    assert softer["conduction"]["max_dx_mm"] == pytest.approx(
        steel["conduction"]["max_dx_mm"] / 2, rel=1e-12)


def test_an_out_of_range_property_is_refused_by_name():
    part = dict(TEMPLATES["structural"], material="x")
    with pytest.raises(PartInputError) as excinfo:
        process_part_from_props(part, PROJECTS[0], k=-45.0, rho=7800.0,
                                cp=0.0)
    problems = [(p["key"], p["code"]) for p in excinfo.value.problems]
    assert problems == [("properties.k", "OUT_OF_RANGE"),
                        ("properties.cp", "OUT_OF_RANGE")]
    assert "k must be > 0 W/m K, got -45.0 (material record)" in str(
        excinfo.value)


def test_record_fields_are_the_keyword_parameters():
    """**MATERIALS[m] works for any record the schema admits, including
    the optional fields the maintainer may fill in later."""
    params = inspect.signature(process_part_from_props).parameters
    keyword_only = {n for n, p in params.items()
                    if p.kind is inspect.Parameter.KEYWORD_ONLY}
    assert keyword_only == set(MATERIAL_RECORD_SCHEMA)


# --- Fluid regions ---

PROJECT = PROJECTS[0]


def _fluid(**overrides):
    part = {"part_id": "FL-001", "component_class": "fluid",
            "convection_zone": "front_end_edges"}
    part.update(overrides)
    return part


class TestFluidRegions:

    def test_a_forced_zone_is_sized_by_its_external_boundary_layer(self):
        r = process_part(_fluid(), PROJECT)
        zone = get_zone("front_end_edges")
        expected = BoundaryLayerCalculator.estimate_mesh(
            U=zone["velocity_ms"], x=0.1, t_fluid=353.15,
            regime="external_forced")
        assert r["boundary_layer"] == expected
        assert r["governing_dx_mm"] == expected["max_dx_surface_mm"]
        assert r["governing_constraint"] == "aero_boundary_layer"
        assert r["all_constraints"] == [
            {"dx_mm": expected["max_dx_surface_mm"],
             "source": "aero_boundary_layer"}]

    def test_no_solid_constraint_runs(self):
        r = process_part(_fluid(), PROJECTS[1])
        for key in ("conduction", "biot", "radiation", "shield",
                    "transient", "solver_advisory", "h_used", "material"):
            assert r[key] is None, key
        assert "lateral" not in r

    def test_a_natural_zone_takes_the_mixed_regime_and_the_wall(self):
        r = process_part(_fluid(convection_zone="exhaust_above",
                                t_surf_K=500.0), PROJECT)
        expected = BoundaryLayerCalculator.estimate_mesh(
            U=0.0, x=0.1, t_fluid=353.15, regime="mixed_unknown",
            t_surf=500.0)
        assert r["boundary_layer"] == expected
        assert r["boundary_layer"]["regime"] == "mixed_unknown"

    def test_velocity_and_regime_stand_in_for_a_zone(self):
        part = {"part_id": "FL-2", "component_class": "fluid",
                "velocity_ms": 12.0, "bl_regime": "external_forced",
                "char_length_mm": 400.0}
        r = process_part(part, PROJECT)
        expected = BoundaryLayerCalculator.estimate_mesh(
            U=12.0, x=0.4, t_fluid=353.15, regime="external_forced")
        assert r["governing_dx_mm"] == expected["max_dx_surface_mm"]

    def test_a_region_without_flow_names_the_alternatives(self):
        with pytest.raises(PartInputError) as excinfo:
            process_part({"part_id": "FL-3", "component_class": "fluid"},
                         PROJECT)
        message = str(excinfo.value)
        assert message.startswith("Missing convection zone: part 'FL-3' "
                                  "has no 'convection_zone' key.")
        assert ("Give 'convection_zone' or both 'velocity_ms' and "
                "'bl_regime'.") in message

    def test_solid_keys_on_a_fluid_region_are_not_read(self):
        bare = process_part(_fluid(), PROJECT)
        with_solid_keys = process_part(
            _fluid(material="air", thickness_mm=-1.0, epsilon=7.0), PROJECT)
        assert _same(with_solid_keys, bare)

    def test_the_properties_path_takes_a_fluid_region_too(self):
        by_props = process_part_from_props(_fluid(), PROJECT,
                                           **MATERIALS["steel_mild"])
        assert _same(by_props, process_part(_fluid(), PROJECT))

    def test_batch_and_summary_table_carry_a_fluid_region(self):
        results = process_batch([_fluid()], PROJECT)
        assert results[0]["error"] is None
        table = summary_table(results)
        assert "FL-001" in table and "aero_boundary_layer" in table


# --- Service temperature limits ---

class TestServiceLimit:

    @staticmethod
    def _codes(result):
        return [w["code"] for w in result["warnings"]]

    def test_a_temperature_above_a_given_limit_warns(self):
        part = dict(TEMPLATES["structural"], material="plastic_pa66")
        record = dict(MATERIALS["plastic_pa66"], t_service_max_K=423.15)
        r = process_part_from_props(part, PROJECT, **record)
        [w] = [w for w in r["warnings"] if w["code"] == "SERVICE_TEMP_EXCEEDED"]
        assert w["severity"] == "warning"
        assert "473.1 K" in w["message"] and "423.1 K" in w["message"]
        assert "'plastic_pa66'" in w["message"]

    def test_a_temperature_below_the_limit_does_not(self):
        part = dict(TEMPLATES["structural"], material="steel_mild")
        record = dict(MATERIALS["steel_mild"], t_service_max_K=900.0)
        r = process_part_from_props(part, PROJECT, **record)
        assert "SERVICE_TEMP_EXCEEDED" not in self._codes(r)

    def test_no_limit_no_warning(self):
        r = process_part(dict(TEMPLATES["structural"], material="plastic_pa66"),
                         PROJECT)
        assert "SERVICE_TEMP_EXCEEDED" not in self._codes(r)

    def test_a_shield_is_held_to_its_solved_temperature(self):
        part = dict(TEMPLATES["shield"], material="steel_mild")
        solved = process_part(part, PROJECT)["shield"]["t_shield_K"]
        below = dict(MATERIALS["steel_mild"], t_service_max_K=solved - 1.0)
        above = dict(MATERIALS["steel_mild"], t_service_max_K=solved + 1.0)
        assert "SERVICE_TEMP_EXCEEDED" in self._codes(
            process_part_from_props(part, PROJECT, **below))
        assert "SERVICE_TEMP_EXCEEDED" not in self._codes(
            process_part_from_props(part, PROJECT, **above))

    def test_the_table_path_reads_a_limit_once_one_is_filled_in(self, monkeypatch):
        filled = dict(MATERIALS["plastic_pa66"], t_service_max_K=423.15,
                      source="a planted value for this test")
        monkeypatch.setitem(batch_mod.MATERIALS, "plastic_pa66", filled)
        part = dict(TEMPLATES["structural"], material="plastic_pa66")
        assert "SERVICE_TEMP_EXCEEDED" in self._codes(process_part(part, PROJECT))
        assert _same(process_part(part, PROJECT),
                     process_part_from_props(part, PROJECT, **filled))


class TestMaterialRecordSchema:

    @pytest.mark.parametrize("name", sorted(MATERIALS))
    def test_every_table_record_conforms(self, name):
        assert validate_material_record(MATERIALS[name]) == []
        assert isinstance(MATERIALS[name].get("description"), str)

    def test_the_optional_fields_may_be_absent(self):
        assert validate_material_record({"k": 1.0, "rho": 1.0, "cp": 1.0}) == []

    @pytest.mark.parametrize("field, value, code", [
        ("t_service_max_K", -5.0, "OUT_OF_RANGE"),
        ("t_service_max_K", "hot", "WRONG_TYPE"),
        ("source", "", "WRONG_TYPE"),
        ("source", 12, "WRONG_TYPE"),
        ("epsilon", 0.9, "UNKNOWN_KEY"),
    ])
    def test_a_bad_field_is_reported(self, field, value, code):
        record = dict(MATERIALS["steel_mild"], **{field: value})
        assert [(p["key"], p["code"]) for p in
                validate_material_record(record)] == [(field, code)]

    def test_a_missing_property_is_reported(self):
        record = dict(MATERIALS["steel_mild"])
        del record["cp"]
        assert [(p["key"], p["code"]) for p in
                validate_material_record(record)] == [("cp", "MISSING_KEY")]
