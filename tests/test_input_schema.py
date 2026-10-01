"""
The input schema is the published form of what process_part() refuses
======================================================================

PART_SCHEMA and PROJECT_SCHEMA declare every key, and validate_part_input()
and validate_project_input() report the problems of a dict.  process_part()
runs the same checks before it computes anything.  These tests hold the two
to one vocabulary, key by key and class by class:

* parity: for every key, every class and a set of planted values (wrong
  types, values outside the range, its boundaries, unknown names), the
  problems PartInputError carries are exactly the problems the validator
  reports, in the same words; a value the validator accepts is sized;
* the census: the keys the sizing code reads, recorded while it runs on a
  part of each class that carries every key, are exactly the keys the
  schema says that class reads;
* the inputs the schema now refuses that the engine used to accept, each
  named.
"""

import math

import pytest

import thermal_mesh_calculators.batch as batch_mod
from thermal_mesh_calculators.batch import (
    CLASS_DEFAULTS,
    MATERIALS,
    PartInputError,
    process_part,
    process_part_from_props,
)
from thermal_mesh_calculators.schema import (
    COMPONENT_CLASSES,
    EXTENSION_PREFIX,
    MATERIAL_RECORD_SCHEMA,
    PART_SCHEMA,
    PROBLEM_CODES,
    PROJECT_SCHEMA,
    validate_material_record,
    validate_part_input,
    validate_project_input,
)
from thermal_mesh_calculators.zones import CONVECTION_ZONES


PROJECT = {"t_fluid_K": 353.15, "t_surr_K": 353.15, "t_exh_K": 1073.15}

BASE = {
    "structural": {"part_id": "B-STR", "material": "steel_mild",
                   "component_class": "structural",
                   "convection_zone": "engine_beside", "thickness_mm": 3.0,
                   "t_surf_K": 473.15, "surface": "painted"},
    "exhaust": {"part_id": "B-EXH", "material": "cast_iron",
                "component_class": "exhaust",
                "convection_zone": "exhaust_internal", "thickness_mm": 6.0,
                "t_surf_K": 1073.15, "surface": "cast_iron_oxidised"},
    "exhaust_adjacent": {"part_id": "B-ADJ", "material": "plastic_pa66_gf30",
                         "component_class": "exhaust_adjacent",
                         "convection_zone": "exhaust_beside",
                         "thickness_mm": 2.5, "t_surf_K": 400.0,
                         "surface": "plastic"},
    "shield": {"part_id": "B-SH", "material": "steel_mild",
               "component_class": "shield",
               "convection_zone_in": "exhaust_beside",
               "convection_zone_out": "engine_beside", "thickness_mm": 0.8,
               "surface_in": "aluminised", "surface_out": "aluminised"},
    "multilayer_shield": {"part_id": "B-ML", "material": "steel_mild",
                          "component_class": "multilayer_shield",
                          "convection_zone_in": "exhaust_beside",
                          "convection_zone_out": "shield_gap_confined",
                          "thickness_mm": 0.6, "surface_in": "aluminised",
                          "surface_out": "aluminised",
                          "surface_g1": "aluminised",
                          "surface_g2": "aluminised", "h_gap": 12.0},
    "fluid": {"part_id": "B-FL", "component_class": "fluid",
              "convection_zone": "front_end_edges"},
}

# A typical value for each number that must be above zero, to plant as a
# value the schema accepts.
TYPICAL = {
    "thickness_mm": 2.5, "t_surf_K": 450.0, "t_fluid_K": 340.0,
    "char_length_mm": 150.0, "t_exh_K": 900.0, "max_dt": 12.0,
    "allowable_flux_error": 400.0, "radius_mm": 30.0, "bl_y_plus": 1.0,
    "t_surr_K": 340.0, "dt": 1.0, "fo_max": 0.5, "safety_factor": 1.5,
    "tau_bc": 30.0, "k": 40.0, "rho": 7800.0, "cp": 470.0,
    "t_service_max_K": 900.0,
}

NAN, INF = float("nan"), float("inf")


def test_the_classes_are_the_keys_of_class_defaults():
    assert set(COMPONENT_CLASSES) == set(CLASS_DEFAULTS)


def test_every_base_part_conforms_and_sizes():
    for cls in COMPONENT_CLASSES:
        assert validate_part_input(BASE[cls]) == [], cls
        assert 0 < process_part(BASE[cls], PROJECT)["governing_dx_mm"] < INF


def planted(key, spec):
    """(value, the schema refuses it) pairs for one key."""
    kind = spec["type"]
    out = []
    if spec.get("nullable"):
        out.append((None, False))
    elif kind != "list":
        out.append((None, True))
    if kind == "number":
        out += [("x", True), (True, True), (NAN, True), (INF, True),
                ([1.0], True)]
        op, lo = spec["range"][0], spec["range"][1]
        if op == ">" and lo == 0.0:
            out += [(0.0, True), (-1.0, True), (TYPICAL[key], False)]
        elif op == ">":
            out += [(lo, True), (lo - 0.5, True), (lo + 0.3, False)]
        elif op == ">=":
            out += [(lo, False), (lo - 0.5, True), (lo + 2.0, False)]
        elif op == "[]":
            hi = spec["range"][2]
            out += [(lo, False), (hi, False), (lo - 0.01, True),
                    (hi + 0.01, True)]
        else:                                            # "(]"
            hi = spec["range"][2]
            out += [(lo, True), (hi, False), (hi + 0.01, True), (0.3, False)]
    elif kind == "integer":
        out += [(0, True), (-3, True), (1.5, True), ("10", True),
                (True, True), (1, False), (50, False)]
    elif kind == "list":
        out += [("x", True), ([], False)]
    elif "vocabulary" in spec:
        good = {"MATERIALS": "cast_iron", "SURFACE_TREATMENTS": "painted",
                "CONVECTION_ZONES": "cabin_tunnel",
                "CLASS_DEFAULTS": "structural"}[spec["vocabulary"]]
        out += [("nope", True), (3, True), (good, False)]
    elif "values" in spec:
        out += [("nope", True), (3, True)] + [(v, False) for v in spec["values"]]
    else:
        out += [("", bool(spec.get("non_empty"))), (3, True), ("P-9", False)]
    return out


def _raised(part, project=PROJECT):
    """The problems process_part() raises with, or None when it sizes."""
    try:
        result = process_part(part, project)
    except PartInputError as exc:
        return exc.problems
    assert 0 < result["governing_dx_mm"] <= INF
    return None


def _project_problems(project):
    return [dict(p, key="project." + p["key"] if p["key"] else "project")
            for p in validate_project_input(project)]


# --- Parity, key by key ---

_PARITY_KEYS = sorted(k for k in PART_SCHEMA
                      if k not in ("component_class", "schema_version",
                                   "bom_errors"))


@pytest.mark.parametrize("key", _PARITY_KEYS)
def test_part_key_parity(key):
    spec = PART_SCHEMA[key]
    for cls in COMPONENT_CLASSES:
        checked = cls in spec["read_by"] or spec.get("checked_for") == "all"
        for value, refused in planted(key, spec):
            part = dict(BASE[cls], **{key: value})
            reported = validate_part_input(part)
            raised = _raised(part)
            where = f"{cls}, {key}={value!r}"
            if checked and refused:
                assert reported, f"the schema accepts {where}"
                assert [p["key"] for p in reported] == [key], where
                assert raised == reported, where
            else:
                assert reported == [], f"the schema refuses {where}: {reported}"
                assert raised is None, f"process_part refuses {where}: {raised}"


@pytest.mark.parametrize("key", sorted(k for k in PART_SCHEMA
                                       if PART_SCHEMA[k].get("required")))
def test_required_key_parity(key):
    for rule in PART_SCHEMA[key]["required"]:
        for cls in rule["classes"]:
            part = dict(BASE[cls])
            for gone in (key,) + tuple(k for alt in rule["unless"] for k in alt):
                part.pop(gone, None)
            reported = validate_part_input(part)
            assert [(p["key"], p["code"]) for p in reported] == [
                (key, "MISSING_KEY")], (cls, reported)
            assert _raised(part) == reported
            for alternative in rule["unless"]:
                given = dict(part)
                for k in alternative:
                    given[k] = {"h_override": 20.0, "h_in_override": 20.0,
                                "h_out_override": 20.0, "velocity_ms": 5.0,
                                "bl_regime": "external_forced",
                                "convection_zone": "exhaust_beside"}[k]
                assert validate_part_input(given) == [], (cls, alternative)
                assert _raised(given) is None, (cls, alternative)


@pytest.mark.parametrize("value, code", [
    ("nope", "UNKNOWN_NAME"), (3, "UNKNOWN_NAME"), (None, "MISSING_KEY"),
])
def test_component_class_problem_is_reported_alone(value, code):
    part = dict(BASE["structural"], component_class=value, thickness_mm=-1.0,
                material="nope")
    if value is None:
        del part["component_class"]
    reported = validate_part_input(part)
    assert [(p["key"], p["code"]) for p in reported] == [("component_class", code)]
    assert _raised(part) == reported


@pytest.mark.parametrize("key", sorted(k for k in PROJECT_SCHEMA
                                       if k != "schema_version"))
def test_project_key_parity(key):
    spec = PROJECT_SCHEMA[key]
    planted_values = planted(key, spec)
    if spec.get("required"):
        planted_values.append(("<absent>", True))
    for value, refused in planted_values:
        project = dict(PROJECT, dt=1.0) if key in ("fo_max", "safety_factor",
                                                   "tau_bc") else dict(PROJECT)
        if value == "<absent>":
            del project[key]
        else:
            project[key] = value
        for cls in COMPONENT_CLASSES:
            reported = _project_problems(project)
            raised = _raised(BASE[cls], project)
            where = f"{cls}, project {key}={value!r}"
            if refused:
                assert [p["key"] for p in reported] == ["project." + key], where
                assert raised == reported, where
            else:
                assert reported == [], where
                assert raised is None, where


@pytest.mark.parametrize("field", sorted(MATERIAL_RECORD_SCHEMA))
def test_material_record_parity(field):
    spec = MATERIAL_RECORD_SCHEMA[field]
    for value, refused in planted(field, dict(spec, nullable=False)):
        record = dict(MATERIALS["steel_mild"], **{field: value})
        reported = [dict(p, key="properties." + p["key"])
                    for p in validate_material_record(record)]
        try:
            process_part_from_props(BASE["structural"], PROJECT, **record)
            raised = None
        except PartInputError as exc:
            raised = exc.problems
        if value is None:
            # None means "not given" to process_part_from_props().
            assert raised is None if field not in ("k", "rho", "cp") else raised
            continue
        assert raised == (reported or None), (field, value)
        assert bool(reported) == refused, (field, value)


# --- The census: the sizing code reads exactly the declared keys ---

class _Recording(dict):
    """A dict that records every key it is asked about."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.read = set()

    def __getitem__(self, key):
        self.read.add(key)
        return super().__getitem__(key)

    def get(self, key, default=None):
        self.read.add(key)
        return super().get(key, default)

    def __contains__(self, key):
        self.read.add(key)
        return super().__contains__(key)


_INTAKE_ONLY = {"schema_version", "bom_row", "bom_errors"}

# Every key a class reads, with values that make every reading path run
# (the boundary layer on, a radius, every shield input given).
_FULL_EXTRAS = {
    "t_fluid_K": 360.0, "char_length_mm": 120.0, "radius_mm": 40.0,
    "velocity_ms": 3.0, "bl_regime": "mixed_unknown", "bl_y_plus": 30.0,
    "bl_growth_ratio": 1.2, "bl_ar_max_prism": 5.0, "bl_fraction": 0.3,
    "max_dt": 12.0, "allowable_flux_error": 400.0,
}
_FULL = {
    "structural": dict(BASE["structural"], h_override=20.0, epsilon=0.8),
    "exhaust": dict(BASE["exhaust"], epsilon=0.8),
    "exhaust_adjacent": dict(BASE["exhaust_adjacent"]),
    "shield": dict(BASE["shield"], convection_zone="engine_beside",
                   eps_in=0.4, eps_out=0.4, h_in_override=20.0,
                   h_out_override=20.0, t_exh_K=900.0, shield_max_iter=60),
    "multilayer_shield": dict(BASE["multilayer_shield"],
                              convection_zone="engine_beside", eps_in=0.4,
                              eps_out=0.4, eps_g1=0.4, eps_g2=0.4, f12=0.9,
                              h_in_override=20.0, h_out_override=20.0,
                              t_exh_K=900.0, shield_max_iter=60),
    "fluid": dict(BASE["fluid"], t_surf_K=420.0),
}
_PROJECT_FULL = dict(PROJECT, dt=1.0, fo_max=0.5, safety_factor=1.0,
                     tau_bc=20.0, transient_scheme="explicit",
                     shield_max_iter=60, max_dt=12.0,
                     allowable_flux_error=400.0, bl_regime="external_forced",
                     bl_y_plus=30.0, bl_growth_ratio=1.2, bl_ar_max_prism=5.0,
                     bl_fraction=0.3)


def _census(cls, part, project):
    part = _Recording(part)
    project = _Recording(project)
    if cls == "fluid":
        batch_mod._size_fluid(part, project)
    else:
        batch_mod._size_solid(part, project, MATERIALS[part["material"]],
                              part["material"])
    return part.read, project.read


@pytest.mark.parametrize("cls", COMPONENT_CLASSES)
def test_census_the_engine_reads_what_the_schema_declares(cls):
    # A given emissivity or h short-circuits the lookup it overrides, so
    # the keys read are the union over a part with every override and the
    # same part without them.
    full = dict(_FULL[cls], **_FULL_EXTRAS)
    plain = dict(BASE[cls], **_FULL_EXTRAS)
    read = set()
    for part in (full, plain):
        assert validate_part_input(part) == []
        read |= _census(cls, part, PROJECT)[0]
    declared = {k for k, s in PART_SCHEMA.items() if cls in s["read_by"]}
    assert read - set(PART_SCHEMA) == set(), "read but not declared"
    assert read - declared == set(), f"{cls} reads keys not declared for it"
    assert declared - _INTAKE_ONLY - read == set(), (
        f"declared for {cls} but never read")


def test_census_of_the_project():
    # The base parts give none of the project's keys, so every project
    # value is the one read.
    read = set()
    for cls in COMPONENT_CLASSES:
        read |= _census(cls, BASE[cls], _PROJECT_FULL)[1]
    assert read - set(PROJECT_SCHEMA) == set(), "read but not declared"
    assert set(PROJECT_SCHEMA) - {"schema_version"} - read == set(), (
        "declared but never read")


# --- Keys the schema does not declare ---

def test_an_unknown_key_is_refused_with_the_key_named():
    part = dict(BASE["structural"], thicknes_mm=3.0)
    reported = validate_part_input(part)
    assert [(p["key"], p["code"]) for p in reported] == [
        ("thicknes_mm", "UNKNOWN_KEY")]
    assert reported[0]["reason"].startswith(
        "Unknown key 'thicknes_mm' (part 'B-STR'). Keys starting with 'x_' "
        "are carried through unread. Available: ")
    assert _raised(part) == reported


def test_an_extension_key_is_carried_and_never_read():
    part = dict(BASE["structural"], x_plm_id="PLM-0042", x_rev={"r": 3})
    assert validate_part_input(part) == []
    assert process_part(part, PROJECT) == process_part(BASE["structural"], PROJECT)
    assert EXTENSION_PREFIX == "x_"


def test_an_unknown_project_key_is_refused():
    project = dict(PROJECT, t_ambient_K=300.0)
    assert [(p["key"], p["code"]) for p in validate_project_input(project)] == [
        ("t_ambient_K", "UNKNOWN_KEY")]
    assert _raised(BASE["fluid"], project) == _project_problems(project)


# --- Versions, collected problems, the error itself ---

@pytest.mark.parametrize("value, code", [
    (2, "UNSUPPORTED_VERSION"), ("1", "WRONG_TYPE"), (True, "WRONG_TYPE"),
])
def test_an_unsupported_schema_version(value, code):
    part = dict(BASE["structural"], schema_version=value)
    reported = validate_part_input(part)
    assert [(p["key"], p["code"]) for p in reported] == [("schema_version", code)]
    assert _raised(part) == reported
    project = dict(PROJECT, schema_version=value)
    assert [p["code"] for p in validate_project_input(project)] == [code]


def test_schema_version_1_is_accepted_and_changes_nothing():
    part = dict(BASE["structural"], schema_version=1)
    assert process_part(part, dict(PROJECT, schema_version=1)) == process_part(
        BASE["structural"], PROJECT)


def test_every_problem_is_collected():
    part = dict(BASE["structural"], thickness_mm=-3.0, epsilon=1.4,
                surface="glossy", wall="thin")
    project = dict(PROJECT, t_surr_K=0.0)
    with pytest.raises(PartInputError) as excinfo:
        process_part(part, project)
    problems = [(p["key"], p["code"]) for p in excinfo.value.problems]
    assert problems == [("wall", "UNKNOWN_KEY"),
                        ("thickness_mm", "OUT_OF_RANGE"),
                        ("surface", "UNKNOWN_NAME"),
                        ("epsilon", "OUT_OF_RANGE"),
                        ("project.t_surr_K", "OUT_OF_RANGE")]
    lines = str(excinfo.value).split("\n")
    assert lines[0] == "5 problems with part 'B-STR':"
    assert lines[2] == "  - thickness_mm must be > 0 mm, got -3.0"
    assert lines[5] == "  - t_surr_K must be > 0 K, got 0.0 (project)"


def test_every_problem_code_is_registered():
    codes = set()
    for cls in COMPONENT_CLASSES:
        for key, spec in PART_SCHEMA.items():
            for value, _ in planted(key, spec):
                part = dict(BASE[cls], **{key: value})
                codes |= {p["code"] for p in validate_part_input(part)}
    assert codes <= set(PROBLEM_CODES)


def test_a_part_that_is_not_a_dict():
    for part in (["part_id"], "P-1", None):
        assert [p["code"] for p in validate_part_input(part)] == ["NOT_AN_OBJECT"]
        with pytest.raises(PartInputError):
            process_part(part, PROJECT)


# --- Inputs refused now that the engine used to accept ---

@pytest.mark.parametrize("cls, change, code", [
    # 0.6.2 echoed a non-string id; summary_table() then failed on it.
    ("structural", {"part_id": 1001}, "WRONG_TYPE"),
    # 0.6.2 sized a negative wall: a negative Biot number, a negative size.
    ("structural", {"thickness_mm": -3.0}, "OUT_OF_RANGE"),
    # 0.6.2 read max_dt = 0 as "use the class default".
    ("structural", {"max_dt": 0.0}, "OUT_OF_RANGE"),
    ("shield", {"allowable_flux_error": 0.0}, "OUT_OF_RANGE"),
    # 0.6.2 ran any regime other than "mixed_unknown" as the forced one.
    ("structural", {"bl_regime": "laminar"}, "UNKNOWN_NAME"),
    # 0.6.2 returned a negative curvature limit and a negative BL size.
    ("structural", {"radius_mm": -5.0}, "OUT_OF_RANGE"),
    ("structural", {"bl_regime": "external_forced", "bl_fraction": -0.3},
     "OUT_OF_RANGE"),
    # 0.6.2 ignored a key it did not know, a misspelt override included.
    ("structural", {"h_overide": 20.0}, "UNKNOWN_KEY"),
])
def test_now_refused(cls, change, code):
    part = dict(BASE[cls], **change)
    reported = validate_part_input(part)
    assert [p["code"] for p in reported] == [code]
    assert _raised(part) == reported


def test_a_shield_keeps_its_ignored_surface_temperature():
    # The documented minimum input gives a shield a t_surf_K it ignores.
    part = dict(BASE["shield"], t_surf_K=0)
    assert validate_part_input(part) == []
    assert process_part(part, PROJECT) == process_part(BASE["shield"], PROJECT)


def test_a_name_is_checked_whatever_the_class():
    part = dict(BASE["structural"], surface_in="glossy")
    assert [(p["key"], p["code"]) for p in validate_part_input(part)] == [
        ("surface_in", "UNKNOWN_NAME")]
    zones_part = dict(BASE["fluid"], convection_zone_out="nowhere")
    assert [p["key"] for p in validate_part_input(zones_part)] == [
        "convection_zone_out"]


def test_the_material_label_of_the_properties_path_is_free():
    part = dict(BASE["structural"], material="house grade 7")
    assert [p["key"] for p in validate_part_input(part)] == ["material"]
    assert validate_part_input(part, material_from_table=False) == []
    no_label = dict(BASE["structural"])
    del no_label["material"]
    assert validate_part_input(no_label, material_from_table=False) == []
    assert [p["code"] for p in validate_part_input(
        dict(part, material=7), material_from_table=False)] == ["WRONG_TYPE"]


def test_zone_vocabulary_is_the_eighteen_keys():
    assert len(CONVECTION_ZONES) == 18
    assert all(math.isfinite(z["velocity_ms"]) for z in CONVECTION_ZONES.values())
