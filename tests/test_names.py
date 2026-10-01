"""
Material and zone names, and the component class a part's keys imply
=====================================================================

The alias tables are data: every alias is stored in its normalised form,
none shadows a canonical key, and every target is a canonical key.  A
canonical key resolves to itself; a family name is refused with the
entries it could mean.  infer_component_class() applies its rules in
order, and never infers a fluid region.
"""

import pytest

from thermal_mesh_calculators.batch import (
    CLASS_DEFAULTS,
    MATERIALS,
    PartInputError,
)
from thermal_mesh_calculators.intake import (
    EXHAUST_ADJACENT_ZONES,
    EXHAUST_ZONES,
    MATERIAL_ALIASES,
    SHIELD_KEYS,
    TWO_LAYER_SHIELD_KEYS,
    ZONE_ALIASES,
    _normalise,
    infer_component_class,
    resolve_material,
    resolve_zone,
)
from thermal_mesh_calculators.schema import PART_SCHEMA, SHIELD_CLASSES
from thermal_mesh_calculators.zones import CONVECTION_ZONES


@pytest.mark.parametrize("aliases, table", [
    (MATERIAL_ALIASES, MATERIALS),
    (ZONE_ALIASES, CONVECTION_ZONES),
], ids=["materials", "zones"])
def test_an_alias_table_is_well_formed(aliases, table):
    for alias, target in aliases.items():
        assert _normalise(alias) == alias, alias
        assert alias not in table, f"{alias} shadows a canonical key"
        assert target in table, f"{alias} -> {target} is not a key"


def test_every_canonical_key_resolves_to_itself():
    assert [resolve_material(k) for k in MATERIALS] == list(MATERIALS)
    assert [resolve_zone(k) for k in CONVECTION_ZONES] == list(CONVECTION_ZONES)


def test_every_alias_resolves_to_its_target():
    for alias, target in MATERIAL_ALIASES.items():
        assert resolve_material(alias) == target
    for alias, target in ZONE_ALIASES.items():
        assert resolve_zone(alias) == target


@pytest.mark.parametrize("written, key", [
    ("Mild steel", "steel_mild"),
    ("STEEL_MILD", "steel_mild"),
    ("steel-mild", "steel_mild"),
    ("Aluminum 6061-T6", "aluminium_6061"),
    ("Galvanized steel", "steel_galvanised"),
    ("Gray cast iron", "cast_iron"),
    ("Fiberglass insulation", "insulation_fibreglass"),
    ("SS-316", "steel_stainless_304"),
    ("PA66-GF30", "plastic_pa66_gf30"),
    ("Ti-6Al-4V", "titanium_6al4v"),
    ("PC/ABS", "plastic_pc_abs"),
    ("Inconel 625", "nickel_alloy"),
    ("Viton", "rubber_fkm"),
])
def test_common_spellings_of_a_material(written, key):
    assert resolve_material(written) == key


@pytest.mark.parametrize("written, key", [
    ("Beside Engine", "engine_beside"),
    ("engine-beside", "engine_beside"),
    ("Above exhaust", "exhaust_above"),
    ("Exhaust gas", "exhaust_internal"),
    ("Clutch outlet", "clutch_outlet_downstream"),
    ("engine_bay_dead_zone", "engine_bay_dead_zone"),   # a legacy key
])
def test_common_spellings_of_a_zone(written, key):
    assert resolve_zone(written) == key


def test_a_family_name_is_refused_with_its_candidates():
    with pytest.raises(PartInputError) as excinfo:
        resolve_material("aluminium")
    message = str(excinfo.value)
    assert message.startswith("Unknown material 'aluminium'. 'aluminium' "
                              "could be any of: aluminium_5052, "
                              "aluminium_6061, cast_aluminium, "
                              "cast_aluminium_a380; give the one you mean.")
    assert message.endswith("(and the spellings in MATERIAL_ALIASES)")
    assert ", ".join(sorted(MATERIALS)) in message


def test_a_family_with_one_member_names_it():
    with pytest.raises(PartInputError, match="Did you mean 'titanium_6al4v'"):
        resolve_material("Titanium")


def test_a_zone_family_is_refused_with_its_candidates():
    with pytest.raises(PartInputError, match="could be any of: exhaust_above, "
                       "exhaust_below, exhaust_beside, exhaust_internal"):
        resolve_zone("exhaust")


@pytest.mark.parametrize("name", ["unobtanium", "", 42, None])
def test_an_unknown_name_names_the_allowed_set(name):
    with pytest.raises(PartInputError) as excinfo:
        resolve_material(name)
    assert "could be any of" not in str(excinfo.value)
    assert "Available: " + ", ".join(sorted(MATERIALS)) in str(excinfo.value)
    assert [p["code"] for p in excinfo.value.problems] == ["UNKNOWN_NAME"]


# --- infer_component_class ---

class TestInferComponentClass:

    def test_a_given_class_is_returned(self):
        assert infer_component_class({"component_class": "shield",
                                      "convection_zone": "exhaust_internal"}) == "shield"

    def test_a_given_unknown_class_is_refused_by_name(self):
        with pytest.raises(PartInputError) as excinfo:
            infer_component_class({"part_id": "P", "component_class": "exhuast"})
        assert str(excinfo.value).startswith(
            "Unknown component class 'exhuast' (part 'P', key 'component_class').")

    @pytest.mark.parametrize("key", TWO_LAYER_SHIELD_KEYS)
    def test_a_two_layer_key_makes_a_two_layer_shield(self, key):
        assert infer_component_class({key: 1, "surface_in": "painted",
                                      "convection_zone": "exhaust_internal"}
                                     ) == "multilayer_shield"

    @pytest.mark.parametrize("key", SHIELD_KEYS)
    def test_a_shield_key_makes_a_shield(self, key):
        assert infer_component_class({key: 1, "convection_zone":
                                      "exhaust_internal"}) == "shield"

    @pytest.mark.parametrize("zone", EXHAUST_ZONES)
    def test_an_exhaust_zone_makes_an_exhaust_part(self, zone):
        assert infer_component_class({"convection_zone": zone}) == "exhaust"

    @pytest.mark.parametrize("zone", EXHAUST_ADJACENT_ZONES)
    def test_a_zone_beside_the_exhaust_makes_an_adjacent_part(self, zone):
        assert infer_component_class({"convection_zone": zone}) == "exhaust_adjacent"

    def test_any_other_part_is_structural(self):
        assert infer_component_class({"convection_zone": "cabin_tunnel"}) == "structural"
        assert infer_component_class({}) == "structural"

    def test_the_zone_is_read_through_its_aliases(self):
        assert infer_component_class({"convection_zone": "Exhaust gas"}) == "exhaust"
        assert infer_component_class({"convection_zone": "above exhaust"}) == "exhaust_adjacent"

    def test_a_fluid_region_is_never_inferred(self):
        inferred = {infer_component_class({"convection_zone": zone})
                    for zone in CONVECTION_ZONES}
        assert "fluid" not in inferred
        assert inferred <= set(CLASS_DEFAULTS)

    def test_the_rule_keys_are_shield_keys_of_the_schema(self):
        for key in TWO_LAYER_SHIELD_KEYS:
            assert PART_SCHEMA[key]["read_by"] == ("multilayer_shield",), key
        for key in SHIELD_KEYS:
            assert set(PART_SCHEMA[key]["read_by"]) == set(SHIELD_CLASSES), key
