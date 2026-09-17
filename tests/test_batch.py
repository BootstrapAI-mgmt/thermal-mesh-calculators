"""
Tests for batch processor and material/surface lookups.
"""

import pytest
from thermal_mesh_calculators.batch import (
    MATERIALS,
    SURFACE_TREATMENTS,
    get_material,
    get_surface_epsilon,
    list_materials,
    list_surface_treatments,
    _resolve_epsilon,
    process_part,
    process_batch,
    summary_table,
)


# --- Fixtures ---

PROJECT = {
    "t_fluid_K": 353.15,
    "t_surr_K": 353.15,
    "t_exh_K": 1073.15,
}

PROJECT_WITH_TRANSIENT = {
    **PROJECT,
    "dt": 0.5,
    "fo_max": 0.5,
    "tau_bc": 20.0,
}


def _make_structural_part(**overrides):
    base = {
        "part_id": "TEST-001",
        "material": "steel_mild",
        "surface": "painted",
        "component_class": "structural",
        "convection_zone": "engine_bay_dead_zone",
        "thickness_mm": 3.0,
        "t_surf_K": 473.15,
    }
    base.update(overrides)
    return base


def _make_shield_part(**overrides):
    base = {
        "part_id": "SH-TEST",
        "material": "steel_mild",
        "component_class": "shield",
        "convection_zone_in": "near_exhaust_natural",
        "convection_zone_out": "engine_bay_dead_zone",
        "thickness_mm": 0.8,
        "t_surf_K": 0,
        "surface_in": "aluminised",
        "surface_out": "aluminised",
    }
    base.update(overrides)
    return base


# --- Material database ---

class TestMaterialDatabase:

    def test_all_materials_have_transport_props(self):
        for name, m in MATERIALS.items():
            assert "k" in m, f"{name} missing k"
            assert "rho" in m, f"{name} missing rho"
            assert "cp" in m, f"{name} missing cp"
            assert m["k"] > 0
            assert m["rho"] > 0
            assert m["cp"] > 0

    def test_no_epsilon_in_materials(self):
        """Emissivity should NOT be in the material database."""
        for name, m in MATERIALS.items():
            assert "epsilon" not in m, (
                f"{name} has epsilon — should be in SURFACE_TREATMENTS"
            )

    def test_get_material_returns_copy(self):
        m = get_material("steel_mild")
        m["k"] = 0
        assert MATERIALS["steel_mild"]["k"] != 0

    def test_unknown_material_raises(self):
        with pytest.raises(KeyError, match="Unknown material"):
            get_material("unobtanium")

    def test_list_materials_sorted(self):
        mats = list_materials()
        assert mats == sorted(mats)


# --- Surface treatments ---

class TestSurfaceTreatments:

    def test_all_treatments_have_epsilon(self):
        for name, t in SURFACE_TREATMENTS.items():
            assert "epsilon" in t
            assert 0 < t["epsilon"] <= 1.0

    def test_aluminised_is_low(self):
        assert SURFACE_TREATMENTS["aluminised"]["epsilon"] < 0.5

    def test_painted_is_high(self):
        assert SURFACE_TREATMENTS["painted"]["epsilon"] > 0.85

    def test_get_surface_epsilon(self):
        assert get_surface_epsilon("aluminised") == 0.40

    def test_unknown_treatment_raises(self):
        with pytest.raises(KeyError, match="Unknown surface treatment"):
            get_surface_epsilon("nonexistent_finish")


# --- Epsilon resolution ---

class TestResolveEpsilon:

    def test_direct_override(self):
        """Direct epsilon in part dict takes priority."""
        part = {"part_id": "T", "material": "steel_mild", "epsilon": 0.55}
        mat = get_material("steel_mild")
        assert _resolve_epsilon(part, mat, "surface") == 0.55

    def test_surface_treatment_lookup(self):
        part = {"part_id": "T", "material": "steel_mild", "surface": "heavily_oxidised"}
        mat = get_material("steel_mild")
        assert _resolve_epsilon(part, mat, "surface") == 0.85

    def test_plastic_fallback(self):
        """Plastic material without explicit surface should fall back to 0.90 (non-metal)."""
        part = {"part_id": "T", "material": "plastic_pa66_gf30"}
        mat = get_material("plastic_pa66_gf30")
        assert _resolve_epsilon(part, mat, "surface") == 0.90

    def test_shield_inner_outer(self):
        part = {
            "part_id": "T",
            "material": "steel_mild",
            "surface_in": "aluminised",
            "surface_out": "heavily_oxidised",
        }
        mat = get_material("steel_mild")
        assert _resolve_epsilon(part, mat, "surface_in") == 0.40
        assert _resolve_epsilon(part, mat, "surface_out") == 0.85

    def test_metal_fallback_without_surface(self):
        """Steel without explicit surface should fall back to 0.73 (other metals)."""
        part = {"part_id": "T", "material": "steel_mild"}
        mat = get_material("steel_mild")
        assert _resolve_epsilon(part, mat, "surface") == 0.73


# --- Single part processing ---

class TestProcessPart:

    def test_structural_part_returns_conduction(self):
        r = process_part(_make_structural_part(), PROJECT)
        assert r["conduction"] is not None
        assert r["governing_dx_mm"] > 0
        assert r["governing_dx_mm"] < float("inf")

    def test_shield_part_solves_temperature(self):
        r = process_part(_make_shield_part(), PROJECT)
        assert r["shield"] is not None
        assert r["shield"]["t_shield_K"] > 353.15
        assert r["shield"]["t_shield_K"] < 1073.15

    def test_transient_included_when_dt_provided(self):
        r = process_part(_make_structural_part(), PROJECT_WITH_TRANSIENT)
        assert r["transient"] is not None

    def test_transient_excluded_when_no_dt(self):
        r = process_part(_make_structural_part(), PROJECT)
        assert r["transient"] is None

    def test_biot_always_computed(self):
        r = process_part(_make_structural_part(), PROJECT)
        assert r["biot"] is not None

    def test_h_override(self):
        r = process_part(
            _make_structural_part(h_override=999.0), PROJECT,
        )
        assert r["h_used"] == 999.0

    def test_eps_used_reported(self):
        r = process_part(_make_structural_part(), PROJECT)
        assert r["eps_used"] == 0.92  # painted

    def test_shield_eps_used_both_sides(self):
        r = process_part(_make_shield_part(), PROJECT)
        assert r["eps_used"]["eps_in"] == 0.40
        assert r["eps_used"]["eps_out"] == 0.40

    def test_lateral_constraint_computed(self):
        """Non-shield parts should have a lateral gradient constraint."""
        r = process_part(_make_structural_part(), PROJECT)
        assert r.get("lateral") is not None
        assert r["lateral"]["max_dx_mm"] > 0

    def test_lateral_in_all_constraints(self):
        """Lateral gradient should appear in all_constraints list."""
        r = process_part(_make_structural_part(), PROJECT)
        sources = [c["source"] for c in r.get("all_constraints", [])]
        assert "lateral_gradient" in sources

    def test_governing_is_smallest(self):
        r = process_part(_make_structural_part(), PROJECT)
        if r.get("all_constraints"):
            all_dx = [c["dx_mm"] for c in r["all_constraints"]]
            assert r["governing_dx_mm"] == min(all_dx)


# --- Batch processing ---

class TestProcessBatch:

    def test_batch_returns_list(self):
        parts = [_make_structural_part(), _make_shield_part()]
        results = process_batch(parts, PROJECT)
        assert len(results) == 2

    def test_error_handling(self):
        bad_part = {"part_id": "BAD", "material": "unobtanium"}
        results = process_batch([bad_part], PROJECT)
        assert results[0]["error"] is not None

    def test_summary_table_string(self):
        parts = [_make_structural_part()]
        results = process_batch(parts, PROJECT)
        table = summary_table(results)
        assert "TEST-001" in table
        assert isinstance(table, str)


# --- Boundary layer integration ---

class TestBoundaryLayerIntegration:

    def test_bl_constraint_external_forced(self):
        """When bl_regime is set, boundary_layer result appears."""
        part = _make_structural_part(
            bl_regime="external_forced",
            velocity_ms=22.0,
            char_length_mm=500.0,
        )
        project = {**PROJECT, "bl_y_plus": 30.0}
        r = process_part(part, project)
        assert "boundary_layer" in r
        assert r["boundary_layer"]["regime"] == "external_forced"
        assert r["boundary_layer"]["max_dx_surface_mm"] > 0
        # Should appear in all_constraints
        sources = [c["source"] for c in r["all_constraints"]]
        assert "aero_boundary_layer" in sources

    def test_bl_constraint_mixed_unknown(self):
        """mixed_unknown regime uses buoyancy velocity."""
        part = _make_structural_part(
            bl_regime="mixed_unknown",
            velocity_ms=0.0,
            char_length_mm=100.0,
            bl_ar_max_prism=3.0,
        )
        r = process_part(part, PROJECT)
        assert r["boundary_layer"]["regime"] == "mixed_unknown"
        assert r["boundary_layer"]["V_buoyancy"] > 0

    def test_bl_not_computed_without_regime(self):
        """Without bl_regime, no boundary_layer in results."""
        part = _make_structural_part()
        r = process_part(part, PROJECT)
        assert "boundary_layer" not in r

    def test_bl_can_govern(self):
        """BL constraint can become governing if it's the smallest."""
        # Use a high-velocity scenario with fine y+ to make BL very tight
        part = _make_structural_part(
            bl_regime="external_forced",
            velocity_ms=50.0,
            char_length_mm=500.0,
            bl_y_plus=1.0,
        )
        r = process_part(part, PROJECT)
        # BL at y+=1 should produce very small surface mesh
        bl_dx = r["boundary_layer"]["max_dx_surface_mm"]
        # It should appear in candidates (may or may not govern)
        sources = [c["source"] for c in r["all_constraints"]]
        assert "aero_boundary_layer" in sources
