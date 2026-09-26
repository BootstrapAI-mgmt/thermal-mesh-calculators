"""
Tests for batch processor and material/surface lookups.
"""

import re

import pytest
import thermal_mesh_calculators.batch as batch_mod
from thermal_mesh_calculators.h_estimator import estimate_h
from thermal_mesh_calculators.radiation import RadiationMeshCalculator
from thermal_mesh_calculators.shields import SingleLayerShieldCalculator
from thermal_mesh_calculators.zones import CONVECTION_ZONES
from thermal_mesh_calculators.batch import (
    CLASS_DEFAULTS,
    MATERIALS,
    SURFACE_TREATMENTS,
    get_material,
    get_surface_epsilon,
    list_materials,
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
        # It should appear in candidates (may or may not govern), with its own size
        sources = [c["source"] for c in r["all_constraints"]]
        assert "aero_boundary_layer" in sources
        assert {"dx_mm": bl_dx, "source": "aero_boundary_layer"} in r["all_constraints"]


# --- Part validation (unknown names raise; defaulted surfaces warn) ---

def _make_multilayer_part(**overrides):
    base = {
        "part_id": "ML-TEST",
        "material": "steel_mild",
        "component_class": "multilayer_shield",
        "convection_zone_in": "near_exhaust_natural",
        "convection_zone_out": "shield_gap_confined",
        "thickness_mm": 0.6,
        "t_surf_K": 0,
        "surface_in": "aluminised",
        "surface_out": "aluminised",
        "surface_g1": "aluminised",
        "surface_g2": "aluminised",
        "h_gap": 12.0,
    }
    base.update(overrides)
    return base


def _warning_codes(result):
    return [w["code"] for w in result["warnings"]]


class TestPartValidation:

    @pytest.mark.parametrize("cls", sorted(CLASS_DEFAULTS))
    def test_every_component_class_is_accepted(self, cls):
        if cls == "shield":
            part = _make_shield_part()
        elif cls == "multilayer_shield":
            part = _make_multilayer_part()
        else:
            part = _make_structural_part(component_class=cls)
        r = process_part(part, PROJECT)
        assert r["component_class"] == cls
        assert 0 < r["governing_dx_mm"] < float("inf")

    def test_unknown_class_names_part_key_and_allowed_set(self):
        with pytest.raises(batch_mod.PartInputError) as excinfo:
            process_part(_make_structural_part(component_class="exhuast"), PROJECT)
        message = str(excinfo.value)
        assert message.startswith("Unknown component class 'exhuast'")
        assert "part 'TEST-001'" in message
        assert "key 'component_class'" in message
        assert message.endswith("Available: " + ", ".join(sorted(CLASS_DEFAULTS)))

    def test_part_input_error_is_a_key_error_and_a_value_error(self):
        """Callers that caught KeyError from the old lookups keep working."""
        assert issubclass(batch_mod.PartInputError, KeyError)
        assert issubclass(batch_mod.PartInputError, ValueError)

    def test_missing_component_class_raises(self):
        part = _make_structural_part()
        del part["component_class"]
        with pytest.raises(batch_mod.PartInputError, match="Missing component class"):
            process_part(part, PROJECT)

    def test_unknown_material_raises_with_the_allowed_set(self):
        with pytest.raises(batch_mod.PartInputError) as excinfo:
            process_part(_make_structural_part(material="unobtanium"), PROJECT)
        message = str(excinfo.value)
        assert message.startswith("Unknown material 'unobtanium'")
        assert message.endswith("Available: " + ", ".join(sorted(MATERIALS)))

    @pytest.mark.parametrize(
        "key", ["convection_zone", "convection_zone_in", "convection_zone_out"],
    )
    def test_unknown_zone_raises_on_every_zone_key(self, key):
        part = _make_shield_part(**{key: "exhaust_abvoe"})
        with pytest.raises(batch_mod.PartInputError) as excinfo:
            process_part(part, PROJECT)
        message = str(excinfo.value)
        assert "Unknown convection zone 'exhaust_abvoe'" in message
        assert "key '%s'" % key in message
        assert message.endswith("Available: " + ", ".join(sorted(CONVECTION_ZONES)))

    @pytest.mark.parametrize(
        "key", ["surface", "surface_in", "surface_out", "surface_g1", "surface_g2"],
    )
    def test_unknown_surface_raises_on_every_surface_key(self, key):
        part = _make_multilayer_part(**{key: "aluminized"})
        with pytest.raises(batch_mod.PartInputError) as excinfo:
            process_part(part, PROJECT)
        message = str(excinfo.value)
        assert "Unknown surface treatment 'aluminized'" in message
        assert "key '%s'" % key in message
        assert message.endswith(
            "Available: " + ", ".join(sorted(SURFACE_TREATMENTS)))

    def test_non_shield_needs_a_zone_or_an_h_override(self):
        part = _make_structural_part()
        del part["convection_zone"]
        with pytest.raises(batch_mod.PartInputError, match="Missing convection zone"):
            process_part(part, PROJECT)
        part["h_override"] = 12.0
        r = process_part(part, PROJECT)
        assert r["h_used"] == 12.0

    def test_shield_side_needs_a_zone_or_an_h_override(self):
        part = _make_shield_part()
        del part["convection_zone_out"]
        with pytest.raises(batch_mod.PartInputError) as excinfo:
            process_part(part, PROJECT)
        assert "'convection_zone_out'" in str(excinfo.value)
        assert "h_out_override" in str(excinfo.value)

    def test_shield_with_both_h_overrides_needs_no_zone(self):
        """The overrides are used as given; no zone is looked up."""
        part = _make_shield_part(h_in_override=22.0, h_out_override=9.0)
        del part["convection_zone_in"]
        del part["convection_zone_out"]
        r = process_part(part, PROJECT)
        assert r["h_used"] == {"h_in": 22.0, "h_out": 9.0}

    def test_process_batch_reports_the_allowed_set(self):
        results = process_batch(
            [_make_structural_part(component_class="exhuast")], PROJECT,
        )
        assert results[0]["error"].startswith("Unknown component class 'exhuast'")
        assert "Available: " in results[0]["error"]


class TestDefaultedSurfaceWarnings:

    def test_given_surface_or_epsilon_is_not_reported(self):
        assert "SURFACE_DEFAULTED" not in _warning_codes(
            process_part(_make_structural_part(), PROJECT))
        part = _make_structural_part(epsilon=0.6)
        del part["surface"]
        r = process_part(part, PROJECT)
        assert r["eps_used"] == 0.6
        assert "SURFACE_DEFAULTED" not in _warning_codes(r)

    def test_defaulted_surface_keeps_the_fallback_value_and_says_so(self):
        part = _make_structural_part(material="plastic_pa66_gf30")
        del part["surface"]
        r = process_part(part, PROJECT)
        assert r["eps_used"] == 0.90
        [w] = [w for w in r["warnings"] if w["code"] == "SURFACE_DEFAULTED"]
        assert w["key"] == "surface"
        assert w["severity"] == "caution"
        assert "0.90" in w["message"]
        assert "non-metal fallback" in w["message"]

    def test_each_defaulted_shield_surface_is_reported(self):
        part = _make_shield_part()
        del part["surface_in"]
        del part["surface_out"]
        r = process_part(part, PROJECT)
        keys = sorted(w["key"] for w in r["warnings"]
                      if w["code"] == "SURFACE_DEFAULTED")
        assert keys == ["surface_in", "surface_out"]

    def test_defaulted_gap_surfaces_are_reported(self):
        part = _make_multilayer_part()
        del part["surface_g1"]
        del part["surface_g2"]
        r = process_part(part, PROJECT)
        keys = sorted(w["key"] for w in r["warnings"]
                      if w["code"] == "SURFACE_DEFAULTED")
        assert keys == ["surface_g1", "surface_g2"]
        assert r["eps_used"]["eps_g1"] == 0.73


# --- One fluid temperature for h and q'' ---

class TestFluidTemperature:

    def test_project_fluid_temperature_is_reported(self):
        r = process_part(_make_structural_part(), PROJECT)
        assert r["t_fluid_K"] == PROJECT["t_fluid_K"]
        assert r["t_fluid_source"] == "project"

    def test_h_moves_with_the_fluid_temperature(self):
        """Natural convection: a larger surface-to-fluid difference gives a
        larger h, so h must fall as t_fluid rises toward t_surf."""
        h = [
            process_part(_make_structural_part(),
                         {**PROJECT, "t_fluid_K": tf, "t_surr_K": tf})["h_used"]
            for tf in (300.0, 353.15, 450.0)
        ]
        assert h[0] > h[1] > h[2] > 0

    def test_part_fluid_temperature_overrides_the_project(self):
        part = _make_structural_part(t_fluid_K=423.15)
        r = process_part(part, PROJECT)
        assert r["t_fluid_K"] == 423.15
        assert r["t_fluid_source"] == "part"
        expected_h = estimate_h(
            velocity=0.0, t_surf=473.15, t_fluid=423.15,
            char_length=0.1, orientation="vertical",
        )["h"]
        assert r["h_used"] == pytest.approx(expected_h, rel=1e-12)
        assert r["conduction"]["q_conv"] == pytest.approx(
            expected_h * (473.15 - 423.15), rel=1e-12)

    def test_h_override_keeps_the_same_fluid_temperature_for_q(self):
        r = process_part(_make_structural_part(h_override=40.0, t_fluid_K=400.0),
                         PROJECT)
        assert r["conduction"]["q_conv"] == pytest.approx(40.0 * (473.15 - 400.0))

    def test_shield_uses_the_part_fluid_temperature(self):
        part = _make_shield_part(t_fluid_K=393.15)
        r = process_part(part, PROJECT)
        direct = SingleLayerShieldCalculator.mesh_size(
            k=54.0, max_dt=15.0, t_exh=1073.15, t_fluid=393.15, t_surr=353.15,
            h_in=r["h_used"]["h_in"], h_out=r["h_used"]["h_out"],
            eps_in=0.40, eps_out=0.40,
        )
        assert r["shield"]["t_shield_K"] == pytest.approx(direct["t_shield_K"])


# --- Shield non-convergence and unbounded sizes are reported ---

class TestShieldSolveWarnings:

    def test_converged_shields_raise_no_solver_warning(self):
        for part in (_make_shield_part(), _make_multilayer_part()):
            codes = _warning_codes(process_part(part, PROJECT))
            assert "SHIELD_NOT_CONVERGED" not in codes
            assert "NO_FINITE_SIZE" not in codes

    @pytest.mark.parametrize("maker", [_make_shield_part, _make_multilayer_part])
    def test_unconverged_shield_warns(self, maker):
        r = process_part(maker(shield_max_iter=1), PROJECT)
        [w] = [w for w in r["warnings"] if w["code"] == "SHIELD_NOT_CONVERGED"]
        assert w["severity"] == "warning"
        assert "shield_max_iter" in w["message"]
        assert r["shield"]["converged"] is False
        assert 0 < r["governing_dx_mm"] < float("inf")

    def test_project_level_iteration_limit_and_part_override(self):
        project = {**PROJECT, "shield_max_iter": 1}
        assert "SHIELD_NOT_CONVERGED" in _warning_codes(
            process_part(_make_shield_part(), project))
        assert "SHIELD_NOT_CONVERGED" not in _warning_codes(
            process_part(_make_shield_part(shield_max_iter=50), project))

    @pytest.mark.parametrize("bad", [0, -3, 1.5, "10", True])
    def test_invalid_iteration_limit_raises(self, bad):
        with pytest.raises(ValueError, match="shield_max_iter"):
            process_part(_make_shield_part(shield_max_iter=bad), PROJECT)

    def test_unbounded_governing_size_warns(self):
        r = process_part(_make_shield_part(t_exh_K=353.15), PROJECT)
        assert r["governing_dx_mm"] == float("inf")
        assert "NO_FINITE_SIZE" in _warning_codes(r)


# --- Transient: the upper bound is the candidate, conflicts are warnings ---

# A bound as the TRANSIENT_CONFLICT message prints it, e.g. "dt to <= 20.24".
_PRINTED_BOUND = re.compile(r"\b(dt|safety_factor) to [<>]= (\d[0-9.eE+-]*)")


class TestTransientInBatch:

    def test_explicit_without_tau_adds_no_transient_candidate(self):
        r = process_part(_make_structural_part(), {**PROJECT, "dt": 1.0})
        assert r["transient"]["max_dx_mm"] == float("inf")
        sources = [c["source"] for c in r["all_constraints"]]
        assert "transient" not in sources
        assert "TRANSIENT_CONFLICT" not in _warning_codes(r)

    def test_drive_cycle_bound_is_the_transient_candidate(self):
        r = process_part(_make_structural_part(), PROJECT_WITH_TRANSIENT)
        [c] = [c for c in r["all_constraints"] if c["source"] == "transient"]
        assert c["dx_mm"] == r["transient"]["drive_cycle_max_dx_mm"]

    def test_conflict_is_a_warning_with_a_remedy_that_closes_it(self):
        project = {**PROJECT, "dt": 100.0}
        r = process_part(_make_structural_part(), project)
        [w] = [w for w in r["warnings"] if w["code"] == "TRANSIENT_CONFLICT"]
        assert w["severity"] == "warning"  # explicit: unstable, not inaccurate
        [remedy] = w["remedies"]
        assert remedy["parameter"] == "dt"
        dx_m = r["governing_dx_mm"] / 1000.0
        assert remedy["max_value"] == pytest.approx(
            0.5 * dx_m * dx_m / r["transient"]["alpha"])
        fixed = process_part(_make_structural_part(),
                             {**project, "dt": remedy["max_value"]})
        assert "TRANSIENT_CONFLICT" not in _warning_codes(fixed)
        assert fixed["governing_dx_mm"] == pytest.approx(r["governing_dx_mm"])

    def test_implicit_conflict_is_a_caution(self):
        project = {**PROJECT, "dt": 1.0e4, "fo_max": 5.0,
                   "safety_factor": 2.0}
        r = process_part(_make_structural_part(), project)
        [w] = [w for w in r["warnings"] if w["code"] == "TRANSIENT_CONFLICT"]
        assert w["severity"] == "caution"

    def test_empty_implicit_window_names_the_parameter_not_dt(self):
        project = {**PROJECT, "dt": 1.0, "transient_scheme": "implicit"}
        r = process_part(_make_structural_part(), project)
        [w] = [w for w in r["warnings"] if w["code"] == "TRANSIENT_CONFLICT"]
        assert w["remedies"][0]["parameter"] == "safety_factor"
        assert "safety_factor" in w["message"]

    @pytest.mark.parametrize("extra", [
        {"dt": 100.0},                                  # dt <= 20.2467 s
        {"dt": 100.0, "transient_scheme": "implicit"},  # safety_factor, then dt
    ], ids=["explicit", "implicit"])
    def test_the_printed_remedy_clears_the_warning(self, extra):
        """Applied as printed, the remedy clears the warning: a maximum is
        printed rounded down and a minimum rounded up, while the structured
        remedies keep the exact values. Rounded to nearest, the explicit
        case printed dt <= 20.25 s and the warning stayed."""
        project = {**PROJECT, **extra}
        r = process_part(_make_structural_part(), project)
        [w] = [w for w in r["warnings"] if w["code"] == "TRANSIENT_CONFLICT"]
        printed = {name: float(value.rstrip("."))
                   for name, value in _PRINTED_BOUND.findall(w["message"])}
        fixed = process_part(_make_structural_part(), {**project, **printed})
        assert "TRANSIENT_CONFLICT" not in _warning_codes(fixed), w["message"]
        exact = {x["parameter"]: x.get("max_value", x.get("min_value"))
                 for x in w["remedies"]}
        assert set(printed) == set(exact)
        if "dt" in exact:
            assert printed["dt"] <= exact["dt"]
        if "safety_factor" in exact:
            assert printed["safety_factor"] >= exact["safety_factor"]


# --- Input ranges, film temperature and correlation fallback in batch ---

def _exhaust_gas_part(**overrides):
    return _make_structural_part(
        material="cast_iron", surface="cast_iron_oxidised",
        component_class="exhaust", convection_zone="exhaust_internal",
        thickness_mm=6.0, t_surf_K=1173.15, t_fluid_K=1073.15, **overrides)


class TestInputRangesInBatch:

    def test_emissivity_above_one_raises(self):
        with pytest.raises(ValueError, match=r"epsilon must lie in \[0, 1\]"):
            process_part(_make_structural_part(epsilon=1.7), PROJECT)

    def test_batch_reports_the_guard_as_the_error(self):
        results = process_batch([_make_structural_part(epsilon=1.7)], PROJECT)
        assert results[0]["error"] == "epsilon must lie in [0, 1], got 1.7"

    def test_film_temperature_warning_names_film_and_evaluation(self):
        r = process_part(_exhaust_gas_part(), PROJECT)
        [w] = [w for w in r["warnings"] if w["code"] == "FILM_TEMP_OUT_OF_RANGE"]
        assert w["severity"] == "caution"
        assert "film temperature 1123 K" in w["message"]
        assert "evaluated at 700 K" in w["message"]
        assert r["h_estimation"]["details"]["film_in_range"] is False

    def test_no_film_warning_in_range(self):
        r = process_part(_make_structural_part(), PROJECT)
        assert "FILM_TEMP_OUT_OF_RANGE" not in _warning_codes(r)
        assert r["h_estimation"]["details"]["film_in_range"] is True

    def test_film_warning_covers_the_boundary_layer(self):
        part = _make_structural_part(
            h_override=20.0, t_surf_K=1573.15, bl_regime="external_forced",
            velocity_ms=10.0, char_length_mm=300.0,
        )
        r = process_part(part, PROJECT)
        [w] = [w for w in r["warnings"] if w["code"] == "FILM_TEMP_OUT_OF_RANGE"]
        assert "the boundary-layer sizes" in w["message"]
        assert r["boundary_layer"]["film_in_range"] is False

    def test_failed_correlation_falls_back_loudly(self):
        r = process_part(_make_structural_part(char_length_mm=0.0), PROJECT)
        assert r["h_estimation"]["method"] == "static_lookup"
        [w] = [w for w in r["warnings"] if w["code"] == "H_CORRELATION_FALLBACK"]
        assert "char_length must be > 0" in w["message"]
        assert r["h_used"] == 15.0  # engine_beside h_high, via the alias


# --- The shield path: radiation constraint, gap inputs, defaults reported ---

class TestShieldPath:

    def test_single_layer_radiation_from_the_exhaust_face(self):
        r = process_part(_make_shield_part(), PROJECT)
        s = r["shield"]
        expected = RadiationMeshCalculator.max_mesh_size(
            t_local=s["t_shield_K"], emissivity=0.40,
            allowable_flux_error=500.0,          # shield class default
            spatial_gradient=s["q_boundary"] / 54.0,  # steel_mild k
        )
        assert r["radiation"]["max_dx_mm"] == pytest.approx(expected["max_dx_mm"])
        assert r["governing_dx_mm"] == min(c["dx_mm"] for c in r["all_constraints"])

    def test_two_layer_radiation_from_layer_one(self):
        r = process_part(_make_multilayer_part(f12=0.9), PROJECT)
        s = r["shield"]
        expected = RadiationMeshCalculator.max_mesh_size(
            t_local=s["t1_K"], emissivity=0.40, allowable_flux_error=500.0,
            spatial_gradient=s["q_layer1"] / 54.0,
        )
        assert r["radiation"]["max_dx_mm"] == pytest.approx(expected["max_dx_mm"])

    def test_f12_reaches_the_solver_and_is_reported(self):
        r = process_part(_make_multilayer_part(f12=0.85), PROJECT)
        assert r["f12_used"] == 0.85
        assert r["shield"]["eps_eff"] == pytest.approx(
            1.0 / (1.0 / 0.4 + 1.0 / 0.4 - 2.0 + 1.0 / 0.85))
        assert "SHIELD_INPUT_DEFAULTED" not in _warning_codes(r)

    def test_defaulted_f12_and_h_gap_are_reported(self):
        part = _make_multilayer_part()
        del part["h_gap"]
        r = process_part(part, PROJECT)
        keys = sorted(w["key"] for w in r["warnings"]
                      if w["code"] == "SHIELD_INPUT_DEFAULTED")
        assert keys == ["f12", "h_gap"]
        assert r["f12_used"] == 1.0
        assert r["h_used"]["h_gap"] == 15.0

    def test_defaulted_exhaust_temperature_is_reported(self):
        project = {"t_fluid_K": 353.15, "t_surr_K": 353.15}
        r = process_part(_make_shield_part(), project)
        [w] = [w for w in r["warnings"] if w["code"] == "SHIELD_INPUT_DEFAULTED"]
        assert w["key"] == "t_exh_K"
        assert "1073.15" in w["message"]
        assert "SHIELD_INPUT_DEFAULTED" not in _warning_codes(
            process_part(_make_shield_part(), PROJECT))

    def test_f12_outside_the_unit_interval_raises(self):
        with pytest.raises(ValueError, match=r"f12 must lie in \[0, 1\]"):
            process_part(_make_multilayer_part(f12=1.5), PROJECT)
