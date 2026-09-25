"""
Fail-loud gate: one planted input per silent input path
========================================================

Each test below feeds the library an input that version 0.6.2 accepted
without a word, and asserts that the library now either raises or returns a
coded warning (the ``code`` field of a dict in ``process_part(...)
["warnings"]``).

Every test in this file fails against 0.6.2: the planted input is the
defect, and the assertion is that it gets reported.  The file therefore
imports only names that 0.6.2 already had, so that against the old code each
case fails on its own assertion rather than on an import.
"""

import math

import pytest

from thermal_mesh_calculators.batch import CLASS_DEFAULTS, MATERIALS, process_part
from thermal_mesh_calculators.conduction import BoundaryDrivenConductionCalculator
from thermal_mesh_calculators.h_estimator import air_properties, estimate_h
from thermal_mesh_calculators.shields import MultilayerShieldCalculator
from thermal_mesh_calculators.transient import TransientMeshCalculator
from thermal_mesh_calculators.zones import CONVECTION_ZONES


PROJECT = {"t_fluid_K": 353.15, "t_surr_K": 353.15, "t_exh_K": 1073.15}


def _structural(**overrides):
    part = {
        "part_id": "GATE-STR",
        "material": "steel_mild",
        "surface": "painted",
        "component_class": "structural",
        "convection_zone": "engine_bay_dead_zone",
        "thickness_mm": 3.0,
        "t_surf_K": 473.15,
    }
    part.update(overrides)
    return part


def _codes(result):
    return [w["code"] for w in result.get("warnings", [])]


class TestUnknownNamesRaise:
    """A misspelt name is an error that names the allowed set, never a default."""

    def test_misspelt_component_class_raises_and_names_the_allowed_set(self):
        # 0.6.2 sized "exhuast" as a structural part and said nothing.
        with pytest.raises(KeyError) as excinfo:
            process_part(_structural(component_class="exhuast"), PROJECT)
        assert isinstance(excinfo.value, ValueError)
        message = str(excinfo.value)
        assert "'exhuast'" in message
        assert ", ".join(sorted(CLASS_DEFAULTS)) in message

    def test_misspelt_zone_raises_even_when_h_is_overridden(self):
        # 0.6.2 never looked the zone up once h_override was given.
        with pytest.raises(KeyError) as excinfo:
            process_part(
                _structural(convection_zone="engine_bay_dead_zon",
                            h_override=10.0),
                PROJECT,
            )
        message = str(excinfo.value)
        assert "'engine_bay_dead_zon'" in message
        assert ", ".join(sorted(CONVECTION_ZONES)) in message


class TestDefaultedSurfaceIsReported:
    """A surface with no treatment and no emissivity still gets a fallback
    emissivity, but the result now says so."""

    def test_missing_surface_warns(self):
        part = _structural()
        del part["surface"]
        result = process_part(part, PROJECT)
        assert result["eps_used"] == pytest.approx(0.73)
        assert "SURFACE_DEFAULTED" in _codes(result)


class TestOneFluidTemperature:
    """h and q'' come from the same fluid temperature: the documented t_fluid_K."""

    @pytest.mark.parametrize("t_fluid", [300.0, 353.15, 450.0])
    def test_h_follows_the_fluid_temperature(self, t_fluid):
        # 0.6.2 estimated h at the zone's upper air temperature whatever
        # t_fluid_K said, and used t_fluid_K only for q''.
        project = dict(PROJECT, t_fluid_K=t_fluid, t_surr_K=t_fluid)
        result = process_part(_structural(), project)
        expected_h = estimate_h(
            velocity=0.0, t_surf=473.15, t_fluid=t_fluid,
            char_length=0.1, orientation="vertical",
        )["h"]
        assert result["h_used"] == pytest.approx(expected_h, rel=1e-12)
        assert result["conduction"]["q_conv"] == pytest.approx(
            result["h_used"] * (473.15 - t_fluid), rel=1e-12)


def _two_layer_shield(**overrides):
    part = {
        "part_id": "GATE-ML",
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
    part.update(overrides)
    return part


class TestShieldNonConvergence:
    """A shield solve that stops unconverged is reported, never an inf."""

    def test_unconverged_two_layer_solve_keeps_finite_sizes(self):
        # 0.6.2 returned no fluxes when out of iterations, so both layer
        # sizes came back as inf.
        result = MultilayerShieldCalculator.mesh_sizes(
            k_metal=45.0, max_dt=15.0,
            t_exh=1073.15, t_fluid=353.15, t_surr=353.15,
            h_in=30.0, h_out=30.0, h_gap=15.0,
            eps_in=0.4, eps_out=0.4, eps_g1=0.4, eps_g2=0.4,
            max_iter=1,
        )
        assert result["converged"] is False
        assert math.isfinite(result["layer1_max_dx_mm"])
        assert math.isfinite(result["layer2_max_dx_mm"])

    def test_unconverged_shield_solve_warns(self):
        result = process_part(_two_layer_shield(shield_max_iter=1), PROJECT)
        assert "SHIELD_NOT_CONVERGED" in _codes(result)
        assert math.isfinite(result["governing_dx_mm"])

    def test_an_infinite_governing_size_is_reported(self):
        # A shield with no heat source: every candidate size is unbounded.
        part = {
            "part_id": "GATE-COLD",
            "material": "steel_mild",
            "component_class": "shield",
            "convection_zone": "engine_bay_dead_zone",
            "thickness_mm": 0.8,
            "t_surf_K": 0,
            "surface_in": "aluminised",
            "surface_out": "aluminised",
            "t_exh_K": 353.15,
        }
        result = process_part(part, PROJECT)
        assert result["governing_dx_mm"] == float("inf")
        assert "NO_FINITE_SIZE" in _codes(result)


class TestExplicitTransientPath:
    """The explicit defaults (Fo <= 0.5, C = 1) give a feasible size or a
    remedy that changes the result, and the batch carries the conflict."""

    @pytest.mark.parametrize("dt", [0.01, 1.0, 100.0])
    def test_defaults_are_feasible_or_carry_a_working_remedy(self, dt):
        # 0.6.2 reported a conflict at every dt (the Fourier minimum was
        # sqrt(2) times the penetration bound whatever dt was) and advised
        # "reduce dt", which cannot change a dt-independent ratio.
        props = MATERIALS["steel_mild"]
        kw = dict(k=props["k"], rho=props["rho"], cp=props["cp"])
        result = TransientMeshCalculator.combined_transient_limits(dt=dt, **kw)
        if not result["feasible"]:
            changes = {r["parameter"]: r.get("max_value", r.get("min_value"))
                       for r in result["conflict"]["remedies"]}
            result = TransientMeshCalculator.combined_transient_limits(
                **dict(kw, dt=dt, **changes))
        assert result["feasible"]
        assert result["fourier_min_dx_mm"] <= result["recommended_dx_mm"]

    @pytest.mark.parametrize("dt", [0.01, 1.0, 100.0])
    def test_batch_never_caps_at_the_minimum_and_reports_conflicts(self, dt):
        # 0.6.2 used the stability minimum as a size cap (the governing
        # "transient" size at dt = 0.01 and 1 s) and never warned.
        result = process_part(_structural(), dict(PROJECT, dt=dt))
        trans = result["transient"]
        fo_min = trans["fourier_min_dx_mm"]
        if result["governing_constraint"] == "transient":
            assert result["governing_dx_mm"] > fo_min
        conflicts = [w for w in result["warnings"]
                     if w["code"] == "TRANSIENT_CONFLICT"]
        if fo_min <= result["governing_dx_mm"]:
            assert conflicts == []
        else:
            [warning] = conflicts
            [remedy] = warning["remedies"]
            again = process_part(
                _structural(), dict(PROJECT, dt=remedy["max_value"]))
            assert "TRANSIENT_CONFLICT" not in _codes(again)


_MANIFOLD = dict(h=150.0, t_surf=1073.15, t_fluid=353.15, t_surr=353.15,
                 max_dt=10.0)


class TestImpossibleInputsRaise:
    """A physically impossible input raises; it never becomes a size."""

    def test_negative_conductivity_raises(self):
        # 0.6.2: k = -45 with epsilon = 1.7 gave max_dx_mm = -1.920.
        with pytest.raises(ValueError, match="k must be > 0"):
            BoundaryDrivenConductionCalculator.max_mesh_size(
                k=-45.0, epsilon=0.85, **_MANIFOLD)

    def test_emissivity_above_one_raises(self):
        with pytest.raises(ValueError, match=r"epsilon must lie in \[0, 1\]"):
            BoundaryDrivenConductionCalculator.max_mesh_size(
                k=45.0, epsilon=1.7, **_MANIFOLD)


class TestFilmTemperatureOutOfRange:
    """Air properties outside their fitted 250-700 K range are reported."""

    def test_air_properties_say_they_were_clamped(self):
        # 0.6.2 returned k_air(1000 K) == k_air(700 K) and said nothing.
        props = air_properties(1000.0)
        assert props["k_air"] == air_properties(700.0)["k_air"]
        assert props["in_range"] is False

    def test_batch_warns_when_h_rests_on_extrapolated_air(self):
        part = _structural(
            material="cast_iron", surface="cast_iron_oxidised",
            component_class="exhaust", convection_zone="exhaust_internal",
            thickness_mm=6.0, t_surf_K=1173.15, t_fluid_K=1073.15,
        )
        result = process_part(part, PROJECT)
        assert "FILM_TEMP_OUT_OF_RANGE" in _codes(result)


def _stainless_shield(**overrides):
    part = {
        "part_id": "GATE-SH",
        "material": "steel_stainless_409",
        "component_class": "shield",
        "convection_zone": "exhaust_above",
        "thickness_mm": 1.0,
        "t_surf_K": 0,
        "surface_in": "aluminised",
        "surface_out": "aluminised",
    }
    part.update(overrides)
    return part


class TestShieldPathComplete:
    """Shields get a radiation constraint; gap inputs are part keys and a
    defaulted one is reported."""

    def test_single_layer_shield_gets_a_radiation_candidate(self):
        # 0.6.2 gave the shield classes no radiation constraint at all.
        result = process_part(_stainless_shield(), PROJECT)
        assert result["radiation"] is not None
        sources = [c["source"] for c in result["all_constraints"]]
        assert "radiation" in sources

    def test_two_layer_shield_reports_its_defaulted_gap_inputs(self):
        result = process_part(
            _stainless_shield(component_class="multilayer_shield"), PROJECT)
        assert result["radiation"] is not None
        defaulted = sorted(
            w.get("key") for w in result["warnings"]
            if w["code"] in ("SURFACE_DEFAULTED", "SHIELD_INPUT_DEFAULTED"))
        assert defaulted == ["f12", "h_gap", "surface_g1", "surface_g2"]

    def test_f12_is_a_part_key(self):
        # 0.6.2 never passed f12 to the solver: the gap was always treated as
        # infinite parallel plates.
        result = process_part(
            _stainless_shield(component_class="multilayer_shield", f12=0.55,
                              h_gap=12.0, surface_g1="aluminised",
                              surface_g2="aluminised"),
            PROJECT,
        )
        expected = 1.0 / (1.0 / 0.4 + 1.0 / 0.4 - 2.0 + 1.0 / 0.55)
        assert result["shield"]["eps_eff"] == pytest.approx(expected)
        assert _codes(result) == []
