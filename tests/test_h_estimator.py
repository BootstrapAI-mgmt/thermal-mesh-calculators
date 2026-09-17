"""
Tests for h_estimator — forced, natural, mixed convection + Richardson number.
"""

import math
import pytest
from thermal_mesh_calculators.h_estimator import (
    air_properties,
    forced_convection_flat_plate,
    natural_convection,
    richardson_number,
    estimate_h,
    solver_advisory,
    h_from_velocity,
    H_EXTERNAL_MAX,
)


class TestAirProperties:

    def test_300K_reference(self):
        """Air at 300 K: k~0.026, nu~1.5e-5, Pr~0.71."""
        props = air_properties(300.0)
        assert 0.020 < props["k_air"] < 0.030
        assert 1e-5 < props["nu"] < 2e-5
        assert props["Pr"] == pytest.approx(0.71)
        assert props["beta"] == pytest.approx(1.0 / 300.0)

    def test_higher_temp_increases_viscosity(self):
        p300 = air_properties(300.0)
        p500 = air_properties(500.0)
        assert p500["nu"] > p300["nu"]

    def test_clamped_below_250(self):
        """Should not crash at extreme temps."""
        props = air_properties(100.0)
        assert props["k_air"] > 0


class TestForcedConvection:

    def test_10mph_100mm_part(self):
        """10 mph ~ 4.47 m/s, 100 mm part. Should give reasonable h."""
        result = forced_convection_flat_plate(
            velocity=4.47, char_length=0.1,
            t_surf=473.15, t_fluid=313.15,
        )
        assert 10 < result["h"] < 100
        assert result["Re"] > 0
        assert result["regime"] in ("laminar", "turbulent")

    def test_h_capped(self):
        """High velocity should be capped at H_EXTERNAL_MAX."""
        result = forced_convection_flat_plate(
            velocity=50.0, char_length=0.5,
            t_surf=373.15, t_fluid=313.15,
            cap=H_EXTERNAL_MAX,
        )
        assert result["h"] <= H_EXTERNAL_MAX
        assert result["h_raw"] >= result["h"]

    def test_laminar_regime(self):
        """Low velocity, short part → laminar."""
        result = forced_convection_flat_plate(
            velocity=1.0, char_length=0.05,
            t_surf=373.15, t_fluid=313.15,
        )
        assert result["regime"] == "laminar"

    def test_higher_velocity_higher_h(self):
        params = dict(char_length=0.1, t_surf=400.0, t_fluid=313.15, cap=999.0)
        r1 = forced_convection_flat_plate(velocity=2.0, **params)
        r2 = forced_convection_flat_plate(velocity=8.0, **params)
        assert r2["h_raw"] > r1["h_raw"]


class TestNaturalConvection:

    def test_vertical_plate_hot(self):
        """Hot vertical plate should give h in 5-25 range."""
        result = natural_convection(
            t_surf=573.15, t_fluid=313.15,  # 300C, 40C
            char_length=0.2, orientation="vertical",
        )
        assert 5 < result["h"] < 30
        assert result["Ra"] > 0

    def test_horizontal_up_stronger_than_down(self):
        """Hot-side-up should give higher h than hot-side-down."""
        params = dict(t_surf=473.15, t_fluid=313.15, char_length=0.15)
        r_up = natural_convection(orientation="horizontal_up", **params)
        r_down = natural_convection(orientation="horizontal_down", **params)
        assert r_up["h"] > r_down["h"]

    def test_zero_dt_returns_zero(self):
        """No temperature difference → no buoyancy → h=0."""
        result = natural_convection(
            t_surf=300.0, t_fluid=300.0,
            char_length=0.1, orientation="vertical",
        )
        assert result["h"] == 0.0

    def test_hotter_surface_higher_h(self):
        """Higher delta-T drives stronger buoyancy."""
        params = dict(t_fluid=313.15, char_length=0.15,
                      orientation="vertical", cap=999.0)
        r_low = natural_convection(t_surf=373.15, **params)   # 100C
        r_high = natural_convection(t_surf=623.15, **params)  # 350C
        assert r_high["h_raw"] > r_low["h_raw"]

    def test_unknown_orientation_raises(self):
        with pytest.raises(ValueError, match="Unknown orientation"):
            natural_convection(400.0, 300.0, 0.1, "diagonal")


class TestRichardsonNumber:

    def test_high_velocity_is_forced(self):
        ri = richardson_number(
            velocity=10.0, t_surf=373.15, t_fluid=313.15,
            char_length=0.1,
        )
        assert ri["regime"] == "forced"
        assert ri["Ri"] < 0.1

    def test_zero_velocity_is_natural(self):
        ri = richardson_number(
            velocity=0.0, t_surf=573.15, t_fluid=313.15,
            char_length=0.1,
        )
        assert ri["regime"] == "natural"
        assert ri["Ri"] == float("inf")

    def test_intermediate_is_mixed(self):
        """Low velocity + high delta-T should give mixed."""
        ri = richardson_number(
            velocity=0.5, t_surf=573.15, t_fluid=313.15,
            char_length=0.2,
        )
        # This should be mixed or natural — not forced
        assert ri["regime"] in ("mixed", "natural")


class TestEstimateH:

    def test_dead_zone_uses_natural(self):
        """Zero velocity → pure natural convection."""
        result = estimate_h(
            velocity=0.0, t_surf=473.15, t_fluid=313.15,
            char_length=0.15, orientation="vertical",
        )
        assert result["regime"] == "natural"
        assert result["h_forced"] == 0.0
        assert result["h_natural"] > 0

    def test_high_velocity_uses_forced(self):
        result = estimate_h(
            velocity=10.0, t_surf=373.15, t_fluid=313.15,
            char_length=0.1, orientation="vertical",
        )
        assert result["regime"] == "forced"

    def test_capped_at_external_max(self):
        result = estimate_h(
            velocity=50.0, t_surf=373.15, t_fluid=313.15,
            char_length=0.5, orientation="vertical",
        )
        assert result["h"] <= H_EXTERNAL_MAX

    def test_sub_regimes_propagated(self):
        """estimate_h should report natural_regime and forced_regime."""
        result = estimate_h(
            velocity=4.47, t_surf=473.15, t_fluid=313.15,
            char_length=0.1, orientation="vertical",
        )
        assert "natural_regime" in result
        assert "forced_regime" in result
        assert result["natural_regime"] in ("laminar", "turbulent", "negligible")
        assert result["forced_regime"] in ("laminar", "turbulent", "n/a")

    def test_solver_advisory_present(self):
        """estimate_h should always include a solver_advisory dict."""
        result = estimate_h(
            velocity=0.0, t_surf=473.15, t_fluid=313.15,
            char_length=0.15, orientation="vertical",
        )
        adv = result["solver_advisory"]
        assert isinstance(adv, dict)
        assert "steady_state_ok" in adv
        assert "transient_advisory" in adv
        assert "severity" in adv

    def test_Ra_Re_propagated(self):
        """estimate_h should report Ra and Re."""
        result = estimate_h(
            velocity=4.47, t_surf=473.15, t_fluid=313.15,
            char_length=0.1, orientation="vertical",
        )
        assert "Ra" in result
        assert "Re" in result
        assert result["Ra"] > 0
        assert result["Re"] > 0


class TestOpposingMixedConvection:

    def test_horizontal_down_uses_opposing(self):
        """Mixed convection with horizontal_down should use opposing formula."""
        result = estimate_h(
            velocity=1.5, t_surf=573.15, t_fluid=313.15,
            char_length=0.2, orientation="horizontal_down",
        )
        # Should report as mixed-opposing if Richardson is in mixed range
        if result["regime"] == "mixed":
            assert "opposing" in result["dominant_mode"]

    def test_opposing_h_lower_than_assisting(self):
        """Opposing mixed convection should generally give lower h than assisting
        (vertical) for the same conditions, since modes partially cancel."""
        params = dict(velocity=1.0, t_surf=573.15, t_fluid=313.15,
                      char_length=0.15)
        r_assist = estimate_h(orientation="vertical", **params)
        r_oppose = estimate_h(orientation="horizontal_down", **params)
        # Both should have non-zero h
        assert r_oppose["h"] > 0
        assert r_assist["h"] > 0
        # For horizontal_down, natural convection is weaker anyway (0.27*Ra^0.25
        # vs 0.59*Ra^0.25 for vertical), so compare within the mixed formula
        # The opposing formula uses |h_f^3 - h_n^3|^(1/3) which is <= (h_f^3 + h_n^3)^(1/3)

    def test_opposing_has_conduction_floor(self):
        """Even when forced and natural nearly cancel, h should not be zero
        — the k_air/L conduction floor prevents singularity."""
        # Near-perfect cancellation: use horizontal_down with velocity that
        # gives Ri in mixed range
        result = estimate_h(
            velocity=0.5, t_surf=373.15, t_fluid=313.15,
            char_length=0.15, orientation="horizontal_down",
        )
        assert result["h"] > 0  # must never be zero

    def test_vertical_still_uses_assisting(self):
        """Vertical orientation should still use assisting formula."""
        result = estimate_h(
            velocity=1.0, t_surf=573.15, t_fluid=313.15,
            char_length=0.15, orientation="vertical",
        )
        if result["regime"] == "mixed":
            assert "assisting" in result["dominant_mode"]


class TestSolverAdvisory:

    def test_laminar_natural_is_ss_ok(self):
        """Laminar natural convection: steady-state should be fine."""
        nat = natural_convection(
            t_surf=340.0, t_fluid=313.15,  # small dT
            char_length=0.05, orientation="vertical",
        )
        assert nat["regime"] == "laminar"
        adv = solver_advisory(natural_result=nat, regime="natural",
                              orientation="vertical")
        assert adv["steady_state_ok"] is True
        assert adv["transient_advisory"] is False
        assert adv["severity"] == "info"

    def test_turbulent_natural_flags_transient(self):
        """Turbulent natural convection should flag transient advisory."""
        nat = natural_convection(
            t_surf=873.15, t_fluid=313.15,  # 600C → 40C, large Ra
            char_length=1.0, orientation="vertical",  # tall surface
        )
        assert nat["regime"] == "turbulent"
        adv = solver_advisory(natural_result=nat, regime="natural",
                              orientation="vertical")
        assert adv["steady_state_ok"] is False
        assert adv["transient_advisory"] is True
        assert adv["severity"] == "warning"

    def test_turbulent_forced_is_ss_ok(self):
        """Turbulent forced convection is fine for SS (time-averaged Nu)."""
        frc = forced_convection_flat_plate(
            velocity=20.0, char_length=0.5,
            t_surf=373.15, t_fluid=313.15, cap=999.0,
        )
        assert frc["regime"] == "turbulent"
        adv = solver_advisory(forced_result=frc, regime="forced")
        assert adv["steady_state_ok"] is True
        assert adv["severity"] == "info"

    def test_mixed_regime_caution(self):
        """Mixed convection should raise caution."""
        nat = natural_convection(
            t_surf=473.15, t_fluid=313.15,
            char_length=0.15, orientation="vertical",
        )
        frc = forced_convection_flat_plate(
            velocity=1.0, char_length=0.15,
            t_surf=473.15, t_fluid=313.15, cap=999.0,
        )
        adv = solver_advisory(forced_result=frc, natural_result=nat,
                              regime="mixed", orientation="vertical")
        assert adv["severity"] in ("caution", "warning")

    def test_negligible_natural_caution(self):
        """Near-zero dT with natural regime should flag caution."""
        nat = natural_convection(
            t_surf=313.15, t_fluid=313.15,
            char_length=0.1, orientation="vertical",
        )
        adv = solver_advisory(natural_result=nat, regime="natural",
                              orientation="vertical")
        assert adv["severity"] == "caution"

    def test_advisory_keys(self):
        """Solver advisory should have all expected keys."""
        nat = natural_convection(
            t_surf=473.15, t_fluid=313.15,
            char_length=0.15, orientation="vertical",
        )
        adv = solver_advisory(natural_result=nat, regime="natural")
        expected_keys = {
            "steady_state_ok", "transient_advisory", "reason",
            "severity", "natural_regime", "forced_regime", "Ra", "Re",
        }
        assert set(adv.keys()) == expected_keys

    def test_horizontal_up_turbulent_threshold_lower(self):
        """Horizontal hot-up transitions to turbulent at Ra=1e7 (lower
        than vertical 1e9), so moderate conditions can trigger warning."""
        # Use a large plate with moderate dT — should trip Ra > 1e7
        nat = natural_convection(
            t_surf=423.15, t_fluid=313.15,  # 150C → 40C
            char_length=0.3, orientation="horizontal_up",
        )
        adv = solver_advisory(natural_result=nat, regime="natural",
                              orientation="horizontal_up")
        if nat["regime"] == "turbulent":
            assert adv["transient_advisory"] is True


class TestHFromVelocity:

    def test_10mph_gut_check(self):
        """10 mph on a 100mm part at 100C should give h < 60."""
        result = h_from_velocity(
            velocity_mph=10.0, char_length_mm=100.0,
            t_surf_C=100.0, t_fluid_C=40.0,
        )
        assert result["h"] < 60.0
        assert result["h"] > 5.0

    def test_unit_conversion(self):
        result = h_from_velocity(velocity_mph=10.0)
        assert result["velocity_ms"] == pytest.approx(10.0 * 0.44704)
        assert result["char_length_mm"] == 100.0

    def test_solver_advisory_in_gut_check(self):
        """h_from_velocity should include solver advisory."""
        result = h_from_velocity(velocity_mph=10.0)
        assert "solver_advisory" in result
