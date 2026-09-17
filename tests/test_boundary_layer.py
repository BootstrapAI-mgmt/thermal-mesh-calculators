"""
Tests for boundary_layer.py — Aerodynamic Boundary Layer Mesh Calculator
========================================================================
"""
import math
import unittest

from thermal_mesh_calculators.boundary_layer import (
    BoundaryLayerCalculator,
    RE_X_MIN,
    RE_X_CRIT,
    Y1_ABSOLUTE_MIN_M,
    _WALL_GRADIENT_PRISM_THRESHOLD,
    ER_V_STABLE,
    ER_V_MARGINAL,
)
from thermal_mesh_calculators.h_estimator import air_properties


class TestSkinFriction(unittest.TestCase):
    """Unit tests for the _skin_friction static method."""

    def test_laminar_regime(self):
        """Re_x = 1e4 should give laminar C_f = 0.664 * Re^-0.5."""
        cf, regime = BoundaryLayerCalculator._skin_friction(1e4)
        self.assertEqual(regime, "laminar")
        self.assertAlmostEqual(cf, 0.664 * 1e4 ** (-0.5), places=6)

    def test_turbulent_regime(self):
        """Re_x = 1e6 should give turbulent C_f = 0.058 * Re^-0.2."""
        cf, regime = BoundaryLayerCalculator._skin_friction(1e6)
        self.assertEqual(regime, "turbulent")
        self.assertAlmostEqual(cf, 0.058 * 1e6 ** (-0.2), places=6)

    def test_leading_edge_clamp(self):
        """Re_x below RE_X_MIN is clamped — should not produce NaN/inf."""
        cf, regime = BoundaryLayerCalculator._skin_friction(0.0)
        self.assertTrue(math.isfinite(cf))
        self.assertGreater(cf, 0)

    def test_transition_boundary(self):
        """Re_x just below 5e5 is laminar, just above is turbulent."""
        _, regime_low = BoundaryLayerCalculator._skin_friction(4.99e5)
        _, regime_high = BoundaryLayerCalculator._skin_friction(5.01e5)
        self.assertEqual(regime_low, "laminar")
        self.assertEqual(regime_high, "turbulent")


class TestBoundaryLayerThickness(unittest.TestCase):
    """Unit tests for _boundary_layer_thickness."""

    def test_laminar_flat_plate(self):
        """δ = 5.0 * x * Re_x^-0.5 for laminar."""
        re_x = 1e4
        x = 0.1  # m
        delta = BoundaryLayerCalculator._boundary_layer_thickness(re_x, x)
        expected = 5.0 * x * re_x ** (-0.5)
        self.assertAlmostEqual(delta, expected, places=6)

    def test_turbulent_flat_plate(self):
        """δ = 0.37 * x * Re_x^-0.2 for turbulent."""
        re_x = 1e6
        x = 0.5  # m
        delta = BoundaryLayerCalculator._boundary_layer_thickness(re_x, x)
        expected = 0.37 * x * re_x ** (-0.2)
        self.assertAlmostEqual(delta, expected, places=6)

    def test_thickness_increases_with_x(self):
        """Boundary layer should be thicker further downstream."""
        re1 = 22 * 0.1 / 2.1e-5  # ~1e5
        re2 = 22 * 0.5 / 2.1e-5  # ~5.2e5
        d1 = BoundaryLayerCalculator._boundary_layer_thickness(re1, 0.1)
        d2 = BoundaryLayerCalculator._boundary_layer_thickness(re2, 0.5)
        self.assertGreater(d2, d1)


class TestBuoyancyVelocity(unittest.TestCase):
    """Unit tests for _buoyancy_velocity."""

    def test_basic_buoyancy(self):
        """V_buoy = sqrt(g * beta * dT * L) for known inputs."""
        beta = 1.0 / 500.0  # 1/K at ~500 K film temp
        dT = 200.0  # K
        L = 0.1  # m gap height
        v = BoundaryLayerCalculator._buoyancy_velocity(beta, dT, L)
        expected = math.sqrt(9.81 * beta * 200.0 * 0.1)
        self.assertAlmostEqual(v, expected, places=4)

    def test_zero_dt_gives_zero(self):
        """No temperature difference → no buoyancy velocity."""
        v = BoundaryLayerCalculator._buoyancy_velocity(1.0 / 300.0, 0.0, 1.0)
        self.assertAlmostEqual(v, 0.0)


class TestInflationLayerCount(unittest.TestCase):
    """Unit tests for _inflation_layer_count."""

    def test_basic_count(self):
        """Known y1, delta, growth should give expected n."""
        y1 = 0.001  # 1 mm
        delta = 0.02  # 20 mm
        r = 1.2
        n = BoundaryLayerCalculator._inflation_layer_count(y1, delta, r)
        # Manual: n = ceil(ln(1 + 20*0.2) / ln(1.2)) = ceil(ln(5)/ln(1.2))
        #        = ceil(1.609/0.182) = ceil(8.83) = 9
        expected = math.ceil(
            math.log(1.0 + (delta / y1) * (r - 1.0)) / math.log(r)
        )
        self.assertEqual(n, expected)

    def test_minimum_three_layers(self):
        """Even for tiny BL, should return at least 3 layers."""
        n = BoundaryLayerCalculator._inflation_layer_count(0.01, 0.01, 1.2)
        self.assertGreaterEqual(n, 3)

    def test_zero_inputs_returns_three(self):
        """Degenerate inputs should return minimum 3."""
        self.assertEqual(
            BoundaryLayerCalculator._inflation_layer_count(0, 0.01, 1.2), 3
        )
        self.assertEqual(
            BoundaryLayerCalculator._inflation_layer_count(0.01, 0, 1.2), 3
        )


class TestEstimateMeshExternalForced(unittest.TestCase):
    """Integration tests for estimate_mesh in external_forced regime."""

    def test_underbody_80kph(self):
        """Typical underbody: 22 m/s, x=0.5m, 80°C air, y+=30."""
        result = BoundaryLayerCalculator.estimate_mesh(
            U=22.0, x=0.5, t_fluid=353.15,
            y_plus_target=30.0, growth_ratio=1.2,
            regime="external_forced",
        )
        # Should produce finite, positive values
        self.assertGreater(result["y1_mm"], 0)
        self.assertTrue(math.isfinite(result["y1_mm"]))
        self.assertGreater(result["n_layers"], 2)
        self.assertGreater(result["delta_mm"], 0)
        self.assertGreater(result["max_dx_surface_mm"], 0)
        self.assertEqual(result["regime"], "external_forced")

    def test_y_plus_1_finer_than_30(self):
        """Wall-resolved y+=1 should give much smaller y1 than y+=30."""
        r1 = BoundaryLayerCalculator.estimate_mesh(
            U=22.0, x=0.5, t_fluid=353.15,
            y_plus_target=1.0, regime="external_forced",
        )
        r30 = BoundaryLayerCalculator.estimate_mesh(
            U=22.0, x=0.5, t_fluid=353.15,
            y_plus_target=30.0, regime="external_forced",
        )
        self.assertLess(r1["y1_mm"], r30["y1_mm"])
        # y1 should scale roughly linearly with y+
        ratio = r30["y1_mm"] / r1["y1_mm"]
        self.assertAlmostEqual(ratio, 30.0, delta=1.0)

    def test_higher_velocity_thinner_bl(self):
        """Higher velocity should produce thinner y1 (more shear)."""
        r_slow = BoundaryLayerCalculator.estimate_mesh(
            U=10.0, x=0.5, t_fluid=353.15, regime="external_forced",
        )
        r_fast = BoundaryLayerCalculator.estimate_mesh(
            U=40.0, x=0.5, t_fluid=353.15, regime="external_forced",
        )
        self.assertLess(r_fast["y1_mm"], r_slow["y1_mm"])

    def test_surface_constraint_is_bl_or_transition(self):
        """Surface constraint label should be one of the expected values."""
        result = BoundaryLayerCalculator.estimate_mesh(
            U=22.0, x=0.5, t_fluid=353.15, regime="external_forced",
        )
        self.assertIn(
            result["surface_constraint"],
            ("prism_transition", "bl_fraction"),
        )

    def test_wall_function_surface_mesh_order_of_magnitude(self):
        """At y+=30, 22 m/s, surface mesh should be on the order of mm."""
        result = BoundaryLayerCalculator.estimate_mesh(
            U=22.0, x=0.5, t_fluid=353.15,
            y_plus_target=30.0, regime="external_forced",
        )
        # Expect surface mesh in the 1–20 mm range for wall functions
        self.assertGreater(result["max_dx_surface_mm"], 0.5)
        self.assertLess(result["max_dx_surface_mm"], 50.0)


class TestEstimateMeshMixedUnknown(unittest.TestCase):
    """Integration tests for estimate_mesh in mixed_unknown regime."""

    def test_buoyancy_dominated_dead_zone(self):
        """Dead zone (U=0) with hot surface should use buoyancy velocity."""
        result = BoundaryLayerCalculator.estimate_mesh(
            U=0.0, x=0.1, t_fluid=353.15,
            y_plus_target=50.0, growth_ratio=1.2,
            regime="mixed_unknown",
            delta_t_buoyancy=200.0,
            ar_max_prism=3.0,
        )
        self.assertGreater(result["V_buoyancy"], 0)
        self.assertEqual(result["U_effective"], max(0.1, result["V_buoyancy"]))
        self.assertEqual(result["regime"], "mixed_unknown")

    def test_buoyancy_from_t_surf(self):
        """If delta_t_buoyancy=0 but t_surf given, should auto-compute."""
        result = BoundaryLayerCalculator.estimate_mesh(
            U=0.0, x=0.1, t_fluid=353.15,
            regime="mixed_unknown",
            t_surf=553.15,  # 200 K above fluid
        )
        self.assertGreater(result["V_buoyancy"], 0)

    def test_mixed_uses_bl_fraction_only(self):
        """Mixed/unknown should use bl_fraction constraint only (no AR coupling).
        This was changed from the original AR-based constraint based on the
        consensus that prism layers are not needed in low-Re mixed zones."""
        r_mixed = BoundaryLayerCalculator.estimate_mesh(
            U=5.0, x=0.3, t_fluid=353.15,
            y_plus_target=50.0, regime="mixed_unknown",
            ar_max_prism=3.0,
        )
        self.assertEqual(r_mixed["surface_constraint"], "bl_fraction")
        # dx should equal bl_fraction * delta
        expected_dx = 0.3 * r_mixed["delta_mm"]
        self.assertAlmostEqual(
            r_mixed["max_dx_surface_mm"], expected_dx, places=3
        )

    def test_exhaust_cargo_bed_scenario(self):
        """Cargo bed above exhaust: U≈0, dT=200K, gap=0.1m, y+=50."""
        result = BoundaryLayerCalculator.estimate_mesh(
            U=0.0, x=0.1, t_fluid=353.15,
            y_plus_target=50.0, growth_ratio=1.2,
            regime="mixed_unknown",
            delta_t_buoyancy=200.0,
            ar_max_prism=3.0,
        )
        # Buoyancy velocity should be ~0.6–0.8 m/s
        # (exact value depends on film-temperature β)
        self.assertAlmostEqual(result["V_buoyancy"], 0.75, delta=0.15)
        # Surface mesh should be small but not absurdly so
        self.assertGreater(result["max_dx_surface_mm"], 0.1)
        self.assertLess(result["max_dx_surface_mm"], 20.0)


class TestLeadingEdgeClamp(unittest.TestCase):
    """Verify the leading-edge singularity protection."""

    def test_very_small_x(self):
        """x = 0.001 m (1 mm) should not produce zero or inf y1."""
        result = BoundaryLayerCalculator.estimate_mesh(
            U=22.0, x=0.001, t_fluid=353.15,
            regime="external_forced",
        )
        self.assertGreater(result["y1_mm"], 0)
        self.assertTrue(math.isfinite(result["y1_mm"]))
        self.assertGreater(result["max_dx_surface_mm"], 0)
        self.assertTrue(math.isfinite(result["max_dx_surface_mm"]))

    def test_zero_x(self):
        """x = 0 should still produce valid results (Re clamped)."""
        result = BoundaryLayerCalculator.estimate_mesh(
            U=22.0, x=0.0, t_fluid=353.15,
            regime="external_forced",
        )
        self.assertGreater(result["y1_mm"], 0)
        self.assertTrue(math.isfinite(result["y1_mm"]))
        self.assertGreaterEqual(result["Re_x"], RE_X_MIN)

    def test_y1_above_absolute_minimum(self):
        """y1 should never go below Y1_ABSOLUTE_MIN_M."""
        result = BoundaryLayerCalculator.estimate_mesh(
            U=100.0, x=5.0, t_fluid=300.0,
            y_plus_target=0.1,  # extremely aggressive
            regime="external_forced",
        )
        self.assertGreaterEqual(result["y1_m"], Y1_ABSOLUTE_MIN_M)


class TestInflationStackGeometry(unittest.TestCase):
    """Verify geometric consistency of the inflation stack."""

    def test_total_height_formula(self):
        """Total prism height should equal geometric series sum."""
        result = BoundaryLayerCalculator.estimate_mesh(
            U=22.0, x=0.5, t_fluid=353.15,
            y_plus_target=30.0, growth_ratio=1.3,
            regime="external_forced",
        )
        y1 = result["y1_m"]
        r = result["growth_ratio"]
        n = result["n_layers"]
        expected_total = y1 * (r ** n - 1.0) / (r - 1.0) * 1000.0
        self.assertAlmostEqual(
            result["total_prism_mm"], expected_total, places=3
        )

    def test_y_last_formula(self):
        """y_last should equal y1 * r^(n-1)."""
        result = BoundaryLayerCalculator.estimate_mesh(
            U=22.0, x=0.5, t_fluid=353.15,
            y_plus_target=30.0, growth_ratio=1.2,
            regime="external_forced",
        )
        y1 = result["y1_mm"]
        r = result["growth_ratio"]
        n = result["n_layers"]
        expected_last = y1 * r ** (n - 1)
        self.assertAlmostEqual(
            result["y_last_mm"], expected_last, places=4
        )

    def test_total_prism_less_than_delta(self):
        """Total prism stack should not greatly exceed BL thickness.
        It should be approximately equal to delta (by design)."""
        result = BoundaryLayerCalculator.estimate_mesh(
            U=22.0, x=0.5, t_fluid=353.15,
            y_plus_target=30.0, growth_ratio=1.2,
            regime="external_forced",
        )
        # Total prism ≈ delta (may overshoot slightly due to ceiling)
        ratio = result["total_prism_mm"] / result["delta_mm"]
        self.assertLess(ratio, 2.0)  # should not be more than 2x delta


class TestOutputCompleteness(unittest.TestCase):
    """Ensure all expected output keys are present."""

    def test_all_keys_present(self):
        expected_keys = {
            "y1_mm", "y1_m", "n_layers", "delta_mm", "y_last_mm",
            "total_prism_mm", "max_dx_surface_mm", "surface_constraint",
            "prisms_recommended", "prism_reason",
            "wall_gradient_Km", "volume_expansion_ratio", "er_v_status",
            "U_effective", "V_buoyancy", "Re_x", "Cf", "u_tau",
            "regime", "bl_regime", "y_plus_target", "growth_ratio",
        }
        result = BoundaryLayerCalculator.estimate_mesh(
            U=22.0, x=0.5, t_fluid=353.15, regime="external_forced",
        )
        self.assertTrue(expected_keys.issubset(set(result.keys())))


class TestPrismRecommendation(unittest.TestCase):
    """Tests for the prisms_recommended flag and supporting metrics."""

    def test_external_forced_always_recommends_prisms(self):
        """External forced regime should always recommend prism layers."""
        result = BoundaryLayerCalculator.estimate_mesh(
            U=22.0, x=0.5, t_fluid=353.15,
            regime="external_forced", t_surf=553.15,
        )
        self.assertTrue(result["prisms_recommended"])
        self.assertIn("Recommended", result["prism_reason"])

    def test_mixed_low_velocity_no_prisms(self):
        """Mixed/unknown at low velocity should NOT recommend prisms."""
        result = BoundaryLayerCalculator.estimate_mesh(
            U=0.0, x=0.1, t_fluid=353.15,
            regime="mixed_unknown",
            delta_t_buoyancy=200.0, t_surf=553.15,
        )
        self.assertFalse(result["prisms_recommended"])
        self.assertIn("Not recommended", result["prism_reason"])

    def test_mixed_high_gradient_edge_case(self):
        """Mixed/unknown with extremely high dT in tiny gap could
        trigger prisms (high wall gradient + small y1).  Just verify
        no crash and flag is boolean."""
        result = BoundaryLayerCalculator.estimate_mesh(
            U=0.0, x=0.005, t_fluid=300.0,
            regime="mixed_unknown",
            delta_t_buoyancy=800.0, t_surf=1100.0,
        )
        self.assertIsInstance(result["prisms_recommended"], bool)
        self.assertIsInstance(result["prism_reason"], str)

    def test_prism_reason_mentions_y1_vs_dx(self):
        """For low-Re mixed zones where y1 > dx, reason should mention it."""
        result = BoundaryLayerCalculator.estimate_mesh(
            U=0.0, x=0.1, t_fluid=353.15,
            regime="mixed_unknown",
            delta_t_buoyancy=200.0, t_surf=553.15,
        )
        # At buoyancy velocities (~0.7 m/s), y1 at y+=30 is large
        if result["y1_mm"] > result["max_dx_surface_mm"]:
            self.assertIn("y1=", result["prism_reason"])


class TestWallGradient(unittest.TestCase):
    """Tests for wall temperature gradient estimation."""

    def test_gradient_with_t_surf(self):
        """Wall gradient should be positive and finite when t_surf given."""
        result = BoundaryLayerCalculator.estimate_mesh(
            U=22.0, x=0.5, t_fluid=353.15,
            regime="external_forced", t_surf=553.15,
        )
        self.assertGreater(result["wall_gradient_Km"], 0)
        self.assertTrue(math.isfinite(result["wall_gradient_Km"]))

    def test_gradient_zero_without_t_surf(self):
        """Wall gradient should be 0 when no surface temperature given."""
        result = BoundaryLayerCalculator.estimate_mesh(
            U=22.0, x=0.5, t_fluid=353.15,
            regime="external_forced",
        )
        self.assertEqual(result["wall_gradient_Km"], 0.0)

    def test_high_re_gradient_exceeds_threshold(self):
        """External forced at typical underbody conditions with hot surface
        should produce gradient well above the prism threshold."""
        result = BoundaryLayerCalculator.estimate_mesh(
            U=22.0, x=0.5, t_fluid=353.15,
            regime="external_forced", t_surf=773.15,  # 500°C surface
        )
        self.assertGreater(
            result["wall_gradient_Km"],
            _WALL_GRADIENT_PRISM_THRESHOLD,
        )

    def test_low_re_gradient_below_threshold(self):
        """Mixed/unknown at buoyancy velocity with moderate dT should
        produce gradient below the threshold."""
        result = BoundaryLayerCalculator.estimate_mesh(
            U=0.0, x=0.1, t_fluid=353.15,
            regime="mixed_unknown",
            delta_t_buoyancy=100.0, t_surf=453.15,
        )
        self.assertLess(
            result["wall_gradient_Km"],
            _WALL_GRADIENT_PRISM_THRESHOLD,
        )


class TestVolumeExpansionRatio(unittest.TestCase):
    """Tests for volume expansion ratio at prism-to-tet transition."""

    def test_er_v_computed_for_external_forced(self):
        """External forced (prisms recommended) should report ER_v > 0."""
        result = BoundaryLayerCalculator.estimate_mesh(
            U=22.0, x=0.5, t_fluid=353.15,
            regime="external_forced", t_surf=553.15,
        )
        self.assertTrue(result["prisms_recommended"])
        self.assertGreater(result["volume_expansion_ratio"], 0)
        self.assertIn(
            result["er_v_status"], ("stable", "marginal", "unstable")
        )

    def test_er_v_na_when_no_prisms(self):
        """Mixed/unknown (prisms not recommended) should report ER_v=0, n/a."""
        result = BoundaryLayerCalculator.estimate_mesh(
            U=0.0, x=0.1, t_fluid=353.15,
            regime="mixed_unknown",
            delta_t_buoyancy=200.0, t_surf=553.15,
        )
        if not result["prisms_recommended"]:
            self.assertEqual(result["volume_expansion_ratio"], 0.0)
            self.assertEqual(result["er_v_status"], "n/a")

    def test_er_v_stable_threshold(self):
        """ER_v at or below ER_V_STABLE should be labelled 'stable'."""
        result = BoundaryLayerCalculator.estimate_mesh(
            U=22.0, x=0.5, t_fluid=353.15,
            y_plus_target=30.0, growth_ratio=1.2,
            regime="external_forced", t_surf=553.15,
        )
        if result["prisms_recommended"]:
            er_v = result["volume_expansion_ratio"]
            status = result["er_v_status"]
            if er_v <= ER_V_STABLE:
                self.assertEqual(status, "stable")
            elif er_v <= ER_V_MARGINAL:
                self.assertEqual(status, "marginal")
            else:
                self.assertEqual(status, "unstable")

    def test_er_v_formula(self):
        """ER_v should equal (dx_surface / y_last)^3."""
        result = BoundaryLayerCalculator.estimate_mesh(
            U=22.0, x=0.5, t_fluid=353.15,
            y_plus_target=30.0, growth_ratio=1.2,
            regime="external_forced", t_surf=553.15,
        )
        if result["prisms_recommended"] and result["y_last_mm"] > 0:
            expected = (result["max_dx_surface_mm"] / result["y_last_mm"]) ** 3
            self.assertAlmostEqual(
                result["volume_expansion_ratio"], expected, places=4
            )


if __name__ == "__main__":
    unittest.main()
