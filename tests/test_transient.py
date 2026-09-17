"""
Tests for TransientMeshCalculator
==================================

Verification strategy:
    - Analytical penetration depth: sqrt(alpha * dt)
    - Fourier number: dx = sqrt(alpha * dt / Fo)
    - Drive-cycle limit: sqrt(alpha * tau_bc)
    - Combined constraint binding logic
    - Dimensional consistency (mm output)
    - Scaling relationships (sqrt(dt), sqrt(k))
"""

import math
import pytest
from thermal_mesh_calculators.transient import TransientMeshCalculator


# --- Material data for tests ---
STEEL = {"k": 45.0, "rho": 7800.0, "cp": 500.0}    # alpha = 1.154e-5
PLASTIC = {"k": 0.25, "rho": 1400.0, "cp": 1600.0}  # alpha = 1.116e-7


class TestMaterialDiffusivity:
    """Test thermal diffusivity calculation."""

    def test_steel(self):
        alpha = TransientMeshCalculator.material_diffusivity(**STEEL)
        assert alpha == pytest.approx(45.0 / (7800.0 * 500.0), rel=1e-10)

    def test_plastic(self):
        alpha = TransientMeshCalculator.material_diffusivity(**PLASTIC)
        assert alpha == pytest.approx(0.25 / (1400.0 * 1600.0), rel=1e-10)


class TestPenetrationDepth:
    """Test thermal penetration depth constraint."""

    def test_analytical_steel(self):
        """Steel, dt=1s: pen = sqrt(1.154e-5 * 1) = 3.397e-3 m = 3.397 mm."""
        result = TransientMeshCalculator.penetration_depth(
            **STEEL, dt=1.0,
        )
        alpha = 45.0 / (7800.0 * 500.0)
        expected = math.sqrt(alpha * 1.0) * 1000.0
        assert result["max_dx_mm"] == pytest.approx(expected, rel=1e-10)
        assert result["alpha"] == pytest.approx(alpha, rel=1e-10)

    def test_safety_factor(self):
        """Safety factor of 2.0 should double the max dx."""
        r1 = TransientMeshCalculator.penetration_depth(**STEEL, dt=1.0, safety_factor=1.0)
        r2 = TransientMeshCalculator.penetration_depth(**STEEL, dt=1.0, safety_factor=2.0)
        assert r2["max_dx_mm"] == pytest.approx(2.0 * r1["max_dx_mm"], rel=1e-10)

    def test_sqrt_dt_scaling(self):
        """Quadrupling dt should double the penetration depth."""
        r1 = TransientMeshCalculator.penetration_depth(**STEEL, dt=1.0)
        r4 = TransientMeshCalculator.penetration_depth(**STEEL, dt=4.0)
        assert r4["max_dx_mm"] == pytest.approx(2.0 * r1["max_dx_mm"], rel=1e-10)

    def test_plastic_much_smaller(self):
        """Plastic has ~100x lower diffusivity → ~10x smaller penetration."""
        r_steel = TransientMeshCalculator.penetration_depth(**STEEL, dt=1.0)
        r_plastic = TransientMeshCalculator.penetration_depth(**PLASTIC, dt=1.0)
        assert r_plastic["max_dx_mm"] < r_steel["max_dx_mm"] / 5.0


class TestFourierNumberLimit:
    """Test Fourier number constraint on minimum element size."""

    def test_analytical_explicit(self):
        """
        Steel, dt=0.5s, Fo_max=0.5:
        dx_min = sqrt(alpha * dt / Fo_max) = sqrt(1.154e-5 * 0.5 / 0.5)
               = sqrt(1.154e-5) = 3.397e-3 m = 3.397 mm
        """
        result = TransientMeshCalculator.fourier_number_limit(
            **STEEL, dt=0.5, fo_max=0.5,
        )
        alpha = 45.0 / (7800.0 * 500.0)
        expected = math.sqrt(alpha * 0.5 / 0.5) * 1000.0
        assert result["min_dx_mm"] == pytest.approx(expected, rel=1e-10)
        assert result["constraint"] == "stability"

    def test_implicit_accuracy(self):
        """Fo_max=5.0 for implicit solver should give 'accuracy' constraint."""
        result = TransientMeshCalculator.fourier_number_limit(
            **STEEL, dt=0.5, fo_max=5.0,
        )
        assert result["constraint"] == "accuracy"
        # Larger Fo_max → smaller min_dx (more permissive stability)
        r_explicit = TransientMeshCalculator.fourier_number_limit(
            **STEEL, dt=0.5, fo_max=0.5,
        )
        assert result["min_dx_mm"] < r_explicit["min_dx_mm"]

    def test_3d_explicit_limit(self):
        """Fo_max = 1/6 for 3D explicit should give tighter constraint."""
        r_1d = TransientMeshCalculator.fourier_number_limit(
            **STEEL, dt=0.5, fo_max=0.5,
        )
        r_3d = TransientMeshCalculator.fourier_number_limit(
            **STEEL, dt=0.5, fo_max=1.0 / 6.0,
        )
        assert r_3d["min_dx_mm"] > r_1d["min_dx_mm"]  # tighter = larger minimum


class TestDriveCycleLimit:
    """Test drive-cycle aware mesh sizing."""

    def test_analytical(self):
        """Steel, tau_bc=20s: dx = sqrt(1.154e-5 * 20) = 0.01519 m = 15.19 mm."""
        result = TransientMeshCalculator.drive_cycle_limit(
            **STEEL, tau_bc=20.0,
        )
        alpha = 45.0 / (7800.0 * 500.0)
        expected = math.sqrt(alpha * 20.0) * 1000.0
        assert result["max_dx_mm"] == pytest.approx(expected, rel=1e-6)

    def test_longer_segment_coarser(self):
        """Longer drive-cycle segment allows coarser mesh."""
        r1 = TransientMeshCalculator.drive_cycle_limit(**STEEL, tau_bc=10.0)
        r2 = TransientMeshCalculator.drive_cycle_limit(**STEEL, tau_bc=40.0)
        assert r2["max_dx_mm"] == pytest.approx(2.0 * r1["max_dx_mm"], rel=1e-10)


class TestCombinedTransientLimits:
    """Test the combined constraint evaluator."""

    def test_normal_implicit_binding(self):
        """Implicit solver with relaxed Fo should bind on penetration or drive-cycle."""
        result = TransientMeshCalculator.combined_transient_limits(
            **PLASTIC, dt=2.0, fo_max=5.0, safety_factor=2.0, tau_bc=30.0,
        )
        assert result["recommended_dx_mm"] > 0
        assert result["binding_constraint"] in ("penetration_depth", "drive_cycle")
        # Drive cycle and penetration should both be present
        assert result["drive_cycle_max_dx_mm"] is not None

    def test_no_drive_cycle(self):
        """Without tau_bc, drive_cycle_max_dx_mm should be None."""
        result = TransientMeshCalculator.combined_transient_limits(
            **STEEL, dt=0.5, fo_max=0.5,
        )
        assert result["drive_cycle_max_dx_mm"] is None

    def test_fourier_conflict_flagged(self):
        """When Fo stability demands larger dx than penetration allows,
        the binding constraint should flag a conflict (reduce dt)."""
        # Use a deliberately large dt to trigger the conflict
        result = TransientMeshCalculator.combined_transient_limits(
            **STEEL, dt=0.5, fo_max=0.5, safety_factor=1.0,
        )
        # Check if conflict exists
        fo_min = result["fourier_min_dx_mm"]
        pen_max = result["penetration_max_dx_mm"]
        if fo_min > pen_max:
            assert "WARNING" in result["binding_constraint"]

    def test_alpha_consistent(self):
        """Alpha in combined result should match material_diffusivity."""
        result = TransientMeshCalculator.combined_transient_limits(
            **STEEL, dt=1.0,
        )
        expected_alpha = TransientMeshCalculator.material_diffusivity(**STEEL)
        assert result["alpha"] == pytest.approx(expected_alpha, rel=1e-10)

    def test_all_keys_present(self):
        """All expected keys should be in the result dict."""
        result = TransientMeshCalculator.combined_transient_limits(
            **STEEL, dt=1.0, tau_bc=10.0,
        )
        expected_keys = {
            "alpha", "penetration_max_dx_mm", "fourier_min_dx_mm",
            "drive_cycle_max_dx_mm", "recommended_dx_mm", "binding_constraint",
        }
        assert set(result.keys()) == expected_keys
