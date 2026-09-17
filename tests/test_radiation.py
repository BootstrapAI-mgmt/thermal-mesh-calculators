"""
Tests for RadiationMeshCalculator
==================================

Verification strategy:
    - Flux sensitivity dq/dT = 4*eps*sigma*T^3 against hand calculation
    - T^4 scaling: sensitivity at 1000 K vs 300 K should be ~(1000/300)^3 ratio
    - Curvature limit against chord-angle geometry
    - Edge cases (zero gradient, zero emissivity)
"""

import math
import pytest
from thermal_mesh_calculators.constants import STEFAN_BOLTZMANN
from thermal_mesh_calculators.radiation import RadiationMeshCalculator


class TestFluxSensitivity:
    """Test dq/dT = 4 * eps * sigma * T^3."""

    def test_at_300K(self):
        """dq/dT at 300 K, eps=0.9 should be ~5.51 W/m^2K."""
        result = RadiationMeshCalculator.flux_sensitivity(
            t_local=300.0, emissivity=0.9,
        )
        expected = 4.0 * 0.9 * STEFAN_BOLTZMANN * 300.0**3
        assert result == pytest.approx(expected, rel=1e-10)

    def test_at_1000K(self):
        """dq/dT at 1000 K, eps=0.9."""
        result = RadiationMeshCalculator.flux_sensitivity(
            t_local=1000.0, emissivity=0.9,
        )
        expected = 4.0 * 0.9 * STEFAN_BOLTZMANN * 1000.0**3
        assert result == pytest.approx(expected, rel=1e-10)

    def test_t_cubed_scaling(self):
        """Sensitivity ratio 1000K/300K should be (1000/300)^3 ≈ 37.0."""
        s300 = RadiationMeshCalculator.flux_sensitivity(300.0, 0.9)
        s1000 = RadiationMeshCalculator.flux_sensitivity(1000.0, 0.9)
        assert s1000 / s300 == pytest.approx((1000.0 / 300.0) ** 3, rel=1e-10)

    def test_zero_emissivity(self):
        """Zero emissivity should give zero sensitivity."""
        result = RadiationMeshCalculator.flux_sensitivity(800.0, 0.0)
        assert result == 0.0


class TestMaxMeshSize:
    """Test radiation-driven mesh size."""

    def test_analytical_800C(self):
        """
        T=1073 K, eps=0.9, allowable error=500 W/m^2, gradient=5000 K/m

        dq/dT = 4 * 0.9 * 5.67e-8 * 1073^3 = 251.4 W/m^2K
        max_dT = 500 / 251.4 = 1.989 K
        max_dx = 1.989 / 5000 = 3.978e-4 m = 0.398 mm
        """
        result = RadiationMeshCalculator.max_mesh_size(
            t_local=1073.0, emissivity=0.9,
            allowable_flux_error=500.0, spatial_gradient=5000.0,
        )
        dq_dt = 4.0 * 0.9 * STEFAN_BOLTZMANN * 1073.0**3
        expected_dt = 500.0 / dq_dt
        expected_dx = expected_dt / 5000.0 * 1000.0
        assert result["max_dx_mm"] == pytest.approx(expected_dx, rel=1e-6)
        assert result["dq_dt"] == pytest.approx(dq_dt, rel=1e-10)

    def test_zero_gradient_returns_inf(self):
        """No spatial gradient → infinite mesh is fine."""
        result = RadiationMeshCalculator.max_mesh_size(
            t_local=800.0, emissivity=0.9,
            allowable_flux_error=500.0, spatial_gradient=0.0,
        )
        assert result["max_dx_mm"] == float("inf")

    def test_higher_temp_gives_finer_mesh(self):
        """At higher T, T^3 sensitivity increases, requiring finer mesh."""
        params = dict(
            emissivity=0.9, allowable_flux_error=500.0, spatial_gradient=3000.0,
        )
        r_low = RadiationMeshCalculator.max_mesh_size(t_local=500.0, **params)
        r_high = RadiationMeshCalculator.max_mesh_size(t_local=1000.0, **params)
        assert r_high["max_dx_mm"] < r_low["max_dx_mm"]


class TestViewFactorCurvatureLimit:
    """Test chord-angle geometry for curved surfaces."""

    def test_analytical_25mm_15deg(self):
        """
        R=25 mm, theta=15 deg
        L_max = 2 * 25 * sin(15/2 * pi/180)
              = 50 * sin(7.5 deg)
              = 50 * 0.13053 = 6.526 mm
        """
        result = RadiationMeshCalculator.view_factor_curvature_limit(
            radius_mm=25.0, max_facet_angle_deg=15.0,
        )
        expected = 2.0 * 25.0 * math.sin(math.radians(15.0 / 2.0))
        assert result == pytest.approx(expected, rel=1e-10)

    def test_doubles_with_radius(self):
        """Doubling radius should double the curvature limit."""
        r1 = RadiationMeshCalculator.view_factor_curvature_limit(25.0, 15.0)
        r2 = RadiationMeshCalculator.view_factor_curvature_limit(50.0, 15.0)
        assert r2 == pytest.approx(2.0 * r1, rel=1e-10)

    def test_90_deg_gives_sqrt2_radius(self):
        """
        At 90 deg facet angle: L = 2R * sin(45) = 2R * sqrt(2)/2 = R*sqrt(2).
        """
        r = 100.0
        result = RadiationMeshCalculator.view_factor_curvature_limit(r, 90.0)
        assert result == pytest.approx(r * math.sqrt(2), rel=1e-10)
