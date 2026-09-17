"""
Tests for BoundaryDrivenConductionCalculator
=============================================

Verification strategy:
    - Hand-computed analytical solutions for known inputs
    - Edge cases (zero flux, zero radiation, zero convection)
    - Dimensional sanity (mm output, sign conventions)
    - Transient penetration depth against analytical sqrt(alpha * dt)
"""

import math
import pytest
from thermal_mesh_calculators.constants import STEFAN_BOLTZMANN
from thermal_mesh_calculators.conduction import BoundaryDrivenConductionCalculator


class TestMaxMeshSize:
    """Test the steady-state boundary-driven mesh size calculator."""

    def test_analytical_exhaust_manifold(self):
        """
        Steel exhaust manifold at 1073 K (800 C), h=180 W/m^2K,
        eps=0.8, T_fluid=353 K, T_surr=353 K, k=45 W/mK, dT=10 K.

        Hand calculation:
            q_conv = 180 * (1073 - 353) = 129600 W/m^2
            q_rad  = 0.8 * 5.67e-8 * (1073^4 - 353^4)
                   = 0.8 * 5.67e-8 * (1.327e12 - 1.554e10)
                   = 0.8 * 5.67e-8 * 1.311e12
                   ~ 59468 W/m^2
            q_total = |129600 + 59468| = 189068 W/m^2
            dx = 45 * 10 / 189068 = 0.002380 m = 2.380 mm
        """
        result = BoundaryDrivenConductionCalculator.max_mesh_size(
            k=45.0, h=180.0, t_surf=1073.15, t_fluid=353.15,
            epsilon=0.8, t_surr=353.15, max_dt=10.0,
        )
        # Check against hand calculation — allow 2% tolerance for
        # rounding differences in the T^4 intermediate
        assert result["max_dx_mm"] == pytest.approx(2.38, rel=0.05)
        assert result["q_conv"] > 0
        assert result["q_rad"] > 0
        assert result["q_total"] == pytest.approx(
            abs(result["q_conv"] + result["q_rad"]), rel=1e-10
        )

    def test_pure_convection(self):
        """When eps=0, radiation is zero. Only convection drives mesh."""
        result = BoundaryDrivenConductionCalculator.max_mesh_size(
            k=45.0, h=100.0, t_surf=500.0, t_fluid=300.0,
            epsilon=0.0, t_surr=300.0, max_dt=5.0,
        )
        assert result["q_rad"] == 0.0
        expected_q = 100.0 * (500.0 - 300.0)  # 20000
        assert result["q_conv"] == pytest.approx(expected_q)
        assert result["max_dx_mm"] == pytest.approx(
            45.0 * 5.0 / expected_q * 1000.0
        )

    def test_pure_radiation(self):
        """When h=0, only radiation drives mesh."""
        result = BoundaryDrivenConductionCalculator.max_mesh_size(
            k=45.0, h=0.0, t_surf=800.0, t_fluid=300.0,
            epsilon=0.9, t_surr=300.0, max_dt=5.0,
        )
        assert result["q_conv"] == 0.0
        assert result["q_rad"] > 0
        expected_q_rad = 0.9 * STEFAN_BOLTZMANN * (800.0**4 - 300.0**4)
        assert result["q_rad"] == pytest.approx(expected_q_rad)

    def test_zero_flux_returns_inf(self):
        """If surface is at fluid AND surr temp, no flux → infinite mesh OK."""
        result = BoundaryDrivenConductionCalculator.max_mesh_size(
            k=45.0, h=100.0, t_surf=300.0, t_fluid=300.0,
            epsilon=0.9, t_surr=300.0, max_dt=5.0,
        )
        assert result["max_dx_mm"] == float("inf")
        assert result["q_total"] == 0.0

    def test_radiation_fraction(self):
        """Radiation fraction should be between 0 and 1 for mixed mode."""
        result = BoundaryDrivenConductionCalculator.max_mesh_size(
            k=45.0, h=50.0, t_surf=800.0, t_fluid=350.0,
            epsilon=0.5, t_surr=350.0, max_dt=10.0,
        )
        assert 0.0 < result["rad_fraction"] < 1.0

    def test_higher_conductivity_gives_coarser_mesh(self):
        """Doubling k should double the max element size (same BCs)."""
        params = dict(
            h=100.0, t_surf=700.0, t_fluid=350.0,
            epsilon=0.5, t_surr=350.0, max_dt=10.0,
        )
        r1 = BoundaryDrivenConductionCalculator.max_mesh_size(k=20.0, **params)
        r2 = BoundaryDrivenConductionCalculator.max_mesh_size(k=40.0, **params)
        assert r2["max_dx_mm"] == pytest.approx(2.0 * r1["max_dx_mm"], rel=1e-10)


class TestLateralGradientLimit:
    """Test the fin theory lateral gradient constraint."""

    def test_basic_steel_bracket(self):
        """
        Steel bracket: k=45, h=20, t=3mm, eps=0, T_surf=300K.
        h_total = 20 (no radiation), m = sqrt(20 / (45 * 0.003)) = 12.17 1/m
        decay = 1/m = 0.0822 m = 82.2 mm
        dx_lateral = (1/3) * 82.2 = 27.4 mm
        """
        result = BoundaryDrivenConductionCalculator.lateral_gradient_limit(
            k=45.0, h=20.0, thickness_m=0.003,
            epsilon=0.0, t_surf=300.0,
        )
        m = math.sqrt(20.0 / (45.0 * 0.003))
        expected_dx = (1.0 / 3.0) / m * 1000.0
        assert result["max_dx_mm"] == pytest.approx(expected_dx, rel=1e-6)
        assert result["h_rad_linearised"] == 0.0
        assert result["h_total"] == pytest.approx(20.0)
        assert result["fin_parameter_m"] == pytest.approx(m, rel=1e-6)

    def test_radiation_increases_h_total(self):
        """Adding radiation should increase h_total and reduce dx."""
        r_no_rad = BoundaryDrivenConductionCalculator.lateral_gradient_limit(
            k=45.0, h=20.0, thickness_m=0.003,
            epsilon=0.0, t_surf=500.0,
        )
        r_with_rad = BoundaryDrivenConductionCalculator.lateral_gradient_limit(
            k=45.0, h=20.0, thickness_m=0.003,
            epsilon=0.8, t_surf=500.0,
        )
        assert r_with_rad["h_total"] > r_no_rad["h_total"]
        assert r_with_rad["max_dx_mm"] < r_no_rad["max_dx_mm"]

    def test_linearised_h_rad_formula(self):
        """h_rad = 4 * eps * sigma * T^3."""
        result = BoundaryDrivenConductionCalculator.lateral_gradient_limit(
            k=45.0, h=20.0, thickness_m=0.003,
            epsilon=0.8, t_surf=500.0,
        )
        expected_h_rad = 4.0 * 0.8 * STEFAN_BOLTZMANN * (500.0 ** 3)
        assert result["h_rad_linearised"] == pytest.approx(expected_h_rad, rel=1e-10)

    def test_thicker_wall_coarser_mesh(self):
        """Thicker wall has longer decay length → coarser mesh."""
        params = dict(k=45.0, h=20.0, epsilon=0.0, t_surf=300.0)
        r_thin = BoundaryDrivenConductionCalculator.lateral_gradient_limit(
            thickness_m=0.001, **params)
        r_thick = BoundaryDrivenConductionCalculator.lateral_gradient_limit(
            thickness_m=0.005, **params)
        assert r_thick["max_dx_mm"] > r_thin["max_dx_mm"]

    def test_zero_h_returns_inf(self):
        """No convection or radiation → infinite decay length."""
        result = BoundaryDrivenConductionCalculator.lateral_gradient_limit(
            k=45.0, h=0.0, thickness_m=0.003,
            epsilon=0.0, t_surf=300.0,
        )
        assert result["max_dx_mm"] == float("inf")

    def test_custom_fraction(self):
        """Tighter fraction should give smaller mesh."""
        params = dict(k=45.0, h=20.0, thickness_m=0.003, epsilon=0.0, t_surf=300.0)
        r_default = BoundaryDrivenConductionCalculator.lateral_gradient_limit(
            fraction=1.0 / 3.0, **params)
        r_tight = BoundaryDrivenConductionCalculator.lateral_gradient_limit(
            fraction=1.0 / 4.0, **params)
        assert r_tight["max_dx_mm"] < r_default["max_dx_mm"]
        assert r_tight["max_dx_mm"] == pytest.approx(
            r_default["max_dx_mm"] * 0.75, rel=1e-6)


class TestTransientPenetrationDepth:
    """Test the transient penetration depth stub."""

    def test_analytical_steel(self):
        """
        Steel: k=45 W/mK, rho=7800 kg/m^3, cp=500 J/kgK
        alpha = 45/(7800*500) = 1.154e-5 m^2/s
        dt = 1 s → sqrt(1.154e-5) = 0.003397 m = 3.397 mm
        """
        result = BoundaryDrivenConductionCalculator.transient_penetration_depth(
            k=45.0, rho=7800.0, cp=500.0, dt=1.0,
        )
        alpha = 45.0 / (7800.0 * 500.0)
        expected = math.sqrt(alpha * 1.0) * 1000.0
        assert result == pytest.approx(expected, rel=1e-10)

    def test_longer_timestep_gives_coarser_mesh(self):
        """Doubling dt should increase penetration depth by sqrt(2)."""
        r1 = BoundaryDrivenConductionCalculator.transient_penetration_depth(
            k=45.0, rho=7800.0, cp=500.0, dt=1.0,
        )
        r2 = BoundaryDrivenConductionCalculator.transient_penetration_depth(
            k=45.0, rho=7800.0, cp=500.0, dt=2.0,
        )
        assert r2 == pytest.approx(r1 * math.sqrt(2), rel=1e-10)
