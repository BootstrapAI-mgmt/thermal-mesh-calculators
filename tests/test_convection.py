"""
Tests for ConvectionMeshCalculator
===================================

Verification strategy:
    - Biot number against hand calculation
    - Shell vs solid transition at Bi = 0.1 boundary
    - h-gradient element count logarithmic scaling
    - CFD mapping ratio linearity
"""

import math
import pytest
from thermal_mesh_calculators.convection import ConvectionMeshCalculator


class TestBiotNumber:
    """Test Biot number evaluation and shell/solid decision."""

    def test_thin_steel_bracket_is_shell(self):
        """
        2 mm steel bracket, h=20 W/m^2K, k=45 W/mK
        Lc = 0.002/2 = 0.001 m
        Bi = 20 * 0.001 / 45 = 0.000444 → shell
        """
        result = ConvectionMeshCalculator.biot_number(
            h=20.0, k=45.0, thickness=0.002,
        )
        assert result["biot"] == pytest.approx(20.0 * 0.001 / 45.0)
        assert result["mesh_type"] == "2D Shell"
        assert result["min_elements_through_thickness"] == 1

    def test_thick_cast_housing_solid(self):
        """
        50 mm cast aluminium, h=200 W/m^2K, k=150 W/mK
        Lc = 0.050/2 = 0.025 m
        Bi = 200 * 0.025 / 150 = 0.0333 → shell (still < 0.1)
        """
        result = ConvectionMeshCalculator.biot_number(
            h=200.0, k=150.0, thickness=0.050,
        )
        assert result["biot"] == pytest.approx(200.0 * 0.025 / 150.0)
        assert result["mesh_type"] == "2D Shell"

    def test_transition_to_solid(self):
        """Force Bi > 0.1 with low conductivity insulator."""
        # Plastic: k=0.2 W/mK, h=50, thickness=0.01 m
        # Lc = 0.005, Bi = 50 * 0.005 / 0.2 = 1.25
        result = ConvectionMeshCalculator.biot_number(
            h=50.0, k=0.2, thickness=0.01,
        )
        assert result["biot"] == pytest.approx(1.25)
        assert result["mesh_type"] == "3D Solid"
        assert result["min_elements_through_thickness"] >= 2

    def test_element_count_heuristic_clamped(self):
        """Element count should be clamped to [2, 10]."""
        # Very high Bi → should clamp at 10
        result = ConvectionMeshCalculator.biot_number(
            h=500.0, k=0.2, thickness=0.02,
        )
        assert result["min_elements_through_thickness"] <= 10

        # Just above 0.1 threshold → at least 2
        # Bi = h * Lc / k, need Bi slightly > 0.1
        # h=10, thickness=0.002, k=0.1 → Bi = 10*0.001/0.1 = 0.1
        # Bump h to 11 → Bi = 0.11
        result2 = ConvectionMeshCalculator.biot_number(
            h=11.0, k=0.1, thickness=0.002,
        )
        assert result2["mesh_type"] == "3D Solid"
        assert result2["min_elements_through_thickness"] >= 2

    def test_cell_biot_heuristic_20x(self):
        """Cell Biot: n = ceil(20*Bi), not ceil(10*Bi).

        Plastic: k=0.25, h=30, thickness=0.006 m
        Lc = 0.003, Bi = 30 * 0.003 / 0.25 = 0.36
        Old: ceil(10*0.36) = 4
        New: ceil(20*0.36) = 8
        """
        result = ConvectionMeshCalculator.biot_number(
            h=30.0, k=0.25, thickness=0.006,
        )
        bi = 30.0 * 0.003 / 0.25
        assert result["biot"] == pytest.approx(bi)
        assert result["mesh_type"] == "3D Solid"
        expected_n = max(2, min(10, math.ceil(20 * bi)))
        assert expected_n == 8  # confirm our hand calc
        assert result["min_elements_through_thickness"] == expected_n

    def test_near_boundary_floating_point(self):
        """Bi ≈ 0.1 — verify the code handles the boundary correctly.

        h=10, thickness=0.002, k=0.1 → Bi = 10*0.001/0.1 = 0.1 analytically.
        Floating point gives 0.0999... which is < 0.1, landing in shell branch.
        This is acceptable — the 0.1 threshold is a guideline, not exact.
        """
        result = ConvectionMeshCalculator.biot_number(
            h=10.0, k=0.1, thickness=0.002,
        )
        assert result["biot"] == pytest.approx(0.1, abs=1e-10)
        # Floating point falls just below 0.1 → shell branch
        assert result["mesh_type"] == "2D Shell"


class TestHGradientMeshLimit:
    """Test h-gradient resolution calculator."""

    def test_basic_gradient(self):
        """h ratio of 16 over 15 mm should give 6 elements."""
        result = ConvectionMeshCalculator.h_gradient_mesh_limit(
            h_max=800.0, h_min=50.0, gradient_length_mm=15.0,
        )
        assert result["h_ratio"] == pytest.approx(16.0)
        expected_n = max(4, math.ceil(2 * math.log(16.0 + 1)))
        assert result["elements_needed"] == expected_n
        assert result["max_dx_mm"] == pytest.approx(15.0 / expected_n)

    def test_uniform_h_uses_minimum_4(self):
        """When h is uniform (ratio=1), still need at least 4 elements."""
        result = ConvectionMeshCalculator.h_gradient_mesh_limit(
            h_max=100.0, h_min=100.0, gradient_length_mm=20.0,
        )
        assert result["h_ratio"] == pytest.approx(1.0)
        assert result["elements_needed"] == 4
        assert result["max_dx_mm"] == pytest.approx(5.0)

    def test_extreme_ratio(self):
        """Very high h ratio should produce many elements but stay reasonable."""
        result = ConvectionMeshCalculator.h_gradient_mesh_limit(
            h_max=10000.0, h_min=1.0, gradient_length_mm=50.0,
        )
        assert result["elements_needed"] >= 4
        assert result["max_dx_mm"] > 0


class TestCFDMappingLimit:
    """Test CFD mapping ratio calculator."""

    def test_default_ratio(self):
        """Default 4:1 mapping → 2 mm face → 8 mm solid."""
        result = ConvectionMeshCalculator.cfd_mapping_mesh_limit(
            fluid_wall_face_mm=2.0,
        )
        assert result == pytest.approx(8.0)

    def test_custom_ratio(self):
        """Custom 2:1 mapping → 3 mm face → 6 mm solid."""
        result = ConvectionMeshCalculator.cfd_mapping_mesh_limit(
            fluid_wall_face_mm=3.0, mapping_ratio=2.0,
        )
        assert result == pytest.approx(6.0)

    def test_linearity(self):
        """Doubling fluid face size should double the result."""
        r1 = ConvectionMeshCalculator.cfd_mapping_mesh_limit(
            fluid_wall_face_mm=2.0, mapping_ratio=3.0,
        )
        r2 = ConvectionMeshCalculator.cfd_mapping_mesh_limit(
            fluid_wall_face_mm=4.0, mapping_ratio=3.0,
        )
        assert r2 == pytest.approx(2.0 * r1)
