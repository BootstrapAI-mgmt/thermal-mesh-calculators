"""
Tests for the input range guards of the primitive calculators
============================================================

Every guarded parameter gets a planted violation: a call that is valid
except for that one input, which must raise ValueError naming the
parameter.  The valid call itself is checked first, so each case fails
for its planted reason and not for another.  The boundary values the
ranges include (emissivity 0 and 1, h = 0) must still be accepted.
"""

import math

import pytest

from thermal_mesh_calculators.boundary_layer import BoundaryLayerCalculator
from thermal_mesh_calculators.conduction import BoundaryDrivenConductionCalculator
from thermal_mesh_calculators.convection import ConvectionMeshCalculator
from thermal_mesh_calculators.h_estimator import (
    air_properties,
    estimate_h,
    forced_convection_flat_plate,
    natural_convection,
    richardson_number,
)
from thermal_mesh_calculators.radiation import RadiationMeshCalculator
from thermal_mesh_calculators.shields import (
    MultilayerShieldCalculator,
    SingleLayerShieldCalculator,
)
from thermal_mesh_calculators.transient import TransientMeshCalculator
from thermal_mesh_calculators.zones import estimate_spatial_gradient


COND = dict(k=45.0, h=100.0, t_surf=800.0, t_fluid=300.0, epsilon=0.8,
            t_surr=300.0, max_dt=5.0)
LATERAL = dict(k=45.0, h=20.0, thickness_m=0.003, epsilon=0.8, t_surf=500.0)
MATERIAL = dict(k=45.0, rho=7800.0, cp=500.0)
SL_SOLVE = dict(t_exh=1073.15, t_fluid=353.15, t_surr=353.15, eps_in=0.4,
                eps_out=0.4, h_in=30.0, h_out=30.0)
ML_SOLVE = dict(t_exh=1073.15, t_fluid=353.15, t_surr=353.15, h_in=30.0,
                h_out=30.0, h_gap=15.0, eps_in=0.4, eps_out=0.4, eps_g1=0.4,
                eps_g2=0.4)

B = BoundaryDrivenConductionCalculator
T = TransientMeshCalculator

# (label, callable, valid keyword arguments, parameter, planted value)
CASES = [
    ("max_mesh_size", B.max_mesh_size, COND, "k", -45.0),
    ("max_mesh_size", B.max_mesh_size, COND, "k", 0.0),
    ("max_mesh_size", B.max_mesh_size, COND, "k", float("nan")),
    ("max_mesh_size", B.max_mesh_size, COND, "k", float("inf")),
    ("max_mesh_size", B.max_mesh_size, COND, "h", -1.0),
    ("max_mesh_size", B.max_mesh_size, COND, "t_surf", 0.0),
    ("max_mesh_size", B.max_mesh_size, COND, "t_fluid", -10.0),
    ("max_mesh_size", B.max_mesh_size, COND, "t_surr", 0.0),
    ("max_mesh_size", B.max_mesh_size, COND, "epsilon", 1.7),
    ("max_mesh_size", B.max_mesh_size, COND, "epsilon", -0.1),
    ("max_mesh_size", B.max_mesh_size, COND, "max_dt", 0.0),
    ("size_wall", B.size_wall, dict(COND, thickness_m=0.004), "epsilon", 2.0),
    ("size_wall", B.size_wall, dict(COND, thickness_m=0.004), "t_surr", -1.0),
    ("lateral_gradient_limit", B.lateral_gradient_limit, LATERAL, "k", 0.0),
    ("lateral_gradient_limit", B.lateral_gradient_limit, LATERAL, "h", -1.0),
    ("lateral_gradient_limit", B.lateral_gradient_limit, LATERAL, "epsilon", 1.5),
    ("lateral_gradient_limit", B.lateral_gradient_limit, LATERAL, "t_surf", 0.0),
    ("transient_penetration_depth", B.transient_penetration_depth,
     dict(MATERIAL, dt=1.0), "rho", 0.0),
    ("transient_penetration_depth", B.transient_penetration_depth,
     dict(MATERIAL, dt=1.0), "dt", 0.0),
    ("biot_number", ConvectionMeshCalculator.biot_number,
     dict(h=20.0, k=45.0, thickness=0.002), "k", 0.0),
    ("biot_number", ConvectionMeshCalculator.biot_number,
     dict(h=20.0, k=45.0, thickness=0.002), "h", -5.0),
    ("flux_sensitivity", RadiationMeshCalculator.flux_sensitivity,
     dict(t_local=800.0, emissivity=0.8), "t_local", 0.0),
    ("flux_sensitivity", RadiationMeshCalculator.flux_sensitivity,
     dict(t_local=800.0, emissivity=0.8), "emissivity", 1.01),
    ("radiation max_mesh_size", RadiationMeshCalculator.max_mesh_size,
     dict(t_local=800.0, emissivity=0.8, allowable_flux_error=500.0,
          spatial_gradient=1000.0), "t_local", -1.0),
    ("radiation max_mesh_size", RadiationMeshCalculator.max_mesh_size,
     dict(t_local=800.0, emissivity=0.8, allowable_flux_error=500.0,
          spatial_gradient=1000.0), "emissivity", -0.5),
    ("radiation max_mesh_size", RadiationMeshCalculator.max_mesh_size,
     dict(t_local=800.0, emissivity=0.8, allowable_flux_error=500.0,
          spatial_gradient=1000.0), "allowable_flux_error", 0.0),
    ("radiation max_mesh_size", RadiationMeshCalculator.max_mesh_size,
     dict(t_local=800.0, emissivity=0.8, allowable_flux_error=500.0,
          spatial_gradient=1000.0), "spatial_gradient", -10.0),
    ("single solve_temperature", SingleLayerShieldCalculator.solve_temperature,
     SL_SOLVE, "t_exh", 0.0),
    ("single solve_temperature", SingleLayerShieldCalculator.solve_temperature,
     SL_SOLVE, "t_fluid", -1.0),
    ("single solve_temperature", SingleLayerShieldCalculator.solve_temperature,
     SL_SOLVE, "t_surr", 0.0),
    ("single solve_temperature", SingleLayerShieldCalculator.solve_temperature,
     SL_SOLVE, "eps_in", 1.2),
    ("single solve_temperature", SingleLayerShieldCalculator.solve_temperature,
     SL_SOLVE, "eps_out", -0.1),
    ("single solve_temperature", SingleLayerShieldCalculator.solve_temperature,
     SL_SOLVE, "h_in", -1.0),
    ("single solve_temperature", SingleLayerShieldCalculator.solve_temperature,
     SL_SOLVE, "h_out", -1.0),
    ("single mesh_size", SingleLayerShieldCalculator.mesh_size,
     dict(SL_SOLVE, k=45.0, max_dt=15.0), "k", 0.0),
    ("single mesh_size", SingleLayerShieldCalculator.mesh_size,
     dict(SL_SOLVE, k=45.0, max_dt=15.0), "max_dt", -15.0),
    ("multi solve_temperatures", MultilayerShieldCalculator.solve_temperatures,
     ML_SOLVE, "t_exh", 0.0),
    ("multi solve_temperatures", MultilayerShieldCalculator.solve_temperatures,
     ML_SOLVE, "t_fluid", 0.0),
    ("multi solve_temperatures", MultilayerShieldCalculator.solve_temperatures,
     ML_SOLVE, "t_surr", -5.0),
    ("multi solve_temperatures", MultilayerShieldCalculator.solve_temperatures,
     ML_SOLVE, "h_in", -1.0),
    ("multi solve_temperatures", MultilayerShieldCalculator.solve_temperatures,
     ML_SOLVE, "h_out", -1.0),
    ("multi solve_temperatures", MultilayerShieldCalculator.solve_temperatures,
     ML_SOLVE, "h_gap", -1.0),
    ("multi solve_temperatures", MultilayerShieldCalculator.solve_temperatures,
     ML_SOLVE, "eps_in", 1.5),
    ("multi solve_temperatures", MultilayerShieldCalculator.solve_temperatures,
     ML_SOLVE, "eps_out", 1.5),
    ("multi solve_temperatures", MultilayerShieldCalculator.solve_temperatures,
     ML_SOLVE, "eps_g1", 1.5),
    ("multi solve_temperatures", MultilayerShieldCalculator.solve_temperatures,
     ML_SOLVE, "eps_g2", -0.2),
    ("multi solve_temperatures", MultilayerShieldCalculator.solve_temperatures,
     ML_SOLVE, "f12", 1.2),
    ("multi mesh_sizes", MultilayerShieldCalculator.mesh_sizes,
     dict(ML_SOLVE, k_metal=45.0, max_dt=15.0), "k_metal", 0.0),
    ("multi mesh_sizes", MultilayerShieldCalculator.mesh_sizes,
     dict(ML_SOLVE, k_metal=45.0, max_dt=15.0), "max_dt", 0.0),
    ("material_diffusivity", T.material_diffusivity, MATERIAL, "k", 0.0),
    ("material_diffusivity", T.material_diffusivity, MATERIAL, "rho", 0.0),
    ("material_diffusivity", T.material_diffusivity, MATERIAL, "cp", -1.0),
    ("penetration_depth", T.penetration_depth, dict(MATERIAL, dt=1.0), "dt", 0.0),
    ("penetration_depth", T.penetration_depth, dict(MATERIAL, dt=1.0),
     "safety_factor", 0.0),
    ("fourier_number_limit", T.fourier_number_limit, dict(MATERIAL, dt=1.0),
     "fo_max", 0.0),
    ("fourier_number_limit", T.fourier_number_limit, dict(MATERIAL, dt=1.0),
     "dt", -1.0),
    ("drive_cycle_limit", T.drive_cycle_limit, dict(MATERIAL, tau_bc=20.0),
     "tau_bc", 0.0),
    ("combined_transient_limits", T.combined_transient_limits,
     dict(MATERIAL, dt=1.0), "rho", 0.0),
    ("combined_transient_limits", T.combined_transient_limits,
     dict(MATERIAL, dt=1.0), "dt", 0.0),
    ("combined_transient_limits", T.combined_transient_limits,
     dict(MATERIAL, dt=1.0, tau_bc=20.0), "tau_bc", -5.0),
    ("air_properties", air_properties, dict(t_film=400.0), "t_film", 0.0),
    ("forced_convection_flat_plate", forced_convection_flat_plate,
     dict(velocity=4.47, char_length=0.1, t_surf=473.15, t_fluid=313.15),
     "t_surf", 0.0),
    ("forced_convection_flat_plate", forced_convection_flat_plate,
     dict(velocity=4.47, char_length=0.1, t_surf=473.15, t_fluid=313.15),
     "velocity", -1.0),
    ("forced_convection_flat_plate", forced_convection_flat_plate,
     dict(velocity=4.47, char_length=0.1, t_surf=473.15, t_fluid=313.15),
     "char_length", 0.0),
    ("natural_convection", natural_convection,
     dict(t_surf=573.15, t_fluid=313.15, char_length=0.2), "t_fluid", -1.0),
    ("natural_convection", natural_convection,
     dict(t_surf=573.15, t_fluid=313.15, char_length=0.2), "char_length", -0.1),
    ("richardson_number", richardson_number,
     dict(velocity=1.0, t_surf=573.15, t_fluid=313.15, char_length=0.1),
     "t_surf", 0.0),
    ("estimate_h", estimate_h,
     dict(velocity=0.0, t_surf=473.15, t_fluid=313.15, char_length=0.15),
     "t_fluid", 0.0),
    ("estimate_h", estimate_h,
     dict(velocity=0.0, t_surf=473.15, t_fluid=313.15, char_length=0.15),
     "velocity", -2.0),
    ("boundary_layer estimate_mesh", BoundaryLayerCalculator.estimate_mesh,
     dict(U=22.0, x=0.5, t_fluid=353.15), "t_fluid", 0.0),
    ("boundary_layer estimate_mesh", BoundaryLayerCalculator.estimate_mesh,
     dict(U=22.0, x=0.5, t_fluid=353.15), "t_surf", -1.0),
    ("estimate_spatial_gradient", estimate_spatial_gradient,
     dict(k=45.0, q_total=9000.0), "k", 0.0),
]


def _case_id(case):
    label, _, _, name, value = case
    return f"{label}-{name}={value!r}"


@pytest.mark.parametrize("label, func, valid, name, bad", CASES,
                         ids=[_case_id(case) for case in CASES])
def test_planted_violation_raises_naming_the_parameter(label, func, valid, name, bad):
    func(**valid)  # the call is valid without the planted value
    with pytest.raises(ValueError) as excinfo:
        func(**dict(valid, **{name: bad}))
    assert str(excinfo.value).startswith(name + " must ")


@pytest.mark.parametrize("bad", ["45", True, None])
def test_a_non_number_raises_type_error(bad):
    with pytest.raises(TypeError, match="k must be a number"):
        B.max_mesh_size(**dict(COND, k=bad))


class TestRangeBoundariesAccepted:

    @pytest.mark.parametrize("eps", [0.0, 1.0])
    def test_emissivity_limits(self, eps):
        result = B.max_mesh_size(**dict(COND, epsilon=eps))
        assert 0 < result["max_dx_mm"] < float("inf")

    def test_zero_h_is_allowed(self):
        result = B.max_mesh_size(**dict(COND, h=0.0))
        assert result["q_conv"] == 0.0

    @pytest.mark.parametrize("name", ["eps_g1", "eps_g2", "f12"])
    def test_zero_gap_exchange_is_the_limit_not_an_error(self, name):
        """eps_eff -> 0 as a gap emissivity or the view factor -> 0."""
        result = MultilayerShieldCalculator.solve_temperatures(
            **dict(ML_SOLVE, **{name: 0.0}))
        assert result["eps_eff"] == 0.0
        assert result["q_gap_rad"] == 0.0
        assert result["converged"] is True
        tiny = MultilayerShieldCalculator.solve_temperatures(
            **dict(ML_SOLVE, **{name: 1e-9}))
        assert tiny["eps_eff"] == pytest.approx(0.0, abs=1e-8)

    def test_boundary_layer_t_surf_zero_means_not_given(self):
        result = BoundaryLayerCalculator.estimate_mesh(
            U=22.0, x=0.5, t_fluid=353.15, t_surf=0.0)
        assert math.isfinite(result["max_dx_surface_mm"])
