"""
Tests for SingleLayerShieldCalculator and MultilayerShieldCalculator
====================================================================

Verification strategy:
    - Energy balance closure: q_in ≈ q_out at solved temperature
    - Newton-Raphson convergence for normal and edge cases
    - Backward compatibility (h_total legacy API)
    - Analytical limits (eps→0 convection-only, h→0 radiation-only)
    - Multilayer: gap energy balance and Jacobian correctness
    - Mesh size consistency: boundary flux → dx_max relationship
"""

import pytest
from thermal_mesh_calculators.shields import (
    SingleLayerShieldCalculator,
    MultilayerShieldCalculator,
)


# -----------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------

def _energy_balance_residual_single(result, t_exh, t_fluid, t_surr, eps_in, eps_out):
    """Compute energy balance residual for single-layer shield.

    At equilibrium: q_rad_in - q_rad_out - q_conv_in - q_conv_out = 0
    """
    return (
        result["q_rad_in"]
        - result["q_rad_out"]
        - result["q_conv_in"]
        - result["q_conv_out"]
    )


def _energy_balance_residual_multi(result):
    """Check both layer balances for multilayer shield.

    Layer 1: q_rad_in + q_conv_in - q_gap_cond - q_gap_rad = 0
    Layer 2: q_gap_cond + q_gap_rad - q_rad_out - q_conv_out = 0
    """
    r1 = (result["q_rad_in"] + result["q_conv_in"]
          - result["q_gap_cond"] - result["q_gap_rad"])
    r2 = (result["q_gap_cond"] + result["q_gap_rad"]
          - result["q_rad_out"] - result["q_conv_out"])
    return r1, r2


# -----------------------------------------------------------------------
# SingleLayerShieldCalculator
# -----------------------------------------------------------------------

class TestSingleLayerSolveTemperature:
    """Test the scalar Newton-Raphson shield solver."""

    def test_convergence_typical(self):
        """Typical aluminised shield should converge quickly."""
        result = SingleLayerShieldCalculator.solve_temperature(
            t_exh=1073.15, t_fluid=353.15, t_surr=353.15,
            eps_in=0.4, eps_out=0.4, h_in=30.0, h_out=30.0,
        )
        assert result["converged"] is True
        assert result["iterations"] < 20

    def test_energy_balance_closure(self):
        """Residual should be near zero at converged temperature."""
        result = SingleLayerShieldCalculator.solve_temperature(
            t_exh=1073.15, t_fluid=353.15, t_surr=353.15,
            eps_in=0.4, eps_out=0.4, h_in=30.0, h_out=30.0,
        )
        residual = _energy_balance_residual_single(
            result, 1073.15, 353.15, 353.15, 0.4, 0.4,
        )
        # Residual should be within solver tolerance (0.1 K drives ~few W/m^2)
        assert abs(residual) < 50.0  # W/m^2

    def test_shield_between_source_and_sink(self):
        """Shield temperature must be between T_surr and T_exh."""
        result = SingleLayerShieldCalculator.solve_temperature(
            t_exh=1073.15, t_fluid=353.15, t_surr=353.15,
            eps_in=0.4, eps_out=0.4, h_in=30.0, h_out=30.0,
        )
        assert result["t_shield_K"] > 353.15
        assert result["t_shield_K"] < 1073.15

    def test_celsius_conversion(self):
        """t_shield_C should be t_shield_K - 273.15."""
        result = SingleLayerShieldCalculator.solve_temperature(
            t_exh=1073.15, t_fluid=353.15, t_surr=353.15,
            eps_in=0.4, eps_out=0.4, h_in=30.0, h_out=30.0,
        )
        assert result["t_shield_C"] == pytest.approx(
            result["t_shield_K"] - 273.15, abs=1e-10,
        )

    def test_backward_compat_h_total(self):
        """Legacy h_total parameter should split evenly to h_in/h_out."""
        r_new = SingleLayerShieldCalculator.solve_temperature(
            t_exh=1073.15, t_fluid=353.15, t_surr=353.15,
            eps_in=0.4, eps_out=0.4, h_in=30.0, h_out=30.0,
        )
        r_legacy = SingleLayerShieldCalculator.solve_temperature(
            t_exh=1073.15, t_fluid=353.15, t_surr=353.15,
            eps_in=0.4, eps_out=0.4, h_total=60.0,
        )
        assert r_new["t_shield_K"] == pytest.approx(r_legacy["t_shield_K"], abs=0.2)

    def test_missing_htc_raises(self):
        """Must provide either (h_in, h_out) or h_total."""
        with pytest.raises(ValueError, match="Provide either"):
            SingleLayerShieldCalculator.solve_temperature(
                t_exh=1073.15, t_fluid=353.15, t_surr=353.15,
                eps_in=0.4, eps_out=0.4,
            )

    def test_asymmetric_convection(self):
        """Asymmetric h_in/h_out should produce different conv fluxes."""
        result = SingleLayerShieldCalculator.solve_temperature(
            t_exh=1073.15, t_fluid=353.15, t_surr=353.15,
            eps_in=0.4, eps_out=0.4, h_in=15.0, h_out=60.0,
        )
        assert result["converged"] is True
        # Both convective fluxes positive (shield hotter than fluid)
        assert result["q_conv_in"] > 0
        assert result["q_conv_out"] > 0
        # Outer convection should be higher (h_out > h_in, same delta-T)
        assert result["q_conv_out"] > result["q_conv_in"]

    def test_high_emissivity_raises_shield_temp(self):
        """Higher emissivity absorbs more radiation → hotter shield."""
        r_low = SingleLayerShieldCalculator.solve_temperature(
            t_exh=1073.15, t_fluid=353.15, t_surr=353.15,
            eps_in=0.2, eps_out=0.2, h_in=30.0, h_out=30.0,
        )
        r_high = SingleLayerShieldCalculator.solve_temperature(
            t_exh=1073.15, t_fluid=353.15, t_surr=353.15,
            eps_in=0.85, eps_out=0.85, h_in=30.0, h_out=30.0,
        )
        assert r_high["t_shield_K"] > r_low["t_shield_K"]

    def test_pure_convection_eps_zero(self):
        """With eps=0, shield temperature should be entirely convection-driven.

        With eps_in=eps_out=0, energy balance reduces to:
            0 - 0 - h_in*(T-T_f) - h_out*(T-T_f) = 0
        Which gives T = T_fluid (no radiation heating).
        """
        result = SingleLayerShieldCalculator.solve_temperature(
            t_exh=1073.15, t_fluid=353.15, t_surr=353.15,
            eps_in=0.0, eps_out=0.0, h_in=30.0, h_out=30.0,
        )
        # Shield should converge to fluid temperature
        assert result["t_shield_K"] == pytest.approx(353.15, abs=1.0)


class TestSingleLayerMeshSize:
    """Test mesh size computation layered on top of temperature solve."""

    def test_mesh_size_positive(self):
        """Mesh size should be a positive finite number."""
        result = SingleLayerShieldCalculator.mesh_size(
            k=45.0, max_dt=15.0,
            t_exh=1073.15, t_fluid=353.15, t_surr=353.15,
            h_in=30.0, h_out=30.0, eps_in=0.4, eps_out=0.4,
        )
        assert result["max_dx_mm"] > 0
        assert result["max_dx_mm"] < float("inf")

    def test_mesh_size_formula(self):
        """Verify dx = k * dT / q_boundary * 1000."""
        result = SingleLayerShieldCalculator.mesh_size(
            k=45.0, max_dt=15.0,
            t_exh=1073.15, t_fluid=353.15, t_surr=353.15,
            h_in=30.0, h_out=30.0, eps_in=0.4, eps_out=0.4,
        )
        expected = (45.0 * 15.0) / result["q_boundary"] * 1000.0
        assert result["max_dx_mm"] == pytest.approx(expected, rel=1e-10)

    def test_boundary_flux_is_inner_face(self):
        """q_boundary should use inner face (rad_in + conv_in), not total."""
        result = SingleLayerShieldCalculator.mesh_size(
            k=45.0, max_dt=15.0,
            t_exh=1073.15, t_fluid=353.15, t_surr=353.15,
            h_in=30.0, h_out=30.0, eps_in=0.4, eps_out=0.4,
        )
        expected_q = abs(result["q_rad_in"]) + abs(result["q_conv_in"])
        assert result["q_boundary"] == pytest.approx(expected_q, rel=1e-10)


# -----------------------------------------------------------------------
# MultilayerShieldCalculator
# -----------------------------------------------------------------------

class TestMultilayerSolveTemperatures:
    """Test the 2x2 multivariate Newton-Raphson solver."""

    def test_convergence_typical(self):
        """Typical dual-wall aluminised shield should converge."""
        result = MultilayerShieldCalculator.solve_temperatures(
            t_exh=1073.15, t_fluid=353.15, t_surr=353.15,
            h_in=30.0, h_out=30.0, h_gap=15.0,
            eps_in=0.4, eps_out=0.4, eps_g1=0.4, eps_g2=0.4,
        )
        assert result["converged"] is True
        assert result["iterations"] < 30

    def test_energy_balance_both_layers(self):
        """Both layer energy balances should close within tolerance."""
        result = MultilayerShieldCalculator.solve_temperatures(
            t_exh=1073.15, t_fluid=353.15, t_surr=353.15,
            h_in=30.0, h_out=30.0, h_gap=15.0,
            eps_in=0.4, eps_out=0.4, eps_g1=0.4, eps_g2=0.4,
        )
        r1, r2 = _energy_balance_residual_multi(result)
        assert abs(r1) < 1.0  # W/m^2 — tighter than single-layer (tol on F)
        assert abs(r2) < 1.0

    def test_layer1_hotter_than_layer2(self):
        """Inner layer (exhaust-facing) should be hotter than outer."""
        result = MultilayerShieldCalculator.solve_temperatures(
            t_exh=1073.15, t_fluid=353.15, t_surr=353.15,
            h_in=30.0, h_out=30.0, h_gap=15.0,
            eps_in=0.4, eps_out=0.4, eps_g1=0.4, eps_g2=0.4,
        )
        assert result["t1_K"] > result["t2_K"]
        assert result["delta_T_C"] > 0

    def test_temperatures_bounded(self):
        """Both layers must be between T_surr and T_exh."""
        result = MultilayerShieldCalculator.solve_temperatures(
            t_exh=1073.15, t_fluid=353.15, t_surr=353.15,
            h_in=30.0, h_out=30.0, h_gap=15.0,
            eps_in=0.4, eps_out=0.4, eps_g1=0.4, eps_g2=0.4,
        )
        assert result["t1_K"] > 353.15
        assert result["t1_K"] < 1073.15
        assert result["t2_K"] > 353.15
        assert result["t2_K"] < 1073.15

    def test_effective_emissivity(self):
        """eps_eff = 1/(1/e1 + 1/e2 - 1) should be computed correctly."""
        result = MultilayerShieldCalculator.solve_temperatures(
            t_exh=1073.15, t_fluid=353.15, t_surr=353.15,
            h_in=30.0, h_out=30.0, h_gap=15.0,
            eps_in=0.4, eps_out=0.4, eps_g1=0.5, eps_g2=0.6,
        )
        expected_eps = 1.0 / (1.0 / 0.5 + 1.0 / 0.6 - 1.0)
        assert result["eps_eff"] == pytest.approx(expected_eps, rel=1e-10)

    def test_high_gap_conductance_reduces_delta_T(self):
        """Higher h_gap should bring layers closer in temperature."""
        r_low = MultilayerShieldCalculator.solve_temperatures(
            t_exh=1073.15, t_fluid=353.15, t_surr=353.15,
            h_in=30.0, h_out=30.0, h_gap=5.0,
            eps_in=0.4, eps_out=0.4, eps_g1=0.4, eps_g2=0.4,
        )
        r_high = MultilayerShieldCalculator.solve_temperatures(
            t_exh=1073.15, t_fluid=353.15, t_surr=353.15,
            h_in=30.0, h_out=30.0, h_gap=100.0,
            eps_in=0.4, eps_out=0.4, eps_g1=0.4, eps_g2=0.4,
        )
        assert r_high["delta_T_C"] < r_low["delta_T_C"]

    def test_celsius_conversion(self):
        """C outputs should be K - 273.15."""
        result = MultilayerShieldCalculator.solve_temperatures(
            t_exh=1073.15, t_fluid=353.15, t_surr=353.15,
            h_in=30.0, h_out=30.0, h_gap=15.0,
            eps_in=0.4, eps_out=0.4, eps_g1=0.4, eps_g2=0.4,
        )
        assert result["t1_C"] == pytest.approx(result["t1_K"] - 273.15, abs=1e-10)
        assert result["t2_C"] == pytest.approx(result["t2_K"] - 273.15, abs=1e-10)

    def test_f12_default_matches_original(self):
        """With F12=1.0 (default), eps_eff should match the original formula."""
        result = MultilayerShieldCalculator.solve_temperatures(
            t_exh=1073.15, t_fluid=353.15, t_surr=353.15,
            h_in=30.0, h_out=30.0, h_gap=15.0,
            eps_in=0.4, eps_out=0.4, eps_g1=0.5, eps_g2=0.6,
            f12=1.0,
        )
        # 1/(1/0.5 + 1/0.6 - 2 + 1/1) = 1/(2 + 1.667 - 2 + 1) = 1/2.667
        expected_eps = 1.0 / (1.0 / 0.5 + 1.0 / 0.6 - 1.0)
        assert result["eps_eff"] == pytest.approx(expected_eps, rel=1e-10)

    def test_f12_less_than_1_reduces_coupling(self):
        """F12 < 1 should reduce eps_eff and produce cooler inner layer."""
        r_f1 = MultilayerShieldCalculator.solve_temperatures(
            t_exh=1073.15, t_fluid=353.15, t_surr=353.15,
            h_in=30.0, h_out=30.0, h_gap=15.0,
            eps_in=0.4, eps_out=0.4, eps_g1=0.4, eps_g2=0.4,
            f12=1.0,
        )
        r_f085 = MultilayerShieldCalculator.solve_temperatures(
            t_exh=1073.15, t_fluid=353.15, t_surr=353.15,
            h_in=30.0, h_out=30.0, h_gap=15.0,
            eps_in=0.4, eps_out=0.4, eps_g1=0.4, eps_g2=0.4,
            f12=0.85,
        )
        assert r_f085["eps_eff"] < r_f1["eps_eff"]
        # Less radiative coupling across gap → higher delta-T (inner stays hotter,
        # outer stays cooler relative to each other)
        assert r_f085["delta_T_C"] > r_f1["delta_T_C"]

    def test_f12_eps_eff_formula(self):
        """Verify eps_eff = 1/(1/e1 + 1/e2 - 2 + 1/F12)."""
        f12 = 0.85
        result = MultilayerShieldCalculator.solve_temperatures(
            t_exh=1073.15, t_fluid=353.15, t_surr=353.15,
            h_in=30.0, h_out=30.0, h_gap=15.0,
            eps_in=0.4, eps_out=0.4, eps_g1=0.5, eps_g2=0.6,
            f12=f12,
        )
        expected = 1.0 / (1.0 / 0.5 + 1.0 / 0.6 - 2.0 + 1.0 / f12)
        assert result["eps_eff"] == pytest.approx(expected, rel=1e-10)

    def test_oxidised_vs_aluminised(self):
        """High emissivity (oxidised) should give hotter inner layer."""
        r_alu = MultilayerShieldCalculator.solve_temperatures(
            t_exh=1073.15, t_fluid=353.15, t_surr=353.15,
            h_in=30.0, h_out=30.0, h_gap=15.0,
            eps_in=0.4, eps_out=0.4, eps_g1=0.4, eps_g2=0.4,
        )
        r_oxi = MultilayerShieldCalculator.solve_temperatures(
            t_exh=1073.15, t_fluid=353.15, t_surr=353.15,
            h_in=30.0, h_out=30.0, h_gap=15.0,
            eps_in=0.85, eps_out=0.85, eps_g1=0.85, eps_g2=0.85,
        )
        assert r_oxi["t1_K"] > r_alu["t1_K"]


class TestMultilayerMeshSizes:
    """Test per-layer mesh size computation."""

    def test_both_layers_have_mesh_sizes(self):
        """Both layers should produce positive finite mesh sizes."""
        result = MultilayerShieldCalculator.mesh_sizes(
            k_metal=45.0, max_dt=15.0,
            t_exh=1073.15, t_fluid=353.15, t_surr=353.15,
            h_in=30.0, h_out=30.0, h_gap=15.0,
            eps_in=0.4, eps_out=0.4, eps_g1=0.4, eps_g2=0.4,
        )
        assert result["layer1_max_dx_mm"] > 0
        assert result["layer1_max_dx_mm"] < float("inf")
        assert result["layer2_max_dx_mm"] > 0
        assert result["layer2_max_dx_mm"] < float("inf")

    def test_inner_layer_finer_than_outer(self):
        """Inner layer (higher flux) should require finer mesh."""
        result = MultilayerShieldCalculator.mesh_sizes(
            k_metal=45.0, max_dt=15.0,
            t_exh=1073.15, t_fluid=353.15, t_surr=353.15,
            h_in=30.0, h_out=30.0, h_gap=15.0,
            eps_in=0.4, eps_out=0.4, eps_g1=0.4, eps_g2=0.4,
        )
        assert result["layer1_max_dx_mm"] < result["layer2_max_dx_mm"]

    def test_mesh_size_formula_layer1(self):
        """Verify layer1 dx = k * dT / (|q_rad_in| + |q_conv_in|) * 1000."""
        result = MultilayerShieldCalculator.mesh_sizes(
            k_metal=45.0, max_dt=15.0,
            t_exh=1073.15, t_fluid=353.15, t_surr=353.15,
            h_in=30.0, h_out=30.0, h_gap=15.0,
            eps_in=0.4, eps_out=0.4, eps_g1=0.4, eps_g2=0.4,
        )
        q1 = abs(result["q_rad_in"]) + abs(result["q_conv_in"])
        expected = 45.0 * 15.0 / q1 * 1000.0
        assert result["layer1_max_dx_mm"] == pytest.approx(expected, rel=1e-6)


# -----------------------------------------------------------------------
# Non-convergence: flagged, with last-iterate fluxes, never an inf size
# -----------------------------------------------------------------------

_ML_KW = dict(
    t_exh=1073.15, t_fluid=353.15, t_surr=353.15,
    h_in=30.0, h_out=30.0, h_gap=15.0,
    eps_in=0.4, eps_out=0.4, eps_g1=0.4, eps_g2=0.4,
)
_SL_KW = dict(
    k=45.0, max_dt=15.0, t_exh=1073.15, t_fluid=353.15, t_surr=353.15,
    h_in=30.0, h_out=30.0, eps_in=0.4, eps_out=0.4,
)


class TestNonConvergence:

    def test_multilayer_out_of_iterations_reports_last_iterate(self):
        result = MultilayerShieldCalculator.solve_temperatures(max_iter=1, **_ML_KW)
        assert result["converged"] is False
        assert result["iterations"] == 1
        for key in ("q_rad_in", "q_conv_in", "q_gap_cond", "q_gap_rad",
                    "q_rad_out", "q_conv_out"):
            assert key in result
        r1, r2 = _energy_balance_residual_multi(result)
        assert result["residual_W_m2"] == pytest.approx(max(abs(r1), abs(r2)))
        assert result["residual_W_m2"] > 0.1  # the tolerance it missed

    def test_multilayer_converged_residual_within_tolerance(self):
        result = MultilayerShieldCalculator.solve_temperatures(**_ML_KW)
        assert result["converged"] is True
        assert result["residual_W_m2"] < 0.1

    def test_multilayer_unconverged_sizes_are_finite(self):
        result = MultilayerShieldCalculator.mesh_sizes(
            k_metal=45.0, max_dt=15.0, max_iter=1, **_ML_KW)
        assert result["converged"] is False
        assert 0 < result["layer1_max_dx_mm"] < float("inf")
        assert 0 < result["layer2_max_dx_mm"] < float("inf")

    def test_single_layer_mesh_size_forwards_max_iter(self):
        result = SingleLayerShieldCalculator.mesh_size(max_iter=1, **_SL_KW)
        assert result["converged"] is False
        assert result["iterations"] == 1
        assert result["residual_W_m2"] > 0
        assert 0 < result["max_dx_mm"] < float("inf")

    def test_single_layer_mesh_size_forwards_tol(self):
        result = SingleLayerShieldCalculator.mesh_size(tol=1000.0, **_SL_KW)
        assert result["converged"] is True
        assert result["iterations"] == 1

    def test_single_layer_default_still_converges(self):
        result = SingleLayerShieldCalculator.mesh_size(**_SL_KW)
        assert result["converged"] is True
        assert result["residual_W_m2"] < 1.0


class TestLayerFluxesReported:

    def test_mesh_sizes_report_the_driving_fluxes(self):
        result = MultilayerShieldCalculator.mesh_sizes(
            k_metal=45.0, max_dt=15.0, **_ML_KW)
        assert result["q_layer1"] == pytest.approx(
            abs(result["q_rad_in"]) + abs(result["q_conv_in"]))
        assert result["q_layer2"] == pytest.approx(
            abs(result["q_rad_out"]) + abs(result["q_conv_out"]))
        assert result["layer1_max_dx_mm"] == pytest.approx(
            45.0 * 15.0 / result["q_layer1"] * 1000.0)
