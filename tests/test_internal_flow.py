"""
Internal-flow heat transfer in ducts: the analytic benchmark, the advisories
and the estimate_h(duct=...) wiring.  test_examples_regression.py holds every
existing result where it was.

The benchmark values are hand-computed from the correlations as published,
one per correlation at a stated Re and Pr, with each step of the arithmetic
written beside the assertion.  They are not the code's own output.
"""

import pytest

from thermal_mesh_calculators import (
    NU_LAMINAR_CONSTANT_HEAT_FLUX,
    NU_LAMINAR_CONSTANT_TEMPERATURE,
    estimate_h,
    friction_factor_petukhov,
    hydraulic_diameter,
    internal_flow_h,
    internal_flow_nusselt,
    nusselt_dittus_boelter,
    nusselt_gnielinski,
    nusselt_hausen,
)
from thermal_mesh_calculators.h_estimator import H_EXTERNAL_MAX

REL = 1e-9   # the hand values below carry ten or more significant digits


# =============================================================================
#  The analytic benchmark: one hand-computed value per correlation
# =============================================================================

class TestBenchmark:

    def test_laminar_fully_developed_constant_wall_temperature(self):
        # Fully developed laminar flow in a circular tube with the wall at
        # one temperature: Nu = 3.66, whatever Re (below 2300) and Pr.
        out = internal_flow_nusselt(1000.0, 0.71, "constant_temperature")
        assert NU_LAMINAR_CONSTANT_TEMPERATURE == 3.66
        assert out["Nu"] == 3.66
        assert out["regime"] == "laminar"
        assert out["correlation"] == "laminar_fully_developed"

    def test_laminar_fully_developed_constant_heat_flux(self):
        # The same flow with a uniform wall heat flux: Nu = 4.36.
        out = internal_flow_nusselt(1000.0, 0.71, "constant_heat_flux")
        assert NU_LAMINAR_CONSTANT_HEAT_FLUX == 4.36
        assert out["Nu"] == 4.36
        assert out["correlation"] == "laminar_fully_developed"

    def test_hausen_thermal_entry(self):
        # Re = 1000, Pr = 0.71, L/D = 20, wall at one temperature.
        #   Gz = Re Pr D/L = 1000 * 0.71 / 20                 = 35.5
        #   Gz^(2/3) = exp((2/3) ln 35.5)                      = 10.8015372755
        #   0.0668 Gz                                          = 2.3714
        #   1 + 0.04 Gz^(2/3)                                  = 1.4320614910
        #   Nu = 3.66 + 2.3714 / 1.4320614910                  = 5.3159344797
        assert nusselt_hausen(1000.0, 0.71, 20.0) == pytest.approx(
            5.3159344797, rel=REL)
        out = internal_flow_nusselt(1000.0, 0.71, "constant_temperature",
                                    length_over_diameter=20.0)
        assert out["Nu"] == pytest.approx(5.3159344797, rel=REL)
        assert out["correlation"] == "hausen"

    def test_hausen_tends_to_the_fully_developed_value(self):
        # Gz -> 0 for a long duct: Nu -> 3.66 from above.
        assert nusselt_hausen(1000.0, 0.71, 1.0e9) == pytest.approx(
            3.66, rel=1e-6)
        assert nusselt_hausen(1000.0, 0.71, 1.0e9) > 3.66

    def test_petukhov_friction_factor(self):
        # Re = 20000:
        #   ln Re                                              = 9.9034875525
        #   0.790 ln Re - 1.64                                 = 6.1837551665
        #   f = 6.1837551665^-2                                = 0.02615142914593
        assert friction_factor_petukhov(20000.0) == pytest.approx(
            0.02615142914593, rel=REL)

    def test_gnielinski(self):
        # Re = 20000, Pr = 0.71, with the Petukhov f above:
        #   f/8                                                = 0.0032689286
        #   (f/8)(Re - 1000) Pr = 0.0032689286 * 19000 * 0.71  = 44.0978473973
        #   sqrt(f/8)                                          = 0.0571745454
        #   Pr^(2/3) = exp((2/3) ln 0.71)                      = 0.7958640711
        #   1 + 12.7 * 0.0571745454 * (0.7958640711 - 1)       = 0.8517734875
        #   Nu = 44.0978473973 / 0.8517734875                  = 51.7718008879
        assert nusselt_gnielinski(20000.0, 0.71) == pytest.approx(
            51.7718008879, rel=REL)
        out = internal_flow_nusselt(20000.0, 0.71, "constant_temperature")
        assert out["Nu"] == pytest.approx(51.7718008879, rel=REL)
        assert out["regime"] == "turbulent"
        assert out["correlation"] == "gnielinski"
        assert out["friction_factor"] == pytest.approx(0.02615142914593, rel=REL)

    def test_dittus_boelter_heating(self):
        # Re = 20000, Pr = 0.71, the fluid heated (n = 0.4):
        #   Re^0.8 = exp(0.8 ln 20000)                         = 2759.4593229224
        #   Pr^0.4 = exp(0.4 ln 0.71)                          = 0.8719736063
        #   Nu = 0.023 * 2759.4593229224 * 0.8719736063        = 55.3420410343
        assert nusselt_dittus_boelter(20000.0, 0.71, heating=True) == \
            pytest.approx(55.3420410343, rel=REL)
        out = internal_flow_nusselt(20000.0, 0.71, "constant_temperature",
                                    correlation="dittus_boelter")
        assert out["Nu"] == pytest.approx(55.3420410343, rel=REL)
        assert out["correlation"] == "dittus_boelter"
        assert out["friction_factor"] is None

    def test_dittus_boelter_cooling(self):
        # The fluid cooled (n = 0.3):
        #   Pr^0.3 = exp(0.3 ln 0.71)                          = 0.9023551568
        #   Nu = 0.023 * 2759.4593229224 * 0.9023551568        = 57.2702840534
        assert nusselt_dittus_boelter(20000.0, 0.71, heating=False) == \
            pytest.approx(57.2702840534, rel=REL)

    def test_hydraulic_diameter_of_a_rectangle(self):
        # 40 mm x 20 mm: A = 8.0e-4 m^2, P = 0.12 m, Dh = 4A/P = 0.0266667 m
        # (= 2ab/(a+b) = 2 * 40 * 20 / 60 mm).
        assert hydraulic_diameter(8.0e-4, 0.12) == pytest.approx(
            0.02666666666666667, rel=1e-15)

    def test_hydraulic_diameter_of_a_circle_is_its_diameter(self):
        d = 0.03
        area = 3.141592653589793 * d * d / 4.0
        assert hydraulic_diameter(area, 3.141592653589793 * d) == \
            pytest.approx(d, rel=1e-15)

    def test_air_in_a_duct_end_to_end_turbulent(self):
        # D = 0.025 m, V = 10 m/s, Ts = 350 K, Tf = 300 K: T_film = 325 K.
        #   k_air = 0.0241 + 7.0e-5 * 25                       = 0.02585 W/m K
        #   nu = 1.5e-5 * (325/300)^1.7                        = 1.7186476149e-5
        #   Re = 10 * 0.025 / 1.7186476149e-5                  = 14546.321063
        #   f (Petukhov)                                       = 0.0284161312
        #   Nu (Gnielinski, Pr = 0.71)                         = 40.4061006928
        #   h = 40.4061006928 * 0.02585 / 0.025                = 41.7799081164
        out = internal_flow_h(10.0, 0.025, 350.0, 300.0, "constant_temperature")
        assert out["Re"] == pytest.approx(14546.321063, rel=1e-9)
        assert out["Pr"] == 0.71
        assert out["Nu"] == pytest.approx(40.4061006928, rel=REL)
        assert out["h"] == pytest.approx(41.7799081164, rel=REL)
        assert out["h_raw"] == out["h"]
        assert out["t_film_K"] == 325.0 and out["film_in_range"] is True

    def test_air_in_a_duct_end_to_end_laminar_entry(self):
        # D = 0.010 m, V = 1.5 m/s, L = 0.5 m, same temperatures:
        #   Re = 1.5 * 0.010 / 1.7186476149e-5                 = 872.779264
        #   Gz = Re Pr D/L = 872.779264 * 0.71 * 0.010 / 0.5   = 12.39346555
        #   Nu (Hausen)                                        = 4.3418246519
        #   h = 4.3418246519 * 0.02585 / 0.010                 = 11.2236167251
        out = internal_flow_h(1.5, 0.010, 350.0, 300.0, "constant_temperature",
                              length=0.5)
        assert out["Re"] == pytest.approx(872.779264, rel=1e-9)
        assert out["Nu"] == pytest.approx(4.3418246519, rel=REL)
        assert out["h"] == pytest.approx(11.2236167251, rel=REL)
        # Fully developed (no length): h = 3.66 * 0.02585 / 0.010 = 9.4611.
        fd = internal_flow_h(1.5, 0.010, 350.0, 300.0, "constant_temperature")
        assert fd["h"] == pytest.approx(9.4611, rel=REL)
        # Uniform heat flux: h = 4.36 * 0.02585 / 0.010 = 11.2706.
        q = internal_flow_h(1.5, 0.010, 350.0, 300.0, "constant_heat_flux")
        assert q["h"] == pytest.approx(11.2706, rel=REL)


# =============================================================================
#  The transition band, the advisories and the input checks
# =============================================================================

class TestTransitionAndAdvisories:

    def test_transition_band_takes_the_higher_value_and_says_so(self):
        # Re = 5000, Pr = 0.71: Gnielinski gives 16.7220171213 (f = 0.0386194727,
        # f/8 = 0.0048274341, numerator 13.7099127932, denominator 0.8198719505);
        # the fully developed laminar value is 3.66.  No blend: the higher.
        out = internal_flow_nusselt(5000.0, 0.71, "constant_temperature")
        assert out["regime"] == "transitional"
        assert out["Nu_laminar"] == 3.66
        assert out["Nu_turbulent"] == pytest.approx(16.7220171213, rel=REL)
        assert out["Nu"] == out["Nu_turbulent"]
        assert out["correlation"] == "gnielinski"
        assert "transitional" in out["advisory"].lower()
        assert "3.66" in out["advisory"] and "16.7" in out["advisory"]

    @pytest.mark.parametrize("reynolds", [2300.0, 10000.0])
    def test_band_edges(self, reynolds):
        out = internal_flow_nusselt(reynolds, 0.71, "constant_temperature")
        assert out["regime"] == ("laminar" if reynolds == 2300.0
                                 else "turbulent")

    def test_gnielinski_below_its_range_is_named(self):
        out = internal_flow_nusselt(2600.0, 0.71, "constant_temperature")
        assert out["regime"] == "transitional"
        assert "3000" in out["advisory"]

    def test_dittus_boelter_in_the_band_is_named_out_of_range(self):
        out = internal_flow_nusselt(5000.0, 0.71, "constant_temperature",
                                    correlation="dittus_boelter")
        assert out["regime"] == "transitional"
        assert "10000" in out["advisory"]

    def test_short_turbulent_duct_is_named(self):
        out = internal_flow_nusselt(20000.0, 0.71, "constant_temperature",
                                    length_over_diameter=5.0)
        assert out["Nu"] == pytest.approx(51.7718008879, rel=REL)
        assert "10 hydraulic diameters" in out["advisory"]

    def test_heat_flux_entry_is_not_modelled_and_says_so(self):
        out = internal_flow_nusselt(1000.0, 0.71, "constant_heat_flux",
                                    length_over_diameter=20.0)
        assert out["Nu"] == 4.36
        assert "entry" in out["advisory"]

    def test_quiet_when_inside_every_range(self):
        out = internal_flow_nusselt(20000.0, 0.71, "constant_temperature",
                                    length_over_diameter=50.0)
        assert out["advisory"] == ""

    def test_film_temperature_out_of_range_is_reported(self):
        out = internal_flow_h(10.0, 0.025, 900.0, 600.0, "constant_temperature")
        assert out["film_in_range"] is False
        assert "700 K" in out["advisory"]

    def test_internal_flow_is_uncapped_unless_asked(self):
        # A small, fast duct: h well above the external cap.
        out = internal_flow_h(40.0, 0.01, 350.0, 300.0, "constant_temperature")
        assert out["h"] == out["h_raw"] > H_EXTERNAL_MAX
        capped = internal_flow_h(40.0, 0.01, 350.0, 300.0,
                                 "constant_temperature", cap=50.0)
        assert capped["h"] == 50.0 and capped["h_raw"] == out["h_raw"]

    @pytest.mark.parametrize("kwargs, words", [
        ({"wall_condition": "isothermal"}, "constant_temperature"),
        ({"correlation": "colburn"}, "gnielinski"),
        ({"velocity": 0.0}, "velocity"),
        ({"hydraulic_diameter": -0.01}, "hydraulic_diameter"),
        ({"length": 0.0}, "length"),
    ])
    def test_bad_inputs_raise_and_name_the_rule(self, kwargs, words):
        args = dict(velocity=10.0, hydraulic_diameter=0.025, t_surf=350.0,
                    t_fluid=300.0, wall_condition="constant_temperature")
        args.update(kwargs)
        with pytest.raises(ValueError, match=words):
            internal_flow_h(**args)


# =============================================================================
#  estimate_h(duct=...): the only way in, and nothing else moves
# =============================================================================

DUCT = {"length": 0.5, "wall_condition": "constant_temperature"}


class TestEstimateHDuct:

    def test_duct_uses_the_internal_correlations(self):
        out = estimate_h(10.0, 350.0, 300.0, 0.025, duct=DUCT)
        direct = internal_flow_h(10.0, 0.025, 350.0, 300.0,
                                 "constant_temperature", length=0.5)
        assert out["h"] == direct["h"]
        assert out["duct"] == direct
        assert out["regime"] == "forced"
        assert out["forced_regime"] == "turbulent"
        assert out["Re"] == direct["Re"]

    def test_duct_result_keeps_the_estimate_h_keys(self):
        plain = estimate_h(10.0, 350.0, 300.0, 0.025)
        out = estimate_h(10.0, 350.0, 300.0, 0.025, duct=DUCT)
        assert set(out) == set(plain) | {"duct"}

    def test_a_transitional_duct_gets_a_caution(self):
        out = estimate_h(3.0, 350.0, 300.0, 0.025,
                         duct={"wall_condition": "constant_temperature"})
        assert out["forced_regime"] == "transitional"
        assert out["solver_advisory"]["severity"] == "caution"

    @pytest.mark.parametrize("duct, words", [
        ({"length": 0.5}, "wall_condition"),
        ({"wall_condition": "constant_temperature", "diameter": 0.02},
         "diameter"),
        ("round", "dict"),
    ])
    def test_a_bad_duct_raises(self, duct, words):
        with pytest.raises((TypeError, ValueError), match=words):
            estimate_h(10.0, 350.0, 300.0, 0.025, duct=duct)

    def test_the_external_cap_is_refused_with_a_duct(self):
        with pytest.raises(ValueError, match="cap"):
            estimate_h(10.0, 350.0, 300.0, 0.025, cap=50.0, duct=DUCT)
