"""
Tests for convection zone lookup and spatial gradient estimation.
"""

import pytest
from thermal_mesh_calculators.zones import (
    CONVECTION_ZONES,
    get_zone,
    list_zones,
    get_conservative_h,
    get_zone_air_temp,
    estimate_spatial_gradient,
    SPATIAL_GRADIENT_DEFAULTS,
)


class TestZoneLookup:

    def test_all_zones_have_required_keys(self):
        required = {"h_low", "h_high", "regime", "velocity_ms",
                    "t_air_C_low", "t_air_C_high", "orientation",
                    "is_internal", "notes"}
        for name, zone in CONVECTION_ZONES.items():
            missing = required - set(zone.keys())
            assert not missing, f"{name} missing keys: {missing}"
            assert zone["h_high"] >= zone["h_low"], (
                f"{name}: h_high < h_low"
            )

    def test_get_zone_returns_copy(self):
        z = get_zone("cooling_pack_downstream")
        z["h_high"] = 9999
        assert CONVECTION_ZONES["cooling_pack_downstream"]["h_high"] != 9999

    def test_unknown_zone_raises(self):
        with pytest.raises(KeyError, match="Unknown convection zone"):
            get_zone("nonexistent_zone")

    def test_list_zones_sorted(self):
        zones = list_zones()
        assert zones == sorted(zones)

    def test_conservative_h_returns_high(self):
        for name in CONVECTION_ZONES:
            h = get_conservative_h(name)
            assert h == CONVECTION_ZONES[name]["h_high"]

    def test_natural_convection_low_h(self):
        """Natural convection zones should have h_high <= 25."""
        for name, zone in CONVECTION_ZONES.items():
            if zone["regime"] == "natural":
                assert zone["h_high"] <= 25.0, (
                    f"{name} natural conv h_high={zone['h_high']} > 25"
                )

    def test_forced_convection_higher_h(self):
        """Forced convection zones should have h_high > 20."""
        for name, zone in CONVECTION_ZONES.items():
            if zone["regime"] == "forced":
                assert zone["h_high"] > 20.0, (
                    f"{name} forced conv h_high={zone['h_high']} <= 20"
                )

    def test_external_zones_capped_at_100(self):
        """All external zones should have h_high <= 100 W/m^2K."""
        for name, zone in CONVECTION_ZONES.items():
            if not zone["is_internal"]:
                assert zone["h_high"] <= 100.0, (
                    f"{name} external h_high={zone['h_high']} > 100"
                )

    def test_natural_zones_zero_velocity(self):
        """Natural convection zones should have velocity_ms = 0."""
        for name, zone in CONVECTION_ZONES.items():
            if zone["regime"] == "natural":
                assert zone["velocity_ms"] == 0.0, (
                    f"{name} natural zone has velocity={zone['velocity_ms']}"
                )

    def test_legacy_aliases_exist(self):
        """Backward compat aliases should work."""
        assert "engine_bay_dead_zone" in CONVECTION_ZONES
        assert "near_exhaust_natural" in CONVECTION_ZONES


class TestZoneAirTemp:

    def test_returns_kelvin(self):
        t = get_zone_air_temp("cooling_pack_downstream", "high")
        assert t > 273.15  # must be in Kelvin

    def test_high_above_low(self):
        t_low = get_zone_air_temp("engine_above", "low")
        t_high = get_zone_air_temp("engine_above", "high")
        assert t_high >= t_low

    def test_mid_between_bounds(self):
        t_low = get_zone_air_temp("engine_above", "low")
        t_high = get_zone_air_temp("engine_above", "high")
        t_mid = get_zone_air_temp("engine_above", "mid")
        assert t_low <= t_mid <= t_high


class TestSpatialGradientEstimation:

    def test_exhaust_returns_max(self):
        grad = estimate_spatial_gradient(k=45.0, component_class="exhaust")
        assert grad == SPATIAL_GRADIENT_DEFAULTS["exhaust_max"]

    def test_exhaust_adjacent_returns_typical(self):
        grad = estimate_spatial_gradient(k=45.0, component_class="exhaust_adjacent")
        assert grad == SPATIAL_GRADIENT_DEFAULTS["exhaust_typical"]

    def test_structural_derives_from_flux(self):
        grad = estimate_spatial_gradient(
            k=45.0, q_total=9000.0, component_class="structural",
        )
        assert grad == pytest.approx(9000.0 / 45.0)

    def test_structural_zero_flux(self):
        grad = estimate_spatial_gradient(
            k=45.0, q_total=0.0, component_class="structural",
        )
        assert grad == 0.0

    def test_structural_no_flux(self):
        grad = estimate_spatial_gradient(
            k=45.0, q_total=None, component_class="structural",
        )
        assert grad == 0.0
