"""
Tests for BoundaryDrivenConductionCalculator.size_wall
======================================================

Verification strategy:
    - The two defect cases that motivated the method (85.8 mm cell for a
      1 mm aluminium shield layer; 19.7 mm for the steel layer) now clamp
      to one cell + a thermally-thin verdict, with the raw dx retained.
    - The genuinely-resolve case (glass-filled PA66 intake) reports a
      correct cell count and Biot number.
    - The pole near zero net flux fires exactly at the convection/radiation
      crossing and NOT slightly off it (guard is relative, not absolute).
    - The physical-consistency check fires when the implied through-wall
      drop exceeds the surface-to-fluid driving difference.
    - Clamp arithmetic, cap, flag, and input validation.
"""

import math
import pytest
from thermal_mesh_calculators.constants import STEFAN_BOLTZMANN
from thermal_mesh_calculators.conduction import (
    BoundaryDrivenConductionCalculator,
    NEAR_EQUILIBRIUM_FRACTION,
    STEEP_BIOT,
    THIN_BIOT,
)


# Shared shield-case boundary conditions (400 C surface, 100 C fluid/sink,
# lightly-oxidised aluminized surface): q_conv = 25*300 = 7500 W/m^2,
# q_rad = 0.35*sigma*(673.15^4 - 373.15^4) ~ 3920 W/m^2, q ~ 11.4 kW/m^2.
SHIELD_BC = dict(
    h=25.0, t_surf=673.15, t_fluid=373.15,
    epsilon=0.35, t_surr=373.15, max_dt=5.0,
)


class TestShieldDefectCases:
    """The 86x / 20x raw-dx defect cases from the original README."""

    def test_aluminium_layer_clamps_to_one_cell(self):
        """k=196, 1 mm layer: raw dx ~ 85.8 mm (86x the part) -> 1 cell."""
        r = BoundaryDrivenConductionCalculator.size_wall(
            k=196.0, thickness_m=0.001, **SHIELD_BC
        )
        assert r["n_cells"] == 1
        assert r["regime"] == "thermally_thin"
        assert r["dx_exceeds_part"] is True
        assert 80.0 < r["dx_raw_mm"] < 92.0
        assert r["dx_used_mm"] == pytest.approx(1.0)
        assert r["biot"] < 0.001
        assert any("not a cell size" in w for w in r["warnings"])

    def test_steel_layer_clamps_to_one_cell(self):
        """k=45, 1 mm layer: raw dx ~ 19.7 mm (20x the part) -> 1 cell."""
        r = BoundaryDrivenConductionCalculator.size_wall(
            k=45.0, thickness_m=0.001, **SHIELD_BC
        )
        assert r["n_cells"] == 1
        assert r["regime"] == "thermally_thin"
        assert 18.0 < r["dx_raw_mm"] < 21.5

    def test_flux_matches_max_mesh_size(self):
        """size_wall composes max_mesh_size — flux components identical."""
        base = BoundaryDrivenConductionCalculator.max_mesh_size(
            k=196.0, **SHIELD_BC
        )
        r = BoundaryDrivenConductionCalculator.size_wall(
            k=196.0, thickness_m=0.001, **SHIELD_BC
        )
        assert r["q_conv"] == base["q_conv"]
        assert r["q_rad"] == base["q_rad"]
        assert r["q_total"] == base["q_total"]
        assert r["dx_raw_mm"] == pytest.approx(base["max_dx_mm"])


class TestResolveCase:
    """Glass-filled PA66 intake — the case that genuinely needs cells."""

    def test_pa66_intake_resolves(self):
        """
        k=0.4, L=3 mm, Ts=120 C, Tf=25 C, h=25, eps=0.9:
            q_conv = 25*95 = 2375 W/m^2
            q_rad  = 0.9*sigma*(393.15^4 - 298.15^4) ~ 816 W/m^2
            q ~ 3191 W/m^2; dt_wall = q*0.003/0.4 ~ 23.9 K
            h_eff = q/95 ~ 33.6; Bi = 33.6*0.003/0.4 ~ 0.252 -> resolve
            dx_raw = 0.4*5/q ~ 0.627 mm -> N = ceil(3/0.627) = 5
        """
        r = BoundaryDrivenConductionCalculator.size_wall(
            k=0.4, h=25.0, t_surf=393.15, t_fluid=298.15,
            epsilon=0.9, t_surr=298.15, max_dt=5.0, thickness_m=0.003,
        )
        assert r["regime"] == "resolve"
        assert r["n_cells"] == 5
        assert r["biot"] == pytest.approx(0.252, rel=0.05)
        assert 20.0 < r["dt_wall_K"] < 28.0
        assert r["dx_exceeds_part"] is False


class TestPoleGuard:
    """Near-zero net flux is a pole in the relation, not in the physics."""

    @staticmethod
    def _crossing_t_surr(t_surf, t_fluid, h, epsilon):
        """Sink temperature at which radiation exactly cancels convection."""
        q_conv = h * (t_surf - t_fluid)
        t4 = t_surf**4 + q_conv / (epsilon * STEFAN_BOLTZMANN)
        return t4**0.25

    def test_fires_at_the_crossing(self):
        """Warm part in a hotter enclosure: net flux -> 0 -> single cell."""
        t_surr = self._crossing_t_surr(350.0, 300.0, 25.0, 0.9)
        r = BoundaryDrivenConductionCalculator.size_wall(
            k=45.0, h=25.0, t_surf=350.0, t_fluid=300.0,
            epsilon=0.9, t_surr=t_surr, max_dt=5.0, thickness_m=0.004,
        )
        assert r["regime"] == "near_equilibrium"
        assert r["n_cells"] == 1
        assert math.isinf(r["dx_raw_mm"])
        assert any("pole" in w for w in r["warnings"])

    def test_does_not_fire_off_the_crossing(self):
        """30 K past the crossing the net is >2% of the terms — no guard."""
        t_surr = self._crossing_t_surr(350.0, 300.0, 25.0, 0.9)
        r = BoundaryDrivenConductionCalculator.size_wall(
            k=45.0, h=25.0, t_surf=350.0, t_fluid=300.0,
            epsilon=0.9, t_surr=t_surr + 30.0, max_dt=5.0,
            thickness_m=0.004,
        )
        assert r["regime"] != "near_equilibrium"
        assert math.isfinite(r["dx_raw_mm"])

    def test_zero_flux_exact(self):
        """Ts == Tf == Tsurr: scale == 0 -> near-equilibrium path."""
        r = BoundaryDrivenConductionCalculator.size_wall(
            k=45.0, h=25.0, t_surf=353.15, t_fluid=353.15,
            epsilon=0.9, t_surr=353.15, max_dt=5.0, thickness_m=0.004,
        )
        assert r["regime"] == "near_equilibrium"
        assert r["n_cells"] == 1


class TestConsistencyCheck:
    """Implied through-wall drop exceeding the driving difference."""

    def test_insulation_batt_flags_inconsistency(self):
        """
        10 mm ceramic-fibre batt (k=0.08) under the shield BC:
            dt_wall = 11420*0.01/0.08 ~ 1428 K >> 300 K driving
            Bi ~ 4.8 -> steep gradient + PHYSICALLY INCONSISTENT warning.
        """
        r = BoundaryDrivenConductionCalculator.size_wall(
            k=0.08, thickness_m=0.010, **SHIELD_BC
        )
        assert r["regime"] == "steep_gradient"
        assert r["biot"] > STEEP_BIOT
        assert r["dt_wall_K"] > 1000.0
        assert any("PHYSICALLY INCONSISTENT" in w for w in r["warnings"])

    def test_thin_wall_does_not_flag(self):
        """The 4 mm cast-iron manifold wall is consistent (7-ish K drop)."""
        r = BoundaryDrivenConductionCalculator.size_wall(
            k=45.0, h=25.0, t_surf=1073.15, t_fluid=373.15,
            epsilon=0.9, t_surr=298.15, max_dt=5.0, thickness_m=0.004,
        )
        assert not any("PHYSICALLY INCONSISTENT" in w for w in r["warnings"])
        assert r["regime"] == "thermally_thin"


class TestClampArithmetic:
    """N_cells = max(1, ceil(L / dx_raw)), capped at max_cells."""

    def test_manifold_two_cells(self):
        """
        800 C manifold, 4 mm wall, max_dt=5:
            q ~ 84.8 kW/m^2 (h=25, eps=0.9, Tf=100 C, Tsurr=25 C)
            dx_raw = 45*5/84800 ~ 2.65 mm -> N = ceil(4/2.65) = 2
        Thermally thin (Bi ~ 0.011) but the count is still reported.
        """
        r = BoundaryDrivenConductionCalculator.size_wall(
            k=45.0, h=25.0, t_surf=1073.15, t_fluid=373.15,
            epsilon=0.9, t_surr=298.15, max_dt=5.0, thickness_m=0.004,
        )
        assert r["n_cells"] == 2
        assert r["regime"] == "thermally_thin"
        assert r["biot"] == pytest.approx(0.011, abs=0.004)
        assert r["dx_used_mm"] == pytest.approx(2.0)

    def test_max_cells_cap(self):
        """Absurdly tight tolerance hits the cap with a warning."""
        r = BoundaryDrivenConductionCalculator.size_wall(
            k=0.08, h=25.0, t_surf=673.15, t_fluid=373.15,
            epsilon=0.35, t_surr=373.15, max_dt=0.05, thickness_m=0.010,
            max_cells=100,
        )
        assert r["n_cells"] == 100
        assert any("capped at 100" in w for w in r["warnings"])

    def test_radiation_dominance_warning(self):
        """800 C manifold: radiation ~79% of flux -> emissivity warning."""
        r = BoundaryDrivenConductionCalculator.size_wall(
            k=45.0, h=25.0, t_surf=1073.15, t_fluid=373.15,
            epsilon=0.9, t_surr=298.15, max_dt=5.0, thickness_m=0.004,
        )
        assert r["rad_fraction"] > 0.6
        assert any("emissivity" in w for w in r["warnings"])


class TestValidation:
    """Input validation raises rather than returning nonsense."""

    def test_zero_thickness_raises(self):
        with pytest.raises(ValueError):
            BoundaryDrivenConductionCalculator.size_wall(
                k=45.0, thickness_m=0.0, **SHIELD_BC
            )

    def test_negative_thickness_raises(self):
        with pytest.raises(ValueError):
            BoundaryDrivenConductionCalculator.size_wall(
                k=45.0, thickness_m=-0.001, **SHIELD_BC
            )

    def test_zero_k_raises(self):
        with pytest.raises(ValueError):
            BoundaryDrivenConductionCalculator.size_wall(
                k=0.0, thickness_m=0.001, **SHIELD_BC
            )

    def test_zero_max_dt_raises(self):
        with pytest.raises(ValueError):
            BoundaryDrivenConductionCalculator.size_wall(
                k=45.0, h=25.0, t_surf=673.15, t_fluid=373.15,
                epsilon=0.35, t_surr=373.15, max_dt=0.0, thickness_m=0.001,
            )

    def test_regime_thresholds_are_module_constants(self):
        """The thresholds are importable, documented policy — not literals."""
        assert THIN_BIOT == 0.1
        assert STEEP_BIOT == 1.0
        assert 0.0 < NEAR_EQUILIBRIUM_FRACTION < 0.1
