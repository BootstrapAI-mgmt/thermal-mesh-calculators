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

import itertools
import math
import re

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
        """A window no element fits is flagged, with advice that closes it.

        Explicit scheme, dt = 0.5 s, tau_bc = 0.5 s: the stability minimum
        sqrt(alpha*dt/0.5) = 3.397 mm exceeds the drive-cycle bound
        sqrt(alpha*tau_bc) = 2.402 mm.  The conflict is there by
        construction (asserted first, unconditionally), and reducing dt to
        Fo_max * tau_bc = 0.25 s must close it.
        """
        result = TransientMeshCalculator.combined_transient_limits(
            **STEEL, dt=0.5, fo_max=0.5, safety_factor=1.0, tau_bc=0.5,
        )
        assert result["fourier_min_dx_mm"] > result["drive_cycle_max_dx_mm"]
        assert result["feasible"] is False
        assert result["binding_constraint"].startswith("fourier_stability")
        assert "WARNING" in result["binding_constraint"]
        assert result["recommended_dx_mm"] == result["fourier_min_dx_mm"]
        [remedy] = result["conflict"]["remedies"]
        assert remedy == {"parameter": "dt", "max_value": pytest.approx(0.25)}

        fixed = TransientMeshCalculator.combined_transient_limits(
            **STEEL, dt=remedy["max_value"], fo_max=0.5, safety_factor=1.0,
            tau_bc=0.5,
        )
        assert fixed["feasible"] is True
        assert fixed["conflict"] is None
        assert fixed["binding_constraint"] == "drive_cycle"

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
            "alpha", "scheme", "penetration_max_dx_mm", "penetration_applied",
            "fourier_min_dx_mm", "drive_cycle_max_dx_mm", "max_dx_mm",
            "feasible", "recommended_dx_mm", "binding_constraint", "conflict",
            "advice",
        }
        assert set(result.keys()) == expected_keys


# The explicit defaults (Fo <= 0.5, C = 1) against the time steps they were
# once infeasible at, whatever the step.
STEEL_MILD = {"k": 54.0, "rho": 7833.0, "cp": 465.0}


class TestExplicitPath:
    """The explicit path is feasible; its conflicts carry remedies that work."""

    @pytest.mark.parametrize("dt", [0.01, 1.0, 100.0])
    def test_explicit_defaults_are_feasible_at_every_step(self, dt):
        result = TransientMeshCalculator.combined_transient_limits(
            **STEEL_MILD, dt=dt,
        )
        alpha = 54.0 / (7833.0 * 465.0)
        assert result["scheme"] == "explicit"
        assert result["penetration_applied"] is False
        assert result["feasible"] is True
        assert result["conflict"] is None
        # No transient upper bound applies: stability sets only a minimum.
        assert result["max_dx_mm"] == float("inf")
        assert result["binding_constraint"] == "none"
        assert result["fourier_min_dx_mm"] == pytest.approx(
            math.sqrt(alpha * dt / 0.5) * 1000.0, rel=1e-12)

    @pytest.mark.parametrize("dt", [0.01, 1.0, 100.0])
    def test_the_two_sqrt_dt_bounds_keep_their_ratio(self, dt):
        """Fourier minimum / penetration = 1 / (C sqrt(Fo_max)) at any dt,
        which is why "reduce dt" could never resolve that pair."""
        result = TransientMeshCalculator.combined_transient_limits(
            **STEEL_MILD, dt=dt,
        )
        ratio = result["fourier_min_dx_mm"] / result["penetration_max_dx_mm"]
        assert ratio == pytest.approx(1.0 / math.sqrt(0.5), rel=1e-12)

    def test_explicit_upper_bound_is_the_drive_cycle(self):
        result = TransientMeshCalculator.combined_transient_limits(
            **STEEL, dt=0.5, fo_max=0.5, safety_factor=1.0, tau_bc=20.0,
        )
        assert result["feasible"] is True
        assert result["binding_constraint"] == "drive_cycle"
        assert result["recommended_dx_mm"] == result["drive_cycle_max_dx_mm"]
        assert (result["fourier_min_dx_mm"] <= result["recommended_dx_mm"]
                < float("inf"))

    @pytest.mark.parametrize("dt", [0.01, 1.0, 100.0])
    def test_drive_cycle_conflict_remedy_closes_it(self, dt):
        """dt <= Fo_max * tau_bc, independent of the material."""
        tau_bc = dt / 10.0  # the stability minimum is sqrt(20) x the bound
        result = TransientMeshCalculator.combined_transient_limits(
            **STEEL_MILD, dt=dt, tau_bc=tau_bc,
        )
        assert result["feasible"] is False
        [remedy] = result["conflict"]["remedies"]
        assert remedy["parameter"] == "dt"
        assert remedy["max_value"] == pytest.approx(0.5 * tau_bc)
        fixed = TransientMeshCalculator.combined_transient_limits(
            **STEEL_MILD, dt=remedy["max_value"], tau_bc=tau_bc,
        )
        assert fixed["feasible"] is True

    def test_implicit_window_empty_at_every_dt_needs_a_parameter(self):
        """Implicit, Fo_max = 0.5, C = 1: the window is empty whatever dt is.
        The remedy is C >= 1/sqrt(Fo_max); halving dt changes nothing."""
        for dt in (0.01, 1.0, 100.0):
            result = TransientMeshCalculator.combined_transient_limits(
                **STEEL_MILD, dt=dt, scheme="implicit",
            )
            assert result["feasible"] is False
            assert result["binding_constraint"].startswith("fourier_accuracy")
            [remedy] = result["conflict"]["remedies"]
            assert remedy["parameter"] == "safety_factor"
            assert remedy["min_value"] == pytest.approx(math.sqrt(2.0))
            halved = TransientMeshCalculator.combined_transient_limits(
                **STEEL_MILD, dt=dt / 2.0, scheme="implicit",
            )
            assert halved["feasible"] is False
            fixed = TransientMeshCalculator.combined_transient_limits(
                **STEEL_MILD, dt=dt, scheme="implicit",
                safety_factor=remedy["min_value"],
            )
            assert fixed["feasible"] is True
            assert fixed["binding_constraint"] == "penetration_depth"

    def test_scheme_is_inferred_from_fo_max(self):
        explicit = TransientMeshCalculator.combined_transient_limits(
            **STEEL, dt=1.0, fo_max=0.5)
        implicit = TransientMeshCalculator.combined_transient_limits(
            **STEEL, dt=1.0, fo_max=5.0)
        assert explicit["scheme"] == "explicit"
        assert implicit["scheme"] == "implicit"
        assert implicit["penetration_applied"] is True
        assert implicit["max_dx_mm"] == implicit["penetration_max_dx_mm"]

    def test_explicit_scheme_ignores_the_penetration_bound(self):
        result = TransientMeshCalculator.combined_transient_limits(
            **STEEL, dt=1.0, fo_max=0.5, safety_factor=2.0, tau_bc=20.0,
            scheme="explicit",
        )
        assert result["penetration_applied"] is False
        assert result["max_dx_mm"] == result["drive_cycle_max_dx_mm"]

    def test_unknown_scheme_raises(self):
        with pytest.raises(ValueError, match="scheme"):
            TransientMeshCalculator.combined_transient_limits(
                **STEEL, dt=1.0, scheme="crank_nicolson")


# A remedy as printed: "reduce dt to <= 20.24 s", or "raise safety_factor to
# >= 1.415 or fo_max to >= 1.563". Steps are joined by "; then ", the
# alternatives within a step by " or ".
_PRINTED_BOUND = re.compile(
    r"\b(dt|safety_factor|fo_max) to ([<>]=) (\d[0-9.eE+-]*)")


def printed_remedy_choices(text):
    """Every way to follow the remedy printed in text, one dict of parameter
    values per choice of one alternative at each step."""
    steps = []
    for step in text.split("; then "):
        options = [{name: float(value.rstrip("."))}
                   for name, _, value in _PRINTED_BOUND.findall(step)]
        if options:
            steps.append(options)
    return [{name: value for option in choice for name, value in option.items()}
            for choice in itertools.product(*steps)]


class TestPrintedRemedies:
    """A remedy closes the conflict when it is applied as printed.

    Each bound is printed rounded toward the side that satisfies it, a maximum
    down and a minimum up. Rounded to nearest, a maximum of 20.2467 s printed
    as "dt <= 20.25 s" and the minimum sqrt(2) as "safety_factor >= 1.414",
    and applying either printed value left the conflict in place. The
    structured remedies keep the exact values.
    """

    CASES = {
        # dt <= Fo_max * tau_bc = 20.2467 s
        "explicit-drive-cycle": dict(dt=100.0, tau_bc=40.4934),
        # safety_factor >= 1/sqrt(Fo_max) = sqrt(2), or fo_max >= 1/0.8**2
        "implicit-empty-window": dict(dt=1.0, scheme="implicit",
                                      safety_factor=0.8),
        # a parameter first, then dt
        "implicit-and-drive-cycle": dict(dt=100.0, scheme="implicit",
                                         tau_bc=40.4934),
    }

    @pytest.mark.parametrize("inputs", list(CASES.values()), ids=list(CASES))
    def test_every_printed_remedy_closes_the_conflict(self, inputs):
        result = TransientMeshCalculator.combined_transient_limits(
            **STEEL_MILD, **inputs)
        assert result["feasible"] is False
        exact = {r["parameter"]: r.get("max_value", r.get("min_value"))
                 for r in result["conflict"]["remedies"]}
        for text in (result["advice"], result["binding_constraint"]):
            choices = printed_remedy_choices(text)
            assert choices, text
            for choice in choices:
                fixed = TransientMeshCalculator.combined_transient_limits(
                    **STEEL_MILD, **{**inputs, **choice})
                assert fixed["feasible"] is True, (text, choice)
                assert fixed["conflict"] is None
                if "dt" in choice:
                    assert choice["dt"] <= exact["dt"]
                if "safety_factor" in choice:
                    assert choice["safety_factor"] >= exact["safety_factor"]

    def test_the_structured_remedies_stay_exact(self):
        result = TransientMeshCalculator.combined_transient_limits(
            **STEEL_MILD, **self.CASES["implicit-and-drive-cycle"])
        assert result["conflict"]["remedies"] == [
            {"parameter": "safety_factor", "min_value": 1.0 / math.sqrt(0.5)},
            {"parameter": "dt", "max_value": 0.5 * 40.4934},
        ]
