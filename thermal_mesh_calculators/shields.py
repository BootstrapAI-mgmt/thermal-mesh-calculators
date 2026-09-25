"""
Heat Shield Temperature & Mesh Solvers
=======================================

In automotive underhood thermal management, exhaust component surface
temperatures are typically fixed in the solver (measured or specified).
Non-exhaust structural parts can usually be estimated to within +/- 30 K.

Heat shields are the primary floating unknown:  they sit between a
known-hot exhaust source and a known-cool ambient, and their equilibrium
temperature is governed by a nonlinear energy balance involving both
convection and radiation.

This module provides two solvers:

1.  SingleLayerShieldCalculator
    For simple stamped shields (steel or aluminium) where the Biot number
    is very small (Bi << 0.1) and internal conduction can be neglected.
    Uses scalar Newton-Raphson.

2.  MultilayerShieldCalculator
    For dual-wall shields with an air gap or insulating layer between
    them (dimpled, offset, or fibre-filled).  The inner and outer layers
    can differ by hundreds of degrees.  Uses multivariate Newton-Raphson
    with a 2x2 Jacobian (pure Python, no numpy dependency).

Both solvers now use a consistent API with separate h_in / h_out
convection coefficients.  The single-layer solver also accepts the
legacy h_total parameter for backward compatibility.

Both solvers output per-layer equilibrium temperatures AND the
corresponding boundary-driven mesh sizes, so the user gets independent
mesh targets for each layer.
"""

import math
from thermal_mesh_calculators.constants import STEFAN_BOLTZMANN
from thermal_mesh_calculators._guards import (
    require_fraction,
    require_non_negative,
    require_positive,
    require_temperatures,
)


# ---------------------------------------------------------------------------
#  Single-layer shield
# ---------------------------------------------------------------------------

class SingleLayerShieldCalculator:
    """
    Solves the 1D steady-state energy balance for a thin, single-layer
    heat shield and computes the boundary-driven mesh size.

    Energy balance at the shield:

        q_rad_in - q_rad_out - q_conv_in - q_conv_out = 0

    where:
        q_rad_in   = eps_in  * sigma * (T_exh^4   - T_shield^4)
        q_rad_out  = eps_out * sigma * (T_shield^4 - T_surr^4)
        q_conv_in  = h_in  * (T_shield - T_fluid)
        q_conv_out = h_out * (T_shield - T_fluid)

    Newton-Raphson derivative:
        F'(T) = -4 sigma T^3 (eps_in + eps_out) - (h_in + h_out)

    API (v0.2):
        Accepts separate h_in / h_out to match the multilayer solver.
        Legacy h_total parameter is still supported — if provided
        without h_in/h_out, it is split evenly as h_in = h_out = h_total / 2.
    """

    @staticmethod
    def _resolve_htc(h_in, h_out, h_total):
        """Resolve convection coefficients from either new or legacy API."""
        if h_in is not None and h_out is not None:
            return h_in, h_out
        if h_total is not None:
            return h_total / 2.0, h_total / 2.0
        raise ValueError(
            "Provide either (h_in, h_out) or h_total for convection."
        )

    @staticmethod
    def solve_temperature(
        t_exh: float,
        t_fluid: float,
        t_surr: float,
        eps_in: float,
        eps_out: float,
        h_in: float = None,
        h_out: float = None,
        h_total: float = None,
        tol: float = 0.1,
        max_iter: int = 50,
    ) -> dict:
        """
        Find the steady-state shield temperature.

        Parameters
        ----------
        t_exh : float      — fixed exhaust surface temperature (K)
        t_fluid : float    — ambient / underhood air temperature (K)
        t_surr : float     — radiation sink temperature (K)
        eps_in : float     — emissivity of the exhaust-facing surface
        eps_out : float    — emissivity of the ambient-facing surface
        h_in : float       — convective HTC on exhaust-facing side (W/m^2 K)
        h_out : float      — convective HTC on ambient-facing side (W/m^2 K)
        h_total : float    — (legacy) total HTC, split evenly if h_in/h_out
                             not provided (W/m^2 K)
        tol : float        — convergence tolerance (K)
        max_iter : int     — iteration limit

        Returns
        -------
        dict with keys:
            t_shield_K  : float — equilibrium temperature (K)
            t_shield_C  : float — equilibrium temperature (deg C)
            converged   : bool
            iterations  : int
            residual_W_m2 : float — |energy-balance residual| at t_shield_K
                                    (W/m^2)
            q_rad_in    : float — radiation absorbed from exhaust (W/m^2)
            q_rad_out   : float — radiation emitted to surroundings (W/m^2)
            q_conv_in   : float — convective flux, exhaust side (W/m^2)
            q_conv_out  : float — convective flux, ambient side (W/m^2)
            h_in        : float — resolved inner HTC used (W/m^2 K)
            h_out       : float — resolved outer HTC used (W/m^2 K)
        """
        h_i, h_o = SingleLayerShieldCalculator._resolve_htc(h_in, h_out, h_total)
        require_temperatures(t_exh=t_exh, t_fluid=t_fluid, t_surr=t_surr)
        require_fraction("eps_in", eps_in)
        require_fraction("eps_out", eps_out)
        require_non_negative("h_in", h_i, "W/m^2 K")
        require_non_negative("h_out", h_o, "W/m^2 K")
        h_sum = h_i + h_o

        t = (t_exh + t_fluid) / 2.0  # initial guess

        for i in range(max_iter):
            q_rad_in = eps_in * STEFAN_BOLTZMANN * (t_exh**4 - t**4)
            q_rad_out = eps_out * STEFAN_BOLTZMANN * (t**4 - t_surr**4)
            q_conv_in = h_i * (t - t_fluid)
            q_conv_out = h_o * (t - t_fluid)

            f = q_rad_in - q_rad_out - q_conv_in - q_conv_out
            fp = -4.0 * STEFAN_BOLTZMANN * (t**3) * (eps_in + eps_out) - h_sum

            t_new = t - f / fp
            if abs(t_new - t) < tol:
                t = t_new
                q_rad_in = eps_in * STEFAN_BOLTZMANN * (t_exh**4 - t**4)
                q_rad_out = eps_out * STEFAN_BOLTZMANN * (t**4 - t_surr**4)
                q_conv_in = h_i * (t - t_fluid)
                q_conv_out = h_o * (t - t_fluid)
                return {
                    "t_shield_K": t,
                    "t_shield_C": t - 273.15,
                    "converged": True,
                    "iterations": i + 1,
                    "residual_W_m2": abs(
                        q_rad_in - q_rad_out - q_conv_in - q_conv_out),
                    "q_rad_in": q_rad_in,
                    "q_rad_out": q_rad_out,
                    "q_conv_in": q_conv_in,
                    "q_conv_out": q_conv_out,
                    "h_in": h_i,
                    "h_out": h_o,
                }
            t = t_new

        # Not converged: report the fluxes at the last iterate, flagged.
        q_rad_in = eps_in * STEFAN_BOLTZMANN * (t_exh**4 - t**4)
        q_rad_out = eps_out * STEFAN_BOLTZMANN * (t**4 - t_surr**4)
        q_conv_in = h_i * (t - t_fluid)
        q_conv_out = h_o * (t - t_fluid)
        return {
            "t_shield_K": t,
            "t_shield_C": t - 273.15,
            "converged": False,
            "iterations": max_iter,
            "residual_W_m2": abs(q_rad_in - q_rad_out - q_conv_in - q_conv_out),
            "q_rad_in": q_rad_in,
            "q_rad_out": q_rad_out,
            "q_conv_in": q_conv_in,
            "q_conv_out": q_conv_out,
            "h_in": h_i,
            "h_out": h_o,
        }

    @classmethod
    def mesh_size(
        cls,
        k: float,
        max_dt: float,
        t_exh: float,
        t_fluid: float,
        t_surr: float,
        eps_in: float,
        eps_out: float,
        h_in: float = None,
        h_out: float = None,
        h_total: float = None,
        tol: float = 0.1,
        max_iter: int = 50,
    ) -> dict:
        """
        Solve for shield temperature, then compute boundary-driven mesh size.

        Parameters
        ----------
        k : float       — shield thermal conductivity (W/m K)
        max_dt : float  — max allowable delta-T per element (K)
        (remaining params, including tol and max_iter, forwarded to
        solve_temperature)

        Returns
        -------
        dict with keys from solve_temperature plus:
            max_dx_mm     : float — maximum element size (mm)
            q_boundary    : float — driving flux at hottest face (W/m^2)

        When the solve does not converge (``converged`` False) the size is
        computed from the last iterate's fluxes: an estimate, not a result.
        """
        require_positive("k", k, "W/m K")
        require_positive("max_dt", max_dt, "K")
        result = cls.solve_temperature(
            t_exh, t_fluid, t_surr, eps_in, eps_out,
            h_in=h_in, h_out=h_out, h_total=h_total,
            tol=tol, max_iter=max_iter,
        )

        # Driving flux at the exhaust face (radiation in + convection on that side)
        q_boundary = abs(result["q_rad_in"]) + abs(result["q_conv_in"])

        if q_boundary > 0:
            result["max_dx_mm"] = (k * max_dt) / q_boundary * 1000.0
        else:
            result["max_dx_mm"] = float("inf")

        result["q_boundary"] = q_boundary
        return result


# ---------------------------------------------------------------------------
#  Multilayer shield
# ---------------------------------------------------------------------------

class MultilayerShieldCalculator:
    """
    Solves the coupled nonlinear energy balance for a 2-layer heat shield
    assembly with an air gap (or insulating layer) between the layers.

    Layer 1 (inner, exhaust-facing):
        F1 = q_rad_in + q_conv_in - q_gap_cond - q_gap_rad = 0

    Layer 2 (outer, ambient-facing):
        F2 = q_gap_cond + q_gap_rad - q_rad_out - q_conv_out = 0

    Gap heat transfer:
        q_gap_cond = h_gap * (T1 - T2)
        q_gap_rad  = eps_eff * sigma * (T1^4 - T2^4)

    where eps_eff = 1 / (1/eps_g1 + 1/eps_g2 - 2 + 1/F12)

    F12 is the geometric view factor between the gap faces.  For large
    parallel plates F12 → 1 and this reduces to the classical formula.
    For finite or offset geometries F12 < 1 reduces the effective
    emissivity (less radiative coupling across the gap).

    Solved with multivariate Newton-Raphson using the 2x2 Jacobian,
    inverted analytically via Cramer's rule (no numpy needed).
    """

    @staticmethod
    def solve_temperatures(
        t_exh: float,
        t_fluid: float,
        t_surr: float,
        h_in: float,
        h_out: float,
        h_gap: float,
        eps_in: float,
        eps_out: float,
        eps_g1: float,
        eps_g2: float,
        f12: float = 1.0,
        tol: float = 0.1,
        max_iter: int = 100,
    ) -> dict:
        """
        Find steady-state temperatures for both layers.

        Parameters
        ----------
        t_exh : float    — fixed exhaust surface temperature (K)
        t_fluid : float  — ambient air temperature (K)
        t_surr : float   — radiation sink temperature (K)
        h_in : float     — convective HTC on inner face (W/m^2 K)
        h_out : float    — convective HTC on outer face (W/m^2 K)
        h_gap : float    — effective gap conductance (W/m^2 K),
                           includes contact conduction + air conduction
        eps_in : float   — emissivity of exhaust-facing surface
        eps_out : float  — emissivity of ambient-facing surface
        eps_g1 : float   — emissivity of inner gap face
        eps_g2 : float   — emissivity of outer gap face
        f12 : float      — geometric view factor between gap faces (0-1).
                           Default 1.0 (infinite parallel plates).
                           Use ~0.85 for typical automotive offset shields;
                           lower for small or highly non-parallel gaps.
                           A zero gap emissivity or view factor means no
                           radiative exchange across the gap (eps_eff = 0).
        tol : float      — convergence tolerance on both balances (W/m^2)
        max_iter : int   — iteration limit

        Returns
        -------
        dict with both layer temperatures, the six fluxes, convergence info
        (converged, iterations) and residual_W_m2, the larger of the two
        |energy-balance residuals| at the returned temperatures.  When the
        solve does not converge the temperatures and fluxes are those of the
        last iterate: an estimate, flagged by ``converged`` False.
        """
        require_temperatures(t_exh=t_exh, t_fluid=t_fluid, t_surr=t_surr)
        for name, value in (("eps_in", eps_in), ("eps_out", eps_out),
                            ("eps_g1", eps_g1), ("eps_g2", eps_g2),
                            ("f12", f12)):
            require_fraction(name, value)
        for name, value in (("h_in", h_in), ("h_out", h_out),
                            ("h_gap", h_gap)):
            require_non_negative(name, value, "W/m^2 K")

        # Effective gap emissivity including geometric view factor
        # Reduces to 1/(1/e1 + 1/e2 - 1) when F12 = 1; tends to 0 as any of
        # the three does (no radiative exchange across the gap).
        if eps_g1 == 0 or eps_g2 == 0 or f12 == 0:
            eps_eff = 0.0
        else:
            eps_eff = 1.0 / ((1.0 / eps_g1) + (1.0 / eps_g2) - 2.0 + (1.0 / f12))

        def balance(t1, t2):
            """The six fluxes at (t1, t2) and the two layer residuals."""
            q = {
                "q_rad_in": eps_in * STEFAN_BOLTZMANN * (t_exh**4 - t1**4),
                "q_conv_in": h_in * (t_fluid - t1),
                "q_gap_cond": h_gap * (t1 - t2),
                "q_gap_rad": eps_eff * STEFAN_BOLTZMANN * (t1**4 - t2**4),
                "q_rad_out": eps_out * STEFAN_BOLTZMANN * (t2**4 - t_surr**4),
                "q_conv_out": h_out * (t2 - t_fluid),
            }
            f1 = q["q_rad_in"] + q["q_conv_in"] - q["q_gap_cond"] - q["q_gap_rad"]
            f2 = q["q_gap_cond"] + q["q_gap_rad"] - q["q_rad_out"] - q["q_conv_out"]
            return q, f1, f2

        def report(t1, t2, converged, iterations, q, f1, f2):
            out = {
                "t1_K": t1,
                "t1_C": t1 - 273.15,
                "t2_K": t2,
                "t2_C": t2 - 273.15,
                "delta_T_C": (t1 - t2),
                "converged": converged,
                "iterations": iterations,
                "eps_eff": eps_eff,
                "residual_W_m2": max(abs(f1), abs(f2)),
            }
            out.update(q)
            return out

        # Initial guesses
        t1 = t_exh - 100.0
        t2 = t_surr + 100.0

        for i in range(max_iter):
            q, f1, f2 = balance(t1, t2)
            if abs(f1) < tol and abs(f2) < tol:
                return report(t1, t2, True, i + 1, q, f1, f2)

            # Jacobian
            j11 = (-4.0 * eps_in * STEFAN_BOLTZMANN * t1**3
                   - h_in - h_gap
                   - 4.0 * eps_eff * STEFAN_BOLTZMANN * t1**3)
            j12 = h_gap + 4.0 * eps_eff * STEFAN_BOLTZMANN * t2**3
            j21 = h_gap + 4.0 * eps_eff * STEFAN_BOLTZMANN * t1**3
            j22 = (-h_gap
                   - 4.0 * eps_eff * STEFAN_BOLTZMANN * t2**3
                   - 4.0 * eps_out * STEFAN_BOLTZMANN * t2**3
                   - h_out)

            det = j11 * j22 - j12 * j21
            if det == 0:
                raise ValueError("Singular Jacobian — check inputs.")

            dt1 = -(f1 * j22 - f2 * j12) / det
            dt2 = -(j11 * f2 - j21 * f1) / det

            t1 += dt1
            t2 += dt2

        # Out of iterations: report the last iterate with its fluxes, so the
        # per-layer sizes stay finite estimates instead of becoming inf.
        q, f1, f2 = balance(t1, t2)
        return report(t1, t2, abs(f1) < tol and abs(f2) < tol, max_iter,
                      q, f1, f2)

    @classmethod
    def mesh_sizes(
        cls,
        k_metal: float,
        max_dt: float,
        **kwargs,
    ) -> dict:
        """
        Solve both layer temperatures, then compute independent mesh sizes.

        Parameters
        ----------
        k_metal : float  — shield layer conductivity (W/m K)
        max_dt : float   — max allowable delta-T per element (K)
        **kwargs          — forwarded to solve_temperatures()

        Returns
        -------
        dict with keys from solve_temperatures() plus:
            layer1_max_dx_mm : float
            layer2_max_dx_mm : float
            q_layer1         : float — layer 1 driving flux,
                                       |q_rad_in| + |q_conv_in| (W/m^2)
            q_layer2         : float — layer 2 driving flux,
                                       |q_rad_out| + |q_conv_out| (W/m^2)

        When the solve does not converge (``converged`` False) the sizes are
        computed from the last iterate's fluxes: estimates, not results.
        """
        require_positive("k_metal", k_metal, "W/m K")
        require_positive("max_dt", max_dt, "K")
        result = cls.solve_temperatures(**kwargs)

        # Layer 1: driven by exhaust radiation + inner convection
        q1 = abs(result.get("q_rad_in", 0)) + abs(result.get("q_conv_in", 0))
        # Layer 2: driven by external radiation + outer convection
        q2 = abs(result.get("q_rad_out", 0)) + abs(result.get("q_conv_out", 0))

        result["layer1_max_dx_mm"] = (k_metal * max_dt / q1 * 1000.0) if q1 > 0 else float("inf")
        result["layer2_max_dx_mm"] = (k_metal * max_dt / q2 * 1000.0) if q2 > 0 else float("inf")
        result["q_layer1"] = q1
        result["q_layer2"] = q2

        return result
