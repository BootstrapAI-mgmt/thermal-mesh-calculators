"""
Transient Mesh Calculator
=========================

In transient thermal simulations the mesh must be fine enough to resolve
the moving thermal front.  Two independent constraints apply:

1.  Penetration Depth Constraint
    A thermal disturbance propagates a distance of order sqrt(alpha * dt)
    in one time step, where alpha = k / (rho * cp) is the thermal
    diffusivity.  If the element is much larger than this, the solver
    cannot resolve the wave front and numerical oscillations or excessive
    smearing result.

    Practical rule:
        dx_max  <=  C * sqrt(alpha * dt)

    where C ~ 1.5–2.0: a resolution guideline for implicit solvers, whose
    stability does not limit dt.  It is not a stability limit.  Written as
    a Fourier number it reads Fo >= 1/C^2, a lower bound on Fo, which an
    explicit scheme (Fo <= Fo_max) can meet only when C^2 * Fo_max >= 1:
    never at C = 1, whatever dt is.  An explicit scheme needs no such bound
    (see "Combining the constraints" below), so it is applied to implicit
    schemes only.

2.  Fourier Number (Fo) Constraint
    Fo = alpha * dt / dx^2 is the dimensionless ratio of the diffusion
    time scale to the element time scale.

    - Explicit (forward Euler):  Fo <= 0.5 for stability (1D).
      In 3D the limit tightens to Fo <= 1/6.
    - Implicit (backward Euler, Crank-Nicolson):  unconditionally stable,
      but Fo >> 1 smears the solution.  Fo <= 5 is a practical accuracy
      limit for most engineering applications.

    Rearranging for dx:
        dx_min  =  sqrt(alpha * dt / Fo_max)

    For explicit solvers, Fo_max is a hard stability limit.
    For implicit solvers, it is a resolution recommendation.

3.  Drive-Cycle Aware Sizing (Automotive)
    In automotive thermal analysis the boundary conditions change on a
    time scale tau_bc (e.g. 10-30 s between steady-state operating
    points in a drive cycle).  The mesh must resolve the thermal response
    within each quasi-steady segment:

        dx_max  <=  sqrt(alpha * tau_bc)

    This is independent of the solver time-step — it constrains the mesh
    to ensure the spatial resolution can represent the thermal field at
    each operating point transition.

Combining the constraints
    Upper bounds (resolution): the drive-cycle limit, and — for implicit
    schemes only — the penetration depth.  Lower bound: the Fourier limit.

    The explicit (forward Euler, lumped) update in 1D,

        T_i^(n+1) = Fo T_(i-1)^n + (1 - 2 Fo) T_i^n + Fo T_(i+1)^n,

    has non-negative weights summing to one for every Fo <= 1/2: each new
    value is a weighted average of old ones, so no new extremum and no
    oscillation can appear however small Fo is.  An element larger than
    the per-step penetration depth (Fo < 1) is simply the normal explicit
    regime, and the penetration bound is not applied to it.

    The two sqrt(dt) bounds have a ratio independent of dt:

        dx_min(Fourier) / dx_max(penetration) = 1 / (C sqrt(Fo_max))

    so a conflict between them cannot be cured by changing dt; only C or
    Fo_max can.  The drive-cycle bound does not depend on dt, so a Fourier
    minimum above it is cured by dt <= Fo_max * tau_bc.

All the constraints produce independent dx limits.  The binding
constraint (smallest applicable upper bound) governs, provided it is not
below the Fourier minimum.

Units: all internal calculations in SI (m, s, K).  Output mesh sizes
in mm (_mm suffix keys).
"""

import math
from typing import Any, Dict, Optional, cast

from thermal_mesh_calculators._guards import require_positive


def _require_material(k, rho, cp):
    require_positive("k", k, "W/m K")
    require_positive("rho", rho, "kg/m^3")
    require_positive("cp", cp, "J/kg K")

# Relative tolerance on "lower bound <= upper bound" comparisons, so that a
# remedy that puts the two bounds exactly equal counts as feasible.
_BOUND_RTOL = 1e-9

# How the constraint labels read in messages.
_LABEL_TEXT = {
    "fourier_stability": "Fourier stability",
    "fourier_accuracy": "Fourier accuracy",
    "penetration_depth": "penetration-depth",
    "drive_cycle": "drive-cycle",
}


class TransientMeshCalculator:
    """
    Computes transient-specific mesh size constraints based on thermal
    diffusivity, time-step, and (optionally) drive-cycle timing.
    """

    @staticmethod
    def material_diffusivity(k: float, rho: float, cp: float) -> float:
        """
        Thermal diffusivity alpha = k / (rho * cp).

        Parameters
        ----------
        k : float    — thermal conductivity (W/m K)
        rho : float  — density (kg/m^3)
        cp : float   — specific heat capacity (J/kg K)

        Returns
        -------
        float — thermal diffusivity (m^2/s)

        Raises ValueError for k, rho or cp not > 0 (as does every method
        of this class, and for dt, fo_max, safety_factor and tau_bc too).
        """
        _require_material(k, rho, cp)
        return k / (rho * cp)

    @staticmethod
    def penetration_depth(
        k: float,
        rho: float,
        cp: float,
        dt: float,
        safety_factor: float = 1.0,
    ) -> dict:
        """
        Maximum element size from the thermal penetration depth.

        Parameters
        ----------
        k : float           — thermal conductivity (W/m K)
        rho : float         — density (kg/m^3)
        cp : float          — specific heat capacity (J/kg K)
        dt : float          — time-step size (s)
        safety_factor : float
            Multiplier C on sqrt(alpha*dt): a resolution guideline for
            implicit solvers (1.5–2.0).  It is not a stability limit; an
            explicit scheme's stability limit is the Fourier minimum
            (fourier_number_limit), a lower bound on the element size.

        Returns
        -------
        dict with keys:
            max_dx_mm       : float — maximum element size (mm)
            alpha           : float — thermal diffusivity (m^2/s)
            penetration_m   : float — raw penetration depth (m)
        """
        _require_material(k, rho, cp)
        require_positive("dt", dt, "s")
        require_positive("safety_factor", safety_factor)
        alpha = k / (rho * cp)
        pen = math.sqrt(alpha * dt)
        max_dx_m = safety_factor * pen

        return {
            "max_dx_mm": max_dx_m * 1000.0,
            "alpha": alpha,
            "penetration_m": pen,
        }

    @staticmethod
    def fourier_number_limit(
        k: float,
        rho: float,
        cp: float,
        dt: float,
        fo_max: float = 0.5,
    ) -> dict:
        """
        Minimum element size to satisfy the Fourier number constraint.

        For explicit solvers, Fo <= 0.5 (1D) is a hard stability limit.
        For implicit solvers, Fo <= 5 is a practical accuracy guideline.

        Rearranging Fo = alpha * dt / dx^2:
            dx_min = sqrt(alpha * dt / Fo_max)

        Parameters
        ----------
        k : float       — thermal conductivity (W/m K)
        rho : float     — density (kg/m^3)
        cp : float      — specific heat capacity (J/kg K)
        dt : float      — time-step size (s)
        fo_max : float  — maximum allowable Fourier number
                          (0.5 for explicit 1D, 1/6 for explicit 3D,
                           5.0 for implicit guideline)

        Returns
        -------
        dict with keys:
            min_dx_mm   : float — minimum element size for stability (mm)
            alpha       : float — thermal diffusivity (m^2/s)
            fo_max      : float — Fourier number limit used
            constraint  : str   — "stability" if fo_max <= 0.5, else "accuracy"
        """
        _require_material(k, rho, cp)
        require_positive("dt", dt, "s")
        require_positive("fo_max", fo_max)
        alpha = k / (rho * cp)
        min_dx_m = math.sqrt(alpha * dt / fo_max)

        constraint = "stability" if fo_max <= 0.5 else "accuracy"

        return {
            "min_dx_mm": min_dx_m * 1000.0,
            "alpha": alpha,
            "fo_max": fo_max,
            "constraint": constraint,
        }

    @staticmethod
    def drive_cycle_limit(
        k: float,
        rho: float,
        cp: float,
        tau_bc: float,
    ) -> dict:
        """
        Maximum element size to resolve thermal response within a
        drive-cycle operating point transition.

        Parameters
        ----------
        k : float       — thermal conductivity (W/m K)
        rho : float     — density (kg/m^3)
        cp : float      — specific heat capacity (J/kg K)
        tau_bc : float  — boundary condition time scale (s),
                          e.g. duration of a quasi-steady segment

        Returns
        -------
        dict with keys:
            max_dx_mm   : float — maximum element size (mm)
            alpha       : float — thermal diffusivity (m^2/s)
            tau_bc      : float — BC time scale used (s)
        """
        _require_material(k, rho, cp)
        require_positive("tau_bc", tau_bc, "s")
        alpha = k / (rho * cp)
        max_dx_m = math.sqrt(alpha * tau_bc)

        return {
            "max_dx_mm": max_dx_m * 1000.0,
            "alpha": alpha,
            "tau_bc": tau_bc,
        }

    @staticmethod
    def resolve_scheme(scheme: Optional[str] = None, fo_max: float = 0.5) -> str:
        """
        The time-integration scheme: "explicit" or "implicit".

        None infers it from fo_max, as fourier_number_limit labels it:
        fo_max <= 0.5 is an explicit stability limit, a larger fo_max an
        implicit accuracy guideline.
        """
        if scheme is None:
            return "explicit" if fo_max <= 0.5 else "implicit"
        if scheme not in ("explicit", "implicit"):
            raise ValueError(
                f"scheme must be 'explicit', 'implicit' or None, got {scheme!r}"
            )
        return scheme

    @classmethod
    def combined_transient_limits(
        cls,
        k: float,
        rho: float,
        cp: float,
        dt: float,
        fo_max: float = 0.5,
        safety_factor: float = 1.0,
        tau_bc: Optional[float] = None,
        scheme: Optional[str] = None,
    ) -> dict:
        """
        Evaluate all applicable transient mesh constraints and identify
        the binding one.

        The element must be at least the Fourier minimum (explicit
        stability, or implicit accuracy) and at most the applicable upper
        bounds: the drive-cycle limit, and for an implicit scheme the
        penetration depth.  An explicit scheme is not bounded by the
        per-step penetration depth (see the module docstring), so with no
        tau_bc it has no transient upper bound at all.

        Parameters
        ----------
        k : float           — thermal conductivity (W/m K)
        rho : float         — density (kg/m^3)
        cp : float          — specific heat capacity (J/kg K)
        dt : float          — time-step size (s)
        fo_max : float      — Fourier number limit (default 0.5 for explicit)
        safety_factor : float — penetration depth multiplier C (default 1.0;
                               applied to implicit schemes only)
        tau_bc : float      — (optional) drive-cycle segment duration (s)
        scheme : str        — "explicit", "implicit", or None to infer it
                              from fo_max (see resolve_scheme)

        Returns
        -------
        dict with keys:
            alpha               : float — thermal diffusivity (m^2/s)
            scheme              : str   — "explicit" or "implicit"
            penetration_max_dx_mm : float — penetration depth limit (mm)
            penetration_applied : bool  — False for an explicit scheme
            fourier_min_dx_mm   : float — Fourier number lower bound (mm)
            drive_cycle_max_dx_mm : float or None — drive-cycle limit (mm)
            max_dx_mm           : float — tightest applicable upper bound
                                          (mm; inf when none applies)
            feasible            : bool  — fourier_min_dx_mm <= max_dx_mm
            recommended_dx_mm   : float — max_dx_mm when feasible (the
                                          coarsest element every transient
                                          bound allows; inf when none
                                          applies), else fourier_min_dx_mm
            binding_constraint  : str   — "penetration_depth",
                                          "drive_cycle", "none", or on a
                                          conflict "fourier_stability
                                          (WARNING: ...)" /
                                          "fourier_accuracy (WARNING: ...)"
            conflict            : None, or dict with "message" and
                                  "remedies": every parameter change the
                                  conflict needs, each a dict with
                                  "parameter" and "max_value" or
                                  "min_value".  Applying all of them makes
                                  the window feasible.
            advice              : str   — the result in words
        """
        alpha = cls.material_diffusivity(k, rho, cp)
        require_positive("dt", dt, "s")
        require_positive("fo_max", fo_max)
        require_positive("safety_factor", safety_factor)
        if tau_bc is not None:
            require_positive("tau_bc", tau_bc, "s")
        scheme = cls.resolve_scheme(scheme, fo_max)
        fourier_label = ("fourier_stability" if scheme == "explicit"
                         else "fourier_accuracy")

        pen = cls.penetration_depth(k, rho, cp, dt, safety_factor)
        fo = cls.fourier_number_limit(k, rho, cp, dt, fo_max)
        upper_pen = pen["max_dx_mm"]
        lower_fo = fo["min_dx_mm"]

        dc_limit = None
        if tau_bc is not None:
            dc_limit = cls.drive_cycle_limit(k, rho, cp, tau_bc)["max_dx_mm"]

        # Applicable upper bounds.
        penetration_applied = scheme == "implicit"
        uppers = []
        if penetration_applied:
            uppers.append((upper_pen, "penetration_depth"))
        if dc_limit is not None:
            uppers.append((dc_limit, "drive_cycle"))
        if uppers:
            max_dx, max_label = min(uppers)
        else:
            max_dx, max_label = float("inf"), None

        feasible = lower_fo <= max_dx * (1.0 + _BOUND_RTOL)

        remedies = []
        conflict: Optional[Dict[str, Any]] = None
        if not feasible:
            if penetration_applied and lower_fo > upper_pen * (1.0 + _BOUND_RTOL):
                # Both scale with sqrt(dt); their ratio 1/(C sqrt(Fo_max))
                # does not depend on dt, so no dt can close this.
                remedies.append({
                    "parameter": "safety_factor",
                    "min_value": 1.0 / math.sqrt(fo_max),
                })
            if dc_limit is not None and lower_fo > dc_limit * (1.0 + _BOUND_RTOL):
                # sqrt(alpha dt / Fo_max) <= sqrt(alpha tau_bc)
                #   <=>  dt <= Fo_max * tau_bc
                remedies.append({  # dc_limit exists only when tau_bc was given
                    "parameter": "dt",
                    "max_value": fo_max * cast(float, tau_bc),
                })
            conflict = {
                "message": (
                    f"The {_LABEL_TEXT[fourier_label]} minimum "
                    f"({lower_fo:.4g} mm, Fo <= {fo_max:g}) exceeds the "
                    # Infeasible means a finite upper bound, so max_label is set.
                    f"{_LABEL_TEXT[cast(str, max_label)]} bound ({max_dx:.4g} mm): "
                    f"no element satisfies both at dt = {dt:g} s. "
                    + _describe_remedies(remedies, fo_max, safety_factor)
                ),
                "remedies": remedies,
            }

        if feasible:
            recommended = max_dx
            binding = max_label if max_label is not None else "none"
            if max_label is None:
                advice = (
                    f"Explicit scheme: dt = {dt:g} s integrates stably any "
                    f"element of {lower_fo:.4g} mm or larger (Fo <= "
                    f"{fo_max:g}). The per-step penetration depth does not "
                    f"bound an explicit scheme, so no transient upper bound "
                    f"applies; give tau_bc for one, otherwise the steady "
                    f"constraints set the size."
                )
            else:
                advice = (
                    f"Use elements between {lower_fo:.4g} mm "
                    f"({fourier_label}) and {max_dx:.4g} mm ({max_label})."
                )
        else:
            recommended = lower_fo
            binding = (
                f"{fourier_label} (WARNING: conflicts with {max_label}; "
                + _describe_remedies(remedies, fo_max, safety_factor,
                                     short=True)
                + ")"
            )
            # conflict was built above, as for every infeasible window.
            advice = cast(Dict[str, Any], conflict)["message"]

        return {
            "alpha": alpha,
            "scheme": scheme,
            "penetration_max_dx_mm": upper_pen,
            "penetration_applied": penetration_applied,
            "fourier_min_dx_mm": lower_fo,
            "drive_cycle_max_dx_mm": dc_limit,
            "max_dx_mm": max_dx,
            "feasible": feasible,
            "recommended_dx_mm": recommended,
            "binding_constraint": binding,
            "conflict": conflict,
            "advice": advice,
        }


def _bound_text(value: float, bound: str) -> str:
    """
    A remedy's bound in words, to 4 significant digits, rounded toward the
    side that satisfies it: a maximum (bound="max") down, a minimum
    (bound="min") up, so that the value as printed closes the conflict.
    Rounded to nearest, a maximum of 20.2467 s printed as 20.25 s, and
    applying the printed value left the conflict in place.

    The value is first taken to 12 significant digits, so floating-point
    noise does not become a digit (a computed 0.49999999999999994 prints as
    0.5, not 0.4999).  That moves a bound by at most 5e-12 relative, far
    inside the 1e-9 relative tolerance of the bound comparisons.  The
    structured remedies keep the exact values.
    """
    if not (value > 0.0 and math.isfinite(value)):
        return f"{value:.4g}"
    mantissa, exponent = f"{value:.11e}".split("e")
    digits = mantissa.replace(".", "")
    head, tail = int(digits[:4]), digits[4:]
    if bound == "min" and tail.strip("0"):
        head += 1
    return f"{head * 10.0 ** (int(exponent) - 3):.4g}"


def _describe_remedies(remedies: list, fo_max: float, safety_factor: float,
                       short: bool = False) -> str:
    """The remedies of a transient conflict, in words (see _bound_text)."""
    parts = []
    for remedy in remedies:
        if remedy["parameter"] == "dt":
            parts.append(
                f"reduce dt to <= {_bound_text(remedy['max_value'], 'max')} s")
        else:
            text = (f"raise safety_factor to >= "
                    f"{_bound_text(remedy['min_value'], 'min')}"
                    f" or fo_max to >= "
                    f"{_bound_text(1.0 / safety_factor ** 2, 'min')}")
            if not short:
                text += " (reducing dt cannot close this: both bounds scale with sqrt(dt))"
            parts.append(text)
    joined = "; then ".join(parts)
    return joined if short else "Remedy: " + joined + "."
