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

    where C ~ 1.0–2.0 depending on the time-integration scheme.  For
    explicit solvers C = 1 is mandatory (CFL stability); for implicit
    solvers C = 2 is a resolution guideline.

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

All three constraints produce independent dx limits.  The binding
constraint (smallest dx) governs.

Units: all internal calculations in SI (m, s, K).  Output mesh sizes
in mm (_mm suffix keys).
"""

import math


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
        """
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
            Multiplier on sqrt(alpha*dt).  Use 1.0 for explicit solvers
            (stability limit), 1.5–2.0 for implicit solvers (resolution).

        Returns
        -------
        dict with keys:
            max_dx_mm       : float — maximum element size (mm)
            alpha           : float — thermal diffusivity (m^2/s)
            penetration_m   : float — raw penetration depth (m)
        """
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
        alpha = k / (rho * cp)
        max_dx_m = math.sqrt(alpha * tau_bc)

        return {
            "max_dx_mm": max_dx_m * 1000.0,
            "alpha": alpha,
            "tau_bc": tau_bc,
        }

    @classmethod
    def combined_transient_limits(
        cls,
        k: float,
        rho: float,
        cp: float,
        dt: float,
        fo_max: float = 0.5,
        safety_factor: float = 1.0,
        tau_bc: float = None,
    ) -> dict:
        """
        Evaluate all applicable transient mesh constraints and identify
        the binding one.

        Parameters
        ----------
        k : float           — thermal conductivity (W/m K)
        rho : float         — density (kg/m^3)
        cp : float          — specific heat capacity (J/kg K)
        dt : float          — time-step size (s)
        fo_max : float      — Fourier number limit (default 0.5 for explicit)
        safety_factor : float — penetration depth multiplier (default 1.0)
        tau_bc : float      — (optional) drive-cycle segment duration (s)

        Returns
        -------
        dict with keys:
            alpha               : float — thermal diffusivity (m^2/s)
            penetration_max_dx_mm : float — penetration depth limit (mm)
            fourier_min_dx_mm   : float — Fourier number lower bound (mm)
            drive_cycle_max_dx_mm : float or None — drive-cycle limit (mm)
            recommended_dx_mm   : float — recommended element size (mm)
            binding_constraint  : str   — which constraint governs
        """
        alpha = k / (rho * cp)

        pen = cls.penetration_depth(k, rho, cp, dt, safety_factor)
        fo = cls.fourier_number_limit(k, rho, cp, dt, fo_max)

        # Upper bound from penetration
        upper_pen = pen["max_dx_mm"]
        # Lower bound from Fourier stability
        lower_fo = fo["min_dx_mm"]

        # Drive cycle upper bound (optional)
        dc_limit = None
        if tau_bc is not None:
            dc = cls.drive_cycle_limit(k, rho, cp, tau_bc)
            dc_limit = dc["max_dx_mm"]

        # Determine recommended size and binding constraint
        # The element must be:
        #   >= lower_fo   (Fourier stability)
        #   <= upper_pen  (penetration resolution)
        #   <= dc_limit   (drive cycle resolution, if applicable)

        upper_bounds = [upper_pen]
        upper_labels = ["penetration_depth"]
        if dc_limit is not None:
            upper_bounds.append(dc_limit)
            upper_labels.append("drive_cycle")

        # Tightest upper bound
        min_upper = min(upper_bounds)
        min_upper_label = upper_labels[upper_bounds.index(min_upper)]

        if lower_fo > min_upper:
            # Fourier stability demands a larger element than resolution allows.
            # This means the time-step is too large for the desired resolution.
            # Report the conflict.
            recommended = lower_fo
            binding = "fourier_stability (WARNING: conflicts with resolution — reduce dt)"
        elif lower_fo <= min_upper:
            # Normal case: Fourier lower bound is below the upper bound.
            # Recommend the tightest upper bound (finest resolution requirement).
            recommended = min_upper
            binding = min_upper_label

        return {
            "alpha": alpha,
            "penetration_max_dx_mm": upper_pen,
            "fourier_min_dx_mm": lower_fo,
            "drive_cycle_max_dx_mm": dc_limit,
            "recommended_dx_mm": recommended,
            "binding_constraint": binding,
        }
