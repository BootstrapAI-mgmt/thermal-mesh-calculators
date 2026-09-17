"""
Aerodynamic Boundary Layer Mesh Calculator
==========================================

Estimates surface mesh size and inflation layer parameters from
aerodynamic boundary layer physics.  Two flow regimes are supported:

1.  **External forced** — attached, predominantly unidirectional flow
    (underbody panels in freestream, cooling-pack downstream).  Uses
    flat-plate skin-friction correlations (Schlichting) to compute
    friction velocity and first-cell height from a target y+.
    Prism layers are recommended: wall gradients are steep (~10⁷ K/m
    for exhaust), and the flow has strong directional anisotropy that
    prisms exploit.

2.  **Mixed / unknown** — regions where flow direction is variable,
    recirculating, impinging, or buoyancy-dominated (cargo bed above
    exhaust, wake zones, engine bay dead zones).  Uses the maximum
    of the local freestream velocity and a buoyancy velocity scale
    √(g·β·ΔT·L).  The surface mesh constraint (dx from δ fraction)
    remains useful, but **prism layers are generally not recommended**
    because the y⁺ formalism produces absurdly large first-cell
    heights at low velocities, and isotropic tets at the thermal
    mesh size already over-resolve the wall-normal gradient.

Prism layer recommendation logic
---------------------------------
The ``prisms_recommended`` output flag evaluates whether structured
inflation layers provide meaningful accuracy improvement over
isotropic tets at the wall.  The decision is based on two criteria:

    1.  **y₁ vs dx_surface**: If the y⁺-derived first cell height
        exceeds the surface mesh size, tets already contain multiple
        cells within the "first cell" zone and over-resolve the
        wall-normal gradient.  Prisms add no resolution benefit.

    2.  **Wall gradient magnitude**: In low-flow buoyancy zones
        (0.1–1 m/s), the wall temperature gradient is ~10⁴ K/m —
        roughly 500× smaller than external forced regions (~10⁷ K/m).
        The finite-volume cross-diffusion error from tet non-
        orthogonality at the wall (Error ∝ k_geo · ∇φ) is therefore
        ~500× smaller in absolute terms.  At ~10–20% HTC error, this
        is within the ±30% physics uncertainty of h itself in these
        mixed/buoyancy zones.  The tet-to-prism transition simply
        moves the same non-orthogonality error a few cells away
        from the wall, where in low-Re flow the gradient has not
        decayed significantly (near-linear profile across the gap).

Volume expansion ratio (external forced):
    At the prism-to-tet transition, the cell volume jump is
    quantified by ER_v = V_tet / V_prism.  The truncation error
    of 2nd-order schemes degrades as ER_v increases:

        ER_v ≤ 5  : stable, 2nd-order accuracy maintained
        5 < ER_v ≤ 10 : marginal, gradient truncation at interface
        ER_v > 10 : unstable, diagonal dominance loss likely

    The calculator reports ER_v estimated as (dx_surface / y_last)³.

Physics basis
-------------
    y₁  = y⁺ · ν / u_τ           first cell height
    u_τ = √(τ_w / ρ)             friction velocity
    τ_w = ½ · C_f · ρ · U²       wall shear stress
    C_f = 0.058 · Re_x^(−0.2)    turbulent flat plate (Schlichting)
    C_f = 0.664 · Re_x^(−0.5)    laminar flat plate

The surface mesh constraint couples to the inflation stack via:

    External forced:
        dx ≤ min(AR_transition · y_last,  C_BL · δ)
        where y_last = y₁ · r^(n-1) is the outermost prism height

    Mixed / unknown:
        dx ≤ C_BL · δ  (BL fraction constraint only)
        Prism AR constraint suppressed; isotropic tets to wall.

Boundary layer thickness:
    Turbulent:  δ = 0.37 · x · Re_x^(−0.2)
    Laminar:    δ = 5.0  · x · Re_x^(−0.5)

Wall temperature gradient (for prism decision):
    dT/dy|_wall = h · ΔT / k_air
    where h is estimated from the flow regime and ΔT = |T_s - T_f|.

Leading-edge singularity:
    As x → 0, C_f → ∞ and y₁ → 0.  A minimum Reynolds number floor
    (Re_x_min = 1000) prevents this, corresponding to ~1 mm at typical
    underbody conditions.

Notes for automotive workflows
------------------------------
- Wall-function meshes (y⁺ ≈ 30–100) are the industry standard for
  full-vehicle underbody thermal analysis.
- Wall-resolved meshes (y⁺ ≈ 1) are reserved for single-component
  validation studies and should not be used with this calculator's
  surface sizing logic without manual override.
- For low-Re regions (buoyancy-dominated gaps, dead zones), isotropic
  tets sized by thermal constraints are sufficient.  The surface dx
  from this calculator confirms the thermal mesh is compatible with
  the flow solver, but inflation layers are not needed.

Units:  SI throughout (m, s, K, W).  Output sizes in mm.
"""

import math
from thermal_mesh_calculators.h_estimator import air_properties


# ---------------------------------------------------------------------------
#  Constants / defaults
# ---------------------------------------------------------------------------

# Minimum Re_x to prevent leading-edge singularity.
# Re_x = 1000 corresponds to ~1 mm at 22 m/s in air at 350 K.
RE_X_MIN = 1000.0

# Critical Re_x for laminar-to-turbulent transition (flat plate).
RE_X_CRIT = 5.0e5

# Absolute minimum first cell height (m).  Manufacturing tolerance
# floor — no mesher should attempt cells below ~1 µm.
Y1_ABSOLUTE_MIN_M = 1.0e-6

# Gravitational acceleration (m/s²).
_G = 9.81

# Wall gradient threshold for prism recommendation (K/m).
# Below this, tet non-orthogonality error is within the physics
# uncertainty of h (~30%).  Derived from: at dT/dy ~ 2e4 K/m,
# the non-orthogonality error in h is ~10-20%, which is within the
# ±30% uncertainty of h itself in mixed/buoyancy correlations.
_WALL_GRADIENT_PRISM_THRESHOLD = 5.0e4  # K/m

# Volume expansion ratio thresholds (prism-to-tet transition).
ER_V_STABLE = 5.0     # 2nd-order accuracy maintained
ER_V_MARGINAL = 10.0  # gradient truncation, solver may struggle


class BoundaryLayerCalculator:
    """
    Computes inflation layer parameters and surface mesh constraints
    from aerodynamic boundary layer physics.
    """

    @staticmethod
    def _skin_friction(re_x: float) -> tuple:
        """
        Flat-plate skin friction coefficient and regime label.

        Parameters
        ----------
        re_x : float — local Reynolds number (clamped to RE_X_MIN)

        Returns
        -------
        (C_f, regime_str)
        """
        re_x = max(re_x, RE_X_MIN)
        if re_x < RE_X_CRIT:
            return 0.664 * re_x ** (-0.5), "laminar"
        else:
            return 0.058 * re_x ** (-0.2), "turbulent"

    @staticmethod
    def _boundary_layer_thickness(re_x: float, x: float) -> float:
        """
        Boundary layer thickness δ (m) from flat-plate correlations.

        Parameters
        ----------
        re_x : float — local Reynolds number (clamped to RE_X_MIN)
        x : float — distance from leading edge (m)

        Returns
        -------
        float — δ in metres
        """
        re_x = max(re_x, RE_X_MIN)
        if re_x < RE_X_CRIT:
            return 5.0 * x * re_x ** (-0.5)
        else:
            return 0.37 * x * re_x ** (-0.2)

    @staticmethod
    def _buoyancy_velocity(beta: float, delta_t: float,
                           char_length: float) -> float:
        """
        Buoyancy-driven characteristic velocity √(g·β·ΔT·L).

        Parameters
        ----------
        beta : float — volumetric expansion coefficient (1/K)
        delta_t : float — |T_surf − T_fluid| (K)
        char_length : float — characteristic length (m), e.g. gap height

        Returns
        -------
        float — buoyancy velocity (m/s)
        """
        return math.sqrt(_G * beta * abs(delta_t) * char_length)

    @staticmethod
    def _inflation_layer_count(y1: float, delta: float,
                               growth_ratio: float) -> int:
        """
        Number of inflation layers to span the boundary layer.

        n = ⌈ ln(1 + (δ/y₁)·(r−1)) / ln(r) ⌉

        Parameters
        ----------
        y1 : float — first cell height (m)
        delta : float — boundary layer thickness (m)
        growth_ratio : float — geometric growth ratio (> 1)

        Returns
        -------
        int — number of inflation layers (minimum 3)
        """
        if y1 <= 0 or delta <= 0 or growth_ratio <= 1.0:
            return 3
        ratio = delta / y1
        if ratio <= 1.0:
            return 3
        n = math.log(1.0 + ratio * (growth_ratio - 1.0)) / math.log(
            growth_ratio
        )
        return max(3, math.ceil(n))

    @staticmethod
    def estimate_mesh(
        U: float,
        x: float,
        t_fluid: float,
        y_plus_target: float = 30.0,
        growth_ratio: float = 1.2,
        regime: str = "external_forced",
        delta_t_buoyancy: float = 0.0,
        ar_max_prism: float = 5.0,
        bl_fraction: float = 0.3,
        ar_transition: float = 3.0,
        t_surf: float = 0.0,
    ) -> dict:
        """
        Estimate surface mesh size and inflation layer parameters.

        Parameters
        ----------
        U : float
            Freestream or local velocity (m/s).  For ``mixed_unknown``
            regime, this is compared against the buoyancy velocity and
            the larger is used.
        x : float
            Characteristic length (m).  For ``external_forced`` this is
            the distance from the leading edge.  For ``mixed_unknown``
            this is the gap height or part dimension.
        t_fluid : float
            Fluid / ambient temperature (K).
        y_plus_target : float
            Target y⁺ for the first inflation cell.  Default 30
            (wall-function).  Use 1.0 only for wall-resolved studies.
        growth_ratio : float
            Geometric growth ratio for inflation layers (> 1.0).
            Default 1.2.
        regime : str
            ``"external_forced"`` or ``"mixed_unknown"``.
        delta_t_buoyancy : float
            |T_surf − T_fluid| for buoyancy velocity calculation (K).
            Only used when ``regime == "mixed_unknown"``.  If 0 and
            t_surf is provided, computed automatically.
        ar_max_prism : float
            Maximum prism cell aspect ratio.  Only enforced in
            ``mixed_unknown`` regime.  Default 5.0.
        bl_fraction : float
            Fraction of δ used as surface mesh constraint (C in
            dx ≤ C·δ).  Default 0.3.  Range 0.2–0.5.
        ar_transition : float
            Aspect ratio of last prism cell to surface mesh in
            ``external_forced`` regime.  Default 3.0.
        t_surf : float
            Surface temperature (K).  Used to compute delta_t_buoyancy
            if not provided directly.  Also used for film temperature.

        Returns
        -------
        dict with keys:
            y1_mm              : float — first cell height (mm)
            y1_m               : float — first cell height (m)
            n_layers           : int   — inflation layer count
            delta_mm           : float — boundary layer thickness (mm)
            y_last_mm          : float — outermost prism height (mm)
            total_prism_mm     : float — total inflation stack height (mm)
            max_dx_surface_mm  : float — maximum surface element size (mm)
            surface_constraint : str   — which constraint drove dx
            prisms_recommended : bool  — whether inflation layers are beneficial
            prism_reason       : str   — human-readable explanation
            wall_gradient_Km   : float — estimated dT/dy at wall (K/m)
            volume_expansion_ratio : float — ER_v at prism-to-tet transition
            er_v_status        : str   — "stable", "marginal", or "unstable"
            U_effective        : float — velocity used for calculation (m/s)
            Re_x               : float — Reynolds number at x
            Cf                 : float — skin friction coefficient
            u_tau              : float — friction velocity (m/s)
            regime             : str   — flow regime used
            bl_regime          : str   — "laminar" or "turbulent"
            y_plus_target      : float — echoed back
            growth_ratio       : float — echoed back
        """
        # --- Resolve effective velocity and film temperature -----------
        if t_surf > 0:
            t_film = (t_surf + t_fluid) / 2.0
        else:
            t_film = t_fluid

        props = air_properties(t_film)
        nu = props["nu"]
        rho = props["rho"]

        # Buoyancy velocity for mixed_unknown regime
        if regime == "mixed_unknown":
            dt_buoy = delta_t_buoyancy
            if dt_buoy <= 0 and t_surf > 0:
                dt_buoy = abs(t_surf - t_fluid)
            v_buoy = BoundaryLayerCalculator._buoyancy_velocity(
                beta=props["beta"], delta_t=dt_buoy, char_length=x,
            ) if dt_buoy > 0 else 0.0
            U_eff = max(U, v_buoy, 0.1)  # floor at 0.1 m/s
        else:
            U_eff = max(U, 0.1)  # floor to prevent division by zero
            v_buoy = 0.0

        # --- Reynolds number (with leading-edge clamp) -----------------
        Re_x = U_eff * x / nu
        Re_x = max(Re_x, RE_X_MIN)

        # --- Skin friction and friction velocity -----------------------
        Cf, bl_regime = BoundaryLayerCalculator._skin_friction(Re_x)
        tau_w = 0.5 * Cf * rho * U_eff ** 2
        u_tau = math.sqrt(tau_w / rho)

        # --- First cell height -----------------------------------------
        y1 = y_plus_target * nu / u_tau
        y1 = max(y1, Y1_ABSOLUTE_MIN_M)

        # --- Boundary layer thickness ----------------------------------
        delta = BoundaryLayerCalculator._boundary_layer_thickness(Re_x, x)

        # --- Inflation layer count -------------------------------------
        n_layers = BoundaryLayerCalculator._inflation_layer_count(
            y1, delta, growth_ratio,
        )

        # --- Outermost prism height ------------------------------------
        y_last = y1 * growth_ratio ** (n_layers - 1)

        # --- Total prism stack height ----------------------------------
        if growth_ratio > 1.0:
            total_prism = y1 * (growth_ratio ** n_layers - 1.0) / (
                growth_ratio - 1.0
            )
        else:
            total_prism = y1 * n_layers

        # --- Surface mesh constraint -----------------------------------
        if regime == "mixed_unknown":
            # BL fraction constraint only — no prism AR coupling.
            # Isotropic tets to wall; dx from BL fraction ensures the
            # thermal mesh is compatible with the flow solver.
            dx_bl = bl_fraction * delta
            dx_surface = dx_bl
            constraint_label = "bl_fraction"
        else:
            # External forced: constrain by last-layer transition and BL
            dx_transition = ar_transition * y_last
            dx_bl = bl_fraction * delta
            dx_surface = min(dx_transition, dx_bl)
            if dx_transition <= dx_bl:
                constraint_label = "prism_transition"
            else:
                constraint_label = "bl_fraction"

        # --- Wall temperature gradient (for prism decision) -----------
        # dT/dy|_wall = q'' / k_air = h * (T_s - T_f) / k_air
        # Use a simple estimate: h ~ Nu * k_air / x, with Nu from
        # the flow regime.  For simplicity, use the empirical range
        # h ~ 5-15 W/m²K for low-Re buoyancy, h ~ 20-100 for forced.
        k_air = props["k_air"]
        dt_wall = abs(t_surf - t_fluid) if t_surf > 0 else 0.0
        if dt_wall > 0 and k_air > 0:
            # Estimate h from Nu correlation at the computed Re
            if Re_x < RE_X_CRIT:
                Nu = 0.664 * Re_x ** 0.5 * 0.71 ** (1.0 / 3.0)
            else:
                Nu = 0.037 * Re_x ** 0.8 * 0.71 ** (1.0 / 3.0)
            h_est = Nu * k_air / max(x, 0.001)
            wall_gradient = h_est * dt_wall / k_air
        else:
            wall_gradient = 0.0

        # --- Prism recommendation logic --------------------------------
        # Criterion 1: y1 vs dx_surface — if y1 > dx_surface, tets
        # already place multiple cells within the "first cell" zone.
        y1_exceeds_dx = (y1 > dx_surface)

        # Criterion 2: wall gradient magnitude — below threshold,
        # tet non-orthogonality error is within physics uncertainty.
        low_gradient = (wall_gradient < _WALL_GRADIENT_PRISM_THRESHOLD)

        if regime == "mixed_unknown":
            # Mixed/unknown: prisms generally not recommended
            prisms_rec = False
            if y1_exceeds_dx and low_gradient:
                prism_reason = (
                    f"Not recommended: y1={y1*1000:.1f}mm > dx={dx_surface*1000:.1f}mm "
                    f"(tets over-resolve wall region), wall gradient "
                    f"{wall_gradient:.0f} K/m is below {_WALL_GRADIENT_PRISM_THRESHOLD:.0e} K/m "
                    f"threshold (tet non-orthogonality error within physics "
                    f"uncertainty of h). Use isotropic tets to wall."
                )
            elif y1_exceeds_dx:
                prism_reason = (
                    f"Not recommended: y1={y1*1000:.1f}mm > dx={dx_surface*1000:.1f}mm "
                    f"(tets over-resolve wall region). Use isotropic tets to wall."
                )
            elif low_gradient:
                prism_reason = (
                    f"Not recommended: wall gradient {wall_gradient:.0f} K/m is "
                    f"below {_WALL_GRADIENT_PRISM_THRESHOLD:.0e} K/m threshold. "
                    f"Tet non-orthogonality error within physics uncertainty."
                )
            else:
                # Unusual case: mixed_unknown with high gradient and small y1.
                # Could occur with high dT in a very small gap.
                prisms_rec = True
                prism_reason = (
                    f"Recommended: wall gradient {wall_gradient:.0f} K/m exceeds "
                    f"threshold and y1={y1*1000:.2f}mm < dx={dx_surface*1000:.1f}mm. "
                    f"Consider low-AR prisms (AR <= {ar_max_prism:.0f})."
                )
        else:
            # External forced: prisms recommended (high-Re, anisotropic BL)
            prisms_rec = True
            prism_reason = (
                f"Recommended: external forced flow at Re={Re_x:.0f}, "
                f"wall gradient {wall_gradient:.0f} K/m. Prism layers "
                f"exploit directional anisotropy of attached boundary layer."
            )

        # --- Volume expansion ratio (prism-to-tet transition) ----------
        # ER_v ~ (dx_surface / y_last)^3 for the volume jump from
        # the outermost prism cell to the first isotropic tet.
        if y_last > 0 and prisms_rec:
            er_v = (dx_surface / y_last) ** 3
            if er_v <= ER_V_STABLE:
                er_v_status = "stable"
            elif er_v <= ER_V_MARGINAL:
                er_v_status = "marginal"
            else:
                er_v_status = "unstable"
        else:
            er_v = 0.0
            er_v_status = "n/a"

        return {
            "y1_mm": y1 * 1000.0,
            "y1_m": y1,
            "n_layers": n_layers,
            "delta_mm": delta * 1000.0,
            "y_last_mm": y_last * 1000.0,
            "total_prism_mm": total_prism * 1000.0,
            "max_dx_surface_mm": dx_surface * 1000.0,
            "surface_constraint": constraint_label,
            "prisms_recommended": prisms_rec,
            "prism_reason": prism_reason,
            "wall_gradient_Km": wall_gradient,
            "volume_expansion_ratio": er_v,
            "er_v_status": er_v_status,
            "U_effective": U_eff,
            "V_buoyancy": v_buoy,
            "Re_x": Re_x,
            "Cf": Cf,
            "u_tau": u_tau,
            "regime": regime,
            "bl_regime": bl_regime,
            "y_plus_target": y_plus_target,
            "growth_ratio": growth_ratio,
        }
