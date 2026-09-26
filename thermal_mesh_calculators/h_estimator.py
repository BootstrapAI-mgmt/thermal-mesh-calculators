"""
Convective HTC Estimator
=========================

Estimates the convective heat transfer coefficient (h) from first
principles for three regimes:

1.  Forced convection (flat plate correlation)
    For zones with known airflow velocity (cooling pack downstream,
    external aero, clutch outlet downstream).

2.  Natural convection (buoyancy-driven)
    For quiescent / dead zones where the only driver is the temperature
    difference between the surface and surrounding air.  Depends on
    surface orientation (vertical, horizontal hot-side-up, horizontal
    hot-side-down).

3.  Mixed convection (Richardson number regime map)
    For intermediate cases where both forced and buoyant effects are
    significant.  The Richardson number Ri = Gr / Re^2 determines
    which regime dominates.

External h cap:
    In a typical off-road vehicle at thermal analysis speeds (~10 mph),
    the practical maximum external h is about 60 W/m^2K.  The absolute
    cap is 100 W/m^2K — anything higher indicates internal flow or a
    correlation extrapolation error.

Air properties:
    Evaluated at the film temperature T_film = (T_surf + T_fluid) / 2.
    Uses polynomial fits valid from 250–700 K (AIR_PROPERTY_RANGE_K).
    Outside that range the properties are those at the nearer end of it;
    air_properties() says so (in_range False) and estimate_h() reports
    film_in_range.  Pure Python, no dependencies.

Units:  SI throughout (m, s, K, W).
"""

from typing import Optional

from thermal_mesh_calculators._guards import (
    require_non_negative,
    require_positive,
    require_temperatures,
)


# ---------------------------------------------------------------------------
#  Air properties at film temperature (polynomial fits, 250–700 K)
# ---------------------------------------------------------------------------

# The film temperatures (K) the air-property fits below are valid over.
AIR_PROPERTY_RANGE_K = (250.0, 700.0)


def film_in_range(t_film: float) -> bool:
    """True when t_film (K) lies within AIR_PROPERTY_RANGE_K."""
    low, high = AIR_PROPERTY_RANGE_K
    return low <= t_film <= high


def air_properties(t_film: float) -> dict:
    """
    Thermophysical properties of air at atmospheric pressure.

    Parameters
    ----------
    t_film : float — film temperature (K)

    Returns
    -------
    dict with keys:
        k_air  : float — thermal conductivity (W/m K)
        nu     : float — kinematic viscosity (m^2/s)
        Pr     : float — Prandtl number
        beta   : float — volumetric expansion coefficient (1/K)
        rho    : float — density (kg/m^3)
        t_film_K : float — the film temperature given (K)
        t_eval_K : float — the temperature the fits were evaluated at: the
                           film temperature clamped to AIR_PROPERTY_RANGE_K
        in_range : bool  — False when the properties are extrapolated
                           (t_eval_K differs from t_film_K)

    Raises ValueError for t_film not > 0.
    """
    require_positive("t_film", t_film, "K")
    low, high = AIR_PROPERTY_RANGE_K
    T = max(low, min(high, t_film))

    # Thermal conductivity (W/m K) — linear fit
    k_air = 0.0241 + 7.0e-5 * (T - 300.0)

    # Kinematic viscosity (m^2/s) — power-law fit
    nu = 1.5e-5 * (T / 300.0) ** 1.7

    # Prandtl number — nearly constant for air
    Pr = 0.71

    # Ideal gas: beta = 1/T
    beta = 1.0 / T

    # Density (ideal gas at 1 atm)
    rho = 101325.0 / (287.05 * T)

    return {
        "k_air": k_air,
        "nu": nu,
        "Pr": Pr,
        "beta": beta,
        "rho": rho,
        "t_film_K": t_film,
        "t_eval_K": T,
        "in_range": film_in_range(t_film),
    }


# ---------------------------------------------------------------------------
#  External h cap
# ---------------------------------------------------------------------------

H_EXTERNAL_MAX = 100.0   # W/m^2 K absolute cap for external surfaces
H_EXTERNAL_TYPICAL = 60.0  # W/m^2 K typical max at ~10 mph


def _cap_h(h: float, cap: float = H_EXTERNAL_MAX) -> float:
    """Clamp h to the external cap."""
    return min(h, cap)


# ---------------------------------------------------------------------------
#  Solver advisory — steady-state vs transient recommendations
# ---------------------------------------------------------------------------

# Rayleigh number thresholds for turbulence onset by orientation
_RA_TURB = {
    "vertical": 1e9,
    "horizontal_up": 1e7,
    "horizontal_down": float("inf"),  # effectively always laminar
}

# Reynolds number threshold for forced turbulence
_RE_TURB = 5e5


def solver_advisory(
    forced_result: Optional[dict] = None,
    natural_result: Optional[dict] = None,
    regime: str = "natural",
    orientation: str = "vertical",
) -> dict:
    """
    Generate solver setting recommendations based on convection regime.

    Turbulent natural convection produces oscillating buoyant plumes that
    make steady-state solutions difficult to converge and potentially
    misleading.  This function flags those conditions and recommends
    solver approach.

    Returns
    -------
    dict with keys:
        steady_state_ok     : bool  — True if SS solve is appropriate
        transient_advisory  : bool  — True if transient solve recommended
        reason              : str   — human-readable explanation
        severity            : str   — "info", "caution", "warning"
        natural_regime      : str   — "laminar", "turbulent", or "negligible"
        forced_regime       : str   — "laminar", "turbulent", or "n/a"
        Ra                  : float — Rayleigh number (if natural conv present)
        Re                  : float — Reynolds number (if forced conv present)
    """
    nat_regime = "n/a"
    frc_regime = "n/a"
    Ra = 0.0
    Re = 0.0

    if natural_result is not None:
        nat_regime = natural_result.get("regime", "negligible")
        Ra = natural_result.get("Ra", 0.0)
    if forced_result is not None:
        frc_regime = forced_result.get("regime", "n/a")
        Re = forced_result.get("Re", 0.0)

    # Default: steady-state is fine
    ss_ok = True
    transient = False
    severity = "info"
    reason = "Laminar convection — steady-state solve is appropriate."

    # --- Check natural convection turbulence ---
    if nat_regime == "turbulent":
        Ra_threshold = _RA_TURB.get(orientation, 1e9)
        if regime in ("natural", "mixed"):
            # Turbulent natural conv is a primary heat transfer mode
            ss_ok = False
            transient = True
            severity = "warning"
            reason = (
                f"Turbulent natural convection (Ra={Ra:.2e}, "
                f"threshold={Ra_threshold:.0e} for {orientation}). "
                f"Buoyant plumes are inherently unsteady — steady-state "
                f"solver may not converge or may give misleading "
                f"time-averaged results. Recommend transient solve "
                f"with time-averaged output, or use relaxation with "
                f"caution in steady-state."
            )
        elif regime == "forced":
            # Turbulent natural conv exists but forced dominates
            severity = "caution"
            reason = (
                f"Forced convection dominates (Ri<0.1), but natural "
                f"convection is turbulent (Ra={Ra:.2e}). If forced flow "
                f"is intermittent or uncertain, buoyancy oscillations "
                f"could affect results. Monitor convergence."
            )

    # --- Check forced convection turbulence ---
    if frc_regime == "turbulent" and regime in ("forced", "mixed"):
        if ss_ok:
            # Turbulent forced is generally fine for SS — just note it
            severity = "info"
            reason = (
                f"Turbulent forced convection (Re={Re:.0f}). "
                f"Steady-state solve is appropriate — turbulent "
                f"forced convection correlations use time-averaged "
                f"Nusselt numbers."
            )

    # --- Check transitional / marginal cases ---
    if regime == "mixed":
        if not transient:
            severity = max(severity, "caution") if severity != "warning" else severity
            if severity == "info":
                severity = "caution"
            reason = (
                "Mixed convection regime (Ri between 0.1 and 10). "
                "Both forced and buoyant effects are significant. "
                "Verify convergence carefully in steady-state, or "
                "consider transient solve for more reliable results."
            )

    # --- Near-zero dT: negligible natural convection ---
    if nat_regime == "negligible" and regime == "natural":
        severity = "caution"
        reason = (
            "Near-zero temperature difference — negligible natural "
            "convection (h ≈ 0). Verify boundary conditions. If the "
            "part truly sees no convection, this may indicate an "
            "adiabatic or radiation-only surface."
        )

    return {
        "steady_state_ok": ss_ok,
        "transient_advisory": transient,
        "reason": reason,
        "severity": severity,
        "natural_regime": nat_regime,
        "forced_regime": frc_regime,
        "Ra": Ra,
        "Re": Re,
    }


# ---------------------------------------------------------------------------
#  Forced convection — flat plate
# ---------------------------------------------------------------------------

def forced_convection_flat_plate(
    velocity: float,
    char_length: float,
    t_surf: float,
    t_fluid: float,
    cap: float = H_EXTERNAL_MAX,
) -> dict:
    """
    Forced convection h from flat plate correlations.

    Laminar (Re < 5e5):  Nu = 0.664 * Re^0.5 * Pr^(1/3)
    Turbulent (Re >= 5e5): Nu = 0.037 * Re^0.8 * Pr^(1/3)

    Parameters
    ----------
    velocity : float      — freestream velocity (m/s)
    char_length : float   — characteristic length (m), typically part
                            dimension in flow direction
    t_surf : float        — surface temperature (K)
    t_fluid : float       — freestream temperature (K)
    cap : float           — maximum h to return (W/m^2 K)

    Returns
    -------
    dict with keys:
        h       : float — convective HTC (W/m^2 K), capped
        h_raw   : float — uncapped value
        Re      : float — Reynolds number
        Nu      : float — Nusselt number
        regime  : str   — "laminar" or "turbulent"
    """
    require_temperatures(t_surf=t_surf, t_fluid=t_fluid)
    require_non_negative("velocity", velocity, "m/s")
    require_positive("char_length", char_length, "m")
    t_film = (t_surf + t_fluid) / 2.0
    props = air_properties(t_film)

    Re = velocity * char_length / props["nu"]
    Pr_term = props["Pr"] ** (1.0 / 3.0)

    if Re < 5e5:
        Nu = 0.664 * (Re ** 0.5) * Pr_term
        regime = "laminar"
    else:
        Nu = 0.037 * (Re ** 0.8) * Pr_term
        regime = "turbulent"

    h_raw = Nu * props["k_air"] / char_length
    h = _cap_h(h_raw, cap)

    return {
        "h": h,
        "h_raw": h_raw,
        "Re": Re,
        "Nu": Nu,
        "regime": regime,
    }


# ---------------------------------------------------------------------------
#  Natural convection — vertical plate and horizontal plate
# ---------------------------------------------------------------------------

def natural_convection(
    t_surf: float,
    t_fluid: float,
    char_length: float,
    orientation: str = "vertical",
    cap: float = H_EXTERNAL_MAX,
) -> dict:
    """
    Natural convection h from buoyancy-driven correlations.

    Vertical plate (Churchill-Chu):
        Laminar  (Ra < 10^9):  Nu = 0.59 * Ra^0.25
        Turbulent (Ra >= 10^9): Nu = 0.10 * Ra^(1/3)

    Horizontal plate, hot side up:
        Laminar  (Ra < 10^7):  Nu = 0.54 * Ra^0.25
        Turbulent (Ra >= 10^7): Nu = 0.15 * Ra^(1/3)

    Horizontal plate, hot side down:
        Nu = 0.27 * Ra^0.25  (laminar only, much weaker)

    Parameters
    ----------
    t_surf : float        — surface temperature (K)
    t_fluid : float       — ambient temperature (K)
    char_length : float   — characteristic length (m)
                            Vertical: height.  Horizontal: area/perimeter.
    orientation : str     — "vertical", "horizontal_up", "horizontal_down"
    cap : float           — maximum h (W/m^2 K)

    Returns
    -------
    dict with keys:
        h       : float — convective HTC (W/m^2 K), capped
        h_raw   : float — uncapped value
        Ra      : float — Rayleigh number
        Gr      : float — Grashof number
        Nu      : float — Nusselt number
        regime  : str   — "laminar" or "turbulent"
    """
    require_temperatures(t_surf=t_surf, t_fluid=t_fluid)
    require_positive("char_length", char_length, "m")
    dt = abs(t_surf - t_fluid)
    if dt < 0.01:
        return {
            "h": 0.0, "h_raw": 0.0,
            "Ra": 0.0, "Gr": 0.0, "Nu": 0.0,
            "regime": "negligible",
        }

    t_film = (t_surf + t_fluid) / 2.0
    props = air_properties(t_film)
    g = 9.81

    Gr = g * props["beta"] * dt * char_length**3 / props["nu"]**2
    Ra = Gr * props["Pr"]

    if orientation == "vertical":
        if Ra < 1e9:
            Nu = 0.59 * (Ra ** 0.25)
            regime = "laminar"
        else:
            Nu = 0.10 * (Ra ** (1.0 / 3.0))
            regime = "turbulent"

    elif orientation == "horizontal_up":
        if Ra < 1e7:
            Nu = 0.54 * (Ra ** 0.25)
            regime = "laminar"
        else:
            Nu = 0.15 * (Ra ** (1.0 / 3.0))
            regime = "turbulent"

    elif orientation == "horizontal_down":
        Nu = 0.27 * (Ra ** 0.25)
        regime = "laminar"

    else:
        raise ValueError(
            f"Unknown orientation '{orientation}'. Use 'vertical', "
            f"'horizontal_up', or 'horizontal_down'."
        )

    h_raw = Nu * props["k_air"] / char_length
    h = _cap_h(h_raw, cap)

    return {
        "h": h,
        "h_raw": h_raw,
        "Ra": Ra,
        "Gr": Gr,
        "Nu": Nu,
        "regime": regime,
    }


# ---------------------------------------------------------------------------
#  Richardson number — regime identification
# ---------------------------------------------------------------------------

def richardson_number(
    velocity: float,
    t_surf: float,
    t_fluid: float,
    char_length: float,
) -> dict:
    """
    Richardson number Ri = Gr / Re^2 to identify the dominant regime.

    Ri << 0.1  :  forced convection dominates
    Ri >  10   :  natural convection dominates
    0.1 <= Ri <= 10 :  mixed convection — both matter

    Parameters
    ----------
    velocity : float      — freestream velocity (m/s)
    t_surf : float        — surface temperature (K)
    t_fluid : float       — ambient temperature (K)
    char_length : float   — characteristic length (m)

    Returns
    -------
    dict with keys:
        Ri      : float — Richardson number
        Gr      : float — Grashof number
        Re      : float — Reynolds number
        regime  : str   — "forced", "natural", or "mixed"
    """
    require_temperatures(t_surf=t_surf, t_fluid=t_fluid)
    require_non_negative("velocity", velocity, "m/s")
    require_positive("char_length", char_length, "m")
    dt = abs(t_surf - t_fluid)
    t_film = (t_surf + t_fluid) / 2.0
    props = air_properties(t_film)
    g = 9.81

    Re = velocity * char_length / props["nu"] if velocity > 0 else 0.0
    Gr = g * props["beta"] * dt * char_length**3 / props["nu"]**2

    if Re < 1.0:
        # Essentially zero velocity → natural convection
        return {
            "Ri": float("inf"),
            "Gr": Gr,
            "Re": Re,
            "regime": "natural",
        }

    Ri = Gr / (Re ** 2)

    if Ri < 0.1:
        regime = "forced"
    elif Ri > 10.0:
        regime = "natural"
    else:
        regime = "mixed"

    return {
        "Ri": Ri,
        "Gr": Gr,
        "Re": Re,
        "regime": regime,
    }


# ---------------------------------------------------------------------------
#  Combined h estimator — picks the right correlation automatically
# ---------------------------------------------------------------------------

def estimate_h(
    velocity: float,
    t_surf: float,
    t_fluid: float,
    char_length: float,
    orientation: str = "vertical",
    cap: float = H_EXTERNAL_MAX,
) -> dict:
    """
    Estimate h by computing both forced and natural convection, then
    selecting or combining based on the Richardson number.

    For mixed regime, uses the asymptotic combination:
        h_mixed = (h_forced^3 + h_natural^3)^(1/3)
    which is standard for assisting mixed convection.

    Parameters
    ----------
    velocity : float      — freestream velocity (m/s), 0 for dead zones
    t_surf : float        — surface temperature (K)
    t_fluid : float       — ambient temperature (K)
    char_length : float   — characteristic length (m)
    orientation : str     — surface orientation for natural conv
    cap : float           — external h cap (W/m^2 K)

    Returns
    -------
    dict with keys:
        h               : float — recommended HTC (W/m^2 K), capped
        h_forced        : float — forced convection estimate
        h_natural       : float — natural convection estimate
        regime          : str   — "forced", "natural", or "mixed"
        Ri              : float — Richardson number
        dominant_mode   : str   — description of what's driving h
        t_film_K        : float — film temperature (K) the air properties
                                  were wanted at
        film_in_range   : bool  — False when that lies outside
                                  AIR_PROPERTY_RANGE_K, so h rests on
                                  extrapolated (clamped) air properties

    Raises ValueError for a temperature or char_length not > 0, or a
    negative velocity.
    """
    require_temperatures(t_surf=t_surf, t_fluid=t_fluid)
    require_non_negative("velocity", velocity, "m/s")
    require_positive("char_length", char_length, "m")
    t_film = (t_surf + t_fluid) / 2.0
    nat = natural_convection(
        t_surf, t_fluid, char_length, orientation, cap=999.0,
    )
    h_nat = nat["h_raw"]

    if velocity < 0.01:
        # No meaningful flow → pure natural convection
        h = _cap_h(h_nat, cap)
        adv = solver_advisory(
            natural_result=nat, regime="natural",
            orientation=orientation,
        )
        return {
            "h": h,
            "h_forced": 0.0,
            "h_natural": h_nat,
            "regime": "natural",
            "natural_regime": nat["regime"],
            "forced_regime": "n/a",
            "Ri": float("inf"),
            "Ra": nat["Ra"],
            "Re": 0.0,
            "dominant_mode": f"natural convection ({orientation})",
            "solver_advisory": adv,
            "t_film_K": t_film,
            "film_in_range": film_in_range(t_film),
        }

    frc = forced_convection_flat_plate(
        velocity, char_length, t_surf, t_fluid, cap=999.0,
    )
    h_frc = frc["h_raw"]

    ri = richardson_number(velocity, t_surf, t_fluid, char_length)

    if ri["regime"] == "forced":
        h_raw = h_frc
        dominant = f"forced convection (Re={frc['Re']:.0f})"
    elif ri["regime"] == "natural":
        h_raw = h_nat
        dominant = f"natural convection ({orientation}, Ra={nat['Ra']:.1e})"
    else:
        # Mixed convection — check for opposing vs assisting flow.
        # Opposing: forced flow opposes buoyancy (e.g. downward forced flow
        # over a hot vertical surface where buoyancy drives upward).
        # Detected when the hot surface is below the flow (horizontal_down
        # with forced flow) or when explicitly flagged.  For vertical
        # surfaces, opposing flow occurs when forced flow is downward
        # against the upward buoyant plume.
        #
        # Assisting: h_mixed = (h_f^3 + h_n^3)^(1/3)
        # Opposing:  h_mixed = max(|h_f^3 - h_n^3|^(1/3), k_air/L)
        #
        # The k_air/L lower bound prevents singularity when the two modes
        # nearly cancel, ensuring at least pure-conduction minimum h.
        is_opposing = (orientation == "horizontal_down")

        if is_opposing:
            props = air_properties(t_film)
            h_conduction_min = props["k_air"] / char_length
            h_raw = max(
                abs(h_frc**3 - h_nat**3) ** (1.0 / 3.0),
                h_conduction_min,
            )
            dominant = f"mixed-opposing (Ri={ri['Ri']:.2f})"
        else:
            # Assisting: standard asymptotic combination
            h_raw = (h_frc**3 + h_nat**3) ** (1.0 / 3.0)
            dominant = f"mixed-assisting (Ri={ri['Ri']:.2f})"

    h = _cap_h(h_raw, cap)

    adv = solver_advisory(
        forced_result=frc, natural_result=nat,
        regime=ri["regime"], orientation=orientation,
    )

    return {
        "h": h,
        "h_forced": h_frc,
        "h_natural": h_nat,
        "regime": ri["regime"],
        "natural_regime": nat["regime"],
        "forced_regime": frc["regime"],
        "Ri": ri["Ri"],
        "Ra": nat["Ra"],
        "Re": frc["Re"],
        "dominant_mode": dominant,
        "solver_advisory": adv,
        "t_film_K": t_film,
        "film_in_range": film_in_range(t_film),
    }


# ---------------------------------------------------------------------------
#  Quick velocity-based gut check
# ---------------------------------------------------------------------------

def h_from_velocity(velocity_mph: float, char_length_mm: float = 100.0,
                    t_surf_C: float = 100.0, t_fluid_C: float = 40.0) -> dict:
    """
    Quick gut-check h estimate from vehicle speed and part size.

    Convenience wrapper using common units (mph, mm, Celsius).

    Parameters
    ----------
    velocity_mph : float     — vehicle / airflow speed (mph)
    char_length_mm : float   — part dimension in flow direction (mm)
    t_surf_C : float         — surface temperature (Celsius)
    t_fluid_C : float        — ambient temperature (Celsius)

    Returns
    -------
    dict with h and supporting info
    """
    velocity_ms = velocity_mph * 0.44704  # mph → m/s
    char_length_m = char_length_mm / 1000.0
    t_surf_K = t_surf_C + 273.15
    t_fluid_K = t_fluid_C + 273.15

    result = estimate_h(
        velocity=velocity_ms,
        t_surf=t_surf_K,
        t_fluid=t_fluid_K,
        char_length=char_length_m,
        orientation="vertical",
    )
    result["velocity_mph"] = velocity_mph
    result["velocity_ms"] = velocity_ms
    result["char_length_mm"] = char_length_mm

    return result
