"""
Internal-Flow Heat Transfer in Ducts
====================================

The convective coefficient h on the wetted wall of a duct — a cooling
passage, an air duct, a pipe — from the textbook correlations for flow
inside a tube.  The flat-plate correlations in h_estimator.py describe a
boundary layer that grows freely along an open surface; inside a duct the
boundary layers meet, and h settles to a value set by the duct's hydraulic
diameter rather than by the distance from a leading edge.

Hydraulic diameter
    Dh = 4 A / P, the flow area A over the wetted perimeter P.  For a
    circular tube it is the diameter; for a rectangle a x b it is
    2ab / (a + b).  Re and Nu are both based on it: Re = V Dh / nu,
    h = Nu k / Dh.

Regimes (Re based on Dh and the mean velocity)
    Re <= 2300             laminar
    2300 < Re < 10000      transitional: reported as such, with an
                           advisory; see below
    Re >= 10000            turbulent

Laminar, fully developed (circular tube)
    Nu = 3.66 with the wall at one temperature; Nu = 4.36 with a uniform
    wall heat flux.  Neither depends on Re or Pr.  For a non-circular
    section of the same Dh the fully developed laminar Nu depends on the
    shape of the section, so these are the circular-tube values.

Laminar, thermal entry (Hausen), wall at one temperature
    Nu = 3.66 + 0.0668 Gz / (1 + 0.04 Gz^(2/3)),  Gz = Re Pr Dh / L
    the mean Nu over a heated length L, for a velocity profile that is
    already developed while the temperature profile develops (the thermal
    entry problem).  It tends to 3.66 for a long duct and rises toward the
    inlet.  With air the velocity profile develops alongside the
    temperature profile, which raises the true mean Nu near the inlet
    above this value.  Used when a duct length is given; for a uniform
    heat flux the package does not model the entry and keeps 4.36, which
    the mean over any finite length exceeds.

Turbulent, Gnielinski with the Petukhov friction factor (the default)
    f  = (0.790 ln Re - 1.64)^-2                     (smooth tube)
    Nu = (f/8)(Re - 1000) Pr / (1 + 12.7 (f/8)^(1/2) (Pr^(2/3) - 1))
    stated for 3000 <= Re <= 5e6 and 0.5 <= Pr <= 2000.

Turbulent, Dittus-Boelter (the simpler alternative)
    Nu = 0.023 Re^0.8 Pr^n, n = 0.4 when the fluid is heated (wall hotter
    than the fluid), 0.3 when it is cooled; stated for Re >= 10000,
    0.6 <= Pr <= 160 and L/Dh >= 10.

Both turbulent correlations are for fully developed flow: in a duct
shorter than 10 hydraulic diameters the entrance raises the mean h, and
the result says so.

The transition band
    Between Re = 2300 and 10000 the flow may be laminar, intermittent or
    turbulent, and no correlation is reliable.  The package does not blend
    the two: it evaluates the laminar value and the chosen turbulent
    correlation at the actual Re, reports both, takes the higher — a
    higher h gives the finer, safer mesh — and says so in the advisory.

Air properties
    From h_estimator.air_properties() at the film temperature
    T_film = (T_wall + T_fluid) / 2, with T_fluid the bulk (mean) fluid
    temperature — the same evaluation the package's external correlations
    use, so the same 250-700 K fit range applies and a film temperature
    outside it is reported.  Pr is 0.71 throughout that range.

The external h cap (H_EXTERNAL_MAX) does not apply: the package treats
every internal flow as uncapped.

Units:  SI throughout (m, s, K, W).
"""

import math
from typing import Optional

from thermal_mesh_calculators._guards import require_positive, require_temperatures
from thermal_mesh_calculators.h_estimator import AIR_PROPERTY_RANGE_K, air_properties

# Regime boundaries (Re based on the hydraulic diameter).
LAMINAR_REYNOLDS_MAX = 2300.0
TURBULENT_REYNOLDS_MIN = 10000.0

# Fully developed laminar Nusselt numbers for a circular tube.
NU_LAMINAR_CONSTANT_TEMPERATURE = 3.66
NU_LAMINAR_CONSTANT_HEAT_FLUX = 4.36

WALL_CONDITIONS = ("constant_temperature", "constant_heat_flux")
TURBULENT_CORRELATIONS = ("gnielinski", "dittus_boelter")

# The ranges the turbulent correlations are stated for.
GNIELINSKI_REYNOLDS_RANGE = (3000.0, 5.0e6)
DITTUS_BOELTER_REYNOLDS_MIN = 10000.0

# Below this length (in hydraulic diameters) the fully developed turbulent
# correlations underestimate the mean h.
FULLY_DEVELOPED_LENGTH_DIAMETERS = 10.0


# ---------------------------------------------------------------------------
#  Geometry
# ---------------------------------------------------------------------------

def hydraulic_diameter(area: float, wetted_perimeter: float) -> float:
    """
    Hydraulic diameter of a duct section, Dh = 4 A / P.

    Parameters
    ----------
    area : float
        Flow area of the section (m^2).
    wetted_perimeter : float
        Perimeter of the section wetted by the flow (m).

    Returns
    -------
    float — the hydraulic diameter (m).  For a circle it is the diameter;
    for a rectangle a x b it is 2ab / (a + b).

    Raises
    ------
    ValueError
        area or wetted_perimeter not > 0.
    """
    require_positive("area", area, "m^2")
    require_positive("wetted_perimeter", wetted_perimeter, "m")
    return 4.0 * area / wetted_perimeter


# ---------------------------------------------------------------------------
#  The correlations, at a given Re and Pr
# ---------------------------------------------------------------------------

def friction_factor_petukhov(reynolds: float) -> float:
    """
    Darcy friction factor of a smooth tube, Petukhov's form.

        f = (0.790 ln Re - 1.64)^-2

    Parameters
    ----------
    reynolds : float — Reynolds number (stated for 3000 <= Re <= 5e6)

    Returns
    -------
    float — f (dimensionless).
    """
    require_positive("reynolds", reynolds)
    return (0.790 * math.log(reynolds) - 1.64) ** -2


def nusselt_gnielinski(reynolds: float, prandtl: float) -> float:
    """
    Turbulent Nusselt number, Gnielinski's correlation.

        Nu = (f/8)(Re - 1000) Pr / (1 + 12.7 (f/8)^(1/2) (Pr^(2/3) - 1))

    with f from friction_factor_petukhov().

    Parameters
    ----------
    reynolds : float — Reynolds number (stated for 3000 <= Re <= 5e6)
    prandtl : float  — Prandtl number (stated for 0.5 <= Pr <= 2000)

    Returns
    -------
    float — Nu (dimensionless), fully developed flow.
    """
    require_positive("prandtl", prandtl)
    f8 = friction_factor_petukhov(reynolds) / 8.0
    return (f8 * (reynolds - 1000.0) * prandtl
            / (1.0 + 12.7 * math.sqrt(f8) * (prandtl ** (2.0 / 3.0) - 1.0)))


def nusselt_dittus_boelter(reynolds: float, prandtl: float,
                           heating: bool = True) -> float:
    """
    Turbulent Nusselt number, the Dittus-Boelter correlation.

        Nu = 0.023 Re^0.8 Pr^n,  n = 0.4 heating the fluid, 0.3 cooling it

    Parameters
    ----------
    reynolds : float — Reynolds number (stated for Re >= 10000)
    prandtl : float  — Prandtl number (stated for 0.6 <= Pr <= 160)
    heating : bool   — True when the wall is hotter than the fluid

    Returns
    -------
    float — Nu (dimensionless), fully developed flow.
    """
    require_positive("reynolds", reynolds)
    require_positive("prandtl", prandtl)
    n = 0.4 if heating else 0.3
    return 0.023 * reynolds ** 0.8 * prandtl ** n


def nusselt_hausen(reynolds: float, prandtl: float,
                   length_over_diameter: float) -> float:
    """
    Mean laminar Nusselt number over a heated length, Hausen's correlation
    for the thermal entry with the wall at one temperature.

        Gz = Re Pr / (L/Dh)
        Nu = 3.66 + 0.0668 Gz / (1 + 0.04 Gz^(2/3))

    Parameters
    ----------
    reynolds : float             — Reynolds number (laminar, Re <= 2300)
    prandtl : float              — Prandtl number
    length_over_diameter : float — heated length over hydraulic diameter

    Returns
    -------
    float — the mean Nu over the length; it tends to 3.66 for a long duct.
    """
    require_positive("reynolds", reynolds)
    require_positive("prandtl", prandtl)
    require_positive("length_over_diameter", length_over_diameter)
    graetz = reynolds * prandtl / length_over_diameter
    return (NU_LAMINAR_CONSTANT_TEMPERATURE
            + 0.0668 * graetz / (1.0 + 0.04 * graetz ** (2.0 / 3.0)))


def _check_choice(name: str, value, allowed) -> None:
    if value not in allowed:
        raise ValueError(
            f"{name} must be one of {', '.join(repr(a) for a in allowed)}, "
            f"got {value!r}")


def _laminar(prandtl: float, reynolds: float, wall_condition: str,
             length_over_diameter: Optional[float]):
    """The laminar Nu, the correlation's name and any advisory."""
    if wall_condition == "constant_heat_flux":
        note = ""
        if length_over_diameter is not None:
            note = ("The thermal entry is not modelled for a uniform wall "
                    "heat flux: Nu = 4.36 is the fully developed value, which "
                    "the mean over any finite length exceeds.")
        return NU_LAMINAR_CONSTANT_HEAT_FLUX, "laminar_fully_developed", note
    if length_over_diameter is None:
        return NU_LAMINAR_CONSTANT_TEMPERATURE, "laminar_fully_developed", ""
    return (nusselt_hausen(reynolds, prandtl, length_over_diameter),
            "hausen",
            "Hausen assumes a developed velocity profile; with air the "
            "velocity develops alongside the temperature, so near the inlet "
            "the true mean Nu is higher.")


def _turbulent(reynolds: float, prandtl: float, correlation: str,
               heating: bool):
    """The turbulent Nu, the friction factor used (or None) and the
    correlation's range note for this Re."""
    if correlation == "gnielinski":
        low, high = GNIELINSKI_REYNOLDS_RANGE
        note = ""
        if reynolds < low:
            note = (f"Gnielinski's correlation is stated for Re >= {low:.0f}; "
                    f"Re = {reynolds:.0f} lies below it.")
        elif reynolds > high:
            note = (f"Gnielinski's correlation is stated for Re <= {high:.0e}; "
                    f"Re = {reynolds:.3g} lies above it.")
        return (nusselt_gnielinski(reynolds, prandtl),
                friction_factor_petukhov(reynolds), note)
    note = ""
    if reynolds < DITTUS_BOELTER_REYNOLDS_MIN:
        note = (f"The Dittus-Boelter correlation is stated for Re >= "
                f"{DITTUS_BOELTER_REYNOLDS_MIN:.0f}; Re = {reynolds:.0f} lies "
                f"below it.")
    return nusselt_dittus_boelter(reynolds, prandtl, heating), None, note


def internal_flow_nusselt(
    reynolds: float,
    prandtl: float,
    wall_condition: str,
    length_over_diameter: Optional[float] = None,
    correlation: str = "gnielinski",
    heating: bool = True,
) -> dict:
    """
    Nusselt number for flow in a duct at a given Re and Pr: picks the
    regime, the correlation and the advisory.

    Parameters
    ----------
    reynolds : float
        Reynolds number based on the hydraulic diameter and the mean
        velocity.
    prandtl : float
        Prandtl number.
    wall_condition : str
        "constant_temperature" (the wall at one temperature) or
        "constant_heat_flux" (a uniform wall heat flux).  It sets the
        laminar value; the turbulent correlations do not distinguish.
    length_over_diameter : float, optional
        Heated length over hydraulic diameter.  None (the default) means
        fully developed flow.  Given, the laminar wall-at-one-temperature
        value is Hausen's mean over that length, and a turbulent duct
        shorter than 10 diameters is named in the advisory.
    correlation : str
        The turbulent correlation: "gnielinski" (the default, with the
        Petukhov friction factor) or "dittus_boelter".
    heating : bool
        True when the wall is hotter than the fluid (Dittus-Boelter's
        exponent).

    Returns
    -------
    dict with keys:
        Nu              : float — the Nusselt number used
        regime          : str   — "laminar", "transitional" or "turbulent"
        correlation     : str   — "laminar_fully_developed", "hausen",
                                  "gnielinski" or "dittus_boelter": the one
                                  whose value Nu is
        Nu_laminar      : float or None — the laminar value (laminar and
                                  transitional flow)
        Nu_turbulent    : float or None — the turbulent correlation's value
                                  (transitional and turbulent flow)
        friction_factor : float or None — the Petukhov f when Gnielinski's
                                  correlation was evaluated
        advisory        : str   — what the caller should know about the
                                  value; "" when it lies inside every
                                  stated range

    Raises
    ------
    ValueError
        reynolds, prandtl or length_over_diameter not > 0, or an unknown
        wall_condition or correlation.
    """
    require_positive("reynolds", reynolds)
    require_positive("prandtl", prandtl)
    if length_over_diameter is not None:
        require_positive("length_over_diameter", length_over_diameter)
    _check_choice("wall_condition", wall_condition, WALL_CONDITIONS)
    _check_choice("correlation", correlation, TURBULENT_CORRELATIONS)

    notes = []
    nu_lam = nu_turb = f = None
    if reynolds <= LAMINAR_REYNOLDS_MAX:
        regime = "laminar"
        nu_lam, used, note = _laminar(prandtl, reynolds, wall_condition,
                                      length_over_diameter)
        nu = nu_lam
        notes.append(note)
    else:
        nu_turb, f, range_note = _turbulent(reynolds, prandtl, correlation,
                                            heating)
        if reynolds < TURBULENT_REYNOLDS_MIN:
            regime = "transitional"
            nu_lam, lam_name, _ = _laminar(prandtl, reynolds, wall_condition,
                                           length_over_diameter)
            if nu_turb >= nu_lam:
                nu, used = nu_turb, correlation
            else:
                nu, used = nu_lam, lam_name
            notes.append(
                f"Transitional flow (Re = {reynolds:.0f}, between "
                f"{LAMINAR_REYNOLDS_MAX:.0f} and {TURBULENT_REYNOLDS_MIN:.0f}): "
                f"the flow may be laminar, intermittent or turbulent. Laminar "
                f"Nu = {nu_lam:.4g}, {correlation} Nu = {nu_turb:.4g}; the "
                f"higher ({used}) is used, not a blend, because a higher h "
                f"gives the finer mesh.")
            notes.append(range_note)
        else:
            regime = "turbulent"
            nu, used = nu_turb, correlation
            notes.append(range_note)
        if (length_over_diameter is not None
                and length_over_diameter < FULLY_DEVELOPED_LENGTH_DIAMETERS):
            notes.append(
                f"The duct is {length_over_diameter:.3g} hydraulic diameters "
                f"long, shorter than 10 hydraulic diameters: the turbulent "
                f"correlations are for fully developed flow, and the entrance "
                f"raises the mean h above this value.")

    return {
        "Nu": nu,
        "regime": regime,
        "correlation": used,
        "Nu_laminar": nu_lam,
        "Nu_turbulent": nu_turb,
        "friction_factor": f,
        "advisory": " ".join(n for n in notes if n),
    }


# ---------------------------------------------------------------------------
#  h for air in a duct
# ---------------------------------------------------------------------------

def internal_flow_h(
    velocity: float,
    hydraulic_diameter: float,
    t_surf: float,
    t_fluid: float,
    wall_condition: str,
    length: Optional[float] = None,
    correlation: str = "gnielinski",
    cap: Optional[float] = None,
) -> dict:
    """
    Convective h on the wall of a duct carrying air.

    Parameters
    ----------
    velocity : float
        Mean velocity of the air in the duct (m/s), > 0.
    hydraulic_diameter : float
        Hydraulic diameter of the section (m), 4A/P (see
        hydraulic_diameter()).
    t_surf : float
        Wall temperature (K).
    t_fluid : float
        Bulk (mean) air temperature (K).
    wall_condition : str
        "constant_temperature" or "constant_heat_flux" (see
        internal_flow_nusselt()).
    length : float, optional
        Heated length of the duct (m).  None (the default) means fully
        developed flow.
    correlation : str
        The turbulent correlation: "gnielinski" (the default) or
        "dittus_boelter".
    cap : float, optional
        A ceiling on h (W/m^2 K).  None (the default) leaves h uncapped,
        as the package treats every internal flow.

    Returns
    -------
    dict with keys:
        h, h_raw        : float — the coefficient (W/m^2 K); h_raw is
                                  uncapped and h equals it unless cap is
                                  given
        Nu, Re, Pr      : float — Nusselt, Reynolds and Prandtl numbers
        regime, correlation, Nu_laminar, Nu_turbulent, friction_factor,
        advisory        : as internal_flow_nusselt(); the advisory also
                          names a film temperature outside the air-property
                          fits
        wall_condition  : str
        hydraulic_diameter_m : float
        length_m        : float or None
        velocity_ms     : float
        t_film_K        : float — film temperature the air properties were
                                  wanted at (K)
        film_in_range   : bool  — False when that lies outside
                                  AIR_PROPERTY_RANGE_K (the properties are
                                  then those at the nearer end of it)

    Raises
    ------
    ValueError
        velocity, hydraulic_diameter, length, cap or a temperature not
        > 0, or an unknown wall_condition or correlation.
    """
    require_positive("velocity", velocity, "m/s")
    require_positive("hydraulic_diameter", hydraulic_diameter, "m")
    require_temperatures(t_surf=t_surf, t_fluid=t_fluid)
    if length is not None:
        require_positive("length", length, "m")
    if cap is not None:
        require_positive("cap", cap, "W/m^2 K")

    t_film = (t_surf + t_fluid) / 2.0
    props = air_properties(t_film)
    reynolds = velocity * hydraulic_diameter / props["nu"]
    out = internal_flow_nusselt(
        reynolds, props["Pr"], wall_condition,
        length_over_diameter=(None if length is None
                              else length / hydraulic_diameter),
        correlation=correlation,
        heating=t_surf >= t_fluid,
    )
    advisory = out["advisory"]
    if not props["in_range"]:
        low, high = AIR_PROPERTY_RANGE_K
        film_note = (f"The air-property fits cover film temperatures of "
                     f"{low:.0f}-{high:.0f} K; this h (film temperature "
                     f"{t_film:.0f} K) used air properties evaluated at "
                     f"{props['t_eval_K']:.0f} K. Treat it as extrapolated.")
        advisory = (advisory + " " + film_note).strip()

    h_raw = out["Nu"] * props["k_air"] / hydraulic_diameter
    return {
        "h": h_raw if cap is None else min(h_raw, cap),
        "h_raw": h_raw,
        "Nu": out["Nu"],
        "Re": reynolds,
        "Pr": props["Pr"],
        "regime": out["regime"],
        "correlation": out["correlation"],
        "Nu_laminar": out["Nu_laminar"],
        "Nu_turbulent": out["Nu_turbulent"],
        "friction_factor": out["friction_factor"],
        "advisory": advisory,
        "wall_condition": wall_condition,
        "hydraulic_diameter_m": hydraulic_diameter,
        "length_m": length,
        "velocity_ms": velocity,
        "t_film_K": t_film,
        "film_in_range": props["in_range"],
    }


__all__ = [
    "LAMINAR_REYNOLDS_MAX",
    "TURBULENT_REYNOLDS_MIN",
    "NU_LAMINAR_CONSTANT_TEMPERATURE",
    "NU_LAMINAR_CONSTANT_HEAT_FLUX",
    "WALL_CONDITIONS",
    "TURBULENT_CORRELATIONS",
    "GNIELINSKI_REYNOLDS_RANGE",
    "DITTUS_BOELTER_REYNOLDS_MIN",
    "FULLY_DEVELOPED_LENGTH_DIAMETERS",
    "hydraulic_diameter",
    "friction_factor_petukhov",
    "nusselt_gnielinski",
    "nusselt_dittus_boelter",
    "nusselt_hausen",
    "internal_flow_nusselt",
    "internal_flow_h",
]
