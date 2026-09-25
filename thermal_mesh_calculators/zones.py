"""
Convection Zone Lookup & Spatial Gradient Estimator
=====================================================

Automation of mesh sizing for hundreds of parts requires that the two
remaining "circular" inputs — convective HTC (h) and spatial temperature
gradient (dT/dx) — be estimatable without prior CFD results.

Convection Zones
----------------
In an automotive underhood environment, airflow patterns are dominated
by a small number of distinct regimes that can be characterised by
vehicle region and local flow condition.

Each zone now carries:
    - h_low / h_high  — static h range (legacy, still usable)
    - velocity_ms     — expected local air velocity for forced zones
    - t_air_C         — typical air temperature range in the zone
    - orientation     — dominant surface orientation for natural conv
    - is_internal     — True for internal flows (exhaust, clutch)
                        that bypass the external h cap

The h_estimator module can compute h from these zone parameters using
proper correlations (flat plate forced, buoyancy natural, Richardson
mixed).  The static h_high values remain as a quick-lookup fallback.

External h cap:
    For external surfaces at typical thermal analysis speeds (~10 mph),
    the practical max is ~60 W/m^2K.  Absolute cap 100 W/m^2K.
    Internal flows (exhaust gas, clutch outlet) are not capped.

Zone Temperature Map (Side-by-Side Off-Road Vehicle)
-----------------------------------------------------
Based on thermal mapping at steady-state max-power conditions:

    Front end:
        Cooling pack outlet:       80–100 C
        Edges (aero):              40–50 C

    Cabin:
        Dead zone / tunnel:        40–50 C

    Engine bay:
        Edges:                     40–50 C
        Beside engine:             60–70 C
        Below engine:              50–60 C
        Above engine:              80–100 C
        Beside exhaust:            90–150 C  (with shielding)
        Below exhaust:             90–150 C  (with shielding)
        Above exhaust:             150–350 C (with shielding)
        Without shielding:         125–500 C (distance-dependent)

    Clutch outlet:                 90–110 C

    Cargo / shield:
        Shield gap:                varies (solved by shield calculator)
        Cargo bed exterior:        40–80 C

Units:  h in W/m^2 K,  velocity in m/s,  temperature in C (zone map)
        and K (internal calculations),  gradient in K/m
"""

from thermal_mesh_calculators._guards import require_positive


# ---------------------------------------------------------------------------
#  Convection zone definitions
# ---------------------------------------------------------------------------

# Each zone now includes:
#   h_low, h_high   — static h range for quick lookup (W/m^2 K)
#   regime          — "forced", "natural", or "mixed"
#   velocity_ms     — expected air velocity (m/s); 0 for dead zones
#   t_air_C_low     — lower bound of typical air temp in zone (C)
#   t_air_C_high    — upper bound of typical air temp in zone (C)
#   orientation     — surface orientation for natural conv calc
#   is_internal     — True if internal flow (not subject to external h cap)
#   notes           — engineering rationale

CONVECTION_ZONES = {
    # --- Front End ---
    "cooling_pack_downstream": {
        "h_low": 40.0,
        "h_high": 100.0,           # capped from 120 to external max
        "regime": "forced",
        "velocity_ms": 3.0,        # ~2-4 m/s through pack at 10 mph
        "t_air_C_low": 80.0,
        "t_air_C_high": 100.0,     # HX outlet temp
        "orientation": "vertical",
        "is_internal": False,
        "notes": (
            "Airflow forced through cooling pack and across components "
            "immediately behind it.  Air temp is HX outlet (80-100 C)."
        ),
    },
    "front_end_edges": {
        "h_low": 30.0,
        "h_high": 80.0,
        "regime": "forced",
        "velocity_ms": 4.5,        # 10 mph = 4.47 m/s, full aero
        "t_air_C_low": 40.0,
        "t_air_C_high": 50.0,      # ambient
        "orientation": "vertical",
        "is_internal": False,
        "notes": (
            "External aerodynamic flow at vehicle edges.  Air is near "
            "ambient (40-50 C).  Max velocity at analysis speed."
        ),
    },
    "front_end_dead_zone": {
        "h_low": 5.0,
        "h_high": 15.0,
        "regime": "natural",
        "velocity_ms": 0.0,
        "t_air_C_low": 40.0,
        "t_air_C_high": 60.0,
        "orientation": "vertical",
        "is_internal": False,
        "notes": (
            "Recirculation regions between flow paths.  Natural "
            "convection only."
        ),
    },

    # --- Cabin ---
    "cabin_dead_zone": {
        "h_low": 3.0,
        "h_high": 10.0,
        "regime": "natural",
        "velocity_ms": 0.0,
        "t_air_C_low": 40.0,
        "t_air_C_high": 50.0,
        "orientation": "vertical",
        "is_internal": False,
        "notes": "Cabin area, mostly stagnant.  Natural convection only.",
    },
    "cabin_tunnel": {
        "h_low": 3.0,
        "h_high": 8.0,
        "regime": "natural",
        "velocity_ms": 0.0,
        "t_air_C_low": 40.0,
        "t_air_C_high": 50.0,
        "orientation": "horizontal_up",     # tunnel floor is hot-side-up
        "is_internal": False,
        "notes": "Lower tunnel beneath cabin floor.  Nearly stagnant.",
    },

    # --- Engine Bay ---
    "engine_bay_edges": {
        "h_low": 20.0,
        "h_high": 60.0,
        "regime": "mixed",
        "velocity_ms": 2.0,        # partial aero penetration
        "t_air_C_low": 40.0,
        "t_air_C_high": 50.0,
        "orientation": "vertical",
        "is_internal": False,
        "notes": (
            "Engine bay perimeter, partial external airflow.  Air near "
            "ambient (40-50 C)."
        ),
    },
    "engine_beside": {
        "h_low": 5.0,
        "h_high": 15.0,
        "regime": "natural",
        "velocity_ms": 0.0,
        "t_air_C_low": 60.0,
        "t_air_C_high": 70.0,
        "orientation": "vertical",
        "is_internal": False,
        "notes": "Beside engine block.  Dead zone, buoyancy only.",
    },
    "engine_below": {
        "h_low": 3.0,
        "h_high": 10.0,
        "regime": "natural",
        "velocity_ms": 0.0,
        "t_air_C_low": 50.0,
        "t_air_C_high": 60.0,
        "orientation": "horizontal_down",   # hot surface facing down
        "is_internal": False,
        "notes": (
            "Below engine.  Hot-side-down natural convection — weakest "
            "mode.  h is low."
        ),
    },
    "engine_above": {
        "h_low": 8.0,
        "h_high": 20.0,
        "regime": "natural",
        "velocity_ms": 0.0,
        "t_air_C_low": 80.0,
        "t_air_C_high": 100.0,
        "orientation": "horizontal_up",     # hot surface facing up → plume
        "is_internal": False,
        "notes": (
            "Above engine.  Hot-side-up natural convection with thermal "
            "plume.  Air temps 80-100 C."
        ),
    },
    "exhaust_beside": {
        "h_low": 8.0,
        "h_high": 20.0,
        "regime": "natural",
        "velocity_ms": 0.0,
        "t_air_C_low": 90.0,
        "t_air_C_high": 150.0,     # with shielding
        "orientation": "vertical",
        "is_internal": False,
        "notes": (
            "Beside exhaust (shielded).  Air temps 90-150 C.  Without "
            "shielding: 125-200 C."
        ),
    },
    "exhaust_below": {
        "h_low": 3.0,
        "h_high": 12.0,
        "regime": "natural",
        "velocity_ms": 0.0,
        "t_air_C_low": 90.0,
        "t_air_C_high": 150.0,
        "orientation": "horizontal_down",
        "is_internal": False,
        "notes": (
            "Below exhaust (shielded).  Hot-side-down → weak natural "
            "convection.  Air 90-150 C."
        ),
    },
    "exhaust_above": {
        "h_low": 10.0,
        "h_high": 25.0,
        "regime": "natural",
        "velocity_ms": 0.0,
        "t_air_C_low": 150.0,
        "t_air_C_high": 350.0,     # with shielding
        "orientation": "horizontal_up",
        "is_internal": False,
        "notes": (
            "Above exhaust (shielded).  Strong buoyancy plume. "
            "Air 150-350 C.  Without shielding: 200-500 C."
        ),
    },
    "clutch_outlet_downstream": {
        "h_low": 25.0,
        "h_high": 60.0,
        "regime": "forced",
        "velocity_ms": 5.0,        # clutch fan outlet velocity
        "t_air_C_low": 90.0,
        "t_air_C_high": 110.0,     # clutch outlet temp
        "orientation": "vertical",
        "is_internal": False,
        "notes": (
            "Downstream of clutch outlet.  Forced convection from "
            "clutch fan.  Air at 90-110 C."
        ),
    },
    "exhaust_internal": {
        "h_low": 100.0,
        "h_high": 300.0,
        "regime": "forced",
        "velocity_ms": 30.0,       # typical exhaust gas velocity
        "t_air_C_low": 500.0,
        "t_air_C_high": 900.0,
        "orientation": "vertical",
        "is_internal": True,        # NOT subject to external h cap
        "notes": (
            "Internal exhaust gas flow.  High velocity turbulent flow. "
            "Not capped — internal flow."
        ),
    },

    # --- Cargo / Shield Area ---
    "shield_gap_confined": {
        "h_low": 5.0,
        "h_high": 15.0,
        "regime": "natural",
        "velocity_ms": 0.0,
        "t_air_C_low": 80.0,
        "t_air_C_high": 200.0,     # varies widely with shield temp
        "orientation": "vertical",
        "is_internal": False,
        "notes": (
            "Air gap between shield and cargo bed.  Confined natural "
            "convection.  Gap radiation typically dominates."
        ),
    },
    "cargo_bed_exterior": {
        "h_low": 10.0,
        "h_high": 40.0,
        "regime": "mixed",
        "velocity_ms": 2.0,
        "t_air_C_low": 40.0,
        "t_air_C_high": 80.0,
        "orientation": "horizontal_up",
        "is_internal": False,
        "notes": (
            "Exterior surface of cargo bed.  Some vehicle-speed "
            "flow, partially sheltered."
        ),
    },
}

# Legacy alias for backward compatibility
CONVECTION_ZONES["engine_bay_dead_zone"] = CONVECTION_ZONES["engine_beside"]
CONVECTION_ZONES["near_exhaust_natural"] = CONVECTION_ZONES["exhaust_beside"]


# ---------------------------------------------------------------------------
#  Spatial gradient defaults (K/m)
# ---------------------------------------------------------------------------

SPATIAL_GRADIENT_DEFAULTS = {
    "exhaust_max": 167.0,       # 50 C over 0.3 m — worst case
    "exhaust_typical": 83.0,    # 25 C over 0.3 m — most components
}


# ---------------------------------------------------------------------------
#  Zone lookup helpers
# ---------------------------------------------------------------------------

def get_zone(zone_name: str) -> dict:
    """
    Look up a convection zone by name.

    Parameters
    ----------
    zone_name : str
        One of the keys in CONVECTION_ZONES.

    Returns
    -------
    dict with zone parameters.

    Raises
    ------
    KeyError if zone_name not found.
    """
    if zone_name not in CONVECTION_ZONES:
        available = ", ".join(sorted(CONVECTION_ZONES.keys()))
        raise KeyError(
            f"Unknown convection zone '{zone_name}'. "
            f"Available zones: {available}"
        )
    return dict(CONVECTION_ZONES[zone_name])


def list_zones() -> list:
    """Return a list of all available zone names."""
    return sorted(CONVECTION_ZONES.keys())


def get_conservative_h(zone_name: str) -> float:
    """
    Return the conservative (upper bound) h for mesh sizing.

    This is the default for automation: using h_high produces finer
    meshes (safe side).  Using h_low would produce coarser meshes
    (risk of under-resolution).

    Parameters
    ----------
    zone_name : str

    Returns
    -------
    float — h_high for the zone (W/m^2 K)
    """
    return get_zone(zone_name)["h_high"]


def get_zone_air_temp(zone_name: str, bound: str = "high") -> float:
    """
    Return the typical air temperature for a zone in Kelvin.

    Parameters
    ----------
    zone_name : str
    bound : str — "low", "high", or "mid"

    Returns
    -------
    float — air temperature (K)
    """
    zone = get_zone(zone_name)
    if bound == "low":
        return zone["t_air_C_low"] + 273.15
    elif bound == "high":
        return zone["t_air_C_high"] + 273.15
    else:
        return (zone["t_air_C_low"] + zone["t_air_C_high"]) / 2.0 + 273.15


def estimate_spatial_gradient(
    k: float,
    q_total: float = None,
    component_class: str = "structural",
) -> float:
    """
    Estimate the spatial temperature gradient for radiation mesh sizing.

    For exhaust-adjacent surfaces, returns a conservative empirical
    bound.  For non-exhaust surfaces, derives the gradient from the
    conduction boundary flux (self-consistent, no circular dependency).

    Parameters
    ----------
    k : float
        Thermal conductivity of the part (W/m K).
    q_total : float or None
        Total boundary flux from the conduction calculator (W/m^2).
        Required for non-exhaust components.  Ignored for exhaust.
    component_class : str
        One of: "exhaust", "exhaust_adjacent", "structural", "shield".

    Returns
    -------
    float — estimated spatial gradient (K/m)

    Raises ValueError for k not > 0.
    """
    require_positive("k", k, "W/m K")
    if component_class == "exhaust":
        return SPATIAL_GRADIENT_DEFAULTS["exhaust_max"]

    if component_class == "exhaust_adjacent":
        return SPATIAL_GRADIENT_DEFAULTS["exhaust_typical"]

    # structural, shield, or anything else:
    # Gradient is conduction-driven → dT/dx ≈ q_total / k
    if q_total is None or q_total == 0:
        return 0.0

    return q_total / k
