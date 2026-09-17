"""
Batch Mesh Sizing Processor
============================

Accepts a list of component definitions (dicts) and runs all applicable
calculators for each component, returning a consolidated results table.

This is the automation entry point for processing hundreds to thousands
of parts from a BOM or component spreadsheet.

Minimum per-part input (6 fields):
    part_id           — unique identifier (string)
    material          — material name (key into material database)
    component_class   — one of: "exhaust", "structural", "shield",
                        "multilayer_shield"
    convection_zone   — key into CONVECTION_ZONES (zones.py)
    thickness_mm      — wall thickness in mm
    t_surf_K          — surface temperature estimate in Kelvin
                        (ignored for shield classes — solved iteratively)

Optional per-part overrides:
    h_override        — use specific h instead of zone lookup (W/m^2 K)
    h_in_override     — shield inner h override
    h_out_override    — shield outer h override
    eps_in            — exhaust-facing emissivity (shields)
    eps_out           — ambient-facing emissivity (shields)
    t_exh_K           — exhaust source temperature (shields)

Project-level defaults (set once, applied to all parts):
    t_fluid_K         — ambient / underhood air temperature
    t_surr_K          — radiation sink temperature
    max_dt            — accuracy target (K per element)
    dt                — solver time-step (s) for transient
    fo_max            — Fourier number limit
    tau_bc            — drive-cycle segment duration (s)
    allowable_flux_error — radiation flux error tolerance (W/m^2)
"""

from thermal_mesh_calculators.conduction import BoundaryDrivenConductionCalculator
from thermal_mesh_calculators.convection import ConvectionMeshCalculator
from thermal_mesh_calculators.radiation import RadiationMeshCalculator
from thermal_mesh_calculators.shields import (
    SingleLayerShieldCalculator,
    MultilayerShieldCalculator,
)
from thermal_mesh_calculators.transient import TransientMeshCalculator
from thermal_mesh_calculators.zones import (
    get_zone,
    get_conservative_h,
    get_zone_air_temp,
    estimate_spatial_gradient,
)
from thermal_mesh_calculators.h_estimator import (
    estimate_h,
    H_EXTERNAL_MAX,
)
from thermal_mesh_calculators.boundary_layer import BoundaryLayerCalculator


# ---------------------------------------------------------------------------
#  Material database (inline for portability — no external files needed)
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
#  Material database — thermal transport properties only (k, rho, cp)
#
#  Emissivity is deliberately NOT stored here.  Emissivity depends on
#  surface treatment (painted, aluminised, oxidised, bare, anodised),
#  not base material.  The same mild steel can be eps=0.25 (bare),
#  eps=0.40 (aluminised), eps=0.85 (oxidised), or eps=0.92 (painted).
#
#  Surface treatment is specified per-part via the SURFACE_TREATMENTS
#  lookup below, or overridden directly with eps / eps_in / eps_out.
# ---------------------------------------------------------------------------

MATERIALS = {
    # =======================================================================
    #  STEELS
    #  Sources: Engineers Edge, Engineering ToolBox, AZoM
    # =======================================================================
    "steel_mild": {
        "k": 54.0, "rho": 7833.0, "cp": 465.0,
        "description": "Mild / carbon steel (0.1–0.5 %C)",
    },
    "steel_high_strength": {
        "k": 46.0, "rho": 7850.0, "cp": 480.0,
        "description": "High-strength low-alloy steel (HSLA, DP600)",
    },
    "steel_galvanised": {
        "k": 52.0, "rho": 7850.0, "cp": 470.0,
        "description": "Galvanised / zinc-coated carbon steel",
    },
    "steel_stainless_304": {
        "k": 16.0, "rho": 8027.0, "cp": 502.0,
        "description": "Austenitic stainless 304 / 316",
    },
    "steel_stainless_409": {
        "k": 25.0, "rho": 7720.0, "cp": 490.0,
        "description": "Ferritic stainless 409 / 430 (exhaust shields, clamps)",
    },

    # =======================================================================
    #  ALUMINIUM ALLOYS
    #  Sources: Engineers Edge, AZoM, MatWeb
    # =======================================================================
    "aluminium_6061": {
        "k": 167.0, "rho": 2710.0, "cp": 896.0,
        "description": "Wrought aluminium 6061-T6",
    },
    "aluminium_5052": {
        "k": 138.0, "rho": 2680.0, "cp": 880.0,
        "description": "Wrought aluminium 5052-H32 (sheet, brackets)",
    },
    "cast_aluminium": {
        "k": 150.0, "rho": 2700.0, "cp": 900.0,
        "description": "Cast aluminium A356-T6 / A380",
    },
    "cast_aluminium_a380": {
        "k": 96.0, "rho": 2740.0, "cp": 963.0,
        "description": "Die-cast aluminium A380 (housings, brackets)",
    },

    # =======================================================================
    #  CAST IRON
    #  Sources: Engineers Edge, Engineering ToolBox
    # =======================================================================
    "cast_iron": {
        "k": 52.0, "rho": 7200.0, "cp": 460.0,
        "description": "Grey cast iron (exhaust manifolds, blocks)",
    },
    "cast_iron_ductile": {
        "k": 36.0, "rho": 7100.0, "cp": 490.0,
        "description": "Ductile / nodular cast iron (SG iron)",
    },

    # =======================================================================
    #  OTHER METALS
    #  Sources: Engineers Edge, AZoM
    # =======================================================================
    "copper": {
        "k": 386.0, "rho": 8954.0, "cp": 380.0,
        "description": "Pure copper (bus bars, heat sinks)",
    },
    "brass": {
        "k": 119.0, "rho": 8800.0, "cp": 380.0,
        "description": "Brass 70Cu-30Zn (fittings, connectors)",
    },
    "magnesium_az91": {
        "k": 72.0, "rho": 1810.0, "cp": 1050.0,
        "description": "Magnesium AZ91D (die-cast housings, covers)",
    },
    "titanium_6al4v": {
        "k": 7.2, "rho": 4430.0, "cp": 560.0,
        "description": "Titanium Ti-6Al-4V (exhaust, fasteners)",
    },
    "zinc": {
        "k": 112.0, "rho": 7144.0, "cp": 384.0,
        "description": "Zinc alloy (die-cast, hardware)",
    },
    "nickel_alloy": {
        "k": 12.0, "rho": 8190.0, "cp": 435.0,
        "description": "Nickel alloy Inconel 625 (turbo, exhaust)",
    },

    # =======================================================================
    #  PLASTICS — GLASS-FILLED
    #  Sources: Professional Plastics, supplier TDS, MatWeb
    # =======================================================================
    "plastic_pa66_gf30": {
        "k": 0.25, "rho": 1400.0, "cp": 1600.0,
        "description": "PA66 GF30 (glass-filled nylon 66)",
    },
    "plastic_pa6_gf30": {
        "k": 0.27, "rho": 1360.0, "cp": 1600.0,
        "description": "PA6 GF30 (glass-filled nylon 6)",
    },
    "plastic_pp_gf30": {
        "k": 0.22, "rho": 1200.0, "cp": 1700.0,
        "description": "PP GF30 (glass-filled polypropylene)",
    },
    "plastic_pbt_gf30": {
        "k": 0.30, "rho": 1530.0, "cp": 1100.0,
        "description": "PBT GF30 (glass-filled, connectors, housings)",
    },
    "plastic_pps_gf40": {
        "k": 0.30, "rho": 1650.0, "cp": 1050.0,
        "description": "PPS GF40 (high-temp, exhaust-adjacent brackets)",
    },

    # =======================================================================
    #  PLASTICS — UNFILLED
    #  Sources: Professional Plastics, Engineering ToolBox, MatWeb
    # =======================================================================
    "plastic_pa66": {
        "k": 0.25, "rho": 1140.0, "cp": 1670.0,
        "description": "PA66 unfilled (nylon 66, clips, fasteners)",
    },
    "plastic_pp": {
        "k": 0.22, "rho": 905.0, "cp": 1920.0,
        "description": "PP unfilled (polypropylene, interior trim, ducts)",
    },
    "plastic_hdpe": {
        "k": 0.45, "rho": 960.0, "cp": 1900.0,
        "description": "HDPE (fuel tanks, liners, shields)",
    },
    "plastic_abs": {
        "k": 0.17, "rho": 1050.0, "cp": 1400.0,
        "description": "ABS (interior panels, bezels, housings)",
    },
    "plastic_pc": {
        "k": 0.20, "rho": 1200.0, "cp": 1250.0,
        "description": "Polycarbonate (lenses, covers, interior)",
    },
    "plastic_pc_abs": {
        "k": 0.19, "rho": 1120.0, "cp": 1350.0,
        "description": "PC/ABS blend (interior, IP components)",
    },
    "plastic_acetal_pom": {
        "k": 0.31, "rho": 1410.0, "cp": 1470.0,
        "description": "Acetal / POM (gears, clips, fuel system)",
    },
    "plastic_pet": {
        "k": 0.29, "rho": 1380.0, "cp": 1170.0,
        "description": "PET (connectors, housings, shields)",
    },

    # =======================================================================
    #  RUBBERS / ELASTOMERS
    #  Sources: Engineering ToolBox, supplier TDS
    # =======================================================================
    "rubber_epdm": {
        "k": 0.25, "rho": 1100.0, "cp": 2000.0,
        "description": "EPDM rubber (seals, hoses, weatherstrip)",
    },
    "rubber_silicone": {
        "k": 0.20, "rho": 1250.0, "cp": 1100.0,
        "description": "Silicone rubber (high-temp seals, boots)",
    },
    "rubber_nbr": {
        "k": 0.25, "rho": 1200.0, "cp": 1960.0,
        "description": "NBR / nitrile rubber (fuel hoses, O-rings)",
    },
    "rubber_natural": {
        "k": 0.15, "rho": 920.0, "cp": 1880.0,
        "description": "Natural rubber (mounts, bushings)",
    },
    "rubber_cr": {
        "k": 0.19, "rho": 1240.0, "cp": 1050.0,
        "description": "Neoprene / chloroprene (CV boots, belts)",
    },
    "rubber_fkm": {
        "k": 0.23, "rho": 1850.0, "cp": 1100.0,
        "description": "FKM / Viton (high-temp fuel seals, O-rings)",
    },

    # =======================================================================
    #  COMPOSITES / SPECIALTY
    #  Sources: MatWeb, supplier TDS
    # =======================================================================
    "composite_smc": {
        "k": 0.50, "rho": 1850.0, "cp": 1200.0,
        "description": "SMC sheet moulding compound (body panels, shields)",
    },
    "composite_cfrp": {
        "k": 1.0, "rho": 1550.0, "cp": 900.0,
        "description": "CFRP (carbon-fibre reinforced, in-plane avg)",
    },
    "glass_soda_lime": {
        "k": 1.0, "rho": 2500.0, "cp": 840.0,
        "description": "Soda-lime glass (windshields, windows)",
    },
    "ceramic_alumina": {
        "k": 25.0, "rho": 3900.0, "cp": 900.0,
        "description": "Alumina ceramic (sensor housings, insulators)",
    },
    "ceramic_cordierite": {
        "k": 2.0, "rho": 2500.0, "cp": 900.0,
        "description": "Cordierite (catalytic converter substrate)",
    },
    "insulation_fibreglass": {
        "k": 0.04, "rho": 24.0, "cp": 700.0,
        "description": "Fibreglass insulation (thermal barriers)",
    },
    "insulation_ceramic_blanket": {
        "k": 0.06, "rho": 128.0, "cp": 1050.0,
        "description": "Ceramic fibre blanket (exhaust wrap, turbo)",
    },
}


# ---------------------------------------------------------------------------
#  Surface treatment database — emissivity only
#
#  Keyed by treatment name, not material.  A part definition specifies
#  a base material (for k, rho, cp) and a surface treatment (for eps).
#  For shields, inner and outer surfaces can have different treatments.
# ---------------------------------------------------------------------------

SURFACE_TREATMENTS = {
    # --- Bare / polished metals ---
    "bare_metal": {
        "epsilon": 0.25,
        "description": "Bare / machined metal surface (generic)",
    },
    "polished_aluminium": {
        "epsilon": 0.05,
        "description": "Polished aluminium (ε 0.04–0.06)",
    },
    "polished_steel": {
        "epsilon": 0.07,
        "description": "Polished / ground steel (ε 0.05–0.10)",
    },
    "polished_copper": {
        "epsilon": 0.04,
        "description": "Polished copper (ε 0.02–0.05)",
    },
    "cast_aluminium_bare": {
        "epsilon": 0.25,
        "description": "Cast aluminium, as-cast / bare",
    },

    # --- Oxide layers ---
    "lightly_oxidised": {
        "epsilon": 0.50,
        "description": "Light oxide layer (mild service conditions)",
    },
    "heavily_oxidised": {
        "epsilon": 0.85,
        "description": "Heavy oxide (exhaust manifold, turbo housing)",
    },
    "cast_iron_oxidised": {
        "epsilon": 0.80,
        "description": "Cast iron, oxidised (manifold, housing, ε 0.60–0.80)",
    },
    "aluminium_oxidised": {
        "epsilon": 0.25,
        "description": "Aluminium heavily oxidised (ε 0.20–0.31)",
    },
    "copper_oxidised": {
        "epsilon": 0.78,
        "description": "Copper, thick oxide layer (ε 0.70–0.80)",
    },
    "stainless_weathered": {
        "epsilon": 0.85,
        "description": "Stainless steel, weathered / in-service (ε 0.80–0.90)",
    },
    "rusted_steel": {
        "epsilon": 0.70,
        "description": "Rusted carbon steel (ε 0.60–0.80)",
    },

    # --- Coatings ---
    "aluminised": {
        "epsilon": 0.40,
        "description": "Aluminised coating (heat shield, exhaust wrap)",
    },
    "painted": {
        "epsilon": 0.92,
        "description": "Painted surface (most colours ≈ 0.90–0.95)",
    },
    "painted_gloss_white": {
        "epsilon": 0.90,
        "description": "White gloss paint / powder coat",
    },
    "painted_matte_black": {
        "epsilon": 0.95,
        "description": "Matte black paint (near-blackbody)",
    },
    "anodised": {
        "epsilon": 0.80,
        "description": "Anodised aluminium (ε 0.70–0.85)",
    },
    "galvanised_new": {
        "epsilon": 0.23,
        "description": "New galvanised / zinc coating (ε 0.20–0.28)",
    },
    "galvanised_weathered": {
        "epsilon": 0.88,
        "description": "Aged / weathered galvanised surface (ε 0.80–0.90)",
    },
    "chrome_plated": {
        "epsilon": 0.06,
        "description": "Chrome plated surface (ε 0.04–0.08)",
    },
    "nickel_plated": {
        "epsilon": 0.12,
        "description": "Nickel plated surface (ε 0.05–0.15)",
    },
    "zinc_plated": {
        "epsilon": 0.20,
        "description": "Zinc plated / electrogalvanised (ε 0.15–0.25)",
    },
    "e_coat": {
        "epsilon": 0.92,
        "description": "E-coat / cathodic electrocoat (ε 0.90–0.95)",
    },
    "powder_coat": {
        "epsilon": 0.92,
        "description": "Powder coat finish (ε 0.90–0.95)",
    },
    "ceramic_coating": {
        "epsilon": 0.60,
        "description": "Ceramic thermal barrier coating (ε 0.50–0.70)",
    },

    # --- Non-metals ---
    "plastic": {
        "epsilon": 0.92,
        "description": "Typical plastic / polymer (PA, PP, PE, ABS)",
    },
    "rubber": {
        "epsilon": 0.90,
        "description": "Rubber / elastomer surface (ε 0.86–0.94)",
    },
    "glass": {
        "epsilon": 0.90,
        "description": "Glass / glazing (ε 0.85–0.94)",
    },
    "composite_smc": {
        "epsilon": 0.90,
        "description": "SMC / composite surface (ε 0.85–0.92)",
    },
    "fabric_woven": {
        "epsilon": 0.90,
        "description": "Woven fabric / textile (seats, headliner)",
    },
}


def get_material(name: str) -> dict:
    """Look up material transport properties by name (k, rho, cp)."""
    if name not in MATERIALS:
        available = ", ".join(sorted(MATERIALS.keys()))
        raise KeyError(
            f"Unknown material '{name}'. Available: {available}"
        )
    return dict(MATERIALS[name])


def list_materials() -> list:
    """Return a list of all available material names."""
    return sorted(MATERIALS.keys())


def get_surface_epsilon(treatment: str) -> float:
    """Look up emissivity for a surface treatment."""
    if treatment not in SURFACE_TREATMENTS:
        available = ", ".join(sorted(SURFACE_TREATMENTS.keys()))
        raise KeyError(
            f"Unknown surface treatment '{treatment}'. Available: {available}"
        )
    return SURFACE_TREATMENTS[treatment]["epsilon"]


def list_surface_treatments() -> list:
    """Return a list of all available surface treatment names."""
    return sorted(SURFACE_TREATMENTS.keys())


def _resolve_epsilon(part: dict, mat: dict, key: str = "surface") -> float:
    """
    Resolve emissivity for a part surface.

    Priority order:
        1. Direct epsilon override in part dict (eps, eps_in, eps_out)
        2. Surface treatment lookup (surface, surface_in, surface_out)
        3. Material-level fallback for plastics/rubber (inferred from type)
        4. KeyError if nothing resolves

    Parameters
    ----------
    part : dict — part definition
    mat : dict — material properties
    key : str — which surface: "surface", "surface_in", "surface_out",
                "surface_g1", "surface_g2"
    """
    # Direct override takes precedence
    eps_key_map = {
        "surface": "epsilon",
        "surface_in": "eps_in",
        "surface_out": "eps_out",
        "surface_g1": "eps_g1",
        "surface_g2": "eps_g2",
    }
    direct_key = eps_key_map.get(key, "epsilon")
    if direct_key in part:
        return part[direct_key]

    # Surface treatment lookup
    if key in part:
        return get_surface_epsilon(part[key])

    # Fallback: material-class rules
    #   Non-metals (plastic, rubber, composite): eps = 0.90
    #   Aluminium (any alloy, any form):         eps = 0.30
    #   All other metals (steel, cast iron):      eps = 0.73
    mat_name = part.get("material", "")
    if any(tag in mat_name for tag in ("plastic", "rubber", "composite",
                                       "ceramic", "glass")):
        return 0.90
    if "aluminium" in mat_name or "aluminum" in mat_name:
        return 0.30
    # All other metals
    return 0.73


# ---------------------------------------------------------------------------
#  Zone-based h estimation helper
# ---------------------------------------------------------------------------

def _estimate_h_from_zone(zone_name: str, t_surf: float, t_fluid: float,
                          char_length_m: float = 0.1,
                          is_internal: bool = False) -> dict:
    """
    Estimate h from zone physics using proper convection correlations.

    For zones with velocity > 0, uses forced/mixed convection correlations.
    For dead zones (velocity == 0), uses natural convection correlations
    with the zone's orientation.

    Falls back to static get_conservative_h() if the correlation fails.

    Parameters
    ----------
    zone_name : str       — convection zone name
    t_surf : float        — surface temperature (K)
    t_fluid : float       — fluid/air temperature (K)
    char_length_m : float — characteristic length (m), default 0.1
    is_internal : bool    — if True, do not apply external h cap

    Returns
    -------
    dict with keys:
        h          : float — estimated HTC (W/m^2 K)
        method     : str   — "correlation" or "static_lookup"
        regime     : str   — "forced", "natural", "mixed", or "static"
        details    : dict  — full estimate_h output (if correlation used)
    """
    zone = get_zone(zone_name)
    cap = 999.0 if (is_internal or zone.get("is_internal", False)) else H_EXTERNAL_MAX

    try:
        result = estimate_h(
            velocity=zone["velocity_ms"],
            t_surf=t_surf,
            t_fluid=t_fluid,
            char_length=char_length_m,
            orientation=zone.get("orientation", "vertical"),
            cap=cap,
        )
        return {
            "h": result["h"],
            "method": "correlation",
            "regime": result["regime"],
            "details": result,
        }
    except Exception:
        # Fallback to static zone lookup if correlation fails
        h = get_conservative_h(zone_name)
        return {
            "h": h,
            "method": "static_lookup",
            "regime": "static",
            "details": None,
        }


# ---------------------------------------------------------------------------
#  Defaults per component class
# ---------------------------------------------------------------------------

CLASS_DEFAULTS = {
    "exhaust": {
        "max_dt": 10.0,
        "allowable_flux_error": 500.0,
    },
    "exhaust_adjacent": {
        "max_dt": 10.0,
        "allowable_flux_error": 500.0,
    },
    "structural": {
        "max_dt": 15.0,
        "allowable_flux_error": 200.0,
    },
    "shield": {
        "max_dt": 15.0,
        "allowable_flux_error": 500.0,
    },
    "multilayer_shield": {
        "max_dt": 15.0,
        "allowable_flux_error": 500.0,
    },
}


# ---------------------------------------------------------------------------
#  Automatic warnings from parametric regime studies
# ---------------------------------------------------------------------------

# Thresholds discovered in analysis/offroad_regime_study.py at 10 mph baseline.
# These flag conditions where the default solver strategy may be inappropriate.

_REGIME_THRESHOLDS = {
    # Horizontal hot-up surfaces transition to turbulent natural convection
    # at surprisingly small characteristic lengths (~130-190 mm).
    # Nearly ALL horizontal dead-zone parts exceed this.
    "horizontal_up_turb_L_mm": 150.0,  # conservative midpoint of 130-190 range

    # Vertical surfaces transition to turbulent at much larger lengths.
    "vertical_turb_L_mm": 600.0,

    # Horizontal hot-down is effectively always laminar (no threshold).

    # Biot number threshold for 3D solid mesh requirement.
    # For plastics, Bi can cross 0.1 at h ~ 25 W/m²K; metals almost never cross.
    "plastic_biot_h_threshold": 25.0,

    # 10 mph forced flow suppresses turbulent buoyancy for all tested scenarios.
    # Critical velocity for suppression is typically 3-9 mph.
    "buoyancy_suppression_velocity_mph": 3.0,
}


def _generate_warnings(part: dict, result: dict, project: dict) -> list:
    """
    Generate automatic warnings based on parametric study thresholds.

    Called after process_part() computes results. Examines the part
    geometry, zone, material, and computed results to flag conditions
    that warrant engineer attention.

    Returns
    -------
    list of dict, each with keys:
        code     : str — machine-readable warning code
        severity : str — "info", "caution", "warning"
        message  : str — human-readable explanation
    """
    warnings = []
    mat_name = part.get("material", "")
    cls = part.get("component_class", "")
    zone_name = part.get("convection_zone", "")
    char_len_mm = part.get("char_length_mm", 100.0)

    # --- 1. Horizontal hot-up dead zone turbulence warning ---
    # If the part is in a zone with horizontal_up orientation and no forced flow,
    # check if characteristic length exceeds the turbulence threshold.
    if cls not in ("shield", "multilayer_shield"):
        try:
            zone = get_zone(zone_name)
        except KeyError:
            zone = {}

        orientation = zone.get("orientation", "vertical")
        velocity = zone.get("velocity_ms", 0.0)

        if orientation == "horizontal_up" and velocity < 0.5:
            threshold = _REGIME_THRESHOLDS["horizontal_up_turb_L_mm"]
            if char_len_mm > threshold:
                warnings.append({
                    "code": "TURB_NAT_HORIZ_UP",
                    "severity": "warning",
                    "message": (
                        f"Horizontal hot-up dead zone with L={char_len_mm:.0f} mm "
                        f"> {threshold:.0f} mm threshold. Natural convection is "
                        f"likely turbulent (Ra > 1e7). Recommend transient solve "
                        f"or verify with solver_advisory."
                    ),
                })

        # --- 2. Vertical dead zone with large characteristic length ---
        if orientation == "vertical" and velocity < 0.5:
            threshold = _REGIME_THRESHOLDS["vertical_turb_L_mm"]
            if char_len_mm > threshold:
                warnings.append({
                    "code": "TURB_NAT_VERTICAL",
                    "severity": "caution",
                    "message": (
                        f"Vertical dead zone with L={char_len_mm:.0f} mm "
                        f"> {threshold:.0f} mm threshold. Natural convection "
                        f"may be turbulent (Ra > 1e9). Check solver_advisory."
                    ),
                })

    # --- 3. Plastic Biot number warning ---
    is_plastic = any(tag in mat_name for tag in ("plastic", "rubber", "composite"))
    if is_plastic and result.get("biot"):
        biot_val = result["biot"]["biot"]
        h_used = result.get("h_used", 0)
        if isinstance(h_used, dict):
            h_used = max(h_used.values())
        if biot_val > 0.08:  # approaching Bi=0.1 threshold
            warnings.append({
                "code": "PLASTIC_BIOT_MARGINAL",
                "severity": "caution" if biot_val < 0.1 else "warning",
                "message": (
                    f"Non-metal part (Bi={biot_val:.3f}) approaching or exceeding "
                    f"Bi=0.1 threshold. Through-thickness gradients may require "
                    f"3D solid mesh instead of shell. h={h_used:.1f} W/m²K."
                ),
            })

    # --- 4. Solver advisory propagation ---
    sa = result.get("solver_advisory")
    if sa and sa.get("severity") == "warning":
        warnings.append({
            "code": "SOLVER_TRANSIENT_RECOMMENDED",
            "severity": "warning",
            "message": sa.get("reason", "Transient solve recommended — see solver_advisory."),
        })

    # --- 5. Low forced velocity with high surface temp (buoyancy may matter) ---
    if cls not in ("shield", "multilayer_shield"):
        t_surf = part.get("t_surf_K", 300)
        t_fluid = project.get("t_fluid_K", 300)
        dt_surf = t_surf - t_fluid
        try:
            zone = get_zone(zone_name)
        except KeyError:
            zone = {}
        velocity = zone.get("velocity_ms", 0.0)
        if 0.5 < velocity < 3.0 and dt_surf > 100:
            warnings.append({
                "code": "LOW_VELOCITY_HIGH_DT",
                "severity": "caution",
                "message": (
                    f"Low forced velocity ({velocity:.1f} m/s) with high surface "
                    f"temperature rise (dT={dt_surf:.0f} K). Mixed convection "
                    f"likely — buoyancy effects may be significant. Check Ri."
                ),
            })

    return warnings


# ---------------------------------------------------------------------------
#  Single-part processor
# ---------------------------------------------------------------------------

def process_part(part: dict, project: dict) -> dict:
    """
    Run all applicable mesh sizing calculators for a single part.

    Parameters
    ----------
    part : dict
        Per-part definition.  Required keys:
            part_id, material, component_class, convection_zone,
            thickness_mm, t_surf_K, surface (surface treatment name)
        Optional keys:
            h_override, h_in_override, h_out_override,
            epsilon (direct override), eps_in, eps_out, eps_g1, eps_g2,
            surface_in, surface_out, surface_g1, surface_g2,
            t_exh_K, h_gap, convection_zone_in, convection_zone_out

    project : dict
        Project-level defaults.  Expected keys:
            t_fluid_K, t_surr_K
        Optional keys:
            max_dt, dt, fo_max, tau_bc, safety_factor,
            allowable_flux_error, t_exh_K (global exhaust temp)

    Returns
    -------
    dict with keys:
        part_id             — echoed back
        material            — echoed back
        component_class     — echoed back
        h_used              — the h value used (W/m^2 K)
        conduction          — dict from conduction calculator (or None)
        biot                — dict from Biot number (or None)
        radiation           — dict from radiation calculator (or None)
        shield              — dict from shield solver (or None)
        transient           — dict from transient calculator (or None)
        governing_dx_mm     — the smallest (most restrictive) mesh size
        governing_constraint — which calculator produced it
    """
    pid = part["part_id"]
    mat = get_material(part["material"])
    cls = part["component_class"]
    t_fluid = project["t_fluid_K"]
    t_surr = project["t_surr_K"]

    # Resolve class defaults (explicit None means "use class default")
    class_def = CLASS_DEFAULTS.get(cls, CLASS_DEFAULTS["structural"])
    max_dt = (part.get("max_dt")
              or project.get("max_dt")
              or class_def["max_dt"])
    flux_err = (part.get("allowable_flux_error")
                or project.get("allowable_flux_error")
                or class_def["allowable_flux_error"])

    k = mat["k"]
    rho = mat["rho"]
    cp = mat["cp"]
    thickness_m = part["thickness_mm"] / 1000.0

    # Resolve general emissivity for non-shield classes.
    # Shield classes resolve eps_in/eps_out separately in their branch.
    if cls not in ("shield", "multilayer_shield"):
        eps = _resolve_epsilon(part, mat, "surface")
    else:
        eps = None  # will be set per-surface in shield branches

    result = {
        "part_id": pid,
        "material": part["material"],
        "component_class": cls,
        "conduction": None,
        "biot": None,
        "radiation": None,
        "shield": None,
        "transient": None,
        "solver_advisory": None,
        "governing_dx_mm": float("inf"),
        "governing_constraint": None,
    }

    dx_candidates = []  # (dx_mm, label) pairs
    h = None        # set in non-shield branch
    h_in = None     # set in shield branches
    h_out = None
    eps_in = None   # set in shield branches

    # ------------------------------------------------------------------
    #  SHIELD classes — solve temperature first, then mesh size
    # ------------------------------------------------------------------
    if cls == "shield":
        t_exh = part.get("t_exh_K", project.get("t_exh_K", 1073.15))
        eps_in = _resolve_epsilon(part, mat, "surface_in")
        eps_out = _resolve_epsilon(part, mat, "surface_out")

        # Resolve h for each side — shields use static lookup because
        # t_surf is not yet known (it's the solve output)
        zone_in = part.get("convection_zone_in", part.get("convection_zone"))
        zone_out = part.get("convection_zone_out", part.get("convection_zone"))
        h_in = part.get("h_in_override", get_conservative_h(zone_in))
        h_out = part.get("h_out_override", get_conservative_h(zone_out))

        shield_result = SingleLayerShieldCalculator.mesh_size(
            k=k, max_dt=max_dt,
            t_exh=t_exh, t_fluid=t_fluid, t_surr=t_surr,
            h_in=h_in, h_out=h_out,
            eps_in=eps_in, eps_out=eps_out,
        )
        result["shield"] = shield_result
        result["h_used"] = {"h_in": h_in, "h_out": h_out}
        result["eps_used"] = {"eps_in": eps_in, "eps_out": eps_out}
        dx_candidates.append((shield_result["max_dx_mm"], "shield"))

        # Use solved shield temp for Biot and radiation
        t_surf = shield_result["t_shield_K"]

    elif cls == "multilayer_shield":
        t_exh = part.get("t_exh_K", project.get("t_exh_K", 1073.15))
        eps_in = _resolve_epsilon(part, mat, "surface_in")
        eps_out = _resolve_epsilon(part, mat, "surface_out")
        eps_g1 = _resolve_epsilon(part, mat, "surface_g1")
        eps_g2 = _resolve_epsilon(part, mat, "surface_g2")
        h_gap = part.get("h_gap", 15.0)

        zone_in = part.get("convection_zone_in", part.get("convection_zone"))
        zone_out = part.get("convection_zone_out", part.get("convection_zone"))
        h_in = part.get("h_in_override", get_conservative_h(zone_in))
        h_out = part.get("h_out_override", get_conservative_h(zone_out))

        shield_result = MultilayerShieldCalculator.mesh_sizes(
            k_metal=k, max_dt=max_dt,
            t_exh=t_exh, t_fluid=t_fluid, t_surr=t_surr,
            h_in=h_in, h_out=h_out, h_gap=h_gap,
            eps_in=eps_in, eps_out=eps_out,
            eps_g1=eps_g1, eps_g2=eps_g2,
        )
        result["shield"] = shield_result
        result["h_used"] = {"h_in": h_in, "h_out": h_out, "h_gap": h_gap}
        result["eps_used"] = {
            "eps_in": eps_in, "eps_out": eps_out,
            "eps_g1": eps_g1, "eps_g2": eps_g2,
        }
        dx_candidates.append((shield_result["layer1_max_dx_mm"], "shield_layer1"))
        dx_candidates.append((shield_result["layer2_max_dx_mm"], "shield_layer2"))

        t_surf = shield_result["t1_K"]  # use hottest layer for Biot/radiation

    else:
        # ------------------------------------------------------------------
        #  NON-SHIELD classes — use estimated or specified t_surf
        # ------------------------------------------------------------------
        t_surf = part["t_surf_K"]

        # Resolve h — correlation-based when possible, static fallback
        if "h_override" in part:
            h = part["h_override"]
            result["h_used"] = h
            result["h_estimation"] = {
                "method": "override", "regime": "user_specified",
            }
        else:
            char_len = part.get("char_length_mm", 100.0) / 1000.0
            zone_t_fluid = get_zone_air_temp(
                part["convection_zone"], "high",
            )
            # Use zone air temp if project t_fluid not explicitly set
            t_fluid_for_h = t_fluid if t_fluid != t_surr else zone_t_fluid
            h_est = _estimate_h_from_zone(
                zone_name=part["convection_zone"],
                t_surf=t_surf,
                t_fluid=zone_t_fluid,
                char_length_m=char_len,
            )
            h = h_est["h"]
            result["h_used"] = h
            result["h_estimation"] = h_est

            # Propagate solver advisory from h estimation
            details = h_est.get("details")
            if details and "solver_advisory" in details:
                result["solver_advisory"] = details["solver_advisory"]

        result["eps_used"] = eps

        # Conduction
        cond = BoundaryDrivenConductionCalculator.max_mesh_size(
            k=k, h=h, t_surf=t_surf, t_fluid=t_fluid,
            epsilon=eps, t_surr=t_surr, max_dt=max_dt,
        )
        result["conduction"] = cond
        dx_candidates.append((cond["max_dx_mm"], "conduction"))

        # Lateral gradient constraint (fin theory)
        lateral = BoundaryDrivenConductionCalculator.lateral_gradient_limit(
            k=k, h=h, thickness_m=thickness_m,
            epsilon=eps, t_surf=t_surf,
        )
        result["lateral"] = lateral
        if lateral["max_dx_mm"] < float("inf"):
            dx_candidates.append((lateral["max_dx_mm"], "lateral_gradient"))

    # ------------------------------------------------------------------
    #  Biot number (all classes)
    # ------------------------------------------------------------------
    h_for_biot = h if cls not in ("shield", "multilayer_shield") else max(h_in, h_out)
    biot = ConvectionMeshCalculator.biot_number(
        h=h_for_biot, k=k, thickness=thickness_m,
    )
    result["biot"] = biot
    # Biot doesn't produce a dx directly, but if 3D solid is required,
    # the through-thickness element size is a constraint
    if biot["mesh_type"] == "3D Solid":
        through_dx = (part["thickness_mm"]
                      / biot["min_elements_through_thickness"])
        dx_candidates.append((through_dx, "biot_through_thickness"))

    # ------------------------------------------------------------------
    #  Radiation mesh size (all classes)
    # ------------------------------------------------------------------
    # For shields, use the exhaust-facing emissivity for radiation sizing
    eps_for_rad = eps if eps is not None else eps_in
    q_total = result["conduction"]["q_total"] if result["conduction"] else None
    grad = estimate_spatial_gradient(k, q_total=q_total, component_class=cls)
    if grad > 0:
        rad = RadiationMeshCalculator.max_mesh_size(
            t_local=t_surf, emissivity=eps_for_rad,
            allowable_flux_error=flux_err,
            spatial_gradient=grad,
        )
        result["radiation"] = rad
        dx_candidates.append((rad["max_dx_mm"], "radiation"))

    # ------------------------------------------------------------------
    #  Transient (if solver time-step provided)
    # ------------------------------------------------------------------
    dt_solver = project.get("dt")
    if dt_solver is not None:
        fo_max = project.get("fo_max", 0.5)
        safety = project.get("safety_factor", 1.0)
        tau_bc = project.get("tau_bc")

        trans = TransientMeshCalculator.combined_transient_limits(
            k=k, rho=rho, cp=cp, dt=dt_solver,
            fo_max=fo_max, safety_factor=safety, tau_bc=tau_bc,
        )
        result["transient"] = trans
        dx_candidates.append((trans["recommended_dx_mm"], "transient"))

    # ------------------------------------------------------------------
    #  View factor curvature limit (if radius provided)
    # ------------------------------------------------------------------
    radius = part.get("radius_mm")
    if radius is not None:
        curv = RadiationMeshCalculator.view_factor_curvature_limit(radius)
        result["curvature_dx_mm"] = curv
        dx_candidates.append((curv, "curvature"))

    # ------------------------------------------------------------------
    #  Aerodynamic boundary layer constraint (if velocity or BL regime
    #  is specified at part or project level)
    # ------------------------------------------------------------------
    # Inputs: velocity_ms (from zone or override), bl_regime
    #   bl_regime = "external_forced" | "mixed_unknown" | None (skip)
    #   bl_y_plus = target y+ (default 30)
    #   bl_growth_ratio = inflation growth ratio (default 1.2)
    #   bl_ar_max_prism = max prism AR for mixed_unknown (default 5)
    #   bl_fraction = C in dx ≤ C·δ (default 0.3)
    bl_regime = part.get("bl_regime", project.get("bl_regime"))
    if bl_regime is not None:
        # Resolve velocity: explicit override > zone lookup > 0
        bl_velocity = part.get("velocity_ms")
        if bl_velocity is None:
            try:
                zone = get_zone(part.get("convection_zone", ""))
                bl_velocity = zone.get("velocity_ms", 0.0)
            except KeyError:
                bl_velocity = 0.0

        bl_char_len = part.get("char_length_mm", 100.0) / 1000.0
        bl_y_plus = part.get("bl_y_plus", project.get("bl_y_plus", 30.0))
        bl_growth = part.get(
            "bl_growth_ratio", project.get("bl_growth_ratio", 1.2)
        )
        bl_ar_max = part.get(
            "bl_ar_max_prism", project.get("bl_ar_max_prism", 5.0)
        )
        bl_frac = part.get(
            "bl_fraction", project.get("bl_fraction", 0.3)
        )

        # For mixed_unknown, compute buoyancy dT from surface & fluid temps
        bl_dt_buoy = 0.0
        if bl_regime == "mixed_unknown":
            bl_dt_buoy = abs(t_surf - t_fluid)

        bl_result = BoundaryLayerCalculator.estimate_mesh(
            U=bl_velocity,
            x=bl_char_len,
            t_fluid=t_fluid,
            y_plus_target=bl_y_plus,
            growth_ratio=bl_growth,
            regime=bl_regime,
            delta_t_buoyancy=bl_dt_buoy,
            ar_max_prism=bl_ar_max,
            bl_fraction=bl_frac,
            t_surf=t_surf,
        )
        result["boundary_layer"] = bl_result
        dx_candidates.append(
            (bl_result["max_dx_surface_mm"], "aero_boundary_layer")
        )

    # ------------------------------------------------------------------
    #  Governing constraint
    # ------------------------------------------------------------------
    if dx_candidates:
        dx_candidates.sort(key=lambda x: x[0])
        result["governing_dx_mm"] = dx_candidates[0][0]
        result["governing_constraint"] = dx_candidates[0][1]
        result["all_constraints"] = [
            {"dx_mm": dx, "source": src} for dx, src in dx_candidates
        ]

    # ------------------------------------------------------------------
    #  Automatic warnings from parametric study thresholds
    # ------------------------------------------------------------------
    result["warnings"] = _generate_warnings(part, result, project)

    return result


# ---------------------------------------------------------------------------
#  Batch processor
# ---------------------------------------------------------------------------

def process_batch(parts: list, project: dict) -> list:
    """
    Process a list of part definitions and return mesh sizing results.

    Parameters
    ----------
    parts : list of dict
        Each dict is a part definition (see process_part).
    project : dict
        Project-level defaults (see process_part).

    Returns
    -------
    list of dict — one result per part, in input order.
    """
    results = []
    for part in parts:
        try:
            r = process_part(part, project)
            r["error"] = None
        except Exception as e:
            r = {
                "part_id": part.get("part_id", "UNKNOWN"),
                "error": str(e),
                "governing_dx_mm": None,
                "governing_constraint": None,
            }
        results.append(r)
    return results


def summary_table(results: list) -> str:
    """
    Format batch results as a human-readable text table.

    Parameters
    ----------
    results : list of dict
        Output from process_batch().

    Returns
    -------
    str — formatted table
    """
    header = (
        f"{'Part ID':<30s} {'Material':<22s} {'Class':<18s} "
        f"{'dx (mm)':>8s}  {'Governing':>22s}  {'Biot':>8s}"
    )
    sep = "-" * len(header)
    lines = [header, sep]

    for r in results:
        if r.get("error"):
            lines.append(
                f"{r['part_id']:<30s} {'ERROR':<22s} {'':<18s} "
                f"{'---':>8s}  {r['error'][:22]:>22s}  {'---':>8s}"
            )
            continue

        dx = r["governing_dx_mm"]
        dx_str = f"{dx:.2f}" if dx < float("inf") else "inf"
        bi_val = r["biot"]["biot"] if r.get("biot") else 0.0
        bi_str = f"{bi_val:.4f}"
        mesh_type = r["biot"]["mesh_type"] if r.get("biot") else "—"

        lines.append(
            f"{r['part_id']:<30s} {r['material']:<22s} "
            f"{r['component_class']:<18s} "
            f"{dx_str:>8s}  {r['governing_constraint']:>22s}  "
            f"{bi_str:>8s}"
        )

    # Append warnings summary
    warn_parts = []
    for r in results:
        for w in r.get("warnings", []):
            if w["severity"] in ("warning", "caution"):
                warn_parts.append((r.get("part_id", "?"), w))

    if warn_parts:
        lines.append("")
        lines.append("=" * len(header))
        lines.append("  WARNINGS")
        lines.append("=" * len(header))
        for pid, w in warn_parts:
            icon = "!!" if w["severity"] == "warning" else " ?"
            lines.append(f"  [{icon}] {pid}: {w['message']}")

    return "\n".join(lines)
