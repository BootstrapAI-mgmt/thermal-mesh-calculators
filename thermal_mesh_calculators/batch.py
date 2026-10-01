"""
Batch Mesh Sizing Processor
============================

Accepts a list of component definitions (dicts) and runs all applicable
calculators for each component, returning a consolidated results table.

This is the automation entry point for processing hundreds to thousands
of parts from a BOM or component spreadsheet; intake.load_bom() reads one
from a CSV or JSON file.

Every key a part or project dict may carry is declared in
schema.PART_SCHEMA and schema.PROJECT_SCHEMA, with its type, range,
unit, the component classes that read it and its default, and every
result conforms to schema.RESULT_SCHEMA.  docs/integration.md walks
through them.

Minimum per-part input (6 fields):
    part_id           — unique identifier (string)
    material          — material name (key into MATERIALS)
    component_class   — one of: "exhaust", "exhaust_adjacent", "structural",
                        "shield", "multilayer_shield", "fluid" (the keys
                        of CLASS_DEFAULTS)
    convection_zone   — key into CONVECTION_ZONES (zones.py)
    thickness_mm      — wall thickness in mm
    t_surf_K          — surface temperature estimate in Kelvin
                        (ignored for shield classes — solved iteratively)

A fluid region (component_class "fluid") gives part_id, the class and its
flow: convection_zone, or both velocity_ms and bl_regime.  It is sized by
the boundary-layer constraint alone and takes no material.
process_part_from_props() sizes a part from material properties given
explicitly (k, rho, cp, ...) instead of a MATERIALS name.

Surface treatment (key into SURFACE_TREATMENTS), or an emissivity directly:
    surface / epsilon            — non-shield classes
    surface_in / eps_in          — shields, exhaust-facing side
    surface_out / eps_out        — shields, ambient-facing side
    surface_g1 / eps_g1,
    surface_g2 / eps_g2          — two-layer shields, the two gap faces
    A surface given neither falls back to a material-class emissivity and
    is reported as a SURFACE_DEFAULTED warning.

Every key is checked before anything is computed: an unknown or missing
material or component class, an unknown convection zone, an unknown
surface treatment, a value of the wrong type or outside its range, a
required key missing and a key no schema declares raise PartInputError,
whose message names the part, the key and the allowed set.  Keys starting
with "x_" are the caller's own: they are carried through and never read.

Optional per-part overrides:
    h_override        — use specific h instead of zone lookup (W/m^2 K)
    h_in_override     — shield inner h override
    h_out_override    — shield outer h override
    t_exh_K           — exhaust source temperature (shields; else the
                        project's; else 1073.15 K)
    h_gap             — two-layer shields: gap conductance (default 15)
    f12               — two-layer shields: view factor between the gap
                        faces (default 1.0, parallel plates)
    A shield input that falls back to its default (t_exh_K, h_gap, f12) is
    reported as a SHIELD_INPUT_DEFAULTED warning naming the key.
    t_fluid_K         — this part's fluid temperature, e.g. exhaust gas
                        inside a pipe (overrides the project's)
    shield_max_iter   — iteration limit of the shield solve (part or
                        project; solver default 50 single-layer, 100
                        two-layer).  A solve that stops unconverged is
                        reported as a SHIELD_NOT_CONVERGED warning.

Project-level defaults (set once, applied to all parts):
    t_fluid_K         — ambient / underhood air temperature.  A part's
                        fluid temperature (its own t_fluid_K, else this)
                        is the one temperature both h and q'' are
                        computed from.
    t_surr_K          — radiation sink temperature
    max_dt            — accuracy target (K per element)
    dt                — solver time-step (s) for transient
    fo_max            — Fourier number limit
    tau_bc            — drive-cycle segment duration (s)
    transient_scheme  — "explicit" or "implicit" (default: inferred from
                        fo_max; see TransientMeshCalculator.resolve_scheme).
                        When the smallest element dt allows exceeds the
                        governing size, the result carries a
                        TRANSIENT_CONFLICT warning with a remedy that
                        closes it.
    allowable_flux_error — radiation flux error tolerance (W/m^2)
"""

from typing import Any, Dict, List, Optional, cast

from thermal_mesh_calculators.conduction import BoundaryDrivenConductionCalculator
from thermal_mesh_calculators.convection import ConvectionMeshCalculator
from thermal_mesh_calculators.radiation import RadiationMeshCalculator
from thermal_mesh_calculators.shields import (
    SingleLayerShieldCalculator,
    MultilayerShieldCalculator,
)
from thermal_mesh_calculators.transient import (
    TransientMeshCalculator,
    _bound_text,
)
from thermal_mesh_calculators.zones import (
    get_zone,
    get_conservative_h,
    estimate_spatial_gradient,
)
from thermal_mesh_calculators.h_estimator import (
    estimate_h,
    AIR_PROPERTY_RANGE_K,
    H_EXTERNAL_MAX,
)
from thermal_mesh_calculators.boundary_layer import BoundaryLayerCalculator
from thermal_mesh_calculators.schema import (
    FLUID_CLASSES,
    RESULT_SCHEMA_VERSION,
    SHIELD_CLASSES,
    validate_material_record,
    validate_part_input,
    validate_project_input,
)


# ---------------------------------------------------------------------------
#  Input errors
# ---------------------------------------------------------------------------

class PartInputError(KeyError, ValueError):
    """
    A part dict, or its project, carries something the library refuses.

    Raised by process_part() before anything is computed: for a missing
    or unknown material or component class, an unknown convection zone, a
    zone the part needs but does not give, an unknown surface treatment,
    a value of the wrong type or out of its range, a key no schema
    declares, and an unsupported schema_version.  The message names the
    part, the key and the allowed set or range.

    ``problems`` lists every problem found, each a dict with key, code
    and reason (schema.PROBLEM_CODES); a message with more than one
    problem lists each on a line of its own.

    It subclasses KeyError, which these lookups raised before the check
    existed (``except KeyError`` keeps working), and ValueError, which is
    what it is: a bad value in an input dict.
    """

    def __init__(self, message: str = "", problems: Optional[list] = None) -> None:
        super().__init__(message)
        self.problems: List[Dict[str, Any]] = list(problems or [])

    def __str__(self) -> str:
        # KeyError.__str__ shows repr() of the message; show the message.
        return str(self.args[0]) if self.args else ""


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

SURFACE_TREATMENTS: Dict[str, Dict[str, Any]] = {
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


# Surface keys, and the direct-emissivity key that overrides each one.
_EPS_KEY_FOR_SURFACE = {
    "surface": "epsilon",
    "surface_in": "eps_in",
    "surface_out": "eps_out",
    "surface_g1": "eps_g1",
    "surface_g2": "eps_g2",
}


def _fallback_epsilon(mat_name: str) -> tuple:
    """
    Material-class emissivity for a surface given no treatment and no eps.

    Returns (epsilon, rule):
        Non-metals (plastic, rubber, composite, ceramic, glass): 0.90
        Aluminium (any alloy, any form):                         0.30
        All other metals (steel, cast iron, ...):                0.73
    """
    if any(tag in mat_name for tag in ("plastic", "rubber", "composite",
                                       "ceramic", "glass")):
        return 0.90, "non-metal"
    if "aluminium" in mat_name or "aluminum" in mat_name:
        return 0.30, "aluminium"
    return 0.73, "metal"


def _resolve_epsilon_with_source(part: dict, key: str = "surface") -> tuple:
    """
    Resolve emissivity for a part surface and say where it came from.

    Priority order:
        1. Direct emissivity in the part dict (epsilon, eps_in, eps_out,
           eps_g1, eps_g2)                       -> source "override"
        2. Surface treatment lookup (surface, surface_in, surface_out,
           surface_g1, surface_g2)               -> source "treatment"
        3. Material-class fallback (see _fallback_epsilon)
                                                 -> source "<rule> fallback"

    An unknown treatment name raises KeyError (get_surface_epsilon).

    Returns
    -------
    (epsilon, source)
    """
    direct_key = _EPS_KEY_FOR_SURFACE.get(key, "epsilon")
    if direct_key in part:
        return part[direct_key], "override"
    if key in part:
        return get_surface_epsilon(part[key]), "treatment"
    eps, rule = _fallback_epsilon(part.get("material", ""))
    return eps, rule + " fallback"


def _resolve_epsilon(part: dict, mat: dict, key: str = "surface") -> float:
    """
    Resolve emissivity for a part surface (see _resolve_epsilon_with_source).

    Parameters
    ----------
    part : dict — part definition
    mat : dict — material properties (unused; kept for compatibility)
    key : str — which surface: "surface", "surface_in", "surface_out",
                "surface_g1", "surface_g2"
    """
    return _resolve_epsilon_with_source(part, key)[0]


def _surface_default_warning(part: dict, key: str, eps: float,
                             source: str) -> dict:
    """The SURFACE_DEFAULTED warning for a surface resolved by fallback."""
    eps_key = _EPS_KEY_FOR_SURFACE[key]
    return {
        "code": "SURFACE_DEFAULTED",
        "severity": "caution",
        "key": key,
        "message": (
            f"No surface treatment ('{key}') and no emissivity "
            f"('{eps_key}') given; assumed eps = {eps:.2f}, the {source} "
            f"for material '{part.get('material', '')}'. Set '{key}' to "
            f"one of list_surface_treatments(), or give '{eps_key}'."
        ),
    }


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

    Falls back to static get_conservative_h() if the correlation fails,
    and records why in "fallback_reason" (process_part() reports it as an
    H_CORRELATION_FALLBACK warning).

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
        fallback_reason : str — the correlation's error (static_lookup only)
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
    except Exception as exc:
        # Fallback to static zone lookup if correlation fails
        h = get_conservative_h(zone_name)
        return {
            "h": h,
            "method": "static_lookup",
            "regime": "static",
            "details": None,
            "fallback_reason": f"{type(exc).__name__}: {exc}",
        }


# ---------------------------------------------------------------------------
#  Defaults per component class
# ---------------------------------------------------------------------------

CLASS_DEFAULTS: Dict[str, Dict[str, Optional[float]]] = {
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
    # A fluid region is sized by its boundary layer alone: no element
    # temperature step and no radiation tolerance apply.
    "fluid": {
        "max_dt": None,
        "allowable_flux_error": None,
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
    # The material's name is read only where there is a Biot number: a
    # fluid region has neither.
    is_plastic = bool(result.get("biot")) and any(
        tag in part.get("material", "")
        for tag in ("plastic", "rubber", "composite"))
    if is_plastic:
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
        t_fluid = result.get("t_fluid_K", project.get("t_fluid_K", 300))
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
#  Input validation
# ---------------------------------------------------------------------------
#
#  What process_part() refuses is what schema.validate_part_input() and
#  schema.validate_project_input() report: the checks are those functions.


def _prefixed(problems: list, prefix: str) -> list:
    """Problems of a project or a material record, their keys prefixed so
    they cannot be read as part keys (t_fluid_K is both)."""
    return [dict(p, key=(prefix + "." + p["key"]) if p["key"] else prefix)
            for p in problems]


def _input_error(part, problems: list) -> PartInputError:
    """PartInputError listing every problem; one problem is its own
    message."""
    pid = part.get("part_id", "?") if isinstance(part, dict) else "?"
    if len(problems) == 1:
        message = problems[0]["reason"]
    else:
        message = (f"{len(problems)} problems with part {pid!r}:"
                   + "".join("\n  - " + p["reason"] for p in problems))
    row = part.get("bom_row") if isinstance(part, dict) else None
    if isinstance(part, dict) and part.get("bom_errors") and row is not None:
        message = f"BOM row {row} did not load: {message}"
    return PartInputError(message, problems)


def _check_inputs(part, project, material_from_table: bool = True,
                  record: Optional[dict] = None) -> None:
    """
    Check a part, its project and (process_part_from_props) its material
    record before anything is computed; raise PartInputError listing
    every problem found.
    """
    problems = validate_part_input(part, material_from_table=material_from_table)
    if record is not None:
        problems += _prefixed(validate_material_record(record), "properties")
    problems += _prefixed(validate_project_input(project), "project")
    if problems:
        raise _input_error(part, problems)


def _shield_solver_options(part: dict, project: dict) -> dict:
    """
    Keyword arguments for the shield solvers: max_iter from the part's
    shield_max_iter, else the project's, else the solver's own default.
    """
    max_iter = part.get("shield_max_iter", project.get("shield_max_iter"))
    if max_iter is None:
        return {}
    if isinstance(max_iter, bool) or not isinstance(max_iter, int) or max_iter < 1:
        raise ValueError(
            f"shield_max_iter must be an integer >= 1, got {max_iter!r} "
            f"(part {part.get('part_id', '?')!r})"
        )
    return {"max_iter": max_iter}


def _shield_not_converged_warning(shield_result: dict, layers: str) -> dict:
    """SHIELD_NOT_CONVERGED: the shield solve stopped at its iteration limit."""
    return {
        "code": "SHIELD_NOT_CONVERGED",
        "severity": "warning",
        "message": (
            f"The {layers} shield solve did not converge in "
            f"{shield_result['iterations']} iterations (energy-balance "
            f"residual {shield_result['residual_W_m2']:.3g} W/m^2). The "
            f"shield temperature(s) and element size(s) are estimates from "
            f"the last iterate. Raise 'shield_max_iter' or check the inputs."
        ),
    }


def _transient_conflict_warning(trans: dict, governing_dx_mm: float,
                                fo_max: float, dt: float):
    """
    TRANSIENT_CONFLICT, or None: the time step cannot integrate the
    governing element.

    Two conflicts are checked.  Fourier minimum above the governing size:
    the minimum scales with sqrt(dt) and the governing size does not
    shrink with dt, so dt <= Fo_max * dx_gov^2 / alpha closes it.  An
    implicit window that is empty at every dt (C^2 * Fo_max < 1): only a
    larger safety_factor or fo_max closes that, and it is reported first.
    """
    remedies = []
    if trans["conflict"] is not None:
        remedies = [r for r in trans["conflict"]["remedies"]
                    if r["parameter"] != "dt"]
    fo_min = trans["fourier_min_dx_mm"]
    below_minimum = (governing_dx_mm < float("inf")
                     and fo_min > governing_dx_mm * (1.0 + 1e-9))
    if below_minimum:
        dx_m = governing_dx_mm / 1000.0
        remedies.append({
            "parameter": "dt",
            "max_value": fo_max * dx_m * dx_m / trans["alpha"],
        })
    if not remedies:
        return None

    explicit = trans["scheme"] == "explicit"
    words = []
    for r in remedies:
        # Each bound is printed rounded toward the side that satisfies it, so
        # the printed value closes the conflict; `remedies` keeps it exact.
        if r["parameter"] == "dt":
            words.append(
                f"reduce dt to <= {_bound_text(r['max_value'], 'max')} s")
        else:
            words.append(f"raise safety_factor to >= "
                         f"{_bound_text(r['min_value'], 'min')} "
                         f"(no dt can close this: both bounds scale with "
                         f"sqrt(dt))")
    if below_minimum:
        head = (
            f"Transient ({trans['scheme']}, dt = {dt:g} s): the smallest "
            f"element this time step allows, {fo_min:.4g} mm (Fo <= "
            f"{fo_max:g}, {'stability' if explicit else 'accuracy'} limit), "
            f"exceeds the governing size, {governing_dx_mm:.4g} mm."
        )
    else:
        head = trans["conflict"]["message"].split(" Remedy:")[0]
    return {
        "code": "TRANSIENT_CONFLICT",
        "severity": "warning" if explicit else "caution",
        "message": head + " Remedy: " + "; then ".join(words) + ".",
        "remedies": remedies,
    }


def _film_temperature_warning(entries: list) -> dict:
    """
    FILM_TEMP_OUT_OF_RANGE: an air-property fit was evaluated outside its
    range.  entries is a list of (what, film temperature in K).
    """
    low, high = AIR_PROPERTY_RANGE_K
    clauses = []
    for what, t_film in entries:
        t_eval = min(max(t_film, low), high)
        clauses.append(f"{what} (film temperature {t_film:.0f} K) used air "
                       f"properties evaluated at {t_eval:.0f} K")
    return {
        "code": "FILM_TEMP_OUT_OF_RANGE",
        "severity": "caution",
        "message": (
            f"The air-property fits cover film temperatures of "
            f"{low:.0f}-{high:.0f} K; " + "; ".join(clauses)
            + ". Treat the result as extrapolated."
        ),
    }


def _shield_input(part: dict, project: dict, key: str, default: float,
                  input_warnings: list, what: str,
                  project_key: bool = False) -> float:
    """
    A shield input from the part (or, when project_key, the project);
    when neither gives it, the default, reported as SHIELD_INPUT_DEFAULTED.
    """
    if key in part:
        return part[key]
    if project_key and key in project:
        return project[key]
    input_warnings.append({
        "code": "SHIELD_INPUT_DEFAULTED",
        "severity": "caution",
        "key": key,
        "message": (
            f"No '{key}' given ({what}); assumed {key} = {default:g}. "
            f"Give '{key}' on the part"
            + (" or the project." if project_key else ".")
        ),
    })
    return default


def _shield_side_h(part: dict) -> tuple:
    """(h_in, h_out) for a shield: the overrides, else each side's zone."""
    zone_in = part.get("convection_zone_in", part.get("convection_zone"))
    zone_out = part.get("convection_zone_out", part.get("convection_zone"))
    # The input check has made sure a side without an override has a zone.
    h_in = (part["h_in_override"] if "h_in_override" in part
            else get_conservative_h(cast(str, zone_in)))
    h_out = (part["h_out_override"] if "h_out_override" in part
             else get_conservative_h(cast(str, zone_out)))
    return h_in, h_out


# ---------------------------------------------------------------------------
#  Single-part processor
# ---------------------------------------------------------------------------

def process_part(part: dict, project: dict) -> dict:
    """
    Run all applicable mesh sizing calculators for a single part.

    Parameters
    ----------
    part : dict
        Per-part definition.  Every key it may carry is declared in
        schema.PART_SCHEMA (type, range, unit, the classes that read it,
        default).  Required keys:
            part_id, material, component_class, convection_zone,
            thickness_mm, t_surf_K
        A shield gives convection_zone_in / convection_zone_out (or
        convection_zone, or h_in_override / h_out_override) and no
        t_surf_K; a fluid region (component_class "fluid") gives part_id,
        component_class and convection_zone (or both velocity_ms and
        bl_regime), and no material or thickness.
        Surface keys (a surface given neither a treatment nor an emissivity
        falls back to a material-class value, reported as a
        SURFACE_DEFAULTED warning):
            surface or epsilon (non-shield classes); surface_in / eps_in and
            surface_out / eps_out (shields); surface_g1 / eps_g1 and
            surface_g2 / eps_g2 (two-layer shields, gap faces)
        Optional keys:
            h_override, h_in_override, h_out_override,
            t_exh_K, h_gap, f12, convection_zone_in, convection_zone_out,
            t_fluid_K (this part's fluid temperature, e.g. exhaust gas
            inside a pipe; overrides the project's),
            shield_max_iter (shield solve iteration limit; also a project
            key), char_length_mm, max_dt, allowable_flux_error,
            radius_mm, velocity_ms, bl_regime, bl_y_plus,
            bl_growth_ratio, bl_ar_max_prism, bl_fraction,
            schema_version (the input schema the part was written for;
            1), and keys starting with "x_" (the caller's own, never
            read)

    project : dict
        Project-level defaults (schema.PROJECT_SCHEMA).  Required keys:
            t_fluid_K, t_surr_K
        The part's fluid temperature (its own t_fluid_K, else the
        project's) is the single temperature used for the convective HTC
        (film temperature and driving difference) and for the boundary
        flux q'' alike.
        Optional keys:
            max_dt, dt, fo_max, tau_bc, safety_factor, transient_scheme,
            allowable_flux_error, t_exh_K (global exhaust temp),
            shield_max_iter, bl_regime, bl_y_plus, bl_growth_ratio,
            bl_ar_max_prism, bl_fraction, schema_version

    Returns
    -------
    dict with keys (schema.RESULT_SCHEMA; schema.validate_result()
    checks one):
        schema_version      — the result schema version (1)
        part_id             — echoed back
        material            — echoed back (None for a fluid region)
        component_class     — echoed back
        t_fluid_K           — the fluid temperature used for h and q'' (K)
        t_fluid_source      — "part" or "project": where it came from
        h_used              — the h value used (W/m^2 K; a dict per face
                              for shields; None for a fluid region)
        conduction          — dict from conduction calculator (or None)
        biot                — dict from Biot number (or None)
        radiation           — dict from radiation calculator (or None)
        shield              — dict from shield solver (or None)
        transient           — dict from transient calculator (or None)
        boundary_layer      — dict from the boundary-layer calculator,
                              when that constraint ran
        governing_dx_mm     — the smallest (most restrictive) mesh size
        governing_constraint — which calculator produced it
        all_constraints     — every candidate, {"dx_mm", "source"}
        bom_row             — the part's, when it has one
        warnings            — list of coded warnings (code, severity,
                              message; input warnings also carry "key");
                              the codes are schema.WARNING_CODES

    Raises
    ------
    PartInputError (a KeyError and a ValueError)
        Before anything is computed, for every problem schema.
        validate_part_input() and schema.validate_project_input() find:
        a missing or unknown material, component class, convection zone
        or surface treatment (the message names the part, the key and
        the allowed set), a value of the wrong type or out of range, a
        required key missing, a key no schema declares, an unsupported
        schema_version.  Its ``problems`` attribute lists them (key,
        code, reason); project keys are prefixed "project.".
    """
    _check_inputs(part, project)
    if part["component_class"] in FLUID_CLASSES:
        return _size_fluid(part, project)
    return _size_solid(part, project, MATERIALS[part["material"]],
                       part["material"])


def process_part_from_props(part: dict, project: dict, *,
                            k: float, rho: float, cp: float,
                            description: Optional[str] = None,
                            t_service_max_K: Optional[float] = None,
                            source: Optional[str] = None) -> dict:
    """
    Run all applicable mesh sizing calculators for a part whose material
    properties are given explicitly instead of named.

    For a caller that holds its own property data.  The keyword arguments
    are the fields of a material record (schema.MATERIAL_RECORD_SCHEMA),
    so an entry of MATERIALS can be passed whole:
    ``process_part_from_props(part, project, **MATERIALS["steel_mild"])``
    returns what ``process_part(dict(part, material="steel_mild"),
    project)`` returns.

    Parameters
    ----------
    part : dict
        As for process_part(), except that ``material`` is optional and is
        a label of the caller's own: it is never looked up.  It is echoed
        into the result, and two rules read its spelling: the fallback
        emissivity of a surface given neither a treatment nor an
        emissivity, and the non-metal Biot warning (see
        docs/integration.md).
    project : dict
        As for process_part().
    k : float     — thermal conductivity (W/m K)
    rho : float   — density (kg/m^3)
    cp : float    — specific heat (J/kg K)
    description : str, optional — what the material is
    t_service_max_K : float, optional
        The material's service temperature limit (K).  A part sized above
        it gets a SERVICE_TEMP_EXCEEDED warning.
    source : str, optional — where the values come from

    Returns
    -------
    dict — as process_part().  A fluid region takes no material
    properties; they are checked, and not used.

    Raises
    ------
    PartInputError
        As process_part(), and for a property out of its range (keys
        prefixed "properties.").
    """
    record: Dict[str, Any] = {"k": k, "rho": rho, "cp": cp}
    for name, value in (("description", description),
                        ("t_service_max_K", t_service_max_K),
                        ("source", source)):
        if value is not None:
            record[name] = value
    _check_inputs(part, project, material_from_table=False, record=record)
    if part["component_class"] in FLUID_CLASSES:
        return _size_fluid(part, project)
    return _size_solid(part, project, record, part.get("material"))


# Boundary-layer inputs: the part's, else the project's, else these.
_BL_DEFAULTS = {
    "bl_y_plus": 30.0,
    "bl_growth_ratio": 1.2,
    "bl_ar_max_prism": 5.0,
    "bl_fraction": 0.3,
}


def _bl_input(part: dict, project: dict, key: str) -> float:
    return part.get(key, project.get(key, _BL_DEFAULTS[key]))


def _size_fluid(part: dict, project: dict) -> dict:
    """
    Size a fluid region: the boundary-layer constraint alone.

    No conduction, lateral, Biot, radiation, shield, transient or
    curvature constraint runs: they are solid-side constraints.  The
    velocity is the part's velocity_ms, else the zone's; the regime is
    the part's or the project's bl_regime, else the zone's (a forced zone
    is "external_forced", any other "mixed_unknown"), else
    "external_forced" when the velocity is above zero.  The wall's
    t_surf_K, when given, sets the film temperature and the buoyancy
    velocity.
    """
    if "t_fluid_K" in part:
        t_fluid, t_fluid_source = part["t_fluid_K"], "part"
    else:
        t_fluid, t_fluid_source = project["t_fluid_K"], "project"
    zone = get_zone(part["convection_zone"]) if "convection_zone" in part else {}
    velocity = part.get("velocity_ms")
    if velocity is None:
        velocity = zone.get("velocity_ms", 0.0)
    regime = part.get("bl_regime", project.get("bl_regime"))
    if regime is None:
        if zone:
            regime = ("external_forced" if zone.get("regime") == "forced"
                      else "mixed_unknown")
        else:
            regime = "external_forced" if velocity > 0 else "mixed_unknown"

    bl = BoundaryLayerCalculator.estimate_mesh(
        U=velocity,
        x=part.get("char_length_mm", 100.0) / 1000.0,
        t_fluid=t_fluid,
        y_plus_target=_bl_input(part, project, "bl_y_plus"),
        growth_ratio=_bl_input(part, project, "bl_growth_ratio"),
        regime=regime,
        ar_max_prism=_bl_input(part, project, "bl_ar_max_prism"),
        bl_fraction=_bl_input(part, project, "bl_fraction"),
        t_surf=part.get("t_surf_K", 0.0),
    )
    dx = bl["max_dx_surface_mm"]
    result: Dict[str, Any] = {
        "schema_version": RESULT_SCHEMA_VERSION,
        "part_id": part["part_id"],
        "material": None,
        "component_class": part["component_class"],
        "t_fluid_K": t_fluid,
        "t_fluid_source": t_fluid_source,
        "h_used": None,
        "conduction": None,
        "biot": None,
        "radiation": None,
        "shield": None,
        "transient": None,
        "solver_advisory": None,
        "boundary_layer": bl,
        "governing_dx_mm": dx,
        "governing_constraint": "aero_boundary_layer",
        "all_constraints": [{"dx_mm": dx, "source": "aero_boundary_layer"}],
    }
    if "bom_row" in part:
        result["bom_row"] = part["bom_row"]
    warnings = []
    if not bl["film_in_range"]:
        warnings.append(_film_temperature_warning(
            [("the boundary-layer sizes", bl["t_film_K"])]))
    result["warnings"] = warnings + _generate_warnings(part, result, project)
    return result


def _size_solid(part: dict, project: dict, mat: dict,
                mat_label: Optional[str]) -> dict:
    """
    Size a solid part from its material record ``mat`` (k, rho, cp and,
    when given, t_service_max_K).  ``mat_label`` is the material name
    echoed into the result: the MATERIALS key, or the caller's own label
    (None when it gives none).  The inputs have been checked.
    """
    input_warnings = []
    solve_warnings = []
    film_out_of_range = []   # (what, film temperature) pairs

    def surface_eps(key):
        eps_value, source = _resolve_epsilon_with_source(part, key)
        if source.endswith("fallback"):
            input_warnings.append(
                _surface_default_warning(part, key, eps_value, source))
        return eps_value

    pid = part["part_id"]
    cls = part["component_class"]
    # One fluid temperature per part, for h and q'' alike: the part's own
    # t_fluid_K when given (e.g. exhaust gas inside a pipe), else the
    # project's.
    if "t_fluid_K" in part:
        t_fluid, t_fluid_source = part["t_fluid_K"], "part"
    else:
        t_fluid, t_fluid_source = project["t_fluid_K"], "project"
    t_surr = project["t_surr_K"]

    # Resolve class defaults (explicit None means "use class default"; the
    # input check has refused zero and negative values, and only the fluid
    # class, which never reaches here, has no defaults).
    class_def = CLASS_DEFAULTS.get(cls, CLASS_DEFAULTS["structural"])
    max_dt = cast(float, part.get("max_dt")
                  or project.get("max_dt")
                  or class_def["max_dt"])
    flux_err = cast(float, part.get("allowable_flux_error")
                    or project.get("allowable_flux_error")
                    or class_def["allowable_flux_error"])

    k = mat["k"]
    rho = mat["rho"]
    cp = mat["cp"]
    thickness_m = part["thickness_mm"] / 1000.0

    # Resolve general emissivity for non-shield classes.
    # Shield classes resolve eps_in/eps_out separately in their branch.
    if cls not in ("shield", "multilayer_shield"):
        eps = surface_eps("surface")
    else:
        eps = None  # will be set per-surface in shield branches

    result = {
        "schema_version": RESULT_SCHEMA_VERSION,
        "part_id": pid,
        "material": mat_label,
        "component_class": cls,
        "t_fluid_K": t_fluid,
        "t_fluid_source": t_fluid_source,
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
        t_exh = _shield_input(part, project, "t_exh_K", 1073.15,
                              input_warnings, "exhaust source temperature, K",
                              project_key=True)
        eps_in = surface_eps("surface_in")
        eps_out = surface_eps("surface_out")

        # Resolve h for each side — shields use static lookup because
        # t_surf is not yet known (it's the solve output)
        h_in, h_out = _shield_side_h(part)

        shield_result = SingleLayerShieldCalculator.mesh_size(
            k=k, max_dt=max_dt,
            t_exh=t_exh, t_fluid=t_fluid, t_surr=t_surr,
            h_in=h_in, h_out=h_out,
            eps_in=eps_in, eps_out=eps_out,
            **_shield_solver_options(part, project),
        )
        if not shield_result["converged"]:
            solve_warnings.append(
                _shield_not_converged_warning(shield_result, "single-layer"))
        result["shield"] = shield_result
        result["h_used"] = {"h_in": h_in, "h_out": h_out}
        result["eps_used"] = {"eps_in": eps_in, "eps_out": eps_out}
        dx_candidates.append((shield_result["max_dx_mm"], "shield"))

        # Use solved shield temp for Biot and radiation
        t_surf = shield_result["t_shield_K"]

    elif cls == "multilayer_shield":
        t_exh = _shield_input(part, project, "t_exh_K", 1073.15,
                              input_warnings, "exhaust source temperature, K",
                              project_key=True)
        eps_in = surface_eps("surface_in")
        eps_out = surface_eps("surface_out")
        eps_g1 = surface_eps("surface_g1")
        eps_g2 = surface_eps("surface_g2")
        h_gap = _shield_input(part, project, "h_gap", 15.0, input_warnings,
                              "gap conductance, W/m^2 K")
        f12 = _shield_input(part, project, "f12", 1.0, input_warnings,
                            "gap view factor; 1.0 is parallel plates, "
                            "offset or curved gaps are lower")

        h_in, h_out = _shield_side_h(part)

        shield_result = MultilayerShieldCalculator.mesh_sizes(
            k_metal=k, max_dt=max_dt,
            t_exh=t_exh, t_fluid=t_fluid, t_surr=t_surr,
            h_in=h_in, h_out=h_out, h_gap=h_gap,
            eps_in=eps_in, eps_out=eps_out,
            eps_g1=eps_g1, eps_g2=eps_g2, f12=f12,
            **_shield_solver_options(part, project),
        )
        if not shield_result["converged"]:
            solve_warnings.append(
                _shield_not_converged_warning(shield_result, "two-layer"))
        result["shield"] = shield_result
        result["h_used"] = {"h_in": h_in, "h_out": h_out, "h_gap": h_gap}
        result["eps_used"] = {
            "eps_in": eps_in, "eps_out": eps_out,
            "eps_g1": eps_g1, "eps_g2": eps_g2,
        }
        result["f12_used"] = f12
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
            # h at the same fluid temperature as q'' below.
            h_est = _estimate_h_from_zone(
                zone_name=part["convection_zone"],
                t_surf=t_surf,
                t_fluid=t_fluid,
                char_length_m=char_len,
            )
            h = h_est["h"]
            result["h_used"] = h
            result["h_estimation"] = h_est
            if h_est["method"] == "static_lookup":
                solve_warnings.append({
                    "code": "H_CORRELATION_FALLBACK",
                    "severity": "caution",
                    "message": (
                        f"The h correlation failed "
                        f"({h_est['fallback_reason']}); used the zone's "
                        f"static h_high = {h:.1f} W/m^2 K instead."
                    ),
                })
            elif not h_est["details"]["film_in_range"]:
                film_out_of_range.append(
                    ("the h estimate", h_est["details"]["t_film_K"]))

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
    #  Service temperature limit (when the material record gives one)
    # ------------------------------------------------------------------
    limit = mat.get("t_service_max_K")
    if limit is not None and t_surf > limit:
        if mat_label:
            whose = f"material {mat_label!r}"
        else:
            whose = repr(mat.get("description", "the given properties"))
        what = ("solved shield temperature" if cls in SHIELD_CLASSES
                else "surface temperature")
        solve_warnings.append({
            "code": "SERVICE_TEMP_EXCEEDED",
            "severity": "warning",
            "message": (
                f"The {what}, {t_surf:.1f} K, is above the service limit "
                f"of {whose}, t_service_max_K = {limit:.1f} K. Check the "
                f"material or the temperature."
            ),
        })

    # ------------------------------------------------------------------
    #  Biot number (all classes)
    # ------------------------------------------------------------------
    # The branches above set h (one surface) or h_in and h_out (shields).
    h_for_biot = (cast(float, h) if cls not in ("shield", "multilayer_shield")
                  else max(cast(float, h_in), cast(float, h_out)))
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
    # For shields, use the exhaust-facing emissivity for radiation sizing,
    # and the exhaust-facing driving flux for the conduction-driven
    # gradient dT/dx ~ q/k (zones.estimate_spatial_gradient).
    eps_for_rad = cast(float, eps if eps is not None else eps_in)
    if result["conduction"]:
        q_total = result["conduction"]["q_total"]
    elif cls == "shield":
        q_total = result["shield"]["q_boundary"]
    else:
        q_total = result["shield"]["q_layer1"]
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
    trans = None
    if dt_solver is not None:
        fo_max = project.get("fo_max", 0.5)
        safety = project.get("safety_factor", 1.0)
        tau_bc = project.get("tau_bc")

        trans = TransientMeshCalculator.combined_transient_limits(
            k=k, rho=rho, cp=cp, dt=dt_solver,
            fo_max=fo_max, safety_factor=safety, tau_bc=tau_bc,
            scheme=project.get("transient_scheme"),
        )
        result["transient"] = trans
        # The candidate is the transient upper bound, never the Fourier
        # minimum: that is a lower bound, checked against the governing
        # size below.
        if trans["max_dx_mm"] < float("inf"):
            dx_candidates.append((trans["max_dx_mm"], "transient"))

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
        bl_y_plus = _bl_input(part, project, "bl_y_plus")
        bl_growth = _bl_input(part, project, "bl_growth_ratio")
        bl_ar_max = _bl_input(part, project, "bl_ar_max_prism")
        bl_frac = _bl_input(part, project, "bl_fraction")

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
        if not bl_result["film_in_range"]:
            film_out_of_range.append(
                ("the boundary-layer sizes", bl_result["t_film_K"]))

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
    if film_out_of_range:
        solve_warnings.append(_film_temperature_warning(film_out_of_range))
    if trans is not None:
        conflict = _transient_conflict_warning(  # trans is set only when dt was given
            trans, result["governing_dx_mm"], fo_max, cast(float, dt_solver))
        if conflict is not None:
            solve_warnings.append(conflict)
    if not result["governing_dx_mm"] < float("inf"):
        solve_warnings.append({
            "code": "NO_FINITE_SIZE",
            "severity": "warning",
            "message": (
                "No constraint produced a finite element size (every "
                "candidate is unbounded, e.g. zero net boundary flux). "
                "Size this part from its geometry."
            ),
        })

    # ------------------------------------------------------------------
    #  Automatic warnings from parametric study thresholds
    # ------------------------------------------------------------------
    if "bom_row" in part:
        result["bom_row"] = part["bom_row"]
    result["warnings"] = (input_warnings + solve_warnings
                          + _generate_warnings(part, result, project))

    return result


# ---------------------------------------------------------------------------
#  Batch processor
# ---------------------------------------------------------------------------

def _error_result(part, message: str, problems: list) -> dict:
    """process_batch()'s result for a part that could not be sized."""
    pid = part.get("part_id") if isinstance(part, dict) else None
    result: Dict[str, Any] = {
        "schema_version": RESULT_SCHEMA_VERSION,
        "part_id": pid if isinstance(pid, str) and pid else "UNKNOWN",
        "error": message,
        "problems": [dict(p) for p in problems],
        "governing_dx_mm": None,
        "governing_constraint": None,
    }
    row = part.get("bom_row") if isinstance(part, dict) else None
    if isinstance(row, int) and not isinstance(row, bool):
        result["bom_row"] = row
    return result


def process_batch(parts: list, project: dict) -> list:
    """
    Process a list of part definitions and return mesh sizing results.

    Parameters
    ----------
    parts : list of dict
        Each dict is a part definition (see process_part), or a row
        returned by load_bom().
    project : dict
        Project-level defaults (see process_part).

    Returns
    -------
    list of dict — one result per part, in input order.  A part that
    cannot be sized gets an error result instead of raising (the
    ERROR_RESULT_SCHEMA of schema.py):
        part_id             — the part's, or "UNKNOWN"
        error               — the message
        problems            — list of dicts: key, code, reason; the
                              PartInputError's problems, or one
                              CALCULATION_FAILED entry when a calculator
                              raised
        governing_dx_mm, governing_constraint — None
        bom_row             — the part's, when it has one
        schema_version      — as on every result
    A row load_bom() could not load is such a part: its result lists
    every problem the row has.  A sized part's result has error None.
    """
    results = []
    for part in parts:
        try:
            r = process_part(part, project)
            r["error"] = None
        except PartInputError as e:
            r = _error_result(part, str(e), e.problems or [{
                "key": None, "code": "CALCULATION_FAILED",
                "reason": str(e)}])
        except Exception as e:
            message = str(e) or type(e).__name__
            r = _error_result(part, message, [{
                "key": None, "code": "CALCULATION_FAILED",
                "reason": f"{type(e).__name__}: {message}"}])
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

        lines.append(
            f"{r['part_id']:<30s} {r.get('material') or '-':<22s} "
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
