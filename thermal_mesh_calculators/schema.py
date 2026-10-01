"""
Input and Result Schemas
========================

process_part() reads a part dict and a project dict and returns a result
dict.  This module publishes all of them as data, each with a version and
a validator that needs nothing but the standard library:

    PART_SCHEMA, PROJECT_SCHEMA  — every key a part or project dict may
                                   carry: type, range, unit, the component
                                   classes that read it, when it is
                                   required, and its default
    MATERIAL_RECORD_SCHEMA       — the fields of a material record: an
                                   entry of MATERIALS, or the properties
                                   given to process_part_from_props()
    RESULT_SCHEMA,               — the dicts process_part() and
    ERROR_RESULT_SCHEMA            process_batch() return, and the
    WARNING_CODES                  warning codes a result may carry

    validate_part_input(), validate_project_input(),
    validate_material_record(), validate_result()
                                 — each returns a list of problems, empty
                                   when the dict conforms

The input checks are the ones process_part() runs before it computes
anything.  When they find a problem, process_part() raises PartInputError
carrying the same problem dicts, worded the same way, so a part that
validates is a part process_part() accepts.  A range check is one of the
calculators' own input guards, called with the key's name.

A key is checked only for the component classes that read it (a shield
does not read t_surf_K, so a shield's t_surf_K is not checked), except
the convection-zone and surface-treatment names, which are checked
whatever the class.  A key no schema declares is refused, unless it
starts with "x_": those keys are the caller's own and are never read.

Versions
--------
INPUT_SCHEMA_VERSION and RESULT_SCHEMA_VERSION are integers.  A part or
project dict may name the input schema it was written for in a
``schema_version`` key (absent means 1); every result carries the result
schema version it conforms to.  docs/integration.md says how keys are
added, renamed and removed across releases.

Problems
--------
A problem is a dict with three keys:
    key     — the key it concerns (None for the dict as a whole)
    code    — one of PROBLEM_CODES
    reason  — a sentence naming the key, the rule and the value given
"""

import math
from typing import Any, Dict, List, Optional

from thermal_mesh_calculators._guards import (
    _number,
    require_fraction,
    require_non_negative,
    require_positive,
)


INPUT_SCHEMA_VERSION = 1
RESULT_SCHEMA_VERSION = 1

# ---------------------------------------------------------------------------
#  Component classes
# ---------------------------------------------------------------------------

NON_SHIELD_SOLID_CLASSES = ("exhaust", "exhaust_adjacent", "structural")
SHIELD_CLASSES = ("shield", "multilayer_shield")
SOLID_CLASSES = NON_SHIELD_SOLID_CLASSES + SHIELD_CLASSES
FLUID_CLASSES = ("fluid",)
COMPONENT_CLASSES = SOLID_CLASSES + FLUID_CLASSES

_TWO_LAYER = ("multilayer_shield",)

BL_REGIMES = ("external_forced", "mixed_unknown")
TRANSIENT_SCHEMES = ("explicit", "implicit")

# Keys the caller may add for its own use; they are never read.
EXTENSION_PREFIX = "x_"

PROBLEM_CODES = {
    "NOT_AN_OBJECT": "The part, project or record is not a dict (a JSON "
                     "element that is not an object, say).",
    "UNKNOWN_KEY": "A key the schema does not declare (and that does not "
                   "start with 'x_').",
    "RESERVED_KEY": "A key that only load_bom() sets (bom_row, bom_errors), "
                    "found in a BOM file.",
    "MISSING_KEY": "A key the part's class requires is absent.",
    "WRONG_TYPE": "A value of the wrong type: a string where a number is "
                  "expected, a bool, a float where an integer is expected.",
    "OUT_OF_RANGE": "A number outside its range, or not finite.",
    "UNKNOWN_NAME": "A name that is not in its vocabulary (a material, "
                    "zone, surface treatment, class or option).",
    "UNSUPPORTED_VERSION": "A schema_version this library does not read.",
    "DUPLICATE_PART_ID": "A part_id an earlier BOM row already uses.",
    "MALFORMED_ROW": "A CSV row with more values than the header has "
                     "columns.",
    "INCONSISTENT": "Result values that contradict each other (the "
                    "governing size is not the smallest candidate, say).",
    "CALCULATION_FAILED": "A calculator raised on inputs the schema "
                          "accepts; process_batch() reports it per part.",
}


# ---------------------------------------------------------------------------
#  Schema entries
# ---------------------------------------------------------------------------
#
#  Each schema maps a key to a dict:
#      type        "number", "integer", "string" or "list"
#      range       for numbers: (">", lo), (">=", lo), ("[]", lo, hi) or
#                  ("(]", lo, hi); a number must also be finite
#      unit        the unit of a number
#      nullable    True when None is accepted and means "use the default"
#      vocabulary  for names: the package table whose keys are the names
#                  ("MATERIALS", "SURFACE_TREATMENTS", "CONVECTION_ZONES",
#                  "CLASS_DEFAULTS")
#      values      for options: the allowed strings
#      non_empty   for strings: the empty string is refused
#      what        how a message names the key ("wall thickness")
#      read_by     (part keys) the component classes that read the key
#      checked_for (part keys) "all" when the key is checked whatever the
#                  class; otherwise only for the classes in read_by
#      required    a tuple of rules (classes, alternatives): the key is
#                  required for a part of one of those classes unless
#                  every key of one alternative is given (present, and
#                  not None); the project schema says True or False
#      default     what applies when the key is absent
#      description one line for the documentation


def _required(classes, *unless):
    return {"classes": tuple(classes), "unless": tuple(tuple(a) for a in unless)}


def _num(what, rng, unit="", **extra):
    entry = {"type": "number", "range": rng, "what": what}
    if unit:
        entry["unit"] = unit
    entry.update(extra)
    return entry


_POS = (">", 0.0)
_NON_NEG = (">=", 0.0)
_FRACTION = ("[]", 0.0, 1.0)

PART_SCHEMA: Dict[str, Dict[str, Any]] = {
    "schema_version": {
        "type": "integer", "what": "input schema version",
        "read_by": COMPONENT_CLASSES, "checked_for": "all",
        "default": "1",
        "description": "The input schema the part was written for.",
    },
    "part_id": {
        "type": "string", "non_empty": True, "what": "part id",
        "read_by": COMPONENT_CLASSES,
        "required": (_required(COMPONENT_CLASSES),),
        "description": "Identifier, echoed into the result.",
    },
    "component_class": {
        "type": "string", "vocabulary": "CLASS_DEFAULTS",
        "what": "component class", "read_by": COMPONENT_CLASSES,
        "required": (_required(COMPONENT_CLASSES),),
        "description": "Which sizing path runs (docs/integration.md).",
    },
    "material": {
        "type": "string", "vocabulary": "MATERIALS", "what": "material",
        "read_by": SOLID_CLASSES,
        "required": (_required(SOLID_CLASSES),),
        "description": "A key of MATERIALS (a label of the caller's own "
                       "with process_part_from_props()).",
    },
    "thickness_mm": _num(
        "wall thickness", _POS, "mm", read_by=SOLID_CLASSES,
        required=(_required(SOLID_CLASSES),),
        description="Wall thickness."),
    "t_surf_K": _num(
        "surface temperature", _POS, "K",
        read_by=NON_SHIELD_SOLID_CLASSES + FLUID_CLASSES,
        required=(_required(NON_SHIELD_SOLID_CLASSES),),
        description="Surface temperature; for a fluid region, its wall's "
                    "(optional).  Shields solve theirs and ignore it."),
    "convection_zone": {
        "type": "string", "vocabulary": "CONVECTION_ZONES",
        "what": "convection zone", "read_by": COMPONENT_CLASSES,
        "checked_for": "all",
        "required": (
            _required(NON_SHIELD_SOLID_CLASSES, ("h_override",)),
            _required(FLUID_CLASSES, ("velocity_ms", "bl_regime")),
        ),
        "description": "A key of CONVECTION_ZONES; a shield's zone for "
                       "a side that names none.",
    },
    "convection_zone_in": {
        "type": "string", "vocabulary": "CONVECTION_ZONES",
        "what": "convection zone", "read_by": SHIELD_CLASSES,
        "checked_for": "all",
        "required": (_required(SHIELD_CLASSES, ("convection_zone",),
                               ("h_in_override",)),),
        "description": "Shields: the exhaust-facing side's zone.",
    },
    "convection_zone_out": {
        "type": "string", "vocabulary": "CONVECTION_ZONES",
        "what": "convection zone", "read_by": SHIELD_CLASSES,
        "checked_for": "all",
        "required": (_required(SHIELD_CLASSES, ("convection_zone",),
                               ("h_out_override",)),),
        "description": "Shields: the ambient-facing side's zone.",
    },
    "surface": {
        "type": "string", "vocabulary": "SURFACE_TREATMENTS",
        "what": "surface treatment", "read_by": NON_SHIELD_SOLID_CLASSES,
        "checked_for": "all",
        "default": "the material-class emissivity (SURFACE_DEFAULTED)",
        "description": "A key of SURFACE_TREATMENTS.",
    },
    "epsilon": _num(
        "emissivity", _FRACTION, read_by=NON_SHIELD_SOLID_CLASSES,
        description="Emissivity; overrides 'surface'."),
    "surface_in": {
        "type": "string", "vocabulary": "SURFACE_TREATMENTS",
        "what": "surface treatment", "read_by": SHIELD_CLASSES,
        "checked_for": "all",
        "default": "the material-class emissivity (SURFACE_DEFAULTED)",
        "description": "Shields: the exhaust-facing surface.",
    },
    "eps_in": _num(
        "emissivity", _FRACTION, read_by=SHIELD_CLASSES,
        description="Shields: overrides 'surface_in'."),
    "surface_out": {
        "type": "string", "vocabulary": "SURFACE_TREATMENTS",
        "what": "surface treatment", "read_by": SHIELD_CLASSES,
        "checked_for": "all",
        "default": "the material-class emissivity (SURFACE_DEFAULTED)",
        "description": "Shields: the ambient-facing surface.",
    },
    "eps_out": _num(
        "emissivity", _FRACTION, read_by=SHIELD_CLASSES,
        description="Shields: overrides 'surface_out'."),
    "surface_g1": {
        "type": "string", "vocabulary": "SURFACE_TREATMENTS",
        "what": "surface treatment", "read_by": _TWO_LAYER,
        "checked_for": "all",
        "default": "the material-class emissivity (SURFACE_DEFAULTED)",
        "description": "Two-layer shields: layer 1's gap face.",
    },
    "eps_g1": _num(
        "emissivity", _FRACTION, read_by=_TWO_LAYER,
        description="Two-layer shields: overrides 'surface_g1'."),
    "surface_g2": {
        "type": "string", "vocabulary": "SURFACE_TREATMENTS",
        "what": "surface treatment", "read_by": _TWO_LAYER,
        "checked_for": "all",
        "default": "the material-class emissivity (SURFACE_DEFAULTED)",
        "description": "Two-layer shields: layer 2's gap face.",
    },
    "eps_g2": _num(
        "emissivity", _FRACTION, read_by=_TWO_LAYER,
        description="Two-layer shields: overrides 'surface_g2'."),
    "h_override": _num(
        "convection coefficient", _NON_NEG, "W/m^2 K",
        read_by=NON_SHIELD_SOLID_CLASSES,
        description="h to use instead of the zone estimate."),
    "h_in_override": _num(
        "convection coefficient", _NON_NEG, "W/m^2 K",
        read_by=SHIELD_CLASSES,
        description="Shields: h of the exhaust-facing side."),
    "h_out_override": _num(
        "convection coefficient", _NON_NEG, "W/m^2 K",
        read_by=SHIELD_CLASSES,
        description="Shields: h of the ambient-facing side."),
    "t_fluid_K": _num(
        "fluid temperature", _POS, "K", read_by=COMPONENT_CLASSES,
        default="the project's t_fluid_K",
        description="This part's fluid temperature (exhaust gas inside a "
                    "pipe, say)."),
    "char_length_mm": _num(
        "characteristic length", _POS, "mm", read_by=COMPONENT_CLASSES,
        default="100",
        description="Length for the h correlation and the boundary layer."),
    "t_exh_K": _num(
        "exhaust temperature", _POS, "K", read_by=SHIELD_CLASSES,
        default="the project's t_exh_K, else 1073.15 "
                "(SHIELD_INPUT_DEFAULTED)",
        description="Shields: the exhaust source temperature."),
    "h_gap": _num(
        "gap conductance", _NON_NEG, "W/m^2 K", read_by=_TWO_LAYER,
        default="15 (SHIELD_INPUT_DEFAULTED)",
        description="Two-layer shields: the gap conductance."),
    "f12": _num(
        "gap view factor", _FRACTION, read_by=_TWO_LAYER,
        default="1.0, parallel plates (SHIELD_INPUT_DEFAULTED)",
        description="Two-layer shields: the view factor across the gap."),
    "shield_max_iter": {
        "type": "integer", "range": (">=", 1), "nullable": True,
        "what": "shield iteration limit", "read_by": SHIELD_CLASSES,
        "default": "the project's, else 50 single-layer, 100 two-layer",
        "description": "Shields: the solve's iteration limit.",
    },
    "max_dt": _num(
        "temperature step", _POS, "K", nullable=True,
        read_by=SOLID_CLASSES,
        default="the project's, else the class default",
        description="Accuracy target: temperature change per element."),
    "allowable_flux_error": _num(
        "radiation flux error", _POS, "W/m^2", nullable=True,
        read_by=SOLID_CLASSES,
        default="the project's, else the class default",
        description="Radiation linearisation error tolerance."),
    "radius_mm": _num(
        "curvature radius", _POS, "mm", nullable=True,
        read_by=SOLID_CLASSES,
        default="no curvature limit",
        description="Radius of a curved radiating surface."),
    "velocity_ms": _num(
        "velocity", _NON_NEG, "m/s", nullable=True,
        read_by=COMPONENT_CLASSES,
        default="the zone's velocity",
        description="Local velocity for the boundary layer."),
    "bl_regime": {
        "type": "string", "values": BL_REGIMES, "nullable": True,
        "what": "boundary-layer regime", "read_by": COMPONENT_CLASSES,
        "default": "the project's; a solid without one gets no "
                   "boundary-layer constraint; a fluid region takes its "
                   "zone's (forced: external_forced, else mixed_unknown)",
        "description": "Turns on the boundary-layer constraint.",
    },
    "bl_y_plus": _num(
        "target y+", _POS, read_by=COMPONENT_CLASSES,
        default="the project's, else 30",
        description="Target y+ of the first inflation cell."),
    "bl_growth_ratio": _num(
        "inflation growth ratio", (">", 1.0), read_by=COMPONENT_CLASSES,
        default="the project's, else 1.2",
        description="Geometric growth ratio of the inflation layers."),
    "bl_ar_max_prism": _num(
        "prism aspect ratio", (">=", 1.0), read_by=COMPONENT_CLASSES,
        default="the project's, else 5",
        description="Largest prism aspect ratio recommended."),
    "bl_fraction": _num(
        "boundary-layer fraction", ("(]", 0.0, 1.0),
        read_by=COMPONENT_CLASSES,
        default="the project's, else 0.3",
        description="C in dx <= C * delta."),
    "bom_row": {
        "type": "integer", "range": (">=", 1), "what": "BOM row",
        "read_by": COMPONENT_CLASSES, "checked_for": "all",
        "description": "Set by load_bom(): the row the part came from; "
                       "echoed into the result.",
    },
    "bom_errors": {
        "type": "list", "what": "BOM row problems",
        "read_by": COMPONENT_CLASSES, "checked_for": "all",
        "description": "Set by load_bom() on a row that did not load; "
                       "process_part() refuses a part that carries one.",
    },
}

PROJECT_SCHEMA: Dict[str, Dict[str, Any]] = {
    "schema_version": {
        "type": "integer", "what": "input schema version", "default": "1",
        "description": "The input schema the project was written for.",
    },
    "t_fluid_K": _num(
        "fluid temperature", _POS, "K", required=True,
        description="Ambient / underhood air temperature."),
    "t_surr_K": _num(
        "radiation sink temperature", _POS, "K", required=True,
        description="Radiation sink temperature."),
    "max_dt": _num(
        "temperature step", _POS, "K", nullable=True,
        default="the class default",
        description="Accuracy target: temperature change per element."),
    "allowable_flux_error": _num(
        "radiation flux error", _POS, "W/m^2", nullable=True,
        default="the class default",
        description="Radiation linearisation error tolerance."),
    "dt": _num(
        "time step", _POS, "s", nullable=True,
        default="no transient constraint",
        description="Solver time step; turns on the transient constraint."),
    "fo_max": _num(
        "Fourier number limit", _POS, default="0.5",
        description="Fourier number limit."),
    "safety_factor": _num(
        "safety factor", _POS, default="1.0",
        description="Penetration-depth safety factor."),
    "tau_bc": _num(
        "drive-cycle segment", _POS, "s", nullable=True,
        default="no drive-cycle limit",
        description="Drive-cycle segment duration."),
    "transient_scheme": {
        "type": "string", "values": TRANSIENT_SCHEMES, "nullable": True,
        "what": "time-integration scheme",
        "default": "inferred from fo_max (<= 0.5 explicit)",
        "description": "'explicit' or 'implicit'.",
    },
    "t_exh_K": _num(
        "exhaust temperature", _POS, "K", default="1073.15",
        description="Exhaust source temperature for shields."),
    "shield_max_iter": {
        "type": "integer", "range": (">=", 1), "nullable": True,
        "what": "shield iteration limit",
        "default": "50 single-layer, 100 two-layer",
        "description": "Iteration limit of the shield solves.",
    },
    "bl_regime": {
        "type": "string", "values": BL_REGIMES, "nullable": True,
        "what": "boundary-layer regime",
        "default": "no boundary-layer constraint on solids",
        "description": "Turns on the boundary-layer constraint for every "
                       "part.",
    },
    "bl_y_plus": _num("target y+", _POS, default="30",
                      description="Target y+ of the first inflation cell."),
    "bl_growth_ratio": _num(
        "inflation growth ratio", (">", 1.0), default="1.2",
        description="Geometric growth ratio of the inflation layers."),
    "bl_ar_max_prism": _num(
        "prism aspect ratio", (">=", 1.0), default="5",
        description="Largest prism aspect ratio recommended."),
    "bl_fraction": _num(
        "boundary-layer fraction", ("(]", 0.0, 1.0), default="0.3",
        description="C in dx <= C * delta."),
}

MATERIAL_RECORD_SCHEMA: Dict[str, Dict[str, Any]] = {
    "k": _num("thermal conductivity", _POS, "W/m K", required=True,
              description="Thermal conductivity."),
    "rho": _num("density", _POS, "kg/m^3", required=True,
                description="Density."),
    "cp": _num("specific heat", _POS, "J/kg K", required=True,
               description="Specific heat."),
    "description": {
        "type": "string", "what": "description",
        "description": "What the material is, for people.",
    },
    "t_service_max_K": _num(
        "service temperature limit", _POS, "K",
        description="The highest temperature the material may see in "
                    "service.  A part sized above it gets a "
                    "SERVICE_TEMP_EXCEEDED warning."),
    "source": {
        "type": "string", "non_empty": True, "what": "source",
        "description": "Where the record's values come from (a citation).",
    },
}


# ---------------------------------------------------------------------------
#  Warning codes
# ---------------------------------------------------------------------------

WARNING_CODES: Dict[str, Dict[str, Any]] = {
    "SURFACE_DEFAULTED": {
        "severities": ("caution",),
        "meaning": "A surface given neither a treatment nor an emissivity "
                   "took the material-class emissivity.",
    },
    "SHIELD_INPUT_DEFAULTED": {
        "severities": ("caution",),
        "meaning": "A shield input (t_exh_K, h_gap, f12) took its default.",
    },
    "H_CORRELATION_FALLBACK": {
        "severities": ("caution",),
        "meaning": "The h correlation raised; the zone's static h was "
                   "used instead.",
    },
    "FILM_TEMP_OUT_OF_RANGE": {
        "severities": ("caution",),
        "meaning": "Air properties were evaluated outside their fitted "
                   "250-700 K range.",
    },
    "SHIELD_NOT_CONVERGED": {
        "severities": ("warning",),
        "meaning": "The shield solve stopped at its iteration limit.",
    },
    "TRANSIENT_CONFLICT": {
        "severities": ("warning", "caution"),
        "meaning": "The time step cannot integrate the governing element "
                   "(a warning for an explicit scheme, a caution for an "
                   "implicit one); the warning names a remedy.",
    },
    "NO_FINITE_SIZE": {
        "severities": ("warning",),
        "meaning": "No constraint produced a finite element size.",
    },
    "SERVICE_TEMP_EXCEEDED": {
        "severities": ("warning",),
        "meaning": "The part's temperature is above its material's "
                   "service limit (t_service_max_K).",
    },
    "TURB_NAT_HORIZ_UP": {
        "severities": ("warning",),
        "meaning": "A horizontal hot-side-up dead-zone part longer than "
                   "150 mm: natural convection is likely turbulent.",
    },
    "TURB_NAT_VERTICAL": {
        "severities": ("caution",),
        "meaning": "A vertical dead-zone part longer than 600 mm: natural "
                   "convection may be turbulent.",
    },
    "PLASTIC_BIOT_MARGINAL": {
        "severities": ("caution", "warning"),
        "meaning": "A non-metal part's Biot number is near (caution) or "
                   "above (warning) 0.1.",
    },
    "SOLVER_TRANSIENT_RECOMMENDED": {
        "severities": ("warning",),
        "meaning": "The solver advisory recommends a transient solve.",
    },
    "LOW_VELOCITY_HIGH_DT": {
        "severities": ("caution",),
        "meaning": "Low forced velocity with a large surface temperature "
                   "rise: buoyancy may matter.",
    },
}

SEVERITIES = ("info", "caution", "warning")


# ---------------------------------------------------------------------------
#  Result schemas
# ---------------------------------------------------------------------------
#
#  Each result key maps to a dict:
#      types     the value types allowed: "number", "integer", "string",
#                "bool", "dict", "list", "none"
#      required  True when every result carries the key
#      keys      for a dict value: the keys it always carries (it may
#                carry more: the calculators' own outputs, documented in
#                docs/api_reference.md)
#      items     for a list value: the schema of each element (closed)

_SUB_RESULT_KEYS = {
    "conduction": ("max_dx_mm", "q_total"),
    "lateral": ("max_dx_mm",),
    "biot": ("biot", "mesh_type", "min_elements_through_thickness"),
    "radiation": ("max_dx_mm",),
    "shield": ("converged", "iterations", "residual_W_m2"),
    "transient": ("scheme", "feasible", "max_dx_mm", "binding_constraint"),
    "boundary_layer": ("max_dx_surface_mm", "y1_mm", "n_layers",
                       "prisms_recommended", "regime"),
    "h_estimation": ("method", "regime"),
    "solver_advisory": ("severity", "steady_state_ok"),
}

_WARNING_ITEM = {
    "code": {"types": ("string",), "required": True},
    "severity": {"types": ("string",), "required": True},
    "message": {"types": ("string",), "required": True},
    "key": {"types": ("string",), "required": False},
    "remedies": {"types": ("list",), "required": False},
}

_CONSTRAINT_ITEM = {
    "dx_mm": {"types": ("number",), "required": True},
    "source": {"types": ("string",), "required": True},
}

_PROBLEM_ITEM = {
    "key": {"types": ("string", "none"), "required": True},
    "code": {"types": ("string",), "required": True},
    "reason": {"types": ("string",), "required": True},
}


def _sub(name, *types, required=True):
    entry = {"types": types, "required": required}
    if name in _SUB_RESULT_KEYS:
        entry["keys"] = _SUB_RESULT_KEYS[name]
    return entry


RESULT_SCHEMA: Dict[str, Dict[str, Any]] = {
    "schema_version": {"types": ("integer",), "required": True},
    "part_id": {"types": ("string",), "required": True},
    "material": {"types": ("string", "none"), "required": True},
    "component_class": {"types": ("string",), "required": True},
    "t_fluid_K": {"types": ("number",), "required": True},
    "t_fluid_source": {"types": ("string",), "required": True,
                       "values": ("part", "project")},
    "h_used": {"types": ("number", "dict", "none"), "required": True},
    "eps_used": {"types": ("number", "dict"), "required": False},
    "h_estimation": _sub("h_estimation", "dict", required=False),
    "solver_advisory": _sub("solver_advisory", "dict", "none"),
    "conduction": _sub("conduction", "dict", "none"),
    "lateral": _sub("lateral", "dict", required=False),
    "biot": _sub("biot", "dict", "none"),
    "radiation": _sub("radiation", "dict", "none"),
    "shield": _sub("shield", "dict", "none"),
    "transient": _sub("transient", "dict", "none"),
    "f12_used": {"types": ("number",), "required": False},
    "curvature_dx_mm": {"types": ("number",), "required": False},
    "boundary_layer": _sub("boundary_layer", "dict", required=False),
    "governing_dx_mm": {"types": ("number",), "required": True},
    "governing_constraint": {"types": ("string", "none"), "required": True},
    "all_constraints": {"types": ("list",), "required": False,
                        "items": _CONSTRAINT_ITEM},
    "warnings": {"types": ("list",), "required": True,
                 "items": _WARNING_ITEM},
    "error": {"types": ("none",), "required": False},
    "bom_row": {"types": ("integer",), "required": False},
}

ERROR_RESULT_SCHEMA: Dict[str, Dict[str, Any]] = {
    "schema_version": {"types": ("integer",), "required": True},
    "part_id": {"types": ("string",), "required": True},
    "error": {"types": ("string",), "required": True},
    "problems": {"types": ("list",), "required": True,
                 "items": _PROBLEM_ITEM},
    "governing_dx_mm": {"types": ("none",), "required": True},
    "governing_constraint": {"types": ("none",), "required": True},
    "bom_row": {"types": ("integer",), "required": False},
}


# ---------------------------------------------------------------------------
#  Checking one value
# ---------------------------------------------------------------------------

def _problem(key: Optional[str], code: str, reason: str) -> Dict[str, Any]:
    return {"key": key, "code": code, "reason": reason}


def _vocabulary(name: str):
    """The keys of a package table, looked up when a check runs."""
    if name == "CONVECTION_ZONES":
        from thermal_mesh_calculators.zones import CONVECTION_ZONES
        return CONVECTION_ZONES
    # batch imports this module, so it is imported here, not at the top.
    from thermal_mesh_calculators import batch
    return getattr(batch, name)


def _bounds_text(rng) -> str:
    op = rng[0]
    if op in (">", ">="):
        return f"{op} {rng[1]:g}"
    left = "[" if op[0] == "[" else "("
    right = "]" if op[1] == "]" else ")"
    return f"lie in {left}{rng[1]:g}, {rng[2]:g}{right}"


def _number_problem(key: str, spec: dict, value, where: str):
    """A problem for a number, or None.  The three common ranges are
    checked by the calculators' guards, so the words are theirs."""
    rng = spec["range"]
    unit = spec.get("unit", "")
    try:
        if rng == _POS:
            require_positive(key, value, unit)
        elif rng == _NON_NEG:
            require_non_negative(key, value, unit)
        elif rng == _FRACTION:
            require_fraction(key, value)
        else:
            number = _number(key, value)
            op = rng[0]
            if op == ">":
                ok = number > rng[1]
            elif op == ">=":
                ok = number >= rng[1]
            elif op == "(]":
                ok = rng[1] < number <= rng[2]
            else:
                ok = rng[1] <= number <= rng[2]
            if not (math.isfinite(number) and ok):
                text = _bounds_text(rng)
                verb = "be " if op in (">", ">=") else ""
                suffix = " " + unit if unit and op in (">", ">=") else ""
                raise ValueError(f"{key} must {verb}{text}{suffix}, "
                                 f"got {value!r}")
    except TypeError as exc:
        return _problem(key, "WRONG_TYPE", str(exc) + where)
    except ValueError as exc:
        return _problem(key, "OUT_OF_RANGE", str(exc) + where)
    return None


def _integer_problem(key: str, spec: dict, value, where: str):
    rng = spec.get("range", (">=", 1))
    text = f"{key} must be an integer {_bounds_text(rng)}, got {value!r}"
    if isinstance(value, bool) or not isinstance(value, int):
        return _problem(key, "WRONG_TYPE", text + where)
    if not value >= rng[1]:
        return _problem(key, "OUT_OF_RANGE", text + where)
    return None


def _name_message(what: str, container: dict, key: str, names,
                  hint: str = "") -> str:
    """The message PartInputError has always given for a name: it names
    the part, the key and the allowed set."""
    pid = container.get("part_id", "?")
    if key in container:
        head = (f"Unknown {what} {container[key]!r} "
                f"(part {pid!r}, key {key!r}).")
    else:
        head = f"Missing {what}: part {pid!r} has no {key!r} key."
    if hint:
        head = head + " " + hint
    return f"{head} Available: {', '.join(sorted(names))}"


def _string_problem(key: str, spec: dict, value, container: dict,
                    where: str):
    vocab = spec.get("vocabulary")
    if vocab is not None:
        names = _vocabulary(vocab)
        if isinstance(value, str) and value in names:
            return None
        return _problem(key, "UNKNOWN_NAME",
                        _name_message(spec["what"], container, key, names))
    values = spec.get("values")
    if values is not None:
        if isinstance(value, str) and value in values:
            return None
        allowed = ", ".join(repr(v) for v in values)
        if spec.get("nullable"):
            allowed += " or None"
        return _problem(key, "UNKNOWN_NAME",
                        f"{key} must be {allowed}, got {value!r}{where}")
    if not isinstance(value, str) or (spec.get("non_empty") and not value):
        kind = "a non-empty string" if spec.get("non_empty") else "a string"
        return _problem(key, "WRONG_TYPE",
                        f"{key} must be {kind}, got {value!r}{where}")
    return None


def _value_problem(key: str, spec: dict, value, container: dict,
                   where: str):
    """The problem with one value against its schema entry, or None."""
    if value is None and spec.get("nullable"):
        return None
    kind = spec["type"]
    if kind == "number":
        return _number_problem(key, spec, value, where)
    if kind == "integer":
        return _integer_problem(key, spec, value, where)
    if kind == "list":
        if isinstance(value, list):
            return None
        return _problem(key, "WRONG_TYPE",
                        f"{key} must be a list, got {value!r}{where}")
    return _string_problem(key, spec, value, container, where)


def _version_problem(container, where: str):
    """A problem with a dict's schema_version, or None."""
    if "schema_version" not in container:
        return None
    value = container["schema_version"]
    if isinstance(value, bool) or not isinstance(value, int):
        return _problem("schema_version", "WRONG_TYPE",
                        f"schema_version must be an integer, got "
                        f"{value!r}{where}")
    if value != INPUT_SCHEMA_VERSION:
        return _problem(
            "schema_version", "UNSUPPORTED_VERSION",
            f"schema_version {value} is not supported: this library reads "
            f"input schema version {INPUT_SCHEMA_VERSION}{where}")
    return None


def _given(container: dict, key: str) -> bool:
    return key in container and container[key] is not None


def _alternatives_text(key: str, unless) -> str:
    choices = [f"'{key}'"]
    for alt in unless:
        if len(alt) == 1:
            choices.append(f"'{alt[0]}'")
        else:
            choices.append("both " + " and ".join(f"'{k}'" for k in alt))
    if len(choices) == 1:
        return choices[0]
    return ", ".join(choices[:-1]) + " or " + choices[-1]


def _unknown_key_problem(key, schema, where: str):
    known = ", ".join(sorted(schema))
    return _problem(
        key if isinstance(key, str) else repr(key), "UNKNOWN_KEY",
        f"Unknown key {key!r}{where}. Keys starting with "
        f"'{EXTENSION_PREFIX}' are carried through unread. "
        f"Available: {known}")


# ---------------------------------------------------------------------------
#  Input validators
# ---------------------------------------------------------------------------

def _part_where(part: dict) -> str:
    return f" (part {part.get('part_id', '?')!r})"


def _material_label_spec() -> dict:
    """With process_part_from_props() the material is a label of the
    caller's own: any string, and optional."""
    spec = dict(PART_SCHEMA["material"])
    spec.pop("vocabulary")
    spec.pop("required")
    spec["what"] = "material label"
    return spec


def validate_part_input(part, material_from_table: bool = True) -> List[Dict[str, Any]]:
    """
    Check a part dict against PART_SCHEMA.

    Parameters
    ----------
    part : dict
        A part definition (see process_part).
    material_from_table : bool
        True (the default) for process_part(): ``material`` must be a key
        of MATERIALS.  False for process_part_from_props(): ``material``
        is an optional label of the caller's own.

    Returns
    -------
    list of dict — every problem found (key, code, reason), in schema
    order; empty when the part conforms.  A part with an unknown or
    missing component class reports that alone, since its class decides
    which keys are required.  A part carrying ``bom_errors`` reports
    those.
    """
    if not isinstance(part, dict):
        return [_problem(None, "NOT_AN_OBJECT",
                         f"A part must be a dict, got {type(part).__name__}")]
    # A value's problem is worded exactly as the calculators' guards word
    # it; the names, the unknown keys and the missing keys name the part.
    where = ""
    version = _version_problem(part, where)
    if version is not None:
        return [version]
    bom_errors = part.get("bom_errors")
    if isinstance(bom_errors, list) and bom_errors:
        return [dict(p) if isinstance(p, dict) else
                _problem(None, "WRONG_TYPE", f"bom_errors entry {p!r}")
                for p in bom_errors]

    cls_problem = _value_problem(
        "component_class", PART_SCHEMA["component_class"],
        part.get("component_class"), part, where)
    if cls_problem is not None:
        if "component_class" not in part:
            cls_problem["code"] = "MISSING_KEY"
        return [cls_problem]
    cls = part["component_class"]

    problems = [_unknown_key_problem(key, PART_SCHEMA, _part_where(part))
                for key in part
                if key not in PART_SCHEMA
                and not (isinstance(key, str)
                         and key.startswith(EXTENSION_PREFIX))]

    for key, spec in PART_SCHEMA.items():
        if key in ("component_class", "schema_version"):
            continue
        if key == "material" and not material_from_table:
            spec = _material_label_spec()
        if key in part:
            if cls in spec["read_by"] or spec.get("checked_for") == "all":
                found = _value_problem(key, spec, part[key], part, where)
                if found is not None:
                    problems.append(found)
            continue
        for rule in spec.get("required", ()):
            if cls not in rule["classes"]:
                continue
            if any(all(_given(part, k) for k in alt)
                   for alt in rule["unless"]):
                continue
            hint = ""
            if rule["unless"]:
                hint = ("Give " + _alternatives_text(key, rule["unless"])
                        + ".")
            if "vocabulary" in spec:
                reason = _name_message(spec["what"], part, key,
                                       _vocabulary(spec["vocabulary"]),
                                       hint)
            else:
                reason = (f"Missing {spec['what']}: part "
                          f"{part.get('part_id', '?')!r} has no {key!r} "
                          f"key." + (" " + hint if hint else ""))
            problems.append(_problem(key, "MISSING_KEY", reason))
    return problems


def _flat_problems(container, schema: dict, where: str,
                   what: str) -> List[Dict[str, Any]]:
    """Problems of a dict whose keys apply whatever the class: a project
    or a material record."""
    if not isinstance(container, dict):
        return [_problem(None, "NOT_AN_OBJECT",
                         f"A {what} must be a dict, got "
                         f"{type(container).__name__}")]
    problems = []
    if "schema_version" in schema:
        version = _version_problem(container, where)
        if version is not None:
            return [version]
    for key in container:
        if key not in schema and not (
                isinstance(key, str) and key.startswith(EXTENSION_PREFIX)):
            problems.append(_unknown_key_problem(key, schema, where))
    for key, spec in schema.items():
        if key == "schema_version":
            continue
        if key in container:
            found = _value_problem(key, spec, container[key], container,
                                   where)
            if found is not None:
                problems.append(found)
        elif spec.get("required"):
            problems.append(_problem(
                key, "MISSING_KEY",
                f"Missing {spec['what']}: the {what} has no {key!r} key."))
    return problems


def validate_project_input(project) -> List[Dict[str, Any]]:
    """
    Check a project dict against PROJECT_SCHEMA.

    Every key the project gives is checked, whatever the parts: a project
    value applies to every part.

    Returns
    -------
    list of dict — every problem found (key, code, reason); empty when
    the project conforms.
    """
    return _flat_problems(project, PROJECT_SCHEMA, " (project)", "project")


def validate_material_record(record) -> List[Dict[str, Any]]:
    """
    Check a material record against MATERIAL_RECORD_SCHEMA.

    k, rho and cp are required; description, t_service_max_K and source
    are optional, and their absence is not a problem.

    Returns
    -------
    list of dict — every problem found (key, code, reason); empty when
    the record conforms.
    """
    return _flat_problems(record, MATERIAL_RECORD_SCHEMA,
                          " (material record)", "material record")


# ---------------------------------------------------------------------------
#  Result validator
# ---------------------------------------------------------------------------

def _type_name(value) -> str:
    if value is None:
        return "none"
    if isinstance(value, bool):
        return "bool"
    if isinstance(value, int):
        return "integer"
    if isinstance(value, float):
        return "number"
    if isinstance(value, str):
        return "string"
    if isinstance(value, dict):
        return "dict"
    if isinstance(value, list):
        return "list"
    return type(value).__name__


def _type_ok(value, types) -> bool:
    name = _type_name(value)
    # An integer is a number; a bool is neither.
    return name in types or (name == "integer" and "number" in types)


def _closed_dict_problems(value, schema: dict, path: str) -> list:
    problems = []
    if not isinstance(value, dict):
        return [_problem(path, "WRONG_TYPE",
                         f"{path} must be a dict, got {value!r}")]
    for key in value:
        if key not in schema:
            problems.append(_problem(f"{path}.{key}", "UNKNOWN_KEY",
                                     f"{path} carries an undeclared key "
                                     f"{key!r}"))
    for key, spec in schema.items():
        if key not in value:
            if spec["required"]:
                problems.append(_problem(f"{path}.{key}", "MISSING_KEY",
                                         f"{path} has no {key!r} key"))
        elif not _type_ok(value[key], spec["types"]):
            problems.append(_problem(
                f"{path}.{key}", "WRONG_TYPE",
                f"{path}.{key} must be {' or '.join(spec['types'])}, got "
                f"{value[key]!r}"))
    return problems


def _warning_problems(warning, path: str) -> list:
    problems = _closed_dict_problems(warning, _WARNING_ITEM, path)
    if problems:
        return problems
    code = warning["code"]
    if code not in WARNING_CODES:
        return [_problem(f"{path}.code", "UNKNOWN_NAME",
                         f"{path} has the unregistered warning code "
                         f"{code!r}")]
    if warning["severity"] not in WARNING_CODES[code]["severities"]:
        return [_problem(f"{path}.severity", "UNKNOWN_NAME",
                         f"{path} ({code}) has severity "
                         f"{warning['severity']!r}; {code} is "
                         f"{' or '.join(WARNING_CODES[code]['severities'])}")]
    return []


def _keys_problems(result: dict, schema: dict) -> list:
    problems = []
    for key in result:
        if key not in schema:
            problems.append(_problem(key, "UNKNOWN_KEY",
                                     f"undeclared result key {key!r}"))
    for key, spec in schema.items():
        if key not in result:
            if spec["required"]:
                problems.append(_problem(key, "MISSING_KEY",
                                         f"the result has no {key!r} key"))
            continue
        value = result[key]
        if not _type_ok(value, spec["types"]):
            problems.append(_problem(
                key, "WRONG_TYPE",
                f"{key} must be {' or '.join(spec['types'])}, got "
                f"{value!r}"))
            continue
        if "values" in spec and value not in spec["values"]:
            problems.append(_problem(key, "UNKNOWN_NAME",
                                     f"{key} must be one of "
                                     f"{', '.join(spec['values'])}, got "
                                     f"{value!r}"))
        if isinstance(value, dict) and "keys" in spec:
            for sub in spec["keys"]:
                if sub not in value:
                    problems.append(_problem(f"{key}.{sub}", "MISSING_KEY",
                                             f"{key} has no {sub!r} key"))
        if isinstance(value, list) and "items" in spec:
            for i, item in enumerate(value):
                path = f"{key}[{i}]"
                if spec["items"] is _WARNING_ITEM:
                    problems += _warning_problems(item, path)
                else:
                    problems += _closed_dict_problems(item, spec["items"],
                                                      path)
    return problems


def validate_result(result) -> List[Dict[str, Any]]:
    """
    Check a result against the result schema.

    A dict whose ``error`` is a non-empty string is checked against
    ERROR_RESULT_SCHEMA (process_batch()'s result for a part that could
    not be sized); any other against RESULT_SCHEMA.  Every top-level key
    must be declared, every warning code registered in WARNING_CODES
    with one of its severities, and the governing size must be the
    smallest candidate in ``all_constraints``.

    Returns
    -------
    list of dict — every problem found (key, code, reason); empty when
    the result conforms.
    """
    if not isinstance(result, dict):
        return [_problem(None, "NOT_AN_OBJECT",
                         f"A result must be a dict, got "
                         f"{type(result).__name__}")]
    is_error = isinstance(result.get("error"), str) and result["error"] != ""
    schema = ERROR_RESULT_SCHEMA if is_error else RESULT_SCHEMA
    problems = _keys_problems(result, schema)
    version = result.get("schema_version")
    if _type_name(version) == "integer" and version != RESULT_SCHEMA_VERSION:
        problems.append(_problem(
            "schema_version", "UNSUPPORTED_VERSION",
            f"result schema_version {version} is not "
            f"{RESULT_SCHEMA_VERSION}"))
    if is_error:
        for i, item in enumerate(result.get("problems") or ()):
            code = item.get("code") if isinstance(item, dict) else None
            if code is not None and code not in PROBLEM_CODES:
                problems.append(_problem(
                    f"problems[{i}].code", "UNKNOWN_NAME",
                    f"problems[{i}] has the unregistered code {code!r}"))
        return problems

    dx = result.get("governing_dx_mm")
    if (isinstance(dx, (int, float)) and not isinstance(dx, bool)
            and not dx > 0):                          # NaN fails too
        problems.append(_problem("governing_dx_mm", "OUT_OF_RANGE",
                                 f"governing_dx_mm must be > 0 (inf when "
                                 f"no constraint is finite), got {dx!r}"))
    candidates = result.get("all_constraints")
    if (isinstance(candidates, list) and candidates
            and all(isinstance(c, dict) and _type_ok(c.get("dx_mm"),
                                                     ("number",))
                    for c in candidates)):
        smallest = min(candidates, key=lambda c: c["dx_mm"])
        if (dx != smallest["dx_mm"]
                or result.get("governing_constraint") != smallest["source"]):
            problems.append(_problem(
                "governing_dx_mm", "INCONSISTENT",
                f"the governing size {dx!r} "
                f"({result.get('governing_constraint')!r}) is not the "
                f"smallest candidate, {smallest['dx_mm']!r} "
                f"({smallest['source']!r})"))
    return problems


__all__ = [
    "INPUT_SCHEMA_VERSION",
    "RESULT_SCHEMA_VERSION",
    "COMPONENT_CLASSES",
    "SOLID_CLASSES",
    "NON_SHIELD_SOLID_CLASSES",
    "SHIELD_CLASSES",
    "FLUID_CLASSES",
    "BL_REGIMES",
    "TRANSIENT_SCHEMES",
    "EXTENSION_PREFIX",
    "PROBLEM_CODES",
    "PART_SCHEMA",
    "PROJECT_SCHEMA",
    "MATERIAL_RECORD_SCHEMA",
    "RESULT_SCHEMA",
    "ERROR_RESULT_SCHEMA",
    "WARNING_CODES",
    "SEVERITIES",
    "validate_part_input",
    "validate_project_input",
    "validate_material_record",
    "validate_result",
]
