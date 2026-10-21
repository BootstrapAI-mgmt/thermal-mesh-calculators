"""
A Sizing Hook for an External Mesher
====================================

One function, size_for_fea(), for a tool outside this package that meshes a
wall and needs its through-thickness cell size: it gives a dict and gets a
dict back, both plain JSON, without reading anything else of the package.

The request names the wall's material (a MATERIALS key or an alias, or the
properties k, rho and cp themselves, as process_part_from_props() takes
them), its thickness and the boundary conditions the sizing rule needs: the
convective coefficient h, the surface, fluid and radiation-sink
temperatures, the emissivity and the temperature step max_dt.  The keys
are the part schema's names (thickness_mm, t_surf_K, ...), and
schema.FEA_SIZING_REQUEST_SCHEMA declares them.

The answer is BoundaryDrivenConductionCalculator.size_wall() — the same
thickness-aware rule process_part() applies — evaluated on that wall and
re-expressed for a mesher:

    char_size_mm        the through-thickness cell size, thickness / n_cells
    n_cells             the cell count through the thickness
    biot                the wall's Biot number (None when undefined)
    binding_constraint  what set n_cells: "conduction" (the max_dt relation),
                        "wall_thickness" (the whole wall is within max_dt, so
                        one cell spans it), "max_cells" (the cap) or
                        "near_equilibrium" (convection and radiation cancel)
    rationale           one line saying why

with every other quantity size_wall() reports, under its own name, and the
material's resolved properties.  The result carries schema_version
(schema.FEA_SIZING_SCHEMA_VERSION); schema.validate_fea_sizing() checks it.
JSON has no infinity, so a quantity size_wall() reports as infinite (the
raw cell size of a wall near equilibrium, the Biot number of a wall whose
surface is at the fluid temperature) is None.

No file is read or written and nothing is printed.
"""

import math
from typing import Any, Dict, List, Optional

from thermal_mesh_calculators.batch import MATERIALS, PartInputError
from thermal_mesh_calculators.conduction import BoundaryDrivenConductionCalculator
from thermal_mesh_calculators.intake import resolve_material
from thermal_mesh_calculators.schema import (
    FEA_SIZING_SCHEMA_VERSION,
    validate_fea_sizing_request,
)

DEFAULT_MAX_CELLS = 100


def _finite_or_none(value: float) -> Optional[float]:
    return value if math.isfinite(value) else None


def _request_error(problems: List[Dict[str, Any]]) -> PartInputError:
    if len(problems) == 1:
        message = problems[0]["reason"]
    else:
        message = (f"{len(problems)} problems with the fea sizing request:"
                   + "".join("\n  - " + p["reason"] for p in problems))
    return PartInputError(message, problems)


def _rationale(binding: str, wall: dict, thickness_mm: float, max_dt: float,
               max_cells: int, raw_cells: int) -> str:
    n = wall["n_cells"]
    size = wall["dx_used_mm"]
    biot = wall["biot"]
    bi = f"Bi = {biot:.3g}" if math.isfinite(biot) else "Bi undefined"
    if binding == "near_equilibrium":
        return (f"One cell of {thickness_mm:.4g} mm: convection and radiation "
                f"nearly cancel, so the steady balance cannot size this wall; "
                f"size it geometrically or from the transient.")
    if binding == "max_cells":
        return (f"Capped at {max_cells} cells of {size:.4g} mm: max_dt = "
                f"{max_dt:g} K would need {raw_cells} cells through the "
                f"{thickness_mm:.4g} mm wall ({bi}).")
    if binding == "wall_thickness":
        return (f"One cell of {thickness_mm:.4g} mm spans the wall: the drop "
                f"across the whole wall, {wall['dt_wall_K']:.3g} K, is within "
                f"max_dt = {max_dt:g} K ({bi}, {wall['regime']}).")
    return (f"{n} cells of {size:.4g} mm through the {thickness_mm:.4g} mm "
            f"wall keep the drop across each cell within max_dt = "
            f"{max_dt:g} K ({bi}, {wall['regime']}).")


def size_for_fea(request: dict) -> dict:
    """
    Through-thickness sizing of a wall, for an external mesher.

    Parameters
    ----------
    request : dict
        The keys of schema.FEA_SIZING_REQUEST_SCHEMA:
            material      : str   — a MATERIALS key or an alias; or instead
            k, rho, cp    : float — W/m K, kg/m^3, J/kg K (all three)
            thickness_mm  : float — wall thickness (mm)
            h             : float — convective coefficient (W/m^2 K)
            t_surf_K      : float — surface temperature (K)
            t_fluid_K     : float — free-stream / bulk fluid temperature (K)
            t_surr_K      : float — radiation sink temperature (K)
            epsilon       : float — surface emissivity (0-1)
            max_dt        : float — temperature step per cell (K)
            max_cells     : int, optional — cell cap (default 100)
            schema_version: int, optional — absent means 1
        Keys starting with "x_" are carried by the caller and never read.

    Returns
    -------
    dict — plain JSON, the keys of schema.FEA_SIZING_RESULT_SCHEMA:
        schema_version, char_size_mm, n_cells, thickness_mm, biot,
        binding_constraint, rationale, and size_wall()'s regime, dx_raw_mm,
        dx_exceeds_part, dt_wall_K, h_effective, q_conv, q_rad, q_total,
        rad_fraction, recommendation and warnings, plus material
        {name, k, rho, cp, alpha_m2_s}.  char_size_mm is size_wall()'s
        dx_used_mm; an infinite quantity is None.

    Raises
    ------
    PartInputError
        For a request that is not a dict, a key the schema does not
        declare, a value of the wrong type or out of its range, a missing
        key, an unknown material, both a material and properties, or an
        unsupported schema_version — every problem at once, each with key,
        code and reason (``problems``).
    """
    problems = validate_fea_sizing_request(request)
    if problems:
        raise _request_error(problems)

    name: Optional[str] = None
    k: Any
    rho: Any
    cp: Any
    if "material" in request:
        key = resolve_material(request["material"])
        record: Dict[str, Any] = MATERIALS[key]
        name = key
        k, rho, cp = record["k"], record["rho"], record["cp"]
    else:
        k, rho, cp = request["k"], request["rho"], request["cp"]

    thickness_mm = float(request["thickness_mm"])
    thickness_m = thickness_mm / 1000.0
    max_dt = request["max_dt"]
    max_cells = request.get("max_cells", DEFAULT_MAX_CELLS)
    wall = BoundaryDrivenConductionCalculator.size_wall(
        k=k, h=request["h"], t_surf=request["t_surf_K"],
        t_fluid=request["t_fluid_K"], epsilon=request["epsilon"],
        t_surr=request["t_surr_K"], max_dt=max_dt,
        thickness_m=thickness_m, max_cells=max_cells,
    )

    # What set the count: size_wall()'s own clamp, recomputed the same way.
    raw_cells = 1
    if wall["regime"] == "near_equilibrium":
        binding = "near_equilibrium"
    else:
        dx_raw_m = k * max_dt / abs(wall["q_conv"] + wall["q_rad"])
        raw_cells = max(1, math.ceil(thickness_m / dx_raw_m))
        if raw_cells > max_cells:
            binding = "max_cells"
        elif raw_cells == 1:
            binding = "wall_thickness"
        else:
            binding = "conduction"

    return {
        "schema_version": FEA_SIZING_SCHEMA_VERSION,
        "char_size_mm": wall["dx_used_mm"],
        "n_cells": wall["n_cells"],
        "thickness_mm": thickness_mm,
        "biot": _finite_or_none(wall["biot"]),
        "binding_constraint": binding,
        "rationale": _rationale(binding, wall, thickness_mm, max_dt,
                                max_cells, raw_cells),
        "regime": wall["regime"],
        "dx_raw_mm": _finite_or_none(wall["dx_raw_mm"]),
        "dx_exceeds_part": wall["dx_exceeds_part"],
        "dt_wall_K": wall["dt_wall_K"],
        "h_effective": _finite_or_none(wall["h_effective"]),
        "q_conv": wall["q_conv"],
        "q_rad": wall["q_rad"],
        "q_total": wall["q_total"],
        "rad_fraction": wall["rad_fraction"],
        "recommendation": wall["recommendation"],
        "warnings": list(wall["warnings"]),
        "material": {
            "name": name,
            "k": k,
            "rho": rho,
            "cp": cp,
            "alpha_m2_s": k / (rho * cp),
        },
    }


__all__ = ["DEFAULT_MAX_CELLS", "size_for_fea"]
