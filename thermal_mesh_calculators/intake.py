"""
Part Intake: BOM Files, Names and Classes
==========================================

The steps between a caller's own data and process_batch():

    load_bom()               — read a bill of materials from a CSV or JSON
                               file into part dicts, checking every row
                               against schema.PART_SCHEMA and reporting
                               every problem of every row instead of
                               stopping at the first
    resolve_material(),      — turn a common spelling of a material or a
    resolve_zone()             zone ("Mild steel", "PA66-GF30",
                               "beside engine") into its canonical key
    MATERIAL_ALIASES,        — the alias tables those two read
    ZONE_ALIASES
    infer_component_class()  — the component class a part's keys imply,
                               for data that does not carry one

Names are compared after normalisation: lower case, every run of
characters other than letters and digits made one underscore, and the
US spellings aluminum, galvanized, gray, fiber, fiberglass and molding
read as the tables' aluminium, galvanised, grey, fibre, fibreglass and
moulding.  A canonical key always resolves to itself.

An alias names the same material or zone as its key: a synonym, a
designation the entry's own description lists, or a word order.  A
family name ("aluminium", "steel", "nylon") is not an alias, because it
does not say which entry applies: resolve_material() refuses it and
names the entries it could mean.
"""

import csv
import io
import json
import os
import re
from typing import Any, Dict, List, Optional, Tuple

from thermal_mesh_calculators.batch import (
    CLASS_DEFAULTS,
    MATERIALS,
    PartInputError,
)
from thermal_mesh_calculators.schema import (
    INPUT_SCHEMA_VERSION,
    PART_SCHEMA,
    SOLID_CLASSES,
    validate_part_input,
)
from thermal_mesh_calculators.zones import CONVECTION_ZONES


# ---------------------------------------------------------------------------
#  Alias tables
# ---------------------------------------------------------------------------
#
#  Keys are normalised names (see _normalise); values are canonical keys.
#  Every alias below is a name the canonical entry's own description in
#  MATERIALS uses (a synonym or designation), or a common spelling of it.

MATERIAL_ALIASES: Dict[str, str] = {
    # --- Steels ---
    "mild_steel": "steel_mild",
    "carbon_steel": "steel_mild",
    "hsla": "steel_high_strength",
    "hsla_steel": "steel_high_strength",
    "high_strength_steel": "steel_high_strength",
    "high_strength_low_alloy_steel": "steel_high_strength",
    "dp600": "steel_high_strength",
    "galvanised_steel": "steel_galvanised",
    "zinc_coated_steel": "steel_galvanised",
    "stainless_304": "steel_stainless_304",
    "stainless_steel_304": "steel_stainless_304",
    "ss304": "steel_stainless_304",
    "ss_304": "steel_stainless_304",
    "stainless_316": "steel_stainless_304",
    "stainless_steel_316": "steel_stainless_304",
    "ss316": "steel_stainless_304",
    "ss_316": "steel_stainless_304",
    "stainless_409": "steel_stainless_409",
    "stainless_steel_409": "steel_stainless_409",
    "ss409": "steel_stainless_409",
    "ss_409": "steel_stainless_409",
    "stainless_430": "steel_stainless_409",
    "stainless_steel_430": "steel_stainless_409",
    "ss430": "steel_stainless_409",
    "ss_430": "steel_stainless_409",
    # --- Aluminium ---
    "al_6061": "aluminium_6061",
    "al6061": "aluminium_6061",
    "aa6061": "aluminium_6061",
    "6061": "aluminium_6061",
    "6061_t6": "aluminium_6061",
    "aluminium_6061_t6": "aluminium_6061",
    "al_5052": "aluminium_5052",
    "al5052": "aluminium_5052",
    "aa5052": "aluminium_5052",
    "5052": "aluminium_5052",
    "5052_h32": "aluminium_5052",
    "aluminium_5052_h32": "aluminium_5052",
    "a356": "cast_aluminium",
    "a356_t6": "cast_aluminium",
    "aluminium_a356": "cast_aluminium",
    "cast_aluminium_a356": "cast_aluminium",
    "a380": "cast_aluminium_a380",
    "aluminium_a380": "cast_aluminium_a380",
    "die_cast_aluminium_a380": "cast_aluminium_a380",
    # --- Cast iron ---
    "grey_cast_iron": "cast_iron",
    "grey_iron": "cast_iron",
    "ductile_iron": "cast_iron_ductile",
    "ductile_cast_iron": "cast_iron_ductile",
    "nodular_iron": "cast_iron_ductile",
    "nodular_cast_iron": "cast_iron_ductile",
    "sg_iron": "cast_iron_ductile",
    # --- Other metals ---
    "pure_copper": "copper",
    "brass_70cu_30zn": "brass",
    "az91": "magnesium_az91",
    "az91d": "magnesium_az91",
    "mg_az91": "magnesium_az91",
    "magnesium_az91d": "magnesium_az91",
    "ti_6al_4v": "titanium_6al4v",
    "ti6al4v": "titanium_6al4v",
    "titanium_ti_6al_4v": "titanium_6al4v",
    "zinc_alloy": "zinc",
    "inconel_625": "nickel_alloy",
    "inconel625": "nickel_alloy",
    "nickel_alloy_625": "nickel_alloy",
    # --- Plastics, glass-filled ---
    "pa66_gf30": "plastic_pa66_gf30",
    "pa66gf30": "plastic_pa66_gf30",
    "nylon_66_gf30": "plastic_pa66_gf30",
    "nylon66_gf30": "plastic_pa66_gf30",
    "pa6_gf30": "plastic_pa6_gf30",
    "pa6gf30": "plastic_pa6_gf30",
    "nylon_6_gf30": "plastic_pa6_gf30",
    "nylon6_gf30": "plastic_pa6_gf30",
    "pp_gf30": "plastic_pp_gf30",
    "ppgf30": "plastic_pp_gf30",
    "pbt_gf30": "plastic_pbt_gf30",
    "pbtgf30": "plastic_pbt_gf30",
    "pps_gf40": "plastic_pps_gf40",
    "ppsgf40": "plastic_pps_gf40",
    # --- Plastics, unfilled ---
    "pa66": "plastic_pa66",
    "nylon_66": "plastic_pa66",
    "nylon66": "plastic_pa66",
    "pp": "plastic_pp",
    "polypropylene": "plastic_pp",
    "hdpe": "plastic_hdpe",
    "high_density_polyethylene": "plastic_hdpe",
    "abs": "plastic_abs",
    "pc": "plastic_pc",
    "polycarbonate": "plastic_pc",
    "pc_abs": "plastic_pc_abs",
    "pom": "plastic_acetal_pom",
    "acetal": "plastic_acetal_pom",
    "polyoxymethylene": "plastic_acetal_pom",
    "pet": "plastic_pet",
    "polyethylene_terephthalate": "plastic_pet",
    # --- Rubbers ---
    "epdm": "rubber_epdm",
    "epdm_rubber": "rubber_epdm",
    "silicone": "rubber_silicone",
    "silicone_rubber": "rubber_silicone",
    "nbr": "rubber_nbr",
    "nitrile": "rubber_nbr",
    "nitrile_rubber": "rubber_nbr",
    "natural_rubber": "rubber_natural",
    "neoprene": "rubber_cr",
    "chloroprene": "rubber_cr",
    "cr_rubber": "rubber_cr",
    "fkm": "rubber_fkm",
    "viton": "rubber_fkm",
    "fkm_rubber": "rubber_fkm",
    # --- Composites and specialty ---
    "smc": "composite_smc",
    "sheet_moulding_compound": "composite_smc",
    "cfrp": "composite_cfrp",
    "carbon_fibre_reinforced_polymer": "composite_cfrp",
    "carbon_fibre_reinforced_plastic": "composite_cfrp",
    "soda_lime_glass": "glass_soda_lime",
    "alumina": "ceramic_alumina",
    "cordierite": "ceramic_cordierite",
    "fibreglass_insulation": "insulation_fibreglass",
    "glass_wool": "insulation_fibreglass",
    "ceramic_fibre_blanket": "insulation_ceramic_blanket",
    "ceramic_blanket": "insulation_ceramic_blanket",
}

#  Zone aliases: the zone labels of the map in zones.py's docstring, in the
#  order "<where> <what>" as well as the keys' "<what> <where>".
ZONE_ALIASES: Dict[str, str] = {
    "cooling_pack_outlet": "cooling_pack_downstream",
    "behind_cooling_pack": "cooling_pack_downstream",
    "beside_engine": "engine_beside",
    "below_engine": "engine_below",
    "above_engine": "engine_above",
    "beside_exhaust": "exhaust_beside",
    "below_exhaust": "exhaust_below",
    "above_exhaust": "exhaust_above",
    "exhaust_gas": "exhaust_internal",
    "internal_exhaust": "exhaust_internal",
    "clutch_outlet": "clutch_outlet_downstream",
    "shield_gap": "shield_gap_confined",
}

# US spellings of the tables' words.
_SPELLINGS = {
    "aluminum": "aluminium",
    "galvanized": "galvanised",
    "gray": "grey",
    "fiber": "fibre",
    "fiberglass": "fibreglass",
    "molding": "moulding",
}

_TOKEN = re.compile(r"[0-9a-z]+")


def _normalise(name: str) -> str:
    tokens = _TOKEN.findall(name.lower())
    return "_".join(_SPELLINGS.get(t, t) for t in tokens)


def _resolve(name, table, aliases) -> Optional[str]:
    if not isinstance(name, str):
        return None
    if name in table:
        return name
    normal = _normalise(name)
    if normal in table:
        return normal
    return aliases.get(normal)


def _candidates(name, table) -> List[str]:
    """Canonical keys that carry every word of name: what a family name
    such as "aluminium" could mean."""
    words = set(_normalise(name).split("_")) if isinstance(name, str) else set()
    words.discard("")
    if not words:
        return []
    return sorted(k for k in table if words <= set(k.split("_")))


def _unknown_name(what: str, name, table, aliases_name: str) -> PartInputError:
    message = f"Unknown {what} {name!r}."
    candidates = _candidates(name, table)
    if candidates and len(candidates) < len(table):
        message += (f" {name!r} could be any of: "
                    f"{', '.join(candidates)}; give the one you mean."
                    if len(candidates) > 1 else
                    f" Did you mean {candidates[0]!r}? Give the key itself.")
    message += (f" Available: {', '.join(sorted(table))}"
                f" (and the spellings in {aliases_name})")
    return PartInputError(message, [{"key": None, "code": "UNKNOWN_NAME",
                                     "reason": message}])


def resolve_material(name: str) -> str:
    """
    The MATERIALS key a material name refers to.

    Parameters
    ----------
    name : str
        A MATERIALS key, or a spelling MATERIAL_ALIASES lists ("Mild
        steel", "PA66-GF30", "SS 316", "aluminum 6061").

    Returns
    -------
    str — the canonical key.

    Raises
    ------
    PartInputError
        For a name that is neither, naming the allowed set; a family name
        ("aluminium") also gets the entries it could mean.
    """
    key = _resolve(name, MATERIALS, MATERIAL_ALIASES)
    if key is None:
        raise _unknown_name("material", name, MATERIALS, "MATERIAL_ALIASES")
    return key


def resolve_zone(name: str) -> str:
    """
    The CONVECTION_ZONES key a zone name refers to.

    Parameters
    ----------
    name : str
        A CONVECTION_ZONES key (the two legacy keys included), or a
        spelling ZONE_ALIASES lists ("Beside engine", "exhaust gas").

    Returns
    -------
    str — the key.

    Raises
    ------
    PartInputError
        For a name that is neither, naming the allowed set.
    """
    key = _resolve(name, CONVECTION_ZONES, ZONE_ALIASES)
    if key is None:
        raise _unknown_name("convection zone", name, CONVECTION_ZONES,
                            "ZONE_ALIASES")
    return key


# ---------------------------------------------------------------------------
#  Component class inference
# ---------------------------------------------------------------------------

# Keys only a two-layer shield reads, then keys only a shield reads.
TWO_LAYER_SHIELD_KEYS = ("surface_g1", "surface_g2", "eps_g1", "eps_g2",
                         "h_gap", "f12")
SHIELD_KEYS = ("surface_in", "surface_out", "eps_in", "eps_out",
               "h_in_override", "h_out_override", "convection_zone_in",
               "convection_zone_out")
# Zones whose air is the exhaust gas, then zones beside the exhaust.
EXHAUST_ZONES = ("exhaust_internal",)
EXHAUST_ADJACENT_ZONES = ("exhaust_beside", "exhaust_below", "exhaust_above",
                          "near_exhaust_natural")


def infer_component_class(part: dict) -> str:
    """
    The component class a part's keys imply.

    The first rule that applies decides:

        1. component_class given        -> that class (checked)
        2. any of TWO_LAYER_SHIELD_KEYS -> "multilayer_shield"
        3. any of SHIELD_KEYS           -> "shield"
        4. convection_zone in EXHAUST_ZONES          -> "exhaust"
        5. convection_zone in EXHAUST_ADJACENT_ZONES -> "exhaust_adjacent"
        6. otherwise                    -> "structural"

    "fluid" is never inferred: no part key says a region is a fluid, so a
    fluid region names its class.  The zone is read through
    resolve_zone() when it resolves.  The rule is a default for data that
    carries no class: a part that is an exhaust component in a zone
    beside the exhaust should say so.

    Parameters
    ----------
    part : dict — a part definition (see process_part)

    Returns
    -------
    str — a key of CLASS_DEFAULTS.

    Raises
    ------
    PartInputError
        When the part gives a component_class that is not a key of
        CLASS_DEFAULTS.
    """
    if "component_class" in part:
        cls = part["component_class"]
        if isinstance(cls, str) and cls in CLASS_DEFAULTS:
            return cls
        problems = validate_part_input({"part_id": part.get("part_id", "?"),
                                        "component_class": cls})
        raise PartInputError(problems[0]["reason"], problems[:1])
    if any(key in part for key in TWO_LAYER_SHIELD_KEYS):
        return "multilayer_shield"
    if any(key in part for key in SHIELD_KEYS):
        return "shield"
    zone = part.get("convection_zone")
    resolved = _resolve(zone, CONVECTION_ZONES, ZONE_ALIASES)
    zone = resolved if resolved is not None else zone
    if zone in EXHAUST_ZONES:
        return "exhaust"
    if zone in EXHAUST_ADJACENT_ZONES:
        return "exhaust_adjacent"
    return "structural"


# ---------------------------------------------------------------------------
#  BOM loader
# ---------------------------------------------------------------------------

_LOADER_KEYS = ("bom_row", "bom_errors")
_ZONE_KEYS = ("convection_zone", "convection_zone_in", "convection_zone_out")


def _problem(key: Optional[str], code: str, reason: str) -> Dict[str, Any]:
    return {"key": key, "code": code, "reason": reason}


def _read_text(source) -> Tuple[str, str]:
    """The text of a path or an open file, and a name for messages."""
    if hasattr(source, "read"):
        data = source.read()
        name = getattr(source, "name", None)
        name = name if isinstance(name, str) else "<file>"
    else:
        path = os.fspath(source)
        name = os.fsdecode(path)
        with open(path, "rb") as handle:
            data = handle.read()
    if isinstance(data, bytes):
        try:
            return data.decode("utf-8-sig"), name
        except UnicodeDecodeError as exc:
            raise ValueError(f"{name}: not UTF-8 text ({exc}); save it as "
                             f"UTF-8") from None
    if data.startswith("﻿"):
        data = data[1:]
    return data, name


def _bom_format(fmt, name: str) -> str:
    if fmt is None:
        ext = os.path.splitext(name)[1].lower()
        if ext not in (".csv", ".json"):
            raise ValueError(f"cannot tell the format of {name!r} from its "
                             f"extension; pass fmt='csv' or fmt='json'")
        return ext[1:]
    if fmt not in ("csv", "json"):
        raise ValueError(f"fmt must be 'csv' or 'json', got {fmt!r}")
    return fmt


def _csv_cell(key: str, text: str):
    """A CSV cell as the type the schema declares for its key; a cell that
    does not parse stays text, which the schema check then reports."""
    spec = PART_SCHEMA.get(key)
    kind = spec["type"] if spec is not None else "string"
    if kind == "number":
        try:
            return float(text)
        except ValueError:
            return text
    if kind == "integer":
        try:
            return int(text)
        except ValueError:
            return text
    return text


def _csv_records(text: str, name: str) -> list:
    """(row, values, problems) per data row; the row is the line the record
    ends on, the header being line 1."""
    reader = csv.reader(io.StringIO(text, newline=""))
    try:
        header = [cell.strip() for cell in next(reader)]
    except StopIteration:
        return []
    duplicates = sorted({c for c in header if c and header.count(c) > 1})
    if duplicates:
        raise ValueError(f"{name}: duplicate column(s) in the header: "
                         f"{', '.join(duplicates)}")
    records: List[Tuple[int, Optional[Dict[str, Any]], List[Dict[str, Any]]]] = []
    for cells in reader:
        cells = [cell.strip() for cell in cells]
        if not any(cells):
            continue          # a blank line, or a row of empty cells
        problems = []
        if len(cells) > len(header) and any(cells[len(header):]):
            problems.append(_problem(
                None, "MALFORMED_ROW",
                f"the row has {len(cells)} values and the header "
                f"{len(header)} columns"))
        values: Dict[str, Any] = {}
        for column, cell in zip(header, cells):
            if cell != "":
                values[column] = _csv_cell(column, cell)
        records.append((reader.line_num, values, problems))
    return records


def _json_records(text: str, name: str) -> list:
    try:
        data = json.loads(text)
    except ValueError as exc:
        raise ValueError(f"{name}: not valid JSON ({exc})") from None
    if isinstance(data, dict):
        extra = sorted(set(data) - {"parts", "schema_version"})
        if "parts" not in data or extra:
            raise ValueError(
                f"{name}: a JSON BOM is an array of part objects, or an "
                f"object with a 'parts' array (and an optional "
                f"'schema_version')"
                + (f"; found {', '.join(extra)}" if extra else ""))
        version = data.get("schema_version", INPUT_SCHEMA_VERSION)
        if isinstance(version, bool) or version != INPUT_SCHEMA_VERSION:
            raise ValueError(f"{name}: schema_version {version!r} is not "
                             f"supported: this library reads input schema "
                             f"version {INPUT_SCHEMA_VERSION}")
        data = data["parts"]
    if not isinstance(data, list):
        raise ValueError(f"{name}: a JSON BOM is an array of part objects, "
                         f"got {type(data).__name__}")
    records: List[Tuple[int, Optional[Dict[str, Any]], List[Dict[str, Any]]]] = []
    for index, item in enumerate(data, start=1):
        if isinstance(item, dict):
            records.append((index, dict(item), []))
        else:
            records.append((index, None, [_problem(
                None, "NOT_AN_OBJECT",
                f"element {index} is a {type(item).__name__}, not an "
                f"object")]))
    return records


def load_bom(source, fmt: Optional[str] = None,
             resolve_names: bool = True) -> List[Dict[str, Any]]:
    """
    Read a bill of materials into part dicts for process_batch().

    Every row is read and checked against schema.PART_SCHEMA; a row with
    problems is not dropped and does not stop the load.  It is returned
    like any other, carrying ``bom_errors``: every problem found in it,
    each a dict with key, code and reason.  process_part() refuses such a
    part, so process_batch() returns an error result for each bad row
    (with its row and problems) beside the sized results of the good
    rows.

    Parameters
    ----------
    source : str, os.PathLike, or an open file
        A CSV file (one column per part key, one row per part; an empty
        cell is a key not given) or a JSON file (an array of part
        objects, or an object with a "parts" array and an optional
        "schema_version").  UTF-8, with or without a byte-order mark.
    fmt : "csv", "json" or None
        None reads the format from the file name's extension.
    resolve_names : bool
        True (the default) reads the material and convection-zone names
        through resolve_material() and resolve_zone(), so a part gets the
        canonical key (a name that does not resolve is a problem of its
        row), and the component class with the same normalisation (case,
        spaces, hyphens).  False takes every name as written.

    Returns
    -------
    list of dict — one part dict per data row, in file order.  Each
    carries ``bom_row``: the row's line in a CSV file (the header is line
    1; a row whose quoted value holds a line break gets the line it ends
    on), or the element's 1-based position in a JSON array.  A row with
    problems also carries ``bom_errors``.  Rows of empty cells are skipped.

    Problem codes a row can carry: those of validate_part_input(), and
    DUPLICATE_PART_ID (a part_id an earlier row uses), MALFORMED_ROW (more
    CSV values than columns), NOT_AN_OBJECT (a JSON element that is not
    an object), RESERVED_KEY (bom_row or bom_errors written in the file).

    Raises
    ------
    ValueError
        For a file that is not a BOM at all: not UTF-8, not valid JSON, a
        JSON value of the wrong shape, a CSV header naming a column
        twice, an unknown format.
    """
    text, name = _read_text(source)
    fmt = _bom_format(fmt, name)
    records = (_csv_records(text, name) if fmt == "csv"
               else _json_records(text, name))

    parts = []
    first_row: Dict[str, int] = {}
    for row, values, problems in records:
        if values is None:                       # not an object
            parts.append({"bom_row": row, "bom_errors": problems})
            continue
        part = dict(values)
        for key in _LOADER_KEYS:
            if key in part:
                del part[key]
                problems.append(_problem(
                    key, "RESERVED_KEY",
                    f"{key!r} is set by load_bom() and may not appear in "
                    f"a BOM file"))
        part["bom_row"] = row

        named = set()
        if resolve_names:
            cls = part.get("component_class")
            if isinstance(cls, str) and _normalise(cls) in CLASS_DEFAULTS:
                part["component_class"] = _normalise(cls)
            lookups = [(key, resolve_zone) for key in _ZONE_KEYS]
            if part.get("component_class") in SOLID_CLASSES:
                lookups.append(("material", resolve_material))
            for key, resolver in lookups:
                if not isinstance(part.get(key), str):
                    continue
                try:
                    part[key] = resolver(part[key])
                except PartInputError as exc:
                    named.add(key)
                    problems.append(_problem(key, "UNKNOWN_NAME", str(exc)))

        problems += [p for p in validate_part_input(part)
                     if p["key"] not in named]
        pid = part.get("part_id")
        if isinstance(pid, str) and pid:
            if pid in first_row:
                problems.append(_problem(
                    "part_id", "DUPLICATE_PART_ID",
                    f"part_id {pid!r} is already used by row "
                    f"{first_row[pid]}"))
            else:
                first_row[pid] = row
        if problems:
            part["bom_errors"] = problems
        parts.append(part)
    return parts


__all__ = [
    "MATERIAL_ALIASES",
    "ZONE_ALIASES",
    "TWO_LAYER_SHIELD_KEYS",
    "SHIELD_KEYS",
    "EXHAUST_ZONES",
    "EXHAUST_ADJACENT_ZONES",
    "resolve_material",
    "resolve_zone",
    "infer_component_class",
    "load_bom",
]
