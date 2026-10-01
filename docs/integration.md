# Integrating thermal-mesh-calculators

For code that builds part dicts from data of its own (a CAD or CFD
preprocessing step, a PLM export, a spreadsheet) and reads the results back.
It covers what a part and a project may say, how a part's component class is
decided, how to size a part from your own material properties or a fluid
region, how names and BOM files are read, what a result promises, and how the
contract changes across releases.

The test suite checks this page against the package: the two key tables below
against `PART_SCHEMA` and `PROJECT_SCHEMA` (`tests/test_result_schema.py`), and
every example as a doctest. The functions are in
[`api_reference.md`](api_reference.md); the schemas are data in
`thermal_mesh_calculators/schema.py`.

<!-- Each pycon block ends with a blank line: without it, doctest reads the
     closing fence as part of the expected output. -->

```pycon
>>> import io
>>> from thermal_mesh_calculators import (
...     MATERIALS, PartInputError, infer_component_class, load_bom,
...     process_batch, process_part, process_part_from_props,
...     resolve_material, resolve_zone, validate_part_input, validate_result,
... )
>>> project = {"t_fluid_K": 353.15, "t_surr_K": 353.15, "t_exh_K": 1073.15}

```

## 1. The intake contract

A part dict and a project dict are checked against a versioned schema,
`INPUT_SCHEMA_VERSION` (1). A dict may name the version it was written for in
a `schema_version` key; absent, it is 1.

`validate_part_input(part)` and `validate_project_input(project)` return every
problem they find, each a dict with `key`, `code` and `reason`.
`process_part()` runs the same checks before it computes anything and raises
`PartInputError` listing the same problems in the same words (its `problems`
attribute; project keys are prefixed `project.`), so a part that validates is
a part `process_part()` accepts. A value's problem is worded as the
calculators' input guards word it.

```pycon
>>> bad = {"part_id": "BRK-009", "material": "steel_mild",
...        "component_class": "structural", "convection_zone": "engine_beside",
...        "thickness_mm": -3.0, "t_surf_K": 473.15, "surface": "glossy"}
>>> for problem in validate_part_input(bad):
...     print(problem["key"], problem["code"])
thickness_mm OUT_OF_RANGE
surface UNKNOWN_NAME
>>> validate_part_input(bad)[0]["reason"]
'thickness_mm must be > 0 mm, got -3.0'

```

Three rules decide what is checked:

- A key is checked only for the component classes that read it: a shield
  solves its temperature, so a shield's `t_surf_K` is not read and not checked.
  The convection-zone and surface-treatment names are the exception: they are
  checked whatever the class, because an unknown name is a misspelling wherever
  it is written.
- A key no schema declares is refused (`UNKNOWN_KEY`), so a misspelt key, say
  `h_overide`, cannot be ignored without a word. A key that starts with `x_` is
  yours: it is carried through and never read.
- A part whose `component_class` is missing or unknown reports that problem
  alone, since the class decides which keys are required.

The problem codes (`PROBLEM_CODES`):

| Code | Meaning |
|---|---|
| `NOT_AN_OBJECT` | The part, project or record is not a dict (a JSON element that is not an object, say). |
| `UNKNOWN_KEY` | A key the schema does not declare (and that does not start with `x_`). |
| `RESERVED_KEY` | A key only `load_bom()` sets (`bom_row`, `bom_errors`), found in a BOM file. |
| `MISSING_KEY` | A key the part's class requires is absent. |
| `WRONG_TYPE` | A value of the wrong type: a string where a number is expected, a bool, a float where an integer is expected. |
| `OUT_OF_RANGE` | A number outside its range, or not finite. |
| `UNKNOWN_NAME` | A name that is not in its vocabulary (a material, zone, surface treatment, class or option). |
| `UNSUPPORTED_VERSION` | A `schema_version` this library does not read. |
| `DUPLICATE_PART_ID` | A `part_id` an earlier BOM row already uses. |
| `MALFORMED_ROW` | A CSV row with more values than the header has columns. |
| `INCONSISTENT` | Result values that contradict each other (the governing size is not the smallest candidate, say). |
| `CALCULATION_FAILED` | A calculator raised on inputs the schema accepts; `process_batch()` reports it per part. |

### Part keys

"Solids" are `exhaust`, `exhaust_adjacent`, `structural`, `shield` and
`multilayer_shield`; "non-shield solids" the first three; "shields" the last
two. `schema_version`, `bom_row` and `bom_errors` are read by the input checks,
not by the sizing.

<!-- part-keys:begin -->
| Key | Type and range | Read by | Required for | Default |
|---|---|---|---|---|
| `schema_version` | integer | all |  | 1 |
| `part_id` | non-empty string | all | all |  |
| `component_class` | a key of `CLASS_DEFAULTS` | all | all |  |
| `material` | a key of `MATERIALS` | solids | solids |  |
| `thickness_mm` | number > 0 mm | solids | solids |  |
| `t_surf_K` | number > 0 K | non-shield solids, fluid | non-shield solids |  |
| `convection_zone` | a key of `CONVECTION_ZONES` | all (name checked for every class) | non-shield solids unless `h_override` given; fluid unless `velocity_ms` and `bl_regime` given |  |
| `convection_zone_in` | a key of `CONVECTION_ZONES` | shields (name checked for every class) | shields unless `convection_zone` or `h_in_override` given |  |
| `convection_zone_out` | a key of `CONVECTION_ZONES` | shields (name checked for every class) | shields unless `convection_zone` or `h_out_override` given |  |
| `surface` | a key of `SURFACE_TREATMENTS` | non-shield solids (name checked for every class) |  | the material-class emissivity (SURFACE_DEFAULTED) |
| `epsilon` | number in [0, 1] | non-shield solids |  |  |
| `surface_in` | a key of `SURFACE_TREATMENTS` | shields (name checked for every class) |  | the material-class emissivity (SURFACE_DEFAULTED) |
| `eps_in` | number in [0, 1] | shields |  |  |
| `surface_out` | a key of `SURFACE_TREATMENTS` | shields (name checked for every class) |  | the material-class emissivity (SURFACE_DEFAULTED) |
| `eps_out` | number in [0, 1] | shields |  |  |
| `surface_g1` | a key of `SURFACE_TREATMENTS` | multilayer_shield (name checked for every class) |  | the material-class emissivity (SURFACE_DEFAULTED) |
| `eps_g1` | number in [0, 1] | multilayer_shield |  |  |
| `surface_g2` | a key of `SURFACE_TREATMENTS` | multilayer_shield (name checked for every class) |  | the material-class emissivity (SURFACE_DEFAULTED) |
| `eps_g2` | number in [0, 1] | multilayer_shield |  |  |
| `h_override` | number >= 0 W/m^2 K | non-shield solids |  |  |
| `h_in_override` | number >= 0 W/m^2 K | shields |  |  |
| `h_out_override` | number >= 0 W/m^2 K | shields |  |  |
| `t_fluid_K` | number > 0 K | all |  | the project's t_fluid_K |
| `char_length_mm` | number > 0 mm | all |  | 100 |
| `t_exh_K` | number > 0 K | shields |  | the project's t_exh_K, else 1073.15 (SHIELD_INPUT_DEFAULTED) |
| `h_gap` | number >= 0 W/m^2 K | multilayer_shield |  | 15 (SHIELD_INPUT_DEFAULTED) |
| `f12` | number in [0, 1] | multilayer_shield |  | 1.0, parallel plates (SHIELD_INPUT_DEFAULTED) |
| `shield_max_iter` | integer >= 1, or `None` | shields |  | the project's, else 50 single-layer, 100 two-layer |
| `max_dt` | number > 0 K, or `None` | solids |  | the project's, else the class default |
| `allowable_flux_error` | number > 0 W/m^2, or `None` | solids |  | the project's, else the class default |
| `radius_mm` | number > 0 mm, or `None` | solids |  | no curvature limit |
| `velocity_ms` | number >= 0 m/s, or `None` | all |  | the zone's velocity |
| `bl_regime` | `"external_forced"` or `"mixed_unknown"`, or `None` | all |  | the project's; a solid without one gets no boundary-layer constraint; a fluid region takes its zone's (forced: external_forced, else mixed_unknown) |
| `bl_y_plus` | number > 0 | all |  | the project's, else 30 |
| `bl_growth_ratio` | number > 1 | all |  | the project's, else 1.2 |
| `bl_ar_max_prism` | number >= 1 | all |  | the project's, else 5 |
| `bl_fraction` | number in (0, 1] | all |  | the project's, else 0.3 |
| `bom_row` | integer >= 1 | all |  |  |
| `bom_errors` | list | all |  |  |
<!-- part-keys:end -->

`None` is accepted only where the table says so, and means "use the default".
A key starting with `x_` is accepted on any part and never read.

### Project keys

Every key a project gives is checked, whatever the parts: a project value
applies to every part.

<!-- project-keys:begin -->
| Key | Type and range | Required | Default |
|---|---|---|---|
| `schema_version` | integer |  | 1 |
| `t_fluid_K` | number > 0 K | yes |  |
| `t_surr_K` | number > 0 K | yes |  |
| `max_dt` | number > 0 K, or `None` |  | the class default |
| `allowable_flux_error` | number > 0 W/m^2, or `None` |  | the class default |
| `dt` | number > 0 s, or `None` |  | no transient constraint |
| `fo_max` | number > 0 |  | 0.5 |
| `safety_factor` | number > 0 |  | 1.0 |
| `tau_bc` | number > 0 s, or `None` |  | no drive-cycle limit |
| `transient_scheme` | `"explicit"` or `"implicit"`, or `None` |  | inferred from fo_max (<= 0.5 explicit) |
| `t_exh_K` | number > 0 K |  | 1073.15 |
| `shield_max_iter` | integer >= 1, or `None` |  | 50 single-layer, 100 two-layer |
| `bl_regime` | `"external_forced"` or `"mixed_unknown"`, or `None` |  | no boundary-layer constraint on solids |
| `bl_y_plus` | number > 0 |  | 30 |
| `bl_growth_ratio` | number > 1 |  | 1.2 |
| `bl_ar_max_prism` | number >= 1 |  | 5 |
| `bl_fraction` | number in (0, 1] |  | 0.3 |
<!-- project-keys:end -->

## 2. Component classes

`component_class` decides which sizing path runs and which defaults apply. The
classes are the keys of `CLASS_DEFAULTS`:

| Class | What it is | Path | Constraints that can govern | `max_dt`, `allowable_flux_error` defaults | Radiation gradient |
|---|---|---|---|---|---|
| `exhaust` | carries the exhaust gas (manifold, pipe) | solid, surface temperature given | conduction, lateral_gradient, biot_through_thickness, radiation, transient, curvature, aero_boundary_layer | 10 K, 500 W/m^2 | 167 K/m, the exhaust worst case |
| `exhaust_adjacent` | sits beside, below or above the exhaust | as `exhaust` | as `exhaust` | 10 K, 500 W/m^2 | 83 K/m |
| `structural` | anything else with a known surface temperature | as `exhaust` | as `exhaust` | 15 K, 200 W/m^2 | boundary flux / k |
| `shield` | a single-layer heat shield | its temperature is solved | shield, biot_through_thickness, radiation, transient, curvature, aero_boundary_layer | 15 K, 500 W/m^2 | exhaust-facing flux / k |
| `multilayer_shield` | a two-layer shield with an air gap | both layer temperatures are solved | shield_layer1, shield_layer2, and as `shield` | 15 K, 500 W/m^2 | layer 1's flux / k |
| `fluid` | a fluid region | the boundary layer alone (section 4) | aero_boundary_layer | none | none |

`transient` applies when the project gives `dt`, `curvature` when the part
gives `radius_mm`, and `aero_boundary_layer` on a solid when the part or the
project gives `bl_regime`.

### How a part's class is inferred

`process_part()` never guesses a class: a part names one. For data that does
not carry one, `infer_component_class(part)` returns the class the part's keys
imply. The first rule that applies decides:

| Rule | Keys that decide | Class |
|---|---|---|
| 1 | `component_class` given | that class (an unknown one raises `PartInputError`) |
| 2 | any of `surface_g1`, `surface_g2`, `eps_g1`, `eps_g2`, `h_gap`, `f12` | `multilayer_shield` |
| 3 | any of `surface_in`, `surface_out`, `eps_in`, `eps_out`, `h_in_override`, `h_out_override`, `convection_zone_in`, `convection_zone_out` | `shield` |
| 4 | `convection_zone` is `exhaust_internal` | `exhaust` |
| 5 | `convection_zone` is `exhaust_beside`, `exhaust_below`, `exhaust_above` or `near_exhaust_natural` | `exhaust_adjacent` |
| 6 | none of the above | `structural` |

Rules 2 and 3 read keys only a shield reads; rules 4 and 5 read the zone, so an
exhaust component sitting in a zone beside the exhaust (a manifold's outer
face, say) is inferred `exhaust_adjacent`: a part that is an exhaust component
says so. The zone is read through `resolve_zone()` when it resolves. `fluid` is
never inferred, since no part key says a region is a fluid: a fluid region
names its class.

The worked examples, one per class:

A bracket in the engine bay: no shield key and no exhaust zone (rule 6).

```pycon
>>> bracket = {"part_id": "BRK-001", "material": "steel_mild",
...            "convection_zone": "engine_beside", "thickness_mm": 3.0,
...            "t_surf_K": 473.15, "surface": "painted"}
>>> bracket["component_class"] = infer_component_class(bracket)
>>> r = process_part(bracket, project)
>>> print(r["component_class"], r["governing_constraint"], f"{r['governing_dx_mm']:.2f}")
structural lateral_gradient 24.50

```

An exhaust manifold: its zone is the exhaust gas itself (rule 4). Its own
`t_fluid_K` is the gas temperature.

```pycon
>>> manifold = {"part_id": "MAN-001", "material": "cast_iron",
...             "convection_zone": "exhaust_internal", "thickness_mm": 6.0,
...             "t_surf_K": 1073.15, "t_fluid_K": 1023.15,
...             "surface": "cast_iron_oxidised"}
>>> manifold["component_class"] = infer_component_class(manifold)
>>> r = process_part(manifold, project)
>>> print(r["component_class"], r["governing_constraint"], f"{r['governing_dx_mm']:.2f}")
exhaust conduction 8.28

```

A plastic housing above the exhaust (rule 5).

```pycon
>>> housing = {"part_id": "HSG-001", "material": "plastic_pa66_gf30",
...            "convection_zone": "exhaust_above", "thickness_mm": 2.5,
...            "t_surf_K": 400.0, "surface": "plastic"}
>>> housing["component_class"] = infer_component_class(housing)
>>> r = process_part(housing, project)
>>> print(r["component_class"], r["governing_constraint"], f"{r['governing_dx_mm']:.2f}")
exhaust_adjacent lateral_gradient 1.90

```

A single-layer shield: it gives each side's zone and surface (rule 3), and no
surface temperature, which is solved.

```pycon
>>> shield = {"part_id": "SH-001", "material": "steel_stainless_409",
...           "convection_zone_in": "exhaust_beside",
...           "convection_zone_out": "engine_beside", "thickness_mm": 0.8,
...           "surface_in": "aluminised", "surface_out": "aluminised"}
>>> shield["component_class"] = infer_component_class(shield)
>>> r = process_part(shield, project)
>>> print(r["component_class"], r["governing_constraint"], f"{r['governing_dx_mm']:.2f}")
shield radiation 9.94
>>> print(f"{r['shield']['t_shield_K']:.1f} K")
769.1 K

```

A two-layer shield: gap faces and a gap conductance (rule 2). Its gap view
factor `f12` is not given, so the result says it took the default.

```pycon
>>> double = {"part_id": "SH-002", "material": "steel_stainless_409",
...           "convection_zone_in": "exhaust_beside",
...           "convection_zone_out": "shield_gap_confined", "thickness_mm": 0.6,
...           "surface_in": "aluminised", "surface_out": "aluminised",
...           "surface_g1": "aluminised", "surface_g2": "aluminised",
...           "h_gap": 12.0}
>>> double["component_class"] = infer_component_class(double)
>>> r = process_part(double, project)
>>> print(r["component_class"], r["governing_constraint"], f"{r['governing_dx_mm']:.2f}")
multilayer_shield radiation 7.85
>>> [w["code"] for w in r["warnings"]]
['SHIELD_INPUT_DEFAULTED']

```

The air at the front of the vehicle: a fluid region names its class (rule 1).

```pycon
>>> air = {"part_id": "AIR-001", "component_class": "fluid",
...        "convection_zone": "front_end_edges"}
>>> infer_component_class(air)
'fluid'
>>> r = process_part(air, project)
>>> print(r["component_class"], r["governing_constraint"], f"{r['governing_dx_mm']:.2f}")
fluid aero_boundary_layer 0.99

```

## 3. Sizing from your own material properties

`process_part_from_props(part, project, k=..., rho=..., cp=...)` sizes a part
whose material properties you hold yourself. The keyword arguments are the
fields of a material record (section 8), so an entry of `MATERIALS` can be
passed whole, and gives what the table path gives:

```pycon
>>> process_part_from_props(bracket, project, **MATERIALS["steel_mild"]) == process_part(bracket, project)
True

```

The part's `material` is then a label of your own: optional, never looked up,
echoed into the result.

```pycon
>>> own = dict(bracket, material="house steel S-7")
>>> r = process_part_from_props(own, project, k=48.0, rho=7850.0, cp=470.0)
>>> print(r["material"], r["governing_constraint"], f"{r['governing_dx_mm']:.2f}")
house steel S-7 lateral_gradient 23.10

```

Two rules read the material's name, so they read the label's spelling: the
fallback emissivity of a surface given neither a treatment nor an emissivity
(0.90 when the name contains plastic, rubber, composite, ceramic or glass; 0.30
when it contains aluminium or aluminum; 0.73 otherwise, reported as
`SURFACE_DEFAULTED`), and the non-metal Biot warning (`PLASTIC_BIOT_MARGINAL`,
when the name contains plastic, rubber or composite). Give every surface a
treatment or an emissivity, and a label that carries the family word, and
neither rule depends on how you name your materials.

## 4. Fluid regions

A part with `component_class` `"fluid"` is a fluid region. It gives
`part_id`, the class, and its flow: `convection_zone`, or both `velocity_ms`
and `bl_regime`. It is sized by the boundary-layer constraint alone; no
conduction, lateral, Biot, radiation, shield, transient or curvature
constraint runs, and it takes no material (keys only a solid reads are not
read). The velocity is the part's `velocity_ms`, else the zone's; the regime is
the part's or the project's `bl_regime`, else the zone's (a forced zone is
`external_forced`, any other `mixed_unknown`). The wall's `t_surf_K`, when
given, sets the film temperature and the buoyancy velocity.

```pycon
>>> r = process_part(air, project)
>>> bl = r["boundary_layer"]
>>> print(bl["regime"], f"{bl['y1_mm']:.2f} mm", bl["n_layers"], r["conduction"])
external_forced 2.81 mm 3 None

```

## 5. Material and zone names

`resolve_material(name)` and `resolve_zone(name)` return the canonical key a
common spelling refers to, from the alias tables `MATERIAL_ALIASES` and
`ZONE_ALIASES`. Names are compared after normalisation: lower case, every run
of characters other than letters and digits made one underscore, and the
spellings aluminum, galvanized, gray, fiber, fiberglass and molding read as
aluminium, galvanised, grey, fibre, fibreglass and moulding. A canonical key
resolves to itself.

```pycon
>>> resolve_material("PA66-GF30"), resolve_material("Aluminum 6061-T6")
('plastic_pa66_gf30', 'aluminium_6061')
>>> resolve_zone("Beside engine")
'engine_beside'

```

An alias names the same material or zone as its key: a synonym, a designation
the entry's description lists ("Austenitic stainless 304 / 316" covers SS 316),
or a word order. A family name says which entries it could mean, and nothing
more, so it is refused with them:

```pycon
>>> try:
...     resolve_material("aluminium")
... except PartInputError as error:
...     print(str(error).split(" Available:")[0])
Unknown material 'aluminium'. 'aluminium' could be any of: aluminium_5052, aluminium_6061, cast_aluminium, cast_aluminium_a380; give the one you mean.

```

`process_part()` reads only canonical keys; `load_bom()` resolves names on the
way in (section 6).

## 6. Reading a BOM file

`load_bom(path_or_file, fmt=None)` reads a CSV file (one column per part key,
one row per part; an empty cell is a key not given) or a JSON file (an array of
part objects, or an object with a `parts` array and an optional
`schema_version`) into part dicts for `process_batch()`. It reads every row,
checks it against the part schema and resolves the material and zone names;
a bad row is neither dropped nor the end of the load. Every part carries
`bom_row` (its line in a CSV file, the header being line 1, or the line it ends
on when a quoted value holds a line break; its position in a JSON array), and a
bad row carries `bom_errors`, every problem it has.

```pycon
>>> text = (
...     "part_id,material,component_class,convection_zone,thickness_mm,t_surf_K,surface\n"
...     "BRK-001,Mild steel,structural,engine_beside,3.0,473.15,painted\n"
...     "BRK-002,steel_mild,structural,engine_beside,-2,473.15,painted\n"
...     "PNL-001,unobtanium,structural,cabin_tunnel,2.0,hot,painted\n"
... )
>>> parts = load_bom(io.StringIO(text), fmt="csv")
>>> for part in parts:
...     print(part["bom_row"], part["part_id"],
...           [(e["key"], e["code"]) for e in part.get("bom_errors", [])])
2 BRK-001 []
3 BRK-002 [('thickness_mm', 'OUT_OF_RANGE')]
4 PNL-001 [('material', 'UNKNOWN_NAME'), ('t_surf_K', 'WRONG_TYPE')]

```

`process_part()` refuses a part that carries `bom_errors`, so
`process_batch()` returns an error result for each bad row, with its row and
problems, beside the sized good rows:

```pycon
>>> for r in process_batch(parts, project):
...     print(r["bom_row"], r["part_id"],
...           r["error"].split("\n")[0] if r["error"] else f"{r['governing_dx_mm']:.2f} mm")
2 BRK-001 24.50 mm
3 BRK-002 BOM row 3 did not load: thickness_mm must be > 0 mm, got -2.0
4 PNL-001 BOM row 4 did not load: 2 problems with part 'PNL-001':

```

A file that is not a BOM at all (not UTF-8, not valid JSON, a CSV header
naming a column twice) raises `ValueError`. A column for data of your own
starts with `x_`.

## 7. Results

Every result carries `schema_version`, the result schema version it conforms
to (`RESULT_SCHEMA_VERSION`, 1), and `validate_result(result)` checks one
against `RESULT_SCHEMA`: every top-level key declared, every promised key of
the calculator sub-results present (for example `biot["mesh_type"]` and
`boundary_layer["max_dx_surface_mm"]`), every warning code registered in
`WARNING_CODES` with one of its severities, and the governing size the
smallest candidate in `all_constraints`. The test suite runs it on every
result it produces.

```pycon
>>> validate_result(process_part(bracket, project))
[]

```

A part `process_batch()` could not size gets an error result instead
(`ERROR_RESULT_SCHEMA`): `part_id`, `error` (the message), `problems` (each
with `key`, `code` and `reason`), `governing_dx_mm` and
`governing_constraint` set to `None`, `bom_row` when the part has one, and
`schema_version`. The warning codes and their severities are listed in
[`quick_reference.md`](quick_reference.md#warning-codes).

## 8. Material records

A material record (an entry of `MATERIALS`, or the keyword arguments of
`process_part_from_props()`) is checked by `validate_material_record()`
against `MATERIAL_RECORD_SCHEMA`:

| Field | Type | Required | Meaning |
|---|---|---|---|
| `k` | number > 0 W/m K | yes | thermal conductivity |
| `rho` | number > 0 kg/m^3 | yes | density |
| `cp` | number > 0 J/kg K | yes | specific heat |
| `description` | string | no | what the material is |
| `t_service_max_K` | number > 0 K | no | the highest temperature the material may see in service |
| `source` | non-empty string | no | where the record's values come from |

A part sized above a limit its record gives gets a `SERVICE_TEMP_EXCEEDED`
warning: for a solid the surface temperature is compared, for a shield the
solved temperature (the exhaust-facing layer of two).

```pycon
>>> r = process_part_from_props(own, project, k=48.0, rho=7850.0, cp=470.0,
...                             t_service_max_K=450.0)
>>> [w["code"] for w in r["warnings"]]
['SERVICE_TEMP_EXCEEDED']

```

No entry of the built-in table gives a service limit or a source yet. Those
values, and the reuse terms of their sources, are the maintainer's to supply;
until then the warning comes only from records you pass yourself.

## 9. Versions and deprecation

The input schema and the result schema are versioned separately
(`INPUT_SCHEMA_VERSION`, `RESULT_SCHEMA_VERSION`), and each changes under these
rules:

1. A patch release (X.Y.1 after X.Y.0) changes neither schema.
2. A minor release may add an optional input key, a result key, a warning code
   or a problem code without changing the schema version. A reader that ignores
   result keys it does not know is unaffected; a part written with a new key is
   refused by an older release (`UNKNOWN_KEY`).
3. A key is renamed or removed in two steps. First the replacement is added and
   the old key is marked deprecated in the changelog's `Deprecated` section;
   both keep working for at least one minor release. Then a later minor release
   removes the old key, increases that schema's version and lists the removal
   under `Removed` with "(behaviour change)".
4. A key never changes its unit or its meaning: a changed meaning is a new key.

What a pin may rely on: within one minor series (what a pin such as
`>=X.Y,<X.(Y+1)` admits) a result carries the keys that series documents,
with the same `schema_version`. A
reader that checks `result["schema_version"]` against the version it was
written for, and refuses another, cannot misread a result whose contract
changed.
