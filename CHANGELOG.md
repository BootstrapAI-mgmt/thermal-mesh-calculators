# Changelog

All notable changes to **thermal-mesh-calculators** are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and
this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Versions prior to `v0.6.0` were developed in-repo but never tagged or published;
their dates below are the dates of the commit that completed the version. `v0.6.0`
is the first release with a git tag, packaging metadata, and a buildable wheel.
Neither `v0.6.0` nor `v0.6.1` is tagged in this repository: its public history
starts at 0.6.2, and `v0.6.2` is its first tag.

## [Unreleased]

Input validation: inputs the library used to accept without a word now raise
a named error or come back with a coded warning. Entries marked **(behaviour
change)** can alter a result, or raise where 0.6.2 returned one.

### Added
- `PartInputError`, exported from the package. `process_part()` raises it
  before computing anything when a part names an unknown or missing
  `material` or `component_class`, an unknown convection zone (any of
  `convection_zone`, `convection_zone_in`, `convection_zone_out`), a zone it
  needs but does not give, or an unknown surface treatment (any of `surface`,
  `surface_in`, `surface_out`, `surface_g1`, `surface_g2`). The message names
  the part, the key and the allowed set. It subclasses both `KeyError` (what
  the old lookups raised) and `ValueError`.
- `SURFACE_DEFAULTED` warning: a surface given neither a treatment nor an
  emissivity still gets the material-class fallback (0.90 non-metal, 0.30
  aluminium, 0.73 other metals), and the result now says so, one warning per
  defaulted surface key. Input warnings carry a `key` field.
- A part-level `t_fluid_K`, for a part whose fluid is not the project's air
  (exhaust gas inside a pipe, say). Results report the fluid temperature used
  (`t_fluid_K`) and where it came from (`t_fluid_source`: `"part"` or
  `"project"`).
- `SHIELD_NOT_CONVERGED` warning, when a shield solve stops at its iteration
  limit, and a `shield_max_iter` part or project key that sets the limit.
  `SingleLayerShieldCalculator.mesh_size()` now takes `tol` and `max_iter`
  and forwards them to the solve.
- `NO_FINITE_SIZE` warning, when every candidate size of a part is unbounded
  (a shield with no heat source, for example), so an infinite governing size
  is never silent.
- `residual_W_m2` in both shield solvers' results: the energy-balance
  residual at the returned temperature(s).
- `TransientMeshCalculator.combined_transient_limits()` takes a `scheme`
  (`"explicit"` or `"implicit"`, inferred from `fo_max` when omitted) and
  reports `scheme`, `penetration_applied`, `max_dx_mm`, `feasible`, `conflict`
  (with `remedies`) and `advice`. `TransientMeshCalculator.resolve_scheme()`
  is the inference rule. `process_part()` reads a `transient_scheme` project
  key and adds a `TRANSIENT_CONFLICT` warning (a `warning` for an explicit
  scheme, a `caution` for an implicit one) when the smallest element the time
  step allows exceeds the governing size, with a remedy that closes it.
- Range guards in the calculators. A `ValueError` names the parameter, its
  range and the value given for: conductivity, density and specific heat not
  > 0; any temperature not > 0 K; an emissivity or the gap view factor `f12`
  outside [0, 1]; a negative convection coefficient (`h`, `h_in`, `h_out`,
  `h_gap`); `max_dt`, `dt`, `fo_max`, `safety_factor`, `tau_bc` or
  `allowable_flux_error` not > 0; a negative `spatial_gradient` or velocity;
  and a correlation `char_length` not > 0. A value that is not a number
  (including a string or a bool) raises `TypeError`.
- `h_estimator.AIR_PROPERTY_RANGE_K` (250–700 K) and `film_in_range()`.
  `air_properties()` reports `t_film_K`, `t_eval_K` (the clamped temperature
  the fits were evaluated at) and `in_range`. `estimate_h()` and
  `BoundaryLayerCalculator.estimate_mesh()` report `t_film_K` and
  `film_in_range`.
- `FILM_TEMP_OUT_OF_RANGE` warning, when the h estimate or the boundary-layer
  sizes rest on air properties clamped to that range (0.6.2 returned
  k_air(1000 K) = k_air(700 K) and said nothing).
- `H_CORRELATION_FALLBACK` warning, when the h correlation raises and the
  zone's static h is used instead. The reason is recorded as
  `fallback_reason`, so an error inside the correlation is never swallowed.
- `f12`, the view factor between the gap faces of a two-layer shield, as a
  part key; `process_part()` passes it to the solver and reports it as
  `f12_used`. `MultilayerShieldCalculator.mesh_sizes()` reports the layer
  driving fluxes `q_layer1` and `q_layer2`.
- `SHIELD_INPUT_DEFAULTED` warning, naming the key, when a shield input falls
  back to its default: `t_exh_K` (neither part nor project gives it; 1073.15
  K), `h_gap` (15 W/m² K) or `f12` (1.0, parallel plates). A defaulted gap
  surface is reported as `SURFACE_DEFAULTED`.

### Changed
- **(behaviour change)** An unknown `component_class` raises. 0.6.2 sized it
  as `structural` without a warning (a part spelt `"exhuast"` came back sized
  like a bracket).
- **(behaviour change)** An unknown convection zone raises even when
  `h_override` is given; 0.6.2 never looked the zone up on that path.
- **(behaviour change)** One fluid temperature per part for both the
  convective coefficient and the boundary flux. 0.6.2 estimated h at the
  zone's upper air temperature (`t_air_C_high`) while computing q″ from the
  project's `t_fluid_K`, so changing `t_fluid_K` moved q″ but never h. Both now
  use the part's `t_fluid_K` if given, else the project's, and so do the
  shield solve, the boundary-layer film temperature and the automatic
  warnings. For a part in the `exhaust_internal` zone, give the gas
  temperature as the part's `t_fluid_K`.
- **(behaviour change)** The explicit transient path is feasible, and its
  advice is correct. With the explicit defaults (`fo_max` 0.5, safety factor
  1) the Fourier minimum is √2 times the penetration bound at every time
  step, so 0.6.2 reported a conflict whatever `dt` was and advised "reduce
  dt", which cannot change a ratio that does not depend on `dt`. An explicit
  scheme's update is monotone for every Fo ≤ ½, so the per-step penetration
  depth does not bound it: that bound now applies to implicit schemes only,
  and an explicit scheme's upper bounds are the drive-cycle limit and the
  steady constraints. A conflict names a remedy that closes it: `dt ≤ fo_max
  · tau_bc` against the drive cycle, or a larger `safety_factor` or `fo_max`
  when an implicit window is empty at every `dt`. The message prints each
  bound rounded toward the side that satisfies it (a maximum down, a minimum
  up), so the value as printed closes the conflict; `remedies` carries the
  exact values. The derivation is in
  `docs/math_derivations.md` §7.4, and example 10 prints the window and a
  working remedy.
- **(behaviour change)** `process_part()` uses the transient upper bound as a
  size candidate instead of `recommended_dx_mm`. 0.6.2 capped the governing
  size at the Fourier stability minimum, a lower bound, and never said so.
- **(behaviour change)** Shield classes get a radiation constraint. 0.6.2
  gave them none (`radiation` was always `None`). The constraint is evaluated
  at the exhaust-facing surface (the shield temperature, or layer 1 for two
  layers) with its emissivity and the conduction-driven gradient that
  `estimate_spatial_gradient()` gives shields: the exhaust-facing driving
  flux divided by k. It can govern a shield's size.
- **(behaviour change)** Out-of-range inputs raise instead of returning a
  number. With k = −45 and ε = 1.7, `max_mesh_size()` returned −1.920 mm.
  `lateral_gradient_limit()` with k ≤ 0 now raises where it returned `inf`.

### Fixed
- A gap emissivity (`eps_g1`, `eps_g2`) or view factor (`f12`) of zero raised
  `ZeroDivisionError` in `MultilayerShieldCalculator.solve_temperatures()`. It
  is now the limit it tends to: no radiative exchange across the gap
  (`eps_eff` = 0).
- **(behaviour change)** `MultilayerShieldCalculator.solve_temperatures()`
  returned no fluxes when it ran out of iterations, so `mesh_sizes()` sized
  both layers as `inf` and batch runs carried that on without a warning. It
  now returns the last iterate's fluxes and temperatures with
  `converged: False`, as the single-layer solver already did, so the sizes are
  finite estimates.
- A shield given both `h_in_override` and `h_out_override` no longer needs a
  convection zone: 0.6.2 looked the zones up before reading the overrides and
  raised `KeyError` when they were absent.

## [0.6.2] — 2026-09-17

Packaging and publication release. No library code changed.

- **Fixed:** the published 0.6.1 sdist shipped `tests/test_open_closed_map.py` without the
  `scripts/` checker and `OPEN-CLOSED-MAP.yaml` it needs, so `pytest tests/` aborted at
  collection for a downstream packager. `MANIFEST.in` now excludes that test from the sdist;
  the sdist's own suite runs 221 tests green.
- **Changed:** the repository is published under the same name from a fresh public history
  (2026-09-17); earlier development history is not public. `CLAUDE.md` and `README.md`
  counts were re-measured (43 materials, 239 tests in-tree).

## [0.6.1] — 2026-09-16

Metadata-only correction, made before the first upload to PyPI. No library code
changed; `0.6.0` and `0.6.1` are functionally identical.

### Changed
- `requires-python` narrowed from `>=3.6` to `>=3.9`, and the `3.6` / `3.7` / `3.8`
  classifiers removed. Those interpreters are end-of-life and were never tested —
  CI exercises 3.9 and 3.11 — so the old floor advertised support the project did
  not have. A published version's metadata is immutable, so this was corrected
  before upload rather than after.

## [0.6.0] — 2026-08-12

First **packaged** release: the library gains distribution metadata and a tagged,
pinnable artifact. No calculator behaviour changed as part of the packaging work.

### Added
- `pyproject.toml` (PEP 621) declaring the distribution as `thermal-mesh-calculators`,
  built with setuptools. The version is read dynamically from
  `thermal_mesh_calculators.__version__`, so the module and the distribution can
  never disagree.
- `LICENSE` — Apache-2.0, copyright BootstrapAI-mgmt. The project previously had
  no stated licence, which made it legally undistributable.
- `NOTICE` — Apache-2.0 attribution notice.
- `MANIFEST.in` — reproducible source-distribution contents.
- `thermal_mesh_calculators/py.typed` — PEP 561 marker so type checkers honour the
  inline annotations in installed copies.
- This `CHANGELOG.md`.
- `BoundaryDrivenConductionCalculator.size_wall()` — thickness-aware sizing that
  reports through-thickness cell count and Biot number instead of a raw `dx` that
  may exceed the wall itself.

### Changed
- Zero-flux pole guard and a prescribed-`Ts` consistency check in the conduction
  calculator.

### Notes
- Runtime dependencies remain **empty by design** — the package is pure standard
  library so it can run inside CAE pre-processor embedded interpreters. Verified
  by an AST scan: the only non-self top-level import across the package is `math`.

## [0.5.1] — 2026-03-15

### Added
- Prism recommendation logic in `boundary_layer.py` (`prisms_recommended`), backed
  by a low-Re tet-sufficiency analysis: y₁-vs-dx and wall-gradient discriminators.
- Volume expansion ratio `ER_v` reporting for the external-forced regime.

## [0.5.0] — 2026-03-14

### Added
- `boundary_layer.py` — `BoundaryLayerCalculator`, a two-regime aerodynamic
  boundary-layer model (`external_forced` and `mixed_unknown`) producing inflation
  layer parameters and surface mesh constraints.
- The boundary-layer constraint wired into the batch processor as the
  `aero_boundary_layer` `dx` candidate.

## [0.4.0] — 2026-03-14

### Added
- Fin-theory lateral gradient constraint — `lateral_gradient_limit()` in
  `conduction.py`, wired into the batch processor.
- Gap view factor `f12` parameter on `MultilayerShieldCalculator.solve_temperatures()`
  (backward compatible; `F₁₂ = 1` reproduces the previous formula).
- Opposing mixed-convection handling in `estimate_h()`, with a conduction floor to
  prevent a singularity.

### Changed
- Cell Biot heuristic tightened from `ceil(10·Bi)` to `ceil(20·Bi)` so that each
  element satisfies `Bi_cell ≤ 0.1`.

## [0.3.0] — 2026-03-13

### Added
- `h_estimator.py` — convective HTC from correlations (forced, natural, mixed),
  Richardson-number regime classification, and the solver advisory system.
- `batch.py` — BOM → mesh sizes, with inline databases of 44 materials and
  30 surface treatments.
- `zones.py` — convection zone lookup with velocity, orientation and air-temperature
  metadata.
- Parametric regime-crossover studies under `analysis/`.

## [0.2.0] — 2026-03-13

### Added
- `transient.py` — penetration depth, Fourier number, and drive-cycle constraints,
  plus a combined evaluator that identifies the binding constraint.

### Changed
- Shield solver API aligned on separate `h_in` / `h_out` coefficients. The legacy
  `h_total` parameter is still accepted by the single-layer solver and split evenly.

## [0.1.0] — 2026-03-13

### Added
- Initial calculators: boundary-driven conduction, convection (Biot, h-gradient,
  CFD mapping), radiation (T⁴ sensitivity, view-factor curvature), and the
  single-layer and multilayer Newton-Raphson heat-shield solvers.

[0.6.2]: https://github.com/BootstrapAI-mgmt/thermal-mesh-calculators/releases/tag/v0.6.2
