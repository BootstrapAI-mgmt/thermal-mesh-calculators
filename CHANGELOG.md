# Changelog

All notable changes to **thermal-mesh-calculators** are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and
this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Versions prior to `v0.6.0` were developed in-repo but never tagged or published;
their dates below are the dates of the commit that completed the version. `v0.6.0`
is the first release with a git tag, packaging metadata, and a buildable wheel.

## [0.6.2.post1] — 2026-09-25

Re-publication of 0.6.2 to PyPI. No library code changed: this tree is the `v0.6.2` tree plus
this entry, the version string and the README's install lines.

- **Changed:** the index files of 0.6.1 and 0.6.2 were removed from PyPI on 2026-09-25, and a removed
  filename cannot be re-uploaded, so the same code ships again as the post-release `0.6.2.post1`
  (`pip install thermal-mesh-calculators` resolves it; a `<0.7` pin still admits it). The README's
  git-pin example now names `v0.6.2.post1`, the first tag of this repository after `v0.6.2`.

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

[0.6.2.post1]: https://github.com/BootstrapAI-mgmt/thermal-mesh-calculators/releases/tag/v0.6.2.post1
[0.6.0]: https://github.com/BootstrapAI-mgmt/thermal-mesh-calculators/releases/tag/v0.6.0
