# CLAUDE.md — Agent Context for thermal-mesh-calculators

## What This Project Is

Physics-driven thermal mesh sizing calculators for automotive CAE (Computer-Aided Engineering). These tools determine the maximum element size for finite element meshes based on heat transfer physics — not arbitrary rules of thumb.

The target users are thermal simulation engineers doing underhood analysis (exhaust systems, heat shields, structural components) in tools like ANSA, HyperMesh, Abaqus, NASTRAN, StarCCM+, etc.

## Repository Owner

**Org:** BootstrapAI-mgmt
**Visibility:** Public (published 2026-09-17; earlier development history is not public)

## Architecture

```
thermal_mesh_calculators/       # Python package (pure stdlib, no dependencies)
    __init__.py                 # Exports all calculator classes + h_estimator (v0.5.0)
    constants.py                # Stefan-Boltzmann constant
    _guards.py                  # Input range guards shared by the calculators (internal)
    conduction.py               # BoundaryDrivenConductionCalculator
    convection.py               # ConvectionMeshCalculator
    radiation.py                # RadiationMeshCalculator
    shields.py                  # SingleLayerShieldCalculator, MultilayerShieldCalculator
    transient.py                # TransientMeshCalculator (penetration depth, Fourier, drive-cycle)
    zones.py                    # Convection zone lookup with velocity/orientation/t_air metadata
    h_estimator.py              # Convective HTC from correlations (forced, natural, mixed, Richardson)
    boundary_layer.py           # BoundaryLayerCalculator (aero BL mesh: y+, inflation, surface dx)
    batch.py                    # BOM → mesh sizes (43 materials, 30 surface treatments, solver advisory)
examples/
    automotive_examples.py      # 11 worked scenarios with realistic inputs
    batch_example.py            # Batch processing demonstration
analysis/
    regime_crossover_study.py   # Parametric study: conduction vs convection dominance thresholds
    offroad_regime_study.py     # Off-road 10 mph study: laminar/turbulent/mixed regime mapping
    exhaust_underbody_scenario.py # Off-road underbody exhaust: pipe and cargo shields, cargo bed
tests/
    test_conduction.py          # 14 tests (incl. lateral gradient)
    test_convection.py          # 12 tests (incl. cell Biot 20x)
    test_radiation.py           # 10 tests
    test_shields.py             # 33 tests (incl. F₁₂ view factor, non-convergence)
    test_transient.py           # 34 tests (incl. the explicit path's window and remedies)
    test_zones.py               # 18 tests
    test_h_estimator.py         # 41 tests — forced, natural, mixed, opposing, Richardson, advisory, film range
    test_boundary_layer.py      # 40 tests — skin friction, BL thickness, buoyancy, inflation, regimes, prism rec, ER_v
    test_batch.py               # 92 tests — material DB, surface DB, h estimation, lateral, BL integration, input validation, warnings
    test_conduction_size_wall.py # 17 tests — thickness-aware sizing, N_cells + Biot reporting
    test_fail_loud.py           # 22 tests — one planted input per formerly silent path, each failing on 0.6.2
    test_input_guards.py        # 82 tests — one planted violation per guarded input
    test_open_closed_map.py     # 18 tests — OPEN-CLOSED-MAP.yaml schema + checker (not physics)
    test_release_statements.py  # 32 tests — version pins, tag links, "On PyPI since" vs CHANGELOG, post-releases, the documented test counts (not physics)
                                # 465 tests total
```

### Design Decisions

1. **Pure Python, no external dependencies.** This is intentional and non-negotiable — the calculators must run inside CAE pre-processor embedded Python environments (HyperMesh, ANSA, etc.) which have limited package availability. Do not add numpy, scipy, or any pip dependency.

2. **All temperatures in Kelvin internally.** Results include `_C` suffixed keys for convenience, but all physics calculations use Kelvin. This avoids sign errors in radiation (T^4) calculations.

3. **All lengths in meters internally, millimeters in output.** The `max_dx_mm` output keys are in mm (standard for CAE mesh sizing). Internal physics uses SI (meters).

4. **Static methods / class methods preferred.** No persistent state needed — each call is a pure function of its inputs. Classes are used only for logical grouping.

5. **Newton-Raphson solvers use pure Python** (Cramer's rule for 2x2 inversion in the multilayer shield solver). No numpy linalg. This is deliberate for portability.

6. **Dict returns, not dataclasses.** Keeps compatibility with Python 3.6+ environments found in older CAE tools. May migrate to dataclasses in a future version if the minimum Python target is raised.

## Key Physics Concepts

### The Core Insight

The conduction heat flux `q''` is a solver *output*, not an input. The original formulation `dx = k * dT / q''` is circular. We resolve this by substituting the surface energy balance:

```
q''_cond = h(Ts - Tf) + eps*sigma*(Ts^4 - Tsurr^4)
```

This makes mesh sizing a function of known/estimatable boundary conditions.

### Input Temperature Strategy (Automotive)

- **Exhaust components**: Surface temperatures are FIXED in the solver (measured or specified). Zero estimation error.
- **Non-exhaust structural**: Estimated to +/- 30 K from experience. Acceptable error.
- **Heat shields**: FLOATING unknown — must be solved iteratively. This is what `shields.py` handles via Newton-Raphson.

### Three Heat Transfer Modes → Three Mesh Constraints (Steady-State)

Each mode imposes an independent constraint. The binding constraint (smallest max element size) wins:

- **Conduction**: `dx_max = k * dT_max / q_total` — dominated by material conductivity
- **Convection**: Biot number (shell vs solid), h-gradient resolution, CFD mapping ratio
- **Radiation**: T^4 sensitivity (critical above ~500°C), view factor curvature limits

### Transient Constraints (v0.2)

For transient simulations, three additional constraints apply:

- **Penetration depth**: `dx_max <= C * sqrt(alpha * dt)` — thermal wave resolution, for implicit schemes only; an explicit scheme's own stability limit (the Fourier minimum) bounds it instead (docs/math_derivations.md §7.4)
- **Fourier number**: `dx_min = sqrt(alpha * dt / Fo_max)` — explicit solver stability (Fo <= 0.5 for 1D, 1/6 for 3D) or implicit accuracy (Fo <= 5)
- **Drive-cycle**: `dx_max <= sqrt(alpha * tau_bc)` — resolution across operating point transitions

### Convective HTC Estimation (v0.3)

The `h_estimator` module computes h from first principles rather than static lookup tables:

- **Forced convection**: Flat plate correlations — laminar (Re < 5e5) and turbulent (Re > 5e5)
- **Natural convection**: Buoyancy-driven — vertical plates, horizontal hot-up (Ra_crit=1e7), horizontal hot-down (always laminar)
- **Richardson number**: Ri = Gr/Re² — classifies regime as forced (Ri < 0.1), natural (Ri > 10), or mixed
- **Mixed convection (assisting)**: Asymptotic blending h_mixed = (h_forced³ + h_natural³)^(1/3)
- **Mixed convection (opposing, v0.4)**: h_mixed = max(|h_f³ − h_n³|^(1/3), k_air/L) — for horizontal_down with forced flow
- **External h cap**: 100 W/m²K absolute (H_EXTERNAL_MAX), 60 typical at 10 mph (H_EXTERNAL_TYPICAL); internal flows uncapped

### Solver Advisory System (v0.3)

The `solver_advisory()` function evaluates whether steady-state solvers are appropriate based on convection regime physics:

- **Turbulent natural convection** (oscillating buoyant plumes) is the primary SS convergence risk → flags `transient_advisory: True`, `severity: "warning"`
- **Turbulent forced convection** is SS-safe (time-averaged Nu correlations hold) → `severity: "info"`
- **Mixed regime** gets `severity: "caution"` — depends on which mode dominates
- **Negligible buoyancy** (near-zero dT) → `severity: "caution"` (numerical issues possible)

Ra thresholds for turbulence: vertical = 1e9, horizontal hot-up = 1e7, horizontal hot-down = inf (always laminar).

### Shield Solver API (v0.2)

Both single-layer and multilayer shield solvers now accept separate `h_in` / `h_out` convection coefficients. The single-layer solver also supports the legacy `h_total` parameter for backward compatibility (split evenly as `h_in = h_out = h_total / 2`).

### Material & Surface Databases (v0.3)

The `batch.py` module includes inline databases (no external files needed):

- **43 materials**: Steels (mild, HSLA, galvanised, SS304, SS409), Aluminium (6061, 5052, cast A356, A380), Cast iron (grey, ductile), Copper, Brass, Mg AZ91, Ti-6Al-4V, Zinc, Inconel 625, Filled plastics (PA66-GF30, PA6-GF30, PP-GF30, PBT-GF30, PPS-GF40), Unfilled plastics (PA66, PP, HDPE, ABS, PC, PC/ABS, POM, PET), Rubbers (EPDM, silicone, NBR, natural, neoprene, FKM), Composites/specialty (SMC, CFRP, glass, alumina, cordierite, fibreglass insulation, ceramic blanket)
- **30 surface treatments**: Bare/oxidised metals, galvanised, chrome/nickel/zinc plated, e-coat, powder coat, ceramic TBC, painted variants, rusted steel, glass, composite, fabric

Emissivity fallback logic: non-metals → 0.90, aluminium → 0.30, other metals → 0.73. A surface resolved this way is reported as a `SURFACE_DEFAULTED` warning.

## Current State (v0.6.2)

### What's Done
- [x] All seven calculator classes implemented and tested
- [x] Boundary-driven conduction (eliminates q'' as input)
- [x] Biot number with shell/solid recommendation and cell Biot heuristic (20·Bi, v0.4)
- [x] Fin theory lateral gradient constraint (v0.4) — `lateral_gradient_limit()` in conduction.py
- [x] h-gradient resolution and CFD mapping limit
- [x] Radiation T^4 sensitivity and view factor curvature limit
- [x] Single-layer shield solver (scalar Newton-Raphson, v0.2 h_in/h_out API)
- [x] Multilayer shield solver (2x2 NR) with gap view factor F₁₂ parameter (v0.4)
- [x] Transient mesh calculator (penetration depth, Fourier number, drive-cycle)
- [x] Combined transient constraint evaluator with binding-constraint identification
- [x] Convection zone system with velocity, orientation, t_air metadata
- [x] h_estimator: forced, natural, mixed (assisting + opposing v0.4) convection
- [x] Richardson number regime classification (forced / natural / mixed)
- [x] Solver advisory system (SS vs transient recommendations by regime)
- [x] Batch processor with 43-material, 30-surface-treatment databases
- [x] h_estimator wired into batch processor (non-shield parts)
- [x] Lateral gradient constraint wired into batch processor (v0.4)
- [x] Aerodynamic boundary layer calculator (v0.5) — two-regime model
- [x] BL constraint wired into batch processor as `aero_boundary_layer` dx_candidate (v0.5)
- [x] Automatic batch warnings from parametric study thresholds (v0.3)
- [x] Parametric regime crossover studies (analysis/ directory)
- [x] 11 worked automotive examples with verified output
- [x] 415 pytest unit tests with analytical verification and energy balance closure
- [x] Inputs checked (unreleased; CHANGELOG `[Unreleased]`): unknown part names raise `PartInputError` naming the allowed set, physical inputs are range-guarded, and defaulted or extrapolated inputs, shield non-convergence and transient conflicts come back as coded warnings
- [x] Full mathematical derivation documentation (docs/math_derivations.md, Sections 1–14)
- [x] Published to PyPI — `pip install thermal-mesh-calculators`. 0.6.1 (2026-09-16,
      `requires-python` narrowed to `>=3.9` for it) and 0.6.2 (2026-09-17) were
      withdrawn from the index; it serves 0.6.2.post1 (2026-09-26), the 0.6.2 code
      under a new file name, tagged `v0.6.2.post1` from `v0.6.2` rather than `main`

### Key Findings from Parametric Studies

1. **Conduction vs convection mesh dominance**: For metals, conduction always governs (Biot << 0.1). For plastics, Biot can cross 0.1 at h ≈ 25 W/m²K but conduction usually still governs.
2. **Laminar → turbulent natural convection**: Vertical safe below ~600mm char length. Horizontal hot-up transitions at ~130–190mm (nearly all dead-zone parts are turbulent). Horizontal hot-down always laminar.
3. **10 mph forced suppresses turbulent buoyancy**: Critical velocity typically 3–9 mph. At 10 mph baseline, all tested scenarios show Ri < 0.1 (forced-dominant, SS-safe).

### v0.4.0 Physics Improvements (External Review)

Four changes based on external peer review of the mathematical derivations:

1. **Fin theory lateral gradient** — `lateral_gradient_limit()` in conduction.py. Ensures surface mesh resolves in-plane hot spots using `Δx ≤ (1/3) · √(k·t/h_total)` where `h_total = h + 4εσT³`. Wired into batch.py as `lateral_gradient` constraint.

2. **Cell Biot number (20·Bi)** — Changed from `ceil(10·Bi)` to `ceil(20·Bi)` in convection.py. Ensures per-element `Bi_cell ≤ 0.1`, preventing numerical oscillation in explicit solvers.

3. **Gap view factor F₁₂** — Added `f12` parameter (default 1.0) to `MultilayerShieldCalculator.solve_temperatures()`. Generalises gap emissivity formula: `ε_eff = 1/(1/ε_g1 + 1/ε_g2 − 2 + 1/F₁₂)`. Backward compatible (F₁₂=1 gives original formula).

4. **Opposing mixed convection** — `estimate_h()` now detects opposing flow (horizontal_down orientation in mixed regime) and uses `h = max(|h_f³ − h_n³|^(1/3), k_air/L)` with a conduction floor to prevent singularity.

### v0.5.0 Aerodynamic Boundary Layer Constraints

New module `boundary_layer.py` with `BoundaryLayerCalculator` class. Computes inflation layer parameters (first cell height, layer count, growth) and surface mesh constraints from aerodynamic boundary layer physics.

**Two-regime model:**

1. **External forced** — attached, unidirectional flow (underbody panels in freestream). Uses flat-plate skin-friction correlations (Schlichting) with surface constraint: `dx ≤ min(AR_transition · y_last, C_BL · δ)`. The last-layer transition constraint prevents degenerate cells at the prism-to-tet interface.

2. **Mixed / unknown** — recirculating, impinging, or buoyancy-dominated regions (cargo bed above exhaust, wake zones). Uses `max(U, V_buoyancy)` as characteristic velocity. Surface constraint: `dx ≤ C_BL · δ` (BL fraction only). **Prism layers are NOT recommended** — isotropic tets to wall.

**Prism recommendation logic (v0.5.1):**
- `prisms_recommended` flag evaluates whether inflation layers provide meaningful accuracy over isotropic tets
- **External forced**: always recommends prisms (steep wall gradients, directional anisotropy)
- **Mixed/unknown**: generally does NOT recommend prisms, based on two criteria:
  1. *y₁ vs dx*: At low velocities, y⁺-derived y₁ (10–60mm) vastly exceeds thermal surface mesh (1–5mm). Tets already over-resolve the wall region.
  2. *Wall gradient*: dT/dy ~ 10⁴ K/m (low-Re) vs ~10⁷ K/m (high-Re). The 500× smaller gradient makes tet non-orthogonality error (~10–20% of h) fall within the ±30% physics uncertainty of h itself.
- Volume expansion ratio ER_v = (dx/y_last)³ reported for external forced: stable ≤ 5, marginal ≤ 10, unstable > 10

**Key safety features:**
- Leading-edge singularity protection (Re_x_min = 1000, y₁ absolute floor)
- Wall-function default (y⁺ = 30) — wall-resolved y⁺=1 available but not recommended for vehicle-level
- Buoyancy velocity scale `V_buoy = √(g·β·ΔT·L)` for dead zones
- Integrated into batch.py as optional `aero_boundary_layer` constraint (activated by `bl_regime` parameter)

### What's NOT Done (Roadmap)
- [ ] Integration with pre-processor APIs (HyperMesh, ANSA)
- [ ] CLI interface (argparse)
- [ ] Parametric sweep / sensitivity plot generation
- [ ] Composite / multi-material through-thickness conduction (currently single-material only)
- [ ] N-layer shield generalisation (currently limited to 2 layers)
- [ ] Temperature-dependent material properties (piecewise linear tables)

## Coding Conventions

- Type hints on all public method signatures
- Docstrings in NumPy format (Parameters / Returns sections)
- Module-level docstrings explain the physics and governing equations
- No abbreviations in public API names (use `max_mesh_size` not `max_ms`)
- Private helpers prefixed with underscore
- Snake_case throughout

## Running the Examples

```bash
cd thermal-mesh-calculators
python -m examples.automotive_examples
```

No install step needed to run from a clone — the package imports from the repo root.

As of v0.6.0 the project is also a real distributable package (`pyproject.toml`,
PEP 621, setuptools backend):

```bash
python -m build                       # -> dist/*.whl + dist/*.tar.gz
pip install thermal-mesh-calculators  # from PyPI: 0.6.2.post1
```

The distribution version is read dynamically from
`thermal_mesh_calculators.__version__`, so bump that one constant and the wheel
follows. Runtime dependencies must stay empty — see Design Decisions #1.

## Testing

```bash
python -m pytest tests/ -v
python -m pytest --doctest-glob=README.md README.md   # the README's examples
python -m ruff check .                                # the rule set written out in pyproject.toml (E4, E7, E9, F)
python -m mypy                                        # the package's annotations ([tool.mypy])
```

CI runs all four on every push and pull request (ruff and mypy at pinned versions),
and also installs the package with `pip install .` and imports it from outside the
checkout.

465 tests across 14 test files (415 physics and input-validation + 18 for the open/closed map checker + 32 for the release statements and the documented test counts).
Test strategy:
- Hand-computed analytical solutions for known inputs
- Edge cases (zero flux, zero emissivity, pure convection/radiation)
- Energy balance closure for Newton-Raphson solvers (q_in ≈ q_out)
- Dimensional sanity (mm output, K→C conversion)
- Scaling relationships (k proportionality, sqrt(dt) scaling)
- Backward compatibility (legacy `h_total` parameter)
- Convection regime classification (Richardson number boundaries)
- Solver advisory correctness (severity levels by regime)
- Material/surface database coverage and fallback logic
- Fail-loud gate (`tests/test_fail_loud.py`): one planted input per formerly silent path, each failing against 0.6.2
- Input guards (`tests/test_input_guards.py`): one planted violation per guarded parameter
- Documented test counts (`tests/test_release_statements.py`): every count stated in this file and in `.github/copilot-instructions.md` is checked against a collection of the suite, so adding a test means updating those counts

## Common Agent Tasks

### "Add a new calculator"
1. Create a new module in `thermal_mesh_calculators/`
2. Add a class with static/class methods
3. Export it from `__init__.py`
4. Add an example scenario in `examples/automotive_examples.py`
5. Add tests in `tests/`
6. Update `README.md` module descriptions

### "Add a new material"
- Add entry to `MATERIALS` dict in `batch.py`
- Keys: `k` (W/mK), `rho` (kg/m³), `cp` (J/kgK), `description` (str); the diffusivity alpha = k/(rho*cp) is computed where it is needed, not stored
- Group in appropriate category comment block (steels, aluminium, plastics, etc.)
- Run `python -m pytest tests/test_batch.py -v` to verify

### "Add a new surface treatment"
- Add entry to `SURFACE_TREATMENTS` dict in `batch.py`
- Key: treatment name (lowercase, underscored)
- Value: dict with `epsilon` (emissivity, 0–1) and `description` (str)
- Run `python -m pytest tests/test_batch.py -v` to verify

### "Add a new convection zone"
- Add entry to `CONVECTION_ZONES` dict in `zones.py`
- Required keys: `h_low`, `h_high`, `velocity_ms`, `t_air_C_low`, `t_air_C_high`, `orientation`, `regime`, `is_internal`, `notes`
- Run `python -m pytest tests/test_zones.py -v` to verify

### "Add CLI interface"
- Use `argparse` (stdlib) not click/typer (external dependency)
- Accept material properties and boundary conditions as arguments
- Output results as formatted table or JSON

## Origin

This project grew from a discussion about the physical dependencies of thermal mesh sizing for automotive simulations. The conversation evolved through several iterations:
1. Initial physics breakdown (conduction, convection, radiation dependencies)
2. First-pass calculators using q'' directly
3. Correction: q'' is unknown → substitute boundary energy balance
4. Single-layer shield solver (Newton-Raphson for floating shield temp)
5. Multilayer shield solver (multivariate Newton-Raphson for dual-wall shields)
6. v0.2: Aligned single-layer shield API (h_in/h_out), added transient calculator, zones, batch processor
7. v0.3: h_estimator (correlation-based HTC), solver advisory, 43 materials, 30 surfaces, regime studies
8. v0.4: Physics improvements (lateral gradient, cell Biot 20x, F₁₂ view factor, opposing mixed convection)
9. v0.5: Aerodynamic boundary layer calculator (two-regime: external_forced + mixed_unknown)
10. v0.5.1: Prism recommendation logic — low-Re tet sufficiency analysis, wall gradient discriminator, ER_v reporting
11. v0.6.0: Thickness-aware `size_wall()` (N_cells + Biot reporting, zero-flux pole guard, prescribed-Ts consistency check); first **packaged** release — `pyproject.toml`, Apache-2.0 `LICENSE`, `CHANGELOG.md`, git tag `v0.6.0` (not in this public repository, whose first tag is `v0.6.2`), buildable wheel
12. v0.6.1: `requires-python` narrowed to >=3.9; first public PyPI release (2026-09-16).
