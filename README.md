# thermal-mesh-calculators

Physics-driven thermal mesh sizing calculators for automotive CAE — conduction, convection, radiation, and heat shield solvers.

## Motivation

Mesh sizing for thermal analysis is often driven by rules of thumb ("use 5 mm everywhere") rather than the underlying physics. This repository provides calculators that derive maximum element sizes directly from material properties, boundary conditions, and heat transfer governing equations. The result is meshes that are fine where the physics demands it and coarse where it doesn't — saving node count without sacrificing accuracy.

## Key Insight: Eliminating Q as an Input

The conduction heat flux `q''` is a solver *output*, not an input — you don't know it before you run. The original Fourier's law formulation (`dx = k * dT / q''`) is therefore circular.

We resolve this by substituting the surface energy balance. At steady state, `q''_cond = q''_conv + q''_rad`, giving:

```
dx_max = k * dT_max / |h(Ts - Tf) + eps*sigma*(Ts^4 - Tsurr^4)|
```

All terms on the right are either known inputs (material props, BCs) or can be estimated. For exhaust components where surface temperatures are fixed in the solver, there is zero estimation error.

## Modules

### `conduction.py` — Boundary-Driven Mesh Sizing
- Substitutes convection + radiation boundary fluxes for the unknown conduction flux
- Works directly with fixed surface temperatures (exhaust) or estimates (structural)
- Includes transient penetration depth calculator for time-dependent problems

### `convection.py` — Solid-Side Convection Constraints
- **Biot number** evaluation: decides shell (2D) vs. solid (3D) meshing with recommended element count through thickness
- **h-gradient resolution**: ensures surface mesh captures spatial variations in heat transfer coefficient (impingement zones, separation regions)
- **CFD mapping limit**: max solid element size to avoid interpolation loss when mapping wall data from CFD

### `radiation.py` — T^4 Sensitivity & View Factors
- Flux sensitivity calculator (`dq/dT = 4*eps*sigma*T^3`) — shows the 120x increase in sensitivity from 300 K to 1000 K
- Element size bound from linearisation error tolerance
- Curvature-based view factor limit: max facet size on curved surfaces to prevent artificial hot/cold spots

### `shields.py` — Heat Shield Solvers
- **Single-layer**: Newton-Raphson solver for thin shield equilibrium temperature, then boundary-driven mesh size
- **Multilayer**: Multivariate Newton-Raphson (2x2 Jacobian, Cramer's rule) for dual-wall shields with air gap. Solves both layer temperatures simultaneously. Outputs independent mesh sizes per layer.

### Also in the package

`transient.py` (penetration depth, Fourier number, drive-cycle limits) ·
`zones.py` (convection zone lookup with velocity / orientation / air-temp metadata) ·
`h_estimator.py` (HTC from correlations — forced, natural and mixed convection,
Richardson-number regime classification, steady-state solver advisory) ·
`boundary_layer.py` (aerodynamic boundary layer: y+, inflation layers, surface dx) ·
`batch.py` (BOM → mesh sizes, with 43 materials and 30 surface treatments inline).

Full signatures are in [`docs/api_reference.md`](docs/api_reference.md); the
derivations behind every constraint are in
[`docs/math_derivations.md`](docs/math_derivations.md).

## Installation

```bash
pip install thermal-mesh-calculators
```

On PyPI. To pin an immutable release straight from git instead:

```bash
pip install "thermal-mesh-calculators @ git+https://github.com/BootstrapAI-mgmt/thermal-mesh-calculators.git@v0.6.2.post1"
```

The package has **no runtime dependencies** — it is pure standard library, by
design, so it can run inside CAE pre-processor embedded interpreters (HyperMesh,
ANSA, …) where installing numpy or scipy is often not an option. Requires Python
**3.9 or newer**.

> The source deliberately avoids syntax newer than 3.6, so the *wheel* runs on
> older embedded interpreters. The declared floor is 3.9 because the PEP 639
> licence metadata needs `setuptools>=77` to **build**, and that needs 3.9+. The
> two statements are both true and are not in conflict: installing the wheel is
> not building the project.

Working from a clone needs no install step at all — the package imports directly
from the repository root.

## Quick Start

```python
from thermal_mesh_calculators import (
    BoundaryDrivenConductionCalculator,
    ConvectionMeshCalculator,
    RadiationMeshCalculator,
    SingleLayerShieldCalculator,
    MultilayerShieldCalculator,
)

# Exhaust manifold — fixed surface temp at 800 C
cond = BoundaryDrivenConductionCalculator()
result = cond.max_mesh_size(
    k=45.0, h=150.0, t_surf=1073.15, t_fluid=353.15,
    epsilon=0.85, t_surr=353.15, max_dt=10.0,
)
print(f"Max element size: {result['max_dx_mm']:.2f} mm")
# -> 2.63 mm

# v0.6.0 — thickness-aware reporting: cell count + Biot number, never a raw
# length. The raw relation above knows nothing about the part it is sizing;
# size_wall clamps against the actual wall and reports the modelling decision.
sized = cond.size_wall(
    k=45.0, h=150.0, t_surf=1073.15, t_fluid=353.15,
    epsilon=0.85, t_surr=353.15, max_dt=10.0,
    thickness_m=0.004,  # 4 mm cast wall
)
print(f"{sized['n_cells']} cell(s), Bi = {sized['biot']:.3f} -> {sized['regime']}")
# -> 2 cell(s), Bi = 0.021 -> thermally_thin  (use shell conduction)

# Dual-wall heat shield
multi = MultilayerShieldCalculator()
m = multi.mesh_sizes(
    k_metal=45.0, max_dt=15.0,
    t_exh=1073.15, t_fluid=353.15, t_surr=353.15,
    h_in=30.0, h_out=30.0, h_gap=15.0,
    eps_in=0.4, eps_out=0.4, eps_g1=0.4, eps_g2=0.4,
)
print(f"Layer 1: {m['t1_C']:.0f} C -> {m['layer1_max_dx_mm']:.1f} mm")
print(f"Layer 2: {m['t2_C']:.0f} C -> {m['layer2_max_dx_mm']:.1f} mm")
# -> Layer 1: 521 C -> 19.7 mm
# -> Layer 2: 282 C -> 85.9 mm (4.4x coarser)
```

## Example Output

Run `python -m examples.automotive_examples` for a full walkthrough covering:

| Scenario | Raw `dx_max` | Thickness-aware (`size_wall`, v0.6.0) |
|---|---|---|
| Steel exhaust manifold (800 C, ~4 mm wall) | 2.6 mm | 2 cells — **thermally thin** (Bi ≈ 0.02): shell conduction |
| Plastic intake manifold (120 C, ~3 mm wall) | 0.6 mm | **5 cells — resolve** (Bi ≈ 0.25, ΔT_wall ≈ 24 K) |
| Single-layer aluminised shield (~1 mm) | 19.3 mm | 1 cell — thermally thin |
| Multilayer shield, Layer 1 (~1 mm) | 19.7 mm | 1 cell — thermally thin (Bi ≈ 0.001) |
| Multilayer shield, Layer 2 (~1 mm) | 85.9 mm | 1 cell — thermally thin (Bi ≈ 0.0002) |
| Radiation limit at 800 C (2 K/mm gradient) | 1.0 mm | (surface constraint — unchanged) |

> **Why the second column exists.** A raw `dx_max` that exceeds the part is the
> relation reporting *"this wall is thermally thin"*, **not a cell size** — an
> 85.9 mm "element" for a ~1 mm shield layer is a category error in the
> reporting layer, found by stress-testing the tool against real part
> thicknesses. `size_wall` closes it: it requires the wall thickness, clamps to
> `N_cells = max(1, ceil(L / dx_max))`, and leads with the two numbers that
> support a modelling decision — the cell count and the Biot number — keeping
> the raw dx only for traceability. It also guards the pole at zero net flux
> (convection and radiation cancelling is equilibrium, not an infinite cell)
> and flags the physically inconsistent case where the implied through-wall
> drop exceeds the surface-to-fluid driving difference (the prescribed-Ts
> assumption has broken down). This is a **Biot-number triage screen** —
> `N_cells = Bi·(Ts−Tf)/ΔT_max` — for deciding which solids in an assembly get
> through-thickness resolution; it is not an error estimator.

## Discussion Notes

### What drives mesh density in practice

**Conduction**: Material conductivity is the dominant factor. Low-k materials (plastics at 0.25 W/mK) require dramatically finer meshes than metals (steel at 45 W/mK) even at moderate temperatures. The 0.6 mm result for the plastic intake manifold vs. 2.6 mm for the steel exhaust manifold — despite the exhaust being 680 C hotter — illustrates this clearly.

**Radiation**: Generally the least sensitive mode for automotive meshes *below ~400 C*. Above 500 C the T^4 dependence becomes aggressive: flux sensitivity increases from 20 W/m^2K at 200 C to 398 W/m^2K at 1000 C. At high temperatures, radiation can become the binding mesh constraint.

**Convection**: The solid mesh constraint from convection is indirect — it's really about resolving the spatial distribution of h mapped from CFD. Impingement zones with steep h-gradients can require very fine surface meshes (2-3 mm) even when the conduction physics would allow coarser.

### Shield solver assumptions
- Shields are treated as lumped capacitance (Bi << 0.1) — valid for typical stamped sheet metal
- Convection on each shield face is parameterised as a single HTC value. For more complex flow fields, use the solved shield temperature as T_surf in the boundary-driven conduction calculator with local h values
- The multilayer gap conductance `h_gap` is an effective value that should include both air conduction and any contact conduction through dimples or spot welds

## Roadmap

Shipped since this list was first written:

- [x] Transient mesh sizing (Fourier number constraint) — `transient.py`,
      `TransientMeshCalculator.fourier_number_limit()`
- [x] Material property database → batch component analysis — `batch.py`, 43
      materials and 30 surface treatments held inline (no external data file, to
      keep the zero-dependency guarantee)
- [x] Multi-component batch report generator — `process_batch()` +
      `summary_table()`; see `examples/batch_example.py`

Still open:

- [ ] Load a BOM from CSV/JSON — `process_batch()` accepts Python dicts only
- [ ] Parametric sweep / sensitivity plots (mesh size vs. temperature, emissivity, etc.)
- [ ] Integration with pre-processor APIs (HyperMesh, ANSA)
- [ ] Composite / multi-material through-thickness conduction (single-material today)
- [ ] N-layer shield generalisation (two layers today)

## License

Licensed under the [Apache License 2.0](LICENSE). Copyright 2026 BootstrapAI-mgmt.
