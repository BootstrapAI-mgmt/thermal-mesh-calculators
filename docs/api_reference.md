# thermal-mesh-calculators API Reference

A pure Python package for mesh sizing and thermal analysis of heat shields and conduction problems. No external dependencies required.

---

## Table of Contents

1. [constants](#constants)
2. [conduction](#conduction)
3. [convection](#convection)
4. [radiation](#radiation)
5. [shields](#shields)
6. [transient](#transient)
7. [zones](#zones)
8. [h_estimator](#h_estimator)
9. [batch](#batch)
10. [Quick Start Examples](#quick-start-examples)

---

## constants

Module containing physical constants used throughout the package.

### STEFAN_BOLTZMANN
Stefan-Boltzmann constant for radiation calculations.

**Value:** `5.67e-8` (W/m²·K⁴)

**Example:**
```python
from thermal_mesh_calculators.constants import STEFAN_BOLTZMANN

# Use in radiation flux calculation
q_rad = STEFAN_BOLTZMANN * emissivity * (T1**4 - T2**4)
```

---

## conduction

Thermal conduction mesh sizing and transient analysis.

### BoundaryDrivenConductionCalculator

Class for analyzing conduction-dominated heat transfer with boundary-driven conditions.

#### max_mesh_size()

Calculates maximum mesh size considering conduction with convective and radiative boundary conditions.

**Signature:**
```python
@staticmethod
def max_mesh_size(k, h, t_surf, t_fluid, epsilon, t_surr, max_dt) -> dict
```

**Parameters:**
- `k` (float): Thermal conductivity (W/m·K)
- `h` (float): Convective heat transfer coefficient (W/m²·K)
- `t_surf` (float): Surface temperature (K)
- `t_fluid` (float): Fluid temperature (K)
- `epsilon` (float): Emissivity (0–1)
- `t_surr` (float): Surroundings temperature (K)
- `max_dt` (float): Maximum allowable temperature drop across one element (K), the accuracy target

**Returns:**
- `dict` with keys:
  - `max_dx_mm` (float): Maximum element size (mm)
  - `q_conv` (float): Convective heat flux (W/m²)
  - `q_rad` (float): Radiative heat flux (W/m²)
  - `q_total` (float): Total heat flux (W/m²)
  - `rad_fraction` (float): Radiation as fraction of total (0–1)

**Example:**
```python
from thermal_mesh_calculators.conduction import BoundaryDrivenConductionCalculator

result = BoundaryDrivenConductionCalculator.max_mesh_size(
    k=50.0,              # W/m·K (carbon steel)
    h=75.0,              # W/m²·K
    t_surf=873.15,       # K (600°C)
    t_fluid=343.15,      # K (70°C)
    epsilon=0.7,
    t_surr=343.15,       # K (70°C)
    max_dt=10.0          # K, allowable drop per element
)

print(f"Max element size: {result['max_dx_mm']:.2f} mm")
print(f"Heat flux: {result['q_total']:.1f} W/m² ({result['rad_fraction']*100:.1f}% radiation)")
```

#### size_wall()

Thickness-aware sizing: the through-thickness cell count and the Biot number of a
wall, rather than a raw length. `max_mesh_size()` knows nothing about the part it
sizes, and for a thin, conductive wall its `max_dx_mm` can exceed the wall itself:
that is the relation reporting "thermally thin", not a cell size. `size_wall()`
clamps against the wall thickness and leads with the cell count and the Biot
number. The derivation is in [`math_derivations.md` §1.6](math_derivations.md).

**Signature:**
```python
@staticmethod
def size_wall(k, h, t_surf, t_fluid, epsilon, t_surr, max_dt, thickness_m, max_cells=100) -> dict
```

**Parameters:**
- `k`, `h`, `t_surf`, `t_fluid`, `epsilon`, `t_surr`, `max_dt`: as for `max_mesh_size()`
- `thickness_m` (float): Wall thickness (m); required
- `max_cells` (int, default 100): Cap on the reported cell count

**Returns:**
- `dict` with keys:
  - `n_cells` (int): Cells through the wall, `max(1, ceil(thickness / raw dx))`, capped at `max_cells`
  - `biot` (float): `h_effective * thickness / k` (infinite when `t_surf == t_fluid`)
  - `regime` (str): `"thermally_thin"` (Bi < 0.1), `"resolve"` (0.1 ≤ Bi < 1), `"steep_gradient"` (Bi ≥ 1), or `"near_equilibrium"` when convection and radiation cancel (then one cell, and `biot` is 0)
  - `dx_used_mm` (float): `thickness / n_cells` (mm)
  - `dx_raw_mm` (float): The unclamped `max_mesh_size()` length (mm), for traceability
  - `dx_exceeds_part` (bool): Whether the raw length exceeds the wall
  - `dt_wall_K` (float): Temperature drop across the wall implied by the surface flux (K)
  - `h_effective` (float): Net surface flux divided by `t_surf - t_fluid` (W/m²·K), radiation included
  - `q_conv`, `q_rad`, `q_total`, `rad_fraction`: as for `max_mesh_size()`
  - `recommendation` (str): The modelling decision, in words
  - `warnings` (list of str): Near equilibrium; the cell cap reached; `t_surf == t_fluid`; the raw length exceeding the wall; physically inconsistent (the implied drop exceeds `t_surf - t_fluid`, so the surface temperature cannot be prescribed independently of the wall); radiation carrying more than 60 % of the flux

`biot` uses the whole wall thickness, where `ConvectionMeshCalculator.biot_number()`
uses half of it, and includes radiation. It carries the sign of `t_surf - t_fluid`:
for a wall colder than its fluid it is negative and `regime` reads
`"thermally_thin"` whatever the wall, so judge such a wall by `abs(biot)`.

**Raises:**
- `ValueError`: `thickness_m`, `k` or `max_dt` not > 0, or any input `max_mesh_size()` rejects

**Example:**
```python
from thermal_mesh_calculators.conduction import BoundaryDrivenConductionCalculator

sized = BoundaryDrivenConductionCalculator.size_wall(
    k=0.25, h=40.0,                  # glass-filled nylon intake manifold
    t_surf=393.15, t_fluid=353.15,   # K (120°C surface, 80°C air)
    epsilon=0.92, t_surr=353.15,     # K
    max_dt=5.0,                      # K per element
    thickness_m=0.003,               # m (3 mm wall)
)
print(sized["n_cells"], sized["regime"], round(sized["biot"], 2))
# 5 resolve 0.61
```

#### transient_penetration_depth()

Calculates thermal penetration depth for transient analysis.

**Signature:**
```python
@staticmethod
def transient_penetration_depth(k, rho, cp, dt) -> float
```

**Parameters:**
- `k` (float): Thermal conductivity (W/m·K)
- `rho` (float): Density (kg/m³)
- `cp` (float): Specific heat (J/kg·K)
- `dt` (float): Time step (s)

**Returns:**
- `float`: Penetration depth (mm)

**Example:**
```python
depth_mm = BoundaryDrivenConductionCalculator.transient_penetration_depth(
    k=50.0,       # W/m·K
    rho=7850.0,   # kg/m³ (steel)
    cp=500.0,     # J/kg·K
    dt=0.1        # s
)
print(f"Penetration depth: {depth_mm:.3f} mm")
```

---

## convection

Convective heat transfer and Biot number analysis.

### ConvectionMeshCalculator

Class for convection-based mesh sizing criteria.

#### biot_number()

Calculates Biot number and recommended mesh refinement for convection dominance.

**Signature:**
```python
@staticmethod
def biot_number(h, k, thickness) -> dict
```

**Parameters:**
- `h` (float): Convective heat transfer coefficient (W/m²·K)
- `k` (float): Thermal conductivity (W/m·K)
- `thickness` (float): Material thickness (m)

**Returns:**
- `dict` with keys:
  - `biot` (float): Biot number `h · (thickness / 2) / k` (dimensionless)
  - `mesh_type` (str): `"2D Shell"` (Bi < 0.1) or `"3D Solid"`
  - `min_elements_through_thickness` (int): 1 for a shell; `ceil(20 · Bi)` clamped to [2, 10] for a solid
  - `rationale` (str): Explanation of mesh type classification

**Example:**
```python
from thermal_mesh_calculators.convection import ConvectionMeshCalculator

result = ConvectionMeshCalculator.biot_number(
    h=100.0,               # W/m²·K
    k=50.0,                # W/m·K (carbon steel)
    thickness=0.001        # m (1 mm)
)

print(f"Biot number: {result['biot']:.4f}")
print(f"Mesh type: {result['mesh_type']}")
print(f"Min elements: {result['min_elements_through_thickness']}")
```

#### h_gradient_mesh_limit()

Calculates maximum mesh size for resolving convection coefficient gradients.

**Signature:**
```python
@staticmethod
def h_gradient_mesh_limit(h_max, h_min, gradient_length_mm) -> dict
```

**Parameters:**
- `h_max` (float): Maximum h value (W/m²·K)
- `h_min` (float): Minimum h value (W/m²·K)
- `gradient_length_mm` (float): Length over which gradient occurs (mm)

**Returns:**
- `dict` with keys:
  - `max_dx_mm` (float): Maximum element size (mm)
  - `h_ratio` (float): Ratio h_max / h_min
  - `elements_needed` (int): Recommended number of elements

**Example:**
```python
result = ConvectionMeshCalculator.h_gradient_mesh_limit(
    h_max=150.0,           # W/m²·K
    h_min=50.0,            # W/m²·K
    gradient_length_mm=50.0
)

print(f"Max element size: {result['max_dx_mm']:.2f} mm")
print(f"Elements needed: {result['elements_needed']}")
```

#### cfd_mapping_mesh_limit()

Calculates mesh size limit for mapping CFD to FEA.

**Signature:**
```python
@staticmethod
def cfd_mapping_mesh_limit(fluid_wall_face_mm, mapping_ratio=4.0) -> float
```

**Parameters:**
- `fluid_wall_face_mm` (float): CFD mesh face size at wall (mm)
- `mapping_ratio` (float, default=4.0): FEA-to-CFD mesh ratio

**Returns:**
- `float`: Maximum FEA mesh element size (mm)

**Example:**
```python
fea_max_dx = ConvectionMeshCalculator.cfd_mapping_mesh_limit(
    fluid_wall_face_mm=0.5,
    mapping_ratio=4.0
)
print(f"FEA max element: {fea_max_dx:.3f} mm")
```

---

## radiation

Radiative heat transfer and mesh sizing.

### RadiationMeshCalculator

Class for radiation-based mesh sizing criteria.

#### flux_sensitivity()

Calculates sensitivity of radiation flux to local temperature changes.

**Signature:**
```python
@staticmethod
def flux_sensitivity(t_local, emissivity) -> float
```

**Parameters:**
- `t_local` (float): Local temperature (K)
- `emissivity` (float): Surface emissivity (0–1)

**Returns:**
- `float`: Flux sensitivity (W/m²·K)

**Example:**
```python
from thermal_mesh_calculators.radiation import RadiationMeshCalculator

dq_dT = RadiationMeshCalculator.flux_sensitivity(
    t_local=873.15,    # K (600°C)
    emissivity=0.7
)
print(f"Flux sensitivity: {dq_dT:.2f} W/m²·K")
```

#### max_mesh_size()

Calculates maximum mesh size to limit radiation flux error.

**Signature:**
```python
@staticmethod
def max_mesh_size(t_local, emissivity, allowable_flux_error, spatial_gradient) -> dict
```

**Parameters:**
- `t_local` (float): Local temperature (K)
- `emissivity` (float): Surface emissivity (0–1)
- `allowable_flux_error` (float): Allowable error in heat flux (W/m²)
- `spatial_gradient` (float): Temperature gradient in solid (K/m)

**Returns:**
- `dict` with keys:
  - `max_dx_mm` (float): Maximum element size (mm)
  - `dq_dt` (float): Flux sensitivity (W/m²·K)
  - `max_dt_element` (float): Maximum time step for element (s)

**Example:**
```python
result = RadiationMeshCalculator.max_mesh_size(
    t_local=873.15,
    emissivity=0.7,
    allowable_flux_error=100.0,  # W/m²
    spatial_gradient=50000.0      # K/m
)

print(f"Max element size: {result['max_dx_mm']:.2f} mm")
```

#### view_factor_curvature_limit()

Calculates mesh size limit based on surface curvature for view factors.

**Signature:**
```python
@staticmethod
def view_factor_curvature_limit(radius_mm, max_facet_angle_deg=15.0) -> float
```

**Parameters:**
- `radius_mm` (float): Surface radius of curvature (mm)
- `max_facet_angle_deg` (float, default=15.0): Maximum subtended angle (degrees)

**Returns:**
- `float`: Maximum element chord length (mm)

**Example:**
```python
max_dx = RadiationMeshCalculator.view_factor_curvature_limit(
    radius_mm=25.0,
    max_facet_angle_deg=15.0
)
print(f"Max element size for curvature: {max_dx:.2f} mm")
```

---

## shields

Single and multilayer thermal shield analysis.

### SingleLayerShieldCalculator

Class for single-layer heat shield thermal analysis.

#### solve_temperature()

Solves steady-state shield temperature given boundary conditions.

**Signature:**
```python
@staticmethod
def solve_temperature(
    t_exh, t_fluid, t_surr, eps_in, eps_out,
    h_in=None, h_out=None, h_total=None,
    tol=0.1, max_iter=50
) -> dict
```

**Parameters:**
- `t_exh` (float): Exhaust/hot side temperature (K)
- `t_fluid` (float): Cooling fluid temperature (K)
- `t_surr` (float): Surroundings temperature (K)
- `eps_in` (float): Inner surface emissivity (0–1)
- `eps_out` (float): Outer surface emissivity (0–1)
- `h_in` (float, optional): Inner surface convection coefficient (W/m²·K)
- `h_out` (float, optional): Outer surface convection coefficient (W/m²·K)
- `h_total` (float, optional): Combined h if uniform (W/m²·K)
- `tol` (float, default=0.1): Convergence tolerance (K)
- `max_iter` (int, default=50): Maximum iterations

**Returns:**
- `dict` with keys:
  - `t_shield_K` (float): Shield temperature (K)
  - `t_shield_C` (float): Shield temperature (°C)
  - `converged` (bool): Whether solution converged; when False, the other keys describe the last iterate
  - `iterations` (int): Number of iterations
  - `residual_W_m2` (float): Energy-balance residual at `t_shield_K` (W/m²)
  - `q_rad_in` (float): Radiation from hot side (W/m²)
  - `q_rad_out` (float): Radiation to surroundings (W/m²)
  - `q_conv_in` (float): Convection from hot side (W/m²)
  - `q_conv_out` (float): Convection to cooling fluid (W/m²)
  - `h_in` (float): Actual h_in used (W/m²·K)
  - `h_out` (float): Actual h_out used (W/m²·K)

**Example:**
```python
from thermal_mesh_calculators.shields import SingleLayerShieldCalculator

result = SingleLayerShieldCalculator.solve_temperature(
    t_exh=873.15,      # K (600°C exhaust)
    t_fluid=343.15,    # K (70°C coolant)
    t_surr=343.15,     # K (70°C surroundings)
    eps_in=0.8,        # Hot side emissivity
    eps_out=0.8,       # Outer emissivity
    h_in=150.0,        # W/m²·K (high-speed exhaust)
    h_out=100.0        # W/m²·K (coolant)
)

print(f"Shield temperature: {result['t_shield_C']:.1f}°C")
print(f"Converged: {result['converged']}")
print(f"Iterations: {result['iterations']}")
```

#### mesh_size()

Calculates mesh constraints for shield temperature solution with material properties.

**Signature:**
```python
@staticmethod
def mesh_size(
    k, max_dt, t_exh, t_fluid, t_surr, eps_in, eps_out,
    h_in=None, h_out=None, h_total=None, tol=0.1, max_iter=50
) -> dict
```

**Parameters:**
- `k` (float): Material thermal conductivity (W/m·K)
- `max_dt` (float): Maximum temperature difference across one element (K)
- Additional parameters same as `solve_temperature()`, including `tol` and `max_iter`, which are forwarded to it

**Returns:**
- `dict` with all keys from `solve_temperature()` plus:
  - `max_dx_mm` (float): Maximum element size (mm)
  - `q_boundary` (float): Boundary heat flux (W/m²)

**Example:**
```python
result = SingleLayerShieldCalculator.mesh_size(
    k=50.0,            # W/m·K (stainless steel)
    max_dt=0.1,        # s
    t_exh=873.15,
    t_fluid=343.15,
    t_surr=343.15,
    eps_in=0.8,
    eps_out=0.8,
    h_in=150.0,
    h_out=100.0
)

print(f"Max element size: {result['max_dx_mm']:.2f} mm")
```

### MultilayerShieldCalculator

Class for multilayer (2-layer) heat shield analysis.

#### solve_temperatures()

Solves steady-state temperatures for both shield layers.

**Signature:**
```python
@staticmethod
def solve_temperatures(
    t_exh, t_fluid, t_surr, h_in, h_out, h_gap,
    eps_in, eps_out, eps_g1, eps_g2,
    f12=1.0, tol=0.1, max_iter=100
) -> dict
```

**Parameters:**
- `t_exh` (float): Exhaust temperature (K)
- `t_fluid` (float): Cooling fluid temperature (K)
- `t_surr` (float): Surroundings temperature (K)
- `h_in` (float): Hot-side convection (W/m²·K)
- `h_out` (float): Coolant-side convection (W/m²·K)
- `h_gap` (float): Gap radiation conductance (W/m²·K)
- `eps_in` (float): Inner layer hot-side emissivity
- `eps_out` (float): Outer layer cold-side emissivity
- `eps_g1` (float): Inner layer gap-side emissivity
- `eps_g2` (float): Outer layer gap-side emissivity
- `f12` (float, default=1.0): View factor between the gap faces (0–1; 1.0 is infinite parallel plates)
- `tol` (float, default=0.1): Convergence tolerance on both energy balances (W/m²)
- `max_iter` (int, default=100): Maximum iterations

**Returns:**
- `dict` with keys:
  - `t1_K` (float): Inner layer temperature (K)
  - `t1_C` (float): Inner layer temperature (°C)
  - `t2_K` (float): Outer layer temperature (K)
  - `t2_C` (float): Outer layer temperature (°C)
  - `delta_T_C` (float): Temperature difference (°C)
  - `converged` (bool): Convergence status; when False, the temperatures and fluxes are those of the last iterate
  - `iterations` (int): Iteration count
  - `eps_eff` (float): Effective emissivity
  - `residual_W_m2` (float): Larger of the two energy-balance residuals (W/m²)
  - All six flux components (`q_rad_in`, `q_conv_in`, `q_gap_cond`, `q_gap_rad`, `q_rad_out`, `q_conv_out`)

**Example:**
```python
from thermal_mesh_calculators.shields import MultilayerShieldCalculator

result = MultilayerShieldCalculator.solve_temperatures(
    t_exh=873.15,      # K
    t_fluid=343.15,    # K
    t_surr=343.15,     # K
    h_in=150.0,        # W/m²·K
    h_out=100.0,       # W/m²·K
    h_gap=50.0,        # W/m²·K
    eps_in=0.8,
    eps_out=0.8,
    eps_g1=0.8,
    eps_g2=0.8
)

print(f"Layer 1: {result['t1_C']:.1f}°C")
print(f"Layer 2: {result['t2_C']:.1f}°C")
print(f"Delta T: {result['delta_T_C']:.1f}°C")
```

#### mesh_sizes()

Calculates mesh constraints for both layers.

**Signature:**
```python
@staticmethod
def mesh_sizes(k_metal, max_dt, **kwargs) -> dict
```

**Parameters:**
- `k_metal` (float): Thermal conductivity (W/m·K)
- `max_dt` (float): Maximum temperature difference across one element (K)
- `**kwargs`: All parameters for `solve_temperatures()`

**Returns:**
- `dict` with all keys from `solve_temperatures()` plus:
  - `layer1_max_dx_mm` (float): Maximum element size for layer 1 (mm)
  - `layer2_max_dx_mm` (float): Maximum element size for layer 2 (mm)
  - `q_layer1` (float): Layer 1 driving flux, `|q_rad_in| + |q_conv_in|` (W/m²)
  - `q_layer2` (float): Layer 2 driving flux, `|q_rad_out| + |q_conv_out|` (W/m²)

**Example:**
```python
result = MultilayerShieldCalculator.mesh_sizes(
    k_metal=50.0,
    max_dt=0.1,
    t_exh=873.15,
    t_fluid=343.15,
    t_surr=343.15,
    h_in=150.0,
    h_out=100.0,
    h_gap=50.0,
    eps_in=0.8,
    eps_out=0.8,
    eps_g1=0.8,
    eps_g2=0.8
)

print(f"Layer 1 max element: {result['layer1_max_dx_mm']:.2f} mm")
print(f"Layer 2 max element: {result['layer2_max_dx_mm']:.2f} mm")
```

---

## transient

Transient thermal analysis and Fourier number constraints.

### TransientMeshCalculator

Class for transient mesh sizing based on material properties and time scales.

#### material_diffusivity()

Calculates thermal diffusivity for a material.

**Signature:**
```python
@staticmethod
def material_diffusivity(k, rho, cp) -> float
```

**Parameters:**
- `k` (float): Thermal conductivity (W/m·K)
- `rho` (float): Density (kg/m³)
- `cp` (float): Specific heat (J/kg·K)

**Returns:**
- `float`: Thermal diffusivity (m²/s)

**Example:**
```python
from thermal_mesh_calculators.transient import TransientMeshCalculator

alpha = TransientMeshCalculator.material_diffusivity(
    k=50.0,        # W/m·K
    rho=7850.0,    # kg/m³ (steel)
    cp=500.0       # J/kg·K
)
print(f"Thermal diffusivity: {alpha:.2e} m²/s")
```

#### penetration_depth()

Calculates maximum element size to resolve thermal penetration.

**Signature:**
```python
@staticmethod
def penetration_depth(k, rho, cp, dt, safety_factor=1.0) -> dict
```

**Parameters:**
- `k` (float): Thermal conductivity (W/m·K)
- `rho` (float): Density (kg/m³)
- `cp` (float): Specific heat (J/kg·K)
- `dt` (float): Time step (s)
- `safety_factor` (float, default=1.0): Multiplier C on √(α·dt): a resolution guideline for implicit schemes (1.5–2.0). A larger C allows a coarser element; it is not a stability limit

**Returns:**
- `dict` with keys:
  - `max_dx_mm` (float): Maximum element size (mm)
  - `alpha` (float): Thermal diffusivity (m²/s)
  - `penetration_m` (float): Penetration depth (m)

**Example:**
```python
result = TransientMeshCalculator.penetration_depth(
    k=50.0,
    rho=7850.0,
    cp=500.0,
    dt=0.1,
    safety_factor=1.0
)

print(f"Max element: {result['max_dx_mm']:.3f} mm")
print(f"Penetration: {result['penetration_m']*1000:.2f} mm")
```

#### fourier_number_limit()

Calculates minimum element size based on Fourier number limit.

**Signature:**
```python
@staticmethod
def fourier_number_limit(k, rho, cp, dt, fo_max=0.5) -> dict
```

**Parameters:**
- `k` (float): Thermal conductivity (W/m·K)
- `rho` (float): Density (kg/m³)
- `cp` (float): Specific heat (J/kg·K)
- `dt` (float): Time step (s)
- `fo_max` (float, default=0.5): Maximum Fourier number

**Returns:**
- `dict` with keys:
  - `min_dx_mm` (float): Minimum element size (mm)
  - `alpha` (float): Thermal diffusivity (m²/s)
  - `fo_max` (float): Maximum Fourier number used
  - `constraint` (str): Description of constraint

**Example:**
```python
result = TransientMeshCalculator.fourier_number_limit(
    k=50.0,
    rho=7850.0,
    cp=500.0,
    dt=0.1,
    fo_max=0.5
)

print(f"Min element: {result['min_dx_mm']:.3f} mm")
```

#### drive_cycle_limit()

Calculates element size based on boundary condition time scale.

**Signature:**
```python
@staticmethod
def drive_cycle_limit(k, rho, cp, tau_bc) -> dict
```

**Parameters:**
- `k` (float): Thermal conductivity (W/m·K)
- `rho` (float): Density (kg/m³)
- `cp` (float): Specific heat (J/kg·K)
- `tau_bc` (float): Boundary condition time scale (s)

**Returns:**
- `dict` with keys:
  - `max_dx_mm` (float): Maximum element size (mm)
  - `alpha` (float): Thermal diffusivity (m²/s)
  - `tau_bc` (float): Time scale used (s)

**Example:**
```python
result = TransientMeshCalculator.drive_cycle_limit(
    k=50.0,
    rho=7850.0,
    cp=500.0,
    tau_bc=10.0      # 10 second transient
)

print(f"Max element: {result['max_dx_mm']:.2f} mm")
```

#### combined_transient_limits()

Computes all transient constraints and recommends element size.

**Signature:**
```python
@classmethod
def combined_transient_limits(
    k, rho, cp, dt, fo_max=0.5, safety_factor=1.0, tau_bc=None, scheme=None
) -> dict
```

**Parameters:**
- `k` (float): Thermal conductivity (W/m·K)
- `rho` (float): Density (kg/m³)
- `cp` (float): Specific heat (J/kg·K)
- `dt` (float): Time step (s)
- `fo_max` (float, default=0.5): Maximum Fourier number
- `safety_factor` (float, default=1.0): Penetration-depth multiplier C, applied to implicit schemes only
- `tau_bc` (float, optional): Boundary condition time scale (s)
- `scheme` (str, optional): `"explicit"` or `"implicit"`; `None` infers it from `fo_max` (≤ 0.5 is explicit)

**Returns:**
- `dict` with keys:
  - `alpha` (float): Thermal diffusivity (m²/s)
  - `scheme` (str): `"explicit"` or `"implicit"`
  - `penetration_max_dx_mm` (float): Limit from penetration depth
  - `penetration_applied` (bool): False for an explicit scheme (see `docs/math_derivations.md` §7.4)
  - `fourier_min_dx_mm` (float): Limit from Fourier number (a lower bound)
  - `drive_cycle_max_dx_mm` (float): Limit from time scale (if tau_bc given)
  - `max_dx_mm` (float): Tightest applicable upper bound (`inf` when none applies)
  - `feasible` (bool): `fourier_min_dx_mm <= max_dx_mm`
  - `recommended_dx_mm` (float): `max_dx_mm` when feasible, else `fourier_min_dx_mm`
  - `binding_constraint` (str): `"penetration_depth"`, `"drive_cycle"`, `"none"`, or a `"fourier_... (WARNING: ...)"` conflict label
  - `conflict` (dict or None): `message` and `remedies`, every parameter change the conflict needs (`dt` maximum, or `safety_factor` minimum). `remedies` holds the exact values; the message prints each bound rounded toward the side that satisfies it, so the printed value closes the conflict too
  - `advice` (str): The result in words

**Example:**
```python
result = TransientMeshCalculator.combined_transient_limits(
    k=50.0,
    rho=7850.0,
    cp=500.0,
    dt=0.1,
    fo_max=0.5,
    safety_factor=1.0,
    tau_bc=10.0
)

print(f"Recommended element: {result['recommended_dx_mm']:.2f} mm")
print(f"Binding constraint: {result['binding_constraint']}")
```

---

## zones

Pre-defined convection zones for automotive thermal analysis.

### CONVECTION_ZONES

Dictionary of 16 predefined zones with convection characteristics.

**Structure:**
Each zone contains keys:
- `h_low` (float): Low estimate for h (W/m²·K)
- `h_high` (float): High estimate for h (W/m²·K)
- `regime` (str): Flow regime description
- `velocity_ms` (float): Characteristic velocity (m/s)
- `t_air_C_low` (float): Low air temperature (°C)
- `t_air_C_high` (float): High air temperature (°C)
- `orientation` (str): Surface orientation
- `is_internal` (bool): Whether internal or external
- `notes` (str): Additional context

**Zones Include:**
exhaust_pipe, underhood_metal, underhood_shield, catalytic_converter, turbo_housing, turbo_scroll, intake_manifold, muffler, radiator, transmission_cooler, drivetrain_housing, suspension_component, internal_fluid_surface, internal_air_gap, vehicle_exterior, thermal_break_surface

**Example:**
```python
from thermal_mesh_calculators.zones import CONVECTION_ZONES

zone_data = CONVECTION_ZONES['exhaust_pipe']
print(f"h range: {zone_data['h_low']:.0f} to {zone_data['h_high']:.0f} W/m²·K")
print(f"Velocity: {zone_data['velocity_ms']:.1f} m/s")
```

### SPATIAL_GRADIENT_DEFAULTS

Dictionary of default spatial temperature gradients by component class.

### get_zone()

Retrieves a zone by name.

**Signature:**
```python
def get_zone(zone_name: str) -> dict
```

**Parameters:**
- `zone_name` (str): Name of the zone (case-insensitive)

**Returns:**
- `dict`: Zone data

**Example:**
```python
from thermal_mesh_calculators.zones import get_zone

zone = get_zone("exhaust_pipe")
print(f"Zone: {zone}")
```

### list_zones()

Lists all available zone names.

**Signature:**
```python
def list_zones() -> list
```

**Returns:**
- `list`: Names of all zones

**Example:**
```python
from thermal_mesh_calculators.zones import list_zones

zones = list_zones()
print(f"Available zones: {zones}")
```

### get_conservative_h()

Gets the conservative (high) h estimate for a zone.

**Signature:**
```python
def get_conservative_h(zone_name: str) -> float
```

**Parameters:**
- `zone_name` (str): Zone name

**Returns:**
- `float`: Conservative h value (W/m²·K)

**Example:**
```python
from thermal_mesh_calculators.zones import get_conservative_h

h = get_conservative_h("exhaust_pipe")
print(f"Conservative h: {h:.0f} W/m²·K")
```

### get_zone_air_temp()

Gets air/fluid temperature for a zone.

**Signature:**
```python
def get_zone_air_temp(zone_name: str, bound: str = "high") -> float
```

**Parameters:**
- `zone_name` (str): Zone name
- `bound` (str, default="high"): "high" or "low"

**Returns:**
- `float`: Temperature (K)

**Example:**
```python
from thermal_mesh_calculators.zones import get_zone_air_temp

t_high = get_zone_air_temp("exhaust_pipe", bound="high")
t_low = get_zone_air_temp("exhaust_pipe", bound="low")
print(f"Air temp range: {t_low-273.15:.0f}°C to {t_high-273.15:.0f}°C")
```

### estimate_spatial_gradient()

Estimates temperature gradient in a component.

**Signature:**
```python
def estimate_spatial_gradient(k, q_total=None, component_class: str = "structural") -> float
```

**Parameters:**
- `k` (float): Thermal conductivity (W/m·K)
- `q_total` (float, optional): Heat flux (W/m²)
- `component_class` (str, default="structural"): Component classification

**Returns:**
- `float`: Spatial temperature gradient (K/m)

**Example:**
```python
from thermal_mesh_calculators.zones import estimate_spatial_gradient

gradient = estimate_spatial_gradient(
    k=50.0,
    q_total=50000.0,
    component_class="structural"
)
print(f"Spatial gradient: {gradient:.0f} K/m")
```

---

## h_estimator

Heat transfer coefficient estimation from convection correlations.

### Module Constants

**H_EXTERNAL_MAX** (float): Maximum capped h value = 100.0 W/m²·K

**H_EXTERNAL_TYPICAL** (float): Typical capped h value = 60.0 W/m²·K

### air_properties()

Calculates thermophysical properties of air at film temperature.

**Signature:**
```python
def air_properties(t_film: float) -> dict
```

**Parameters:**
- `t_film` (float): Film temperature (K)

**Returns:**
- `dict` with keys:
  - `k_air` (float): Thermal conductivity (W/m·K)
  - `nu` (float): Kinematic viscosity (m²/s)
  - `Pr` (float): Prandtl number
  - `beta` (float): Volumetric expansion coefficient (1/K)
  - `rho` (float): Density (kg/m³)
  - `t_film_K` (float): The film temperature given (K)
  - `t_eval_K` (float): The temperature the fits were evaluated at: `t_film` clamped to `AIR_PROPERTY_RANGE_K` (250–700 K)
  - `in_range` (bool): False when the properties are extrapolated

**Example:**
```python
from thermal_mesh_calculators.h_estimator import air_properties

props = air_properties(t_film=400.0)  # K
print(f"k_air: {props['k_air']:.4f} W/m·K")
print(f"Pr: {props['Pr']:.2f}")
```

### forced_convection_flat_plate()

Estimates h for forced convection over a flat plate.

**Signature:**
```python
def forced_convection_flat_plate(
    velocity, char_length, t_surf, t_fluid,
    cap=H_EXTERNAL_MAX
) -> dict
```

**Parameters:**
- `velocity` (float): Free stream velocity (m/s)
- `char_length` (float): Characteristic length (m)
- `t_surf` (float): Surface temperature (K)
- `t_fluid` (float): Fluid temperature (K)
- `cap` (float, default=100.0): Maximum h cap (W/m²·K)

**Returns:**
- `dict` with keys:
  - `h` (float): Convection coefficient (W/m²·K)
  - `h_raw` (float): Uncapped value
  - `Re` (float): Reynolds number
  - `Nu` (float): Nusselt number
  - `regime` (str): "laminar", "turbulent", or "transition"

**Example:**
```python
from thermal_mesh_calculators.h_estimator import forced_convection_flat_plate

result = forced_convection_flat_plate(
    velocity=30.0,         # m/s (108 km/h)
    char_length=0.1,       # m (100 mm)
    t_surf=450.0,          # K
    t_fluid=300.0          # K
)

print(f"h: {result['h']:.0f} W/m²·K")
print(f"Re: {result['Re']:.0f}")
print(f"Regime: {result['regime']}")
```

### natural_convection()

Estimates h for natural (free) convection.

**Signature:**
```python
def natural_convection(
    t_surf, t_fluid, char_length,
    orientation="vertical", cap=H_EXTERNAL_MAX
) -> dict
```

**Parameters:**
- `t_surf` (float): Surface temperature (K)
- `t_fluid` (float): Ambient temperature (K)
- `char_length` (float): Characteristic length (m)
- `orientation` (str, default="vertical"): "vertical" or "horizontal"
- `cap` (float, default=100.0): Maximum h cap (W/m²·K)

**Returns:**
- `dict` with keys:
  - `h` (float): Convection coefficient (W/m²·K)
  - `h_raw` (float): Uncapped value
  - `Ra` (float): Rayleigh number
  - `Gr` (float): Grashof number
  - `Nu` (float): Nusselt number
  - `regime` (str): "laminar", "turbulent", or "transition"

**Example:**
```python
from thermal_mesh_calculators.h_estimator import natural_convection

result = natural_convection(
    t_surf=450.0,          # K
    t_fluid=300.0,         # K
    char_length=0.05,      # m (50 mm)
    orientation="vertical"
)

print(f"h: {result['h']:.0f} W/m²·K")
print(f"Ra: {result['Ra']:.2e}")
```

### richardson_number()

Calculates Richardson number and flow regimes.

**Signature:**
```python
def richardson_number(velocity, t_surf, t_fluid, char_length) -> dict
```

**Parameters:**
- `velocity` (float): Free stream velocity (m/s)
- `t_surf` (float): Surface temperature (K)
- `t_fluid` (float): Fluid temperature (K)
- `char_length` (float): Characteristic length (m)

**Returns:**
- `dict` with keys:
  - `Ri` (float): Richardson number
  - `Gr` (float): Grashof number
  - `Re` (float): Reynolds number
  - `regime` (str): "forced_dominant", "mixed", or "natural_dominant"

**Example:**
```python
from thermal_mesh_calculators.h_estimator import richardson_number

result = richardson_number(
    velocity=5.0,          # m/s (18 km/h)
    t_surf=450.0,          # K
    t_fluid=300.0,         # K
    char_length=0.1        # m
)

print(f"Ri: {result['Ri']:.3f}")
print(f"Regime: {result['regime']}")
```

### estimate_h()

Comprehensive h estimation considering both forced and natural convection.

**Signature:**
```python
def estimate_h(
    velocity, t_surf, t_fluid, char_length,
    orientation="vertical", cap=H_EXTERNAL_MAX
) -> dict
```

**Parameters:**
- `velocity` (float): Free stream velocity (m/s)
- `t_surf` (float): Surface temperature (K)
- `t_fluid` (float): Fluid temperature (K)
- `char_length` (float): Characteristic length (m)
- `orientation` (str, default="vertical"): Surface orientation
- `cap` (float, default=100.0): Maximum h cap (W/m²·K)

**Returns:**
- `dict` with keys:
  - `h` (float): Recommended h (W/m²·K)
  - `h_forced` (float): Forced convection component
  - `h_natural` (float): Natural convection component
  - `regime` (str): Overall regime
  - `natural_regime` (str): Natural convection regime
  - `forced_regime` (str): Forced convection regime
  - `Ri` (float): Richardson number
  - `Ra` (float): Rayleigh number
  - `Re` (float): Reynolds number
  - `dominant_mode` (str): "forced", "mixed", or "natural"
  - `solver_advisory` (dict): Solver recommendations
  - `t_film_K` (float): Film temperature (K)
  - `film_in_range` (bool): False when the air properties behind h were extrapolated

**Example:**
```python
from thermal_mesh_calculators.h_estimator import estimate_h

result = estimate_h(
    velocity=10.0,         # m/s
    t_surf=450.0,          # K
    t_fluid=300.0,         # K
    char_length=0.05,      # m
    orientation="vertical"
)

print(f"h: {result['h']:.0f} W/m²·K")
print(f"Dominant mode: {result['dominant_mode']}")
print(f"Solver advisory: {result['solver_advisory']}")
```

### solver_advisory()

Recommends solver settings based on convection regime.

**Signature:**
```python
def solver_advisory(
    forced_result=None, natural_result=None,
    regime="natural", orientation="vertical"
) -> dict
```

**Parameters:**
- `forced_result` (dict, optional): Result from forced_convection_flat_plate()
- `natural_result` (dict, optional): Result from natural_convection()
- `regime` (str, default="natural"): Overall regime
- `orientation` (str, default="vertical"): Surface orientation

**Returns:**
- `dict` with keys:
  - `steady_state_ok` (bool): Whether steady-state solve is adequate
  - `transient_advisory` (str): Transient analysis recommendation
  - `reason` (str): Explanation
  - `severity` (str): "low", "medium", or "high"
  - `natural_regime` (str): Natural convection regime
  - `forced_regime` (str): Forced convection regime
  - `Ra` (float): Rayleigh number
  - `Re` (float): Reynolds number

**Example:**
```python
from thermal_mesh_calculators.h_estimator import solver_advisory

advisory = solver_advisory(
    regime="mixed",
    orientation="vertical"
)

print(f"Steady-state OK: {advisory['steady_state_ok']}")
print(f"Advisory: {advisory['transient_advisory']}")
```

### h_from_velocity()

Quick convenience function for h estimation from velocity in mph.

**Signature:**
```python
def h_from_velocity(
    velocity_mph, char_length_mm=100.0,
    t_surf_C=100.0, t_fluid_C=40.0
) -> dict
```

**Parameters:**
- `velocity_mph` (float): Velocity (miles per hour)
- `char_length_mm` (float, default=100.0): Characteristic length (mm)
- `t_surf_C` (float, default=100.0): Surface temperature (°C)
- `t_fluid_C` (float, default=40.0): Fluid temperature (°C)

**Returns:**
- `dict`: Same as estimate_h()

**Example:**
```python
from thermal_mesh_calculators.h_estimator import h_from_velocity

result = h_from_velocity(
    velocity_mph=60.0,      # 60 mph highway
    char_length_mm=100.0,   # 100 mm part
    t_surf_C=150.0,         # 150°C surface
    t_fluid_C=40.0          # 40°C ambient
)

print(f"h at 60 mph: {result['h']:.0f} W/m²·K")
```

---

## batch

Batch processing for complete part analysis and materials database.

### MATERIALS

Dictionary of 44 materials with properties.

Keys include: steel_301, aluminum_6061, titanium_grade2, nickel, copper, stainless_304, stainless_316, stainless_321, superalloy_inconel_718, and 35 others.

**Each entry contains:**
- `k` (float): Thermal conductivity (W/m·K)
- `rho` (float): Density (kg/m³)
- `cp` (float): Specific heat (J/kg·K)
- `common_name` (str): Human-readable name
- `temp_range_C` (tuple): Applicable temperature range

**Example:**
```python
from thermal_mesh_calculators.batch import MATERIALS

steel = MATERIALS['stainless_316']
print(f"Conductivity: {steel['k']:.1f} W/m·K")
```

### SURFACE_TREATMENTS

Dictionary of 30 surface treatments with emissivity data.

Keys include: oxidized_steel, bare_aluminum, ceramic_coating, anodized_aluminum, polished_steel, painted_black, and 24 others.

**Each entry contains:**
- `emissivity` (float): Hemispherical emissivity (0–1)
- `description` (str): Treatment description
- `typical_temps_C` (list): Applicable temperatures

**Example:**
```python
from thermal_mesh_calculators.batch import SURFACE_TREATMENTS

coating = SURFACE_TREATMENTS['ceramic_coating_high_temp']
print(f"Emissivity: {coating['emissivity']:.2f}")
```

### get_material()

Retrieves material properties by name.

**Signature:**
```python
def get_material(name: str) -> dict
```

**Parameters:**
- `name` (str): Material name (case-insensitive)

**Returns:**
- `dict`: Material properties

**Example:**
```python
from thermal_mesh_calculators.batch import get_material

mat = get_material("stainless_316")
print(mat)
```

### list_materials()

Lists all available material names.

**Signature:**
```python
def list_materials() -> list
```

**Returns:**
- `list`: Material names

**Example:**
```python
from thermal_mesh_calculators.batch import list_materials

mats = list_materials()
print(f"Available materials: {len(mats)}")
```

### get_surface_epsilon()

Retrieves emissivity for a surface treatment.

**Signature:**
```python
def get_surface_epsilon(treatment: str) -> float
```

**Parameters:**
- `treatment` (str): Treatment name

**Returns:**
- `float`: Emissivity (0–1)

**Example:**
```python
from thermal_mesh_calculators.batch import get_surface_epsilon

eps = get_surface_epsilon("oxidized_steel")
print(f"Emissivity: {eps:.2f}")
```

### list_surface_treatments()

Lists all available surface treatments.

**Signature:**
```python
def list_surface_treatments() -> list
```

**Returns:**
- `list`: Treatment names

**Example:**
```python
from thermal_mesh_calculators.batch import list_surface_treatments

treatments = list_surface_treatments()
print(f"Available treatments: {len(treatments)}")
```

### process_part()

Runs every applicable calculator for one part and returns the governing (smallest)
element size. The module docstring of `thermal_mesh_calculators.batch` lists every
key a part or project may carry.

**Signature:**
```python
def process_part(part: dict, project: dict) -> dict
```

**Parameters:**

`part` dict, required keys:
- `part_id` (str): Part identifier
- `material` (str): A key of `MATERIALS`
- `component_class` (str): `"exhaust"`, `"exhaust_adjacent"`, `"structural"`, `"shield"` or `"multilayer_shield"` (the keys of `CLASS_DEFAULTS`)
- `convection_zone` (str): A key of `CONVECTION_ZONES`; shields may give `convection_zone_in` and `convection_zone_out` instead
- `thickness_mm` (float): Wall thickness (mm)
- `t_surf_K` (float): Surface temperature (K); ignored for the shield classes, whose temperature is solved

`part` dict, surface keys: `surface` (a key of `SURFACE_TREATMENTS`) or `epsilon` for the non-shield classes; `surface_in` / `eps_in` and `surface_out` / `eps_out` for shields; `surface_g1` / `eps_g1` and `surface_g2` / `eps_g2` for the gap faces of a two-layer shield. A surface given neither falls back to a material-class emissivity and is reported as a `SURFACE_DEFAULTED` warning.

`part` dict, optional keys (among others): `h_override`, `h_in_override`, `h_out_override`, `char_length_mm` (default 100), `t_fluid_K` (this part's fluid temperature), `t_exh_K`, `h_gap`, `f12`, `shield_max_iter`.

`project` dict: `t_fluid_K` and `t_surr_K` (K); optionally `max_dt` (K per element), `dt`, `fo_max`, `tau_bc`, `safety_factor`, `transient_scheme`, `allowable_flux_error`, `t_exh_K`, `shield_max_iter`.

**Returns:**
- `dict` with keys:
  - `part_id`, `material`, `component_class`: echoed back
  - `t_fluid_K` (float), `t_fluid_source` (str): the fluid temperature used for h and q″, and `"part"` or `"project"`
  - `h_used`, `eps_used`: the h (W/m²·K) and emissivity used; for shields, a dict per face
  - `h_estimation` (dict): The h estimate (`h`, `method`, `regime`, `details`); `None` for shields
  - `solver_advisory` (dict): Steady-state vs. transient advice; `None` for shields
  - `conduction`, `biot`, `radiation`, `lateral`, `shield`, `transient` (dict or `None`): Each calculator's result
  - `governing_dx_mm` (float): The most restrictive element size (mm)
  - `governing_constraint` (str): The constraint that produced it
  - `all_constraints` (list): Every candidate as `{"dx_mm": ..., "source": ...}`
  - `warnings` (list): Coded warnings, each with `code`, `severity` and `message` (input warnings also carry `key`)

**Raises:**
- `PartInputError` (a `KeyError` and a `ValueError`): before anything is computed, for a missing or unknown material, component class or convection zone, or an unknown surface treatment; the message names the part, the key and the allowed set

**Example:**
```python
from thermal_mesh_calculators.batch import process_part

part = {
    "part_id": "BRK-001",
    "material": "steel_mild",
    "component_class": "structural",
    "convection_zone": "engine_beside",
    "thickness_mm": 3.0,
    "t_surf_K": 473.15,
    "surface": "painted",
}
project = {"t_fluid_K": 353.15, "t_surr_K": 353.15}

result = process_part(part, project)
print(f"{result['part_id']}: {result['governing_dx_mm']:.1f} mm ({result['governing_constraint']})")
# BRK-001: 24.5 mm (lateral_gradient)
```

### process_batch()

Analyzes multiple parts.

**Signature:**
```python
def process_batch(parts: list, project: dict) -> list
```

**Parameters:**
- `parts` (list): List of part dicts (same format as process_part)
- `project` (dict): Project parameters (same as process_part)

**Returns:**
- `list`: List of result dicts from process_part()

**Example:**
```python
from thermal_mesh_calculators.batch import process_batch

parts = [
    {
        'id': 'shield_1',
        'material': 'stainless_316',
        'component_class': 'structural',
        'thickness_mm': 1.0,
        'characteristic_length_mm': 100.0,
        'zone': 'exhaust_pipe'
    },
    {
        'id': 'shield_2',
        'material': 'aluminum_6061',
        'component_class': 'structural',
        'thickness_mm': 2.0,
        'characteristic_length_mm': 80.0,
        'zone': 'underhood_metal'
    }
]

project = {
    't_exh_K': 873.15,
    't_fluid_K': 343.15,
    't_surr_K': 343.15,
    'epsilon_in': 0.8,
    'epsilon_out': 0.8,
    'max_dt': 0.1
}

results = process_batch(parts, project)

for result in results:
    print(f"{result['part_id']}: {result['governing_dx_mm']:.2f} mm")
```

### summary_table()

Generates a formatted summary table from batch results.

**Signature:**
```python
def summary_table(results: list) -> str
```

**Parameters:**
- `results` (list): List of results from process_batch()

**Returns:**
- `str`: Formatted table suitable for printing or logging

**Example:**
```python
from thermal_mesh_calculators.batch import process_batch, summary_table

results = process_batch(parts, project)
table = summary_table(results)
print(table)
```

---

## Quick Start Examples

### Example 1: Simple Single-Part Conduction Calculation

```python
from thermal_mesh_calculators.conduction import BoundaryDrivenConductionCalculator

k = 50.0              # Stainless steel, W/m·K
h = 75.0              # Convection coefficient, W/m²·K
t_surf = 873.15       # Surface temp, K (600°C)
t_fluid = 343.15      # Fluid temp, K (70°C)
epsilon = 0.7         # Emissivity
t_surr = 343.15       # Surroundings, K
max_dt = 0.1          # Max time step, s

result = BoundaryDrivenConductionCalculator.max_mesh_size(
    k, h, t_surf, t_fluid, epsilon, t_surr, max_dt
)

print(f"Thermal conductivity: {k} W/m·K")
print(f"Max mesh element: {result['max_dx_mm']:.2f} mm")
print(f"Heat flux: {result['q_total']:.0f} W/m²")
print(f"Radiation fraction: {result['rad_fraction']*100:.1f}%")
```

### Example 2: Biot Number Check

```python
from thermal_mesh_calculators.convection import ConvectionMeshCalculator

cases = [
    {"name": "Thin wall (1 mm)", "h": 100, "k": 50, "thickness": 0.001},
    {"name": "Medium wall (5 mm)", "h": 100, "k": 50, "thickness": 0.005},
    {"name": "Thick wall (20 mm)", "h": 100, "k": 50, "thickness": 0.020},
]

for case in cases:
    result = ConvectionMeshCalculator.biot_number(
        case["h"], case["k"], case["thickness"]
    )
    print(f"\n{case['name']}")
    print(f"  Biot = {result['biot']:.4f}")
    print(f"  Mesh type: {result['mesh_type']}")
    print(f"  Min elements: {result['min_elements_through_thickness']}")
```

### Example 3: Shield Temperature Solve

```python
from thermal_mesh_calculators.shields import SingleLayerShieldCalculator

t_exh = 873.15       # Exhaust, K (600°C)
t_fluid = 343.15     # Coolant, K (70°C)
t_surr = 343.15      # Surroundings, K (70°C)
epsilon = 0.8        # Both surfaces
h_in = 150.0         # Hot-side convection, W/m²·K
h_out = 100.0        # Coolant-side convection, W/m²·K

result = SingleLayerShieldCalculator.solve_temperature(
    t_exh, t_fluid, t_surr, epsilon, epsilon,
    h_in=h_in, h_out=h_out
)

print(f"Exhaust temperature: {t_exh - 273.15:.0f}°C")
print(f"Shield temperature: {result['t_shield_C']:.1f}°C")
print(f"Coolant temperature: {t_fluid - 273.15:.0f}°C")
print(f"Thermal relief: {t_exh - 273.15 - result['t_shield_C']:.0f}°C")
print(f"Converged: {result['converged']} ({result['iterations']} iterations)")
```

### Example 4: Batch Processing

```python
from thermal_mesh_calculators.batch import process_batch, summary_table

parts = [
    {
        'id': 'outer_shield',
        'material': 'stainless_316',
        'component_class': 'structural',
        'thickness_mm': 1.0,
        'characteristic_length_mm': 150.0,
        'zone': 'exhaust_pipe'
    },
    {
        'id': 'inner_shield',
        'material': 'stainless_304',
        'component_class': 'structural',
        'thickness_mm': 0.8,
        'characteristic_length_mm': 140.0,
        'zone': 'exhaust_pipe'
    },
    {
        'id': 'bracket',
        'material': 'aluminum_6061',
        'component_class': 'structural',
        'thickness_mm': 3.0,
        'characteristic_length_mm': 50.0,
        'zone': 'underhood_metal'
    }
]

project = {
    't_exh_K': 873.15,
    't_fluid_K': 343.15,
    't_surr_K': 343.15,
    'epsilon_in': 0.8,
    'epsilon_out': 0.8,
    'max_dt': 0.1
}

results = process_batch(parts, project)
print(summary_table(results))

for result in results:
    print(f"\n{result['part_id']}:")
    print(f"  Mesh element: {result['governing_dx_mm']:.2f} mm")
    print(f"  Binding constraint: {result['governing_constraint']}")
```

### Example 5: Convection h Estimation from Velocity

```python
from thermal_mesh_calculators.h_estimator import estimate_h, h_from_velocity

result = estimate_h(
    velocity=30.0,         # m/s (108 km/h)
    t_surf=450.0,          # K (177°C)
    t_fluid=300.0,         # K (27°C)
    char_length=0.1        # m (100 mm)
)

print("Method 1: SI units")
print(f"  Velocity: 30 m/s (108 km/h)")
print(f"  h: {result['h']:.0f} W/m²·K")
print(f"  Re: {result['Re']:.0f}")
print(f"  Regime: {result['dominant_mode']}")

result2 = h_from_velocity(
    velocity_mph=67.0,      # ~108 km/h
    char_length_mm=100.0,
    t_surf_C=177.0,
    t_fluid_C=27.0
)

print("\nMethod 2: Convenience (mph)")
print(f"  Velocity: 67 mph")
print(f"  h: {result2['h']:.0f} W/m²·K")
```

### Example 6: Quick Velocity Gut-Check

```python
from thermal_mesh_calculators.h_estimator import h_from_velocity

speeds_mph = [20, 35, 50, 65, 80]

print("Quick h estimates (100 mm part, 150°C surface, 40°C ambient)")
print("-" * 50)

for speed in speeds_mph:
    result = h_from_velocity(
        velocity_mph=speed,
        char_length_mm=100.0,
        t_surf_C=150.0,
        t_fluid_C=40.0
    )
    print(f"{speed:3d} mph -> h = {result['h']:5.0f} W/m²·K")

print("\nScaling note:")
print("  h scales proportionally to velocity^0.5 for laminar flow")
print("  h scales proportionally to velocity^0.8 for turbulent flow")
```

---

## Summary

The thermal-mesh-calculators package provides comprehensive thermal analysis and mesh sizing for heat shield and conduction applications. Key capabilities include:

- Steady-state analysis: Shield temperatures, flux distributions
- Transient analysis: Fourier number limits, penetration depth, time scales
- Convection modeling: Biot numbers, h estimation, zone databases
- Radiation effects: View factors, flux sensitivity, emissivity
- Batch processing: Multi-part analysis with unified mesh recommendations
- Material and surface databases: 44 materials and 30 surface treatments

All functions are pure Python with no external dependencies.
