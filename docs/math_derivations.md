# Mathematical Derivations for Thermal Mesh Calculators

**Audience:** Thermal simulation engineers
**Purpose:** Comprehensive derivations from governing equations through implementation details

---

## 1. Boundary-Driven Conduction Mesh Sizing

### 1.1 Fourier's Law and Discretization

The fundamental heat conduction equation is Fourier's law:

```
q'' = -k dT/dx
```

where:
- `q''` = heat flux (W/m²)
- `k` = thermal conductivity (W/m·K)
- `dT/dx` = temperature gradient (K/m)

For finite differences across a characteristic length `Δx_max`, taking magnitude:

```
|q''| ≈ k · ΔT_max / Δx_max
```

Solving for the maximum element size:

```
Δx_max = k · ΔT_max / |q''|
```

### 1.2 The Challenge: q'' is Unknown

In a thermal simulation, the heat flux at the surface is not prescribed
independently — it is a **solver output**. The formulation above is circular:
we need q'' to compute mesh size, but we need the mesh to compute q''.

### 1.3 Surface Energy Balance — Eliminating q''

The key insight is to substitute the surface energy balance for q''. At
steady state, the conduction flux through the solid at the boundary must
equal the sum of convective and radiative losses from the surface:

```
q''_cond = q''_conv + q''_rad
```

Expanding each term:

```
q''_conv = h(T_s - T_∞)         [Newton's law of cooling]
q''_rad  = εσ(T_s^4 - T_surr^4)  [Stefan-Boltzmann net radiation]
```

where:
- `h` = convection coefficient (W/m²·K)
- `T_s` = surface temperature (K)
- `T_∞` = fluid ambient temperature (K)
- `ε` = emissivity (0–1)
- `σ` = Stefan-Boltzmann constant = 5.67 × 10⁻⁸ W/m²·K⁴
- `T_surr` = surroundings radiation sink temperature (K)

Therefore:

```
q'' = h(T_s - T_∞) + εσ(T_s^4 - T_surr^4)
```

All terms on the right-hand side are either known inputs (material properties,
boundary conditions) or can be reasonably estimated (surface temperature).

**Implementation:** `conduction.py` lines 94–96:
```python
q_conv = h * (t_surf - t_fluid)
q_rad = epsilon * STEFAN_BOLTZMANN * (t_surf**4 - t_surr**4)
q_total = abs(q_conv + q_rad)
```

### 1.4 Combined Mesh Size Equation

Substituting the energy balance into the discretized Fourier relation:

```
Δx_max = k · ΔT_max / |h(T_s - T_∞) + εσ(T_s^4 - T_surr^4)|
```

This is the **key result**: mesh size depends only on material conductivity,
the accuracy target ΔT_max, and boundary conditions that are known or estimatable.

**Implementation:** `conduction.py` line 107:
```python
max_dx_m = (k * max_dt) / q_total
```

### 1.5 Radiation Fraction

To assess the relative importance of radiation vs. convection:

```
q_conv = h(T_s - T_∞)
q_rad = εσ(T_s^4 - T_surr^4)

f_rad = |q_rad| / (|q_conv| + |q_rad|)
```

Note: the denominator sums **absolute values** individually, not the absolute
value of their sum. This prevents cancellation when convection and radiation
oppose each other (e.g., a cold surface gaining heat by radiation but losing
it by convection).

- `f_rad ≈ 0` → convection dominated → radiation mesh constraint is loose
- `f_rad > 0.5` → radiation important → both constraints matter
- `f_rad ≈ 1` → radiation dominated → finest mesh required

**Physical Insight:** Radiation flux increases as T⁴, making fine
discretization essential at high temperatures.

---

## 2. Biot Number and Lumped Capacitance Criterion

### 2.1 Definition

The Biot number characterizes the ratio of internal conduction resistance
to surface convection resistance:

```
Bi = h · L_c / k
```

where:
- `h` = convective HTC (W/m²·K)
- `L_c = thickness / 2` = half-thickness (shortest path from center to surface)
- `k` = solid thermal conductivity (W/m·K)

The resistance interpretation:

```
R_conduction = L_c / k    [internal, K·m²/W]
R_convection = 1 / h      [surface, K·m²/W]

Bi = R_conduction / R_convection = (L_c / k) / (1 / h) = h · L_c / k
```

**Implementation:** `convection.py` lines 76–77:
```python
lc = thickness / 2.0
bi = (h * lc) / k
```

### 2.2 Physical Interpretation

```
Bi << 1   →  Conduction resistance << convection resistance
            →  Heat moves easily through the solid
            →  Temperature nearly uniform across thickness
            →  Surface convection is the bottleneck

Bi >> 1   →  Conduction resistance >> convection resistance
            →  Steep temperature gradients exist within the material
            →  Internal conduction is the bottleneck
            →  Full spatial discretization required through thickness
```

### 2.3 Lumped Capacitance Criterion

**Criterion:** `Bi < 0.1` allows lumped capacitance model

When Bi < 0.1, internal temperature gradients are less than ~5% of the
surface-to-fluid temperature difference. A single element through the
thickness (shell mesh) is physically valid.

- **Metals at typical h:** Almost always Bi << 0.1. Steel with h = 20 W/m²K
  and 3mm thickness: Bi = 20 × 0.0015 / 54 = 0.0006.
- **Plastics at moderate h:** Can cross Bi = 0.1 at h ≈ 25 W/m²K.
  PA66-GF30 with h = 25 and 3.5mm: Bi = 25 × 0.00175 / 0.25 = 0.175.

### 2.4 Element Count Heuristic

For `Bi ≥ 0.1`, through-thickness elements are needed. The heuristic:

```
n = max(2, min(10, ceil(10 · Bi)))
```

Examples:
- `Bi = 0.05` → `ceil(0.5) = 1` → `max(2, min(10, 1)) = 2` (minimum)
- `Bi = 0.2`  → `ceil(2) = 2`   → `max(2, min(10, 2)) = 2`
- `Bi = 0.5`  → `ceil(5) = 5`   → `max(2, min(10, 5)) = 5`
- `Bi = 1.0`  → `ceil(10) = 10` → `max(2, min(10, 10)) = 10` (maximum)

**Implementation:** `convection.py` line 92:
```python
n_elem = max(2, min(10, math.ceil(10 * bi)))
```

---

## 3. Radiation Flux Sensitivity and Nonlinear Discretization Error

### 3.1 Radiation Flux-Temperature Relationship

The net radiation heat flux emitted by a surface is:

```
q'' = εσT^4
```

The **sensitivity** of radiation flux to temperature change is:

```
dq/dT = 4εσT³
```

### 3.2 Discretization Error Analysis

When we discretize temperature with element spacing `Δx` and the spatial
gradient is `dT/dx`, the temperature variation across one element is:

```
ΔT_element = (dT/dx) · Δx
```

The resulting radiation flux error has two components:

**First-order (linear) error bound:**

```
δq ≈ (dq/dT) · ΔT_element = 4εσT³ · (dT/dx) · Δx
```

This is the dominant term and gives a conservative bound because it
overestimates the error (the T⁴ curve is convex, so the linear
approximation always overestimates the flux at intermediate points).

**Second-order (curvature) error:**

```
δq₂ = ½ · (d²q/dT²) · ΔT² = ½ · 12εσT² · ΔT²
```

This is the true linearization error from the Taylor expansion, but is
smaller than the first-order bound for practical ΔT values.

The code implements the **first-order bound** because it is simpler,
conservative, and adequate for mesh sizing:

```
δq ≤ δq_allow

→ 4εσT³ · (dT/dx) · Δx ≤ δq_allow

→ Δx_max = δq_allow / (4εσT³ · |dT/dx|)
```

**Implementation:** `radiation.py` lines 110–112:
```python
dq_dt = 4.0 * emissivity * STEFAN_BOLTZMANN * (t_local ** 3)
max_dt_element = allowable_flux_error / dq_dt
max_dx_m = max_dt_element / spatial_gradient
```

### 3.3 Numerical Sensitivity at Key Temperatures

The sensitivity coefficient `dq/dT = 4εσT³` at representative temperatures
(ε = 0.8, σ = 5.67 × 10⁻⁸):

| Temperature (K) | 4εσT³ (W/m²·K) | Flux error per 1 K uncertainty |
|-----------------|-----------------|-------------------------------|
| 300             | 4.9             | ~5 W/m²                       |
| 500             | 22.7            | ~23 W/m²                      |
| 800             | 92.9            | ~93 W/m²                      |
| 1000            | 181.4           | ~181 W/m²                     |

**Verification:** At T = 1000 K: `4 × 0.8 × 5.67e-8 × (1000)³ = 4 × 0.8 × 56.7 = 181.4`

**Interpretation:** At 1000 K, a 1 K temperature error produces ~181 W/m²
flux error — a 37× increase over the 300 K value. At 10 K uncertainty,
the flux error is 1814 W/m². This necessitates fine meshing near
high-temperature surfaces.

The ratio between any two temperatures scales as `(T₁/T₂)³`:
- 1000/300 ratio: (10/3)³ = 37.0 → 181.4/4.9 = 37.0 ✓
- 800/500 ratio: (8/5)³ = 4.096 → 92.9/22.7 = 4.09 ✓

---

## 4. View Factor Curvature Limit

### 4.1 Geometric Setup

For radiation calculations on curved surfaces (exhaust pipes, turbo
housings), view factor accuracy degrades when flat mesh facets poorly
approximate the true surface curvature.

A chord of length `L` on a circle of radius `R` subtends an angle:

```
θ = 2 · arcsin(L / (2R))
```

### 4.2 Inversion for Maximum Chord Length

Setting a maximum allowable facet angle `θ_max` and inverting:

```
sin(θ_max / 2) = L_max / (2R)

L_max = 2R · sin(θ_max / 2)
```

**Implementation:** `radiation.py` line 148:
```python
theta_rad = math.radians(max_facet_angle_deg)
return 2.0 * radius_mm * math.sin(theta_rad / 2.0)
```

### 4.3 Practical Application

The code defaults to `θ_max = 15°` (0.262 rad). This balances view factor
accuracy with practical mesh density.

Example: For a 50 mm radius exhaust pipe at θ_max = 15°:
```
L_max = 2 × 50 × sin(7.5°) = 100 × 0.1305 = 13.05 mm
```

Tighter thresholds (e.g., 10°) may be appropriate for concave surfaces
facing a concentrated heat source.

---

## 5. Single-Layer Radiation Shield — Newton-Raphson Iteration

### 5.1 Physical Setup

A thin heat shield sits between a hot exhaust source (temperature `T_exh`,
fixed boundary condition) and a cooler ambient environment (temperature
`T_fluid` for convection, `T_surr` for radiation). The shield temperature
`T_s` is unknown and must be solved.

Key assumption: the shield is thin enough that Bi << 0.1, so the
temperature is uniform through its thickness.

### 5.2 Energy Balance

The shield absorbs radiation from the exhaust and loses heat by radiation
to surroundings and by convection from both faces to the ambient air.

At steady state, net energy into shield = 0:

```
q_rad_in − q_rad_out − q_conv_in − q_conv_out = 0
```

where (using the code's sign convention — convection terms are positive
when the shield is hotter than fluid, i.e., heat leaving the shield):

```
q_rad_in  = ε_in  · σ · (T_exh⁴ − T_s⁴)     [net radiation absorbed from exhaust]
q_rad_out = ε_out · σ · (T_s⁴ − T_surr⁴)     [net radiation emitted to surroundings]
q_conv_in = h_in  · (T_s − T_fluid)            [convection from exhaust-facing surface]
q_conv_out = h_out · (T_s − T_fluid)           [convection from ambient-facing surface]
```

Note: the code uses **separate emissivities** (`ε_in`, `ε_out`) for the
two shield faces, since they may have different surface treatments (e.g.,
aluminised exhaust side, bare ambient side). Both convection terms use the
same `T_fluid` (ambient air).

The residual function:

```
F(T_s) = ε_in·σ·(T_exh⁴ − T_s⁴) − ε_out·σ·(T_s⁴ − T_surr⁴)
         − h_in·(T_s − T_fluid) − h_out·(T_s − T_fluid)
```

Expanding and collecting terms in T_s:

```
F(T_s) = ε_in·σ·T_exh⁴ + ε_out·σ·T_surr⁴ + (h_in + h_out)·T_fluid
         − (ε_in + ε_out)·σ·T_s⁴ − (h_in + h_out)·T_s
```

**Implementation:** `shields.py` lines 129–134:
```python
q_rad_in = eps_in * STEFAN_BOLTZMANN * (t_exh**4 - t**4)
q_rad_out = eps_out * STEFAN_BOLTZMANN * (t**4 - t_surr**4)
q_conv_in = h_i * (t - t_fluid)
q_conv_out = h_o * (t - t_fluid)
f = q_rad_in - q_rad_out - q_conv_in - q_conv_out
```

### 5.3 Derivative for Newton-Raphson

Differentiating F(T_s) with respect to T_s:

```
F'(T_s) = −4σ·T_s³·(ε_in + ε_out) − (h_in + h_out)
```

This derivative is always negative (F is monotonically decreasing in T_s),
which guarantees convergence of Newton-Raphson from any reasonable
initial guess.

**Implementation:** `shields.py` line 135:
```python
fp = -4.0 * STEFAN_BOLTZMANN * (t**3) * (eps_in + eps_out) - h_sum
```

### 5.4 Newton-Raphson Iteration

Update rule:

```
T_s^{n+1} = T_s^n − F(T_s^n) / F'(T_s^n)
```

Initial guess: `T_s^0 = (T_exh + T_fluid) / 2`

Convergence tolerance: `|T_s^{n+1} − T_s^n| < 0.1 K`

Typically converges in 3–5 iterations due to quadratic convergence rate.

### 5.5 Mesh Size from Solved Temperature

Once `T_s` is known, the boundary-driven conduction mesh size follows from
Section 1. The driving flux is taken at the hottest face (exhaust side):

```
q_boundary = |q_rad_in| + |q_conv_in|

Δx_max = k · ΔT_max / q_boundary
```

---

## 6. Multilayer Radiation Shield — 2×2 Newton-Raphson System

### 6.1 Problem Setup

Two shield layers separated by an air gap. Layer 1 (inner) faces the
exhaust; Layer 2 (outer) faces the ambient environment. Each layer has
its own unknown equilibrium temperature (T₁, T₂).

### 6.2 Gap Emissivity (Parallel-Plate Model)

Radiation exchange across the gap between two parallel surfaces uses an
effective emissivity:

```
ε_eff = 1 / (1/ε_g1 + 1/ε_g2 − 1)
```

where `ε_g1` and `ε_g2` are the emissivities of the gap-facing surfaces
of layers 1 and 2 respectively.

Example: ε_g1 = ε_g2 = 0.8:
```
ε_eff = 1 / (1.25 + 1.25 − 1) = 1/1.5 = 0.667
```

**Implementation:** `shields.py` line 279:
```python
eps_eff = 1.0 / ((1.0 / eps_g1) + (1.0 / eps_g2) - 1.0)
```

### 6.3 Residual Equations

The code models the following energy balance for each layer:

**Layer 1 (inner, exhaust-facing):**

```
F₁(T₁, T₂) = q_rad_in + q_conv_in − q_gap_cond − q_gap_rad
```

where:
```
q_rad_in   = ε_in · σ · (T_exh⁴ − T₁⁴)      [net radiation from exhaust]
q_conv_in  = h_in · (T_fluid − T₁)             [convection on inner face]
q_gap_cond = h_gap · (T₁ − T₂)                 [conduction/convection across gap]
q_gap_rad  = ε_eff · σ · (T₁⁴ − T₂⁴)          [radiation across gap]
```

Note: `q_conv_in` uses `(T_fluid − T₁)`, i.e., heat flows into the shield
when the fluid is hotter. For typical automotive shields, the ambient air
is cooler than the shield, so this term is negative (net convective loss).

**Layer 2 (outer, ambient-facing):**

```
F₂(T₁, T₂) = q_gap_cond + q_gap_rad − q_rad_out − q_conv_out
```

where:
```
q_rad_out  = ε_out · σ · (T₂⁴ − T_surr⁴)     [radiation to surroundings]
q_conv_out = h_out · (T₂ − T_fluid)            [convection from outer face]
```

**Implementation:** `shields.py` lines 287–295:
```python
q_rad_in = eps_in * STEFAN_BOLTZMANN * (t_exh**4 - t1**4)
q_conv_in = h_in * (t_fluid - t1)
q_gap_cond = h_gap * (t1 - t2)
q_gap_rad = eps_eff * STEFAN_BOLTZMANN * (t1**4 - t2**4)
q_rad_out = eps_out * STEFAN_BOLTZMANN * (t2**4 - t_surr**4)
q_conv_out = h_out * (t2 - t_fluid)

f1 = q_rad_in + q_conv_in - q_gap_cond - q_gap_rad
f2 = q_gap_cond + q_gap_rad - q_rad_out - q_conv_out
```

### 6.4 Jacobian Matrix

The 2×2 Jacobian for the Newton-Raphson system:

```
J = [ ∂F₁/∂T₁   ∂F₁/∂T₂ ]
    [ ∂F₂/∂T₁   ∂F₂/∂T₂ ]
```

Differentiating each residual term by term:

**∂F₁/∂T₁:**

```
∂(q_rad_in)/∂T₁   = −4·ε_in·σ·T₁³           [from −ε_in·σ·T₁⁴]
∂(q_conv_in)/∂T₁   = −h_in                     [from h_in·(T_fluid − T₁)]
∂(−q_gap_cond)/∂T₁ = −h_gap                    [from −h_gap·(T₁ − T₂)]
∂(−q_gap_rad)/∂T₁  = −4·ε_eff·σ·T₁³           [from −ε_eff·σ·T₁⁴]

∂F₁/∂T₁ = −4·ε_in·σ·T₁³ − h_in − h_gap − 4·ε_eff·σ·T₁³
```

**∂F₁/∂T₂:**

```
∂(−q_gap_cond)/∂T₂ = +h_gap                    [from −h_gap·(T₁ − T₂)]
∂(−q_gap_rad)/∂T₂  = +4·ε_eff·σ·T₂³           [from +ε_eff·σ·T₂⁴]

∂F₁/∂T₂ = h_gap + 4·ε_eff·σ·T₂³
```

**∂F₂/∂T₁:**

```
∂(q_gap_cond)/∂T₁  = +h_gap                    [from h_gap·(T₁ − T₂)]
∂(q_gap_rad)/∂T₁   = +4·ε_eff·σ·T₁³           [from ε_eff·σ·T₁⁴]

∂F₂/∂T₁ = h_gap + 4·ε_eff·σ·T₁³
```

**∂F₂/∂T₂:**

```
∂(q_gap_cond)/∂T₂  = −h_gap                    [from h_gap·(T₁ − T₂)]
∂(q_gap_rad)/∂T₂   = −4·ε_eff·σ·T₂³           [from −ε_eff·σ·T₂⁴]
∂(−q_rad_out)/∂T₂  = −4·ε_out·σ·T₂³           [from −ε_out·σ·T₂⁴]
∂(−q_conv_out)/∂T₂ = −h_out                     [from −h_out·(T₂ − T_fluid)]

∂F₂/∂T₂ = −h_gap − 4·ε_eff·σ·T₂³ − 4·ε_out·σ·T₂³ − h_out
```

**Implementation:** `shields.py` lines 316–324:
```python
j11 = (-4.0 * eps_in * STEFAN_BOLTZMANN * t1**3
       - h_in - h_gap
       - 4.0 * eps_eff * STEFAN_BOLTZMANN * t1**3)
j12 = h_gap + 4.0 * eps_eff * STEFAN_BOLTZMANN * t2**3
j21 = h_gap + 4.0 * eps_eff * STEFAN_BOLTZMANN * t1**3
j22 = (-h_gap
       - 4.0 * eps_eff * STEFAN_BOLTZMANN * t2**3
       - 4.0 * eps_out * STEFAN_BOLTZMANN * t2**3
       - h_out)
```

### 6.5 Newton-Raphson Update (Cramer's Rule)

Rather than inverting the 2×2 Jacobian with a linear algebra library,
the code uses Cramer's rule (pure Python, no numpy dependency):

```
det(J) = J₁₁·J₂₂ − J₁₂·J₂₁
```

The update vector [ΔT₁, ΔT₂] is:

```
ΔT₁ = −(F₁·J₂₂ − F₂·J₁₂) / det(J)
ΔT₂ = −(J₁₁·F₂ − J₂₁·F₁) / det(J)
```

Update:

```
T₁^{n+1} = T₁^n + ΔT₁
T₂^{n+1} = T₂^n + ΔT₂
```

**Implementation:** `shields.py` lines 326–331:
```python
det = j11 * j22 - j12 * j21
dt1 = -(f1 * j22 - f2 * j12) / det
dt2 = -(j11 * f2 - j21 * f1) / det
```

Convergence: typically 4–8 iterations to reach |F₁| < 0.1 and |F₂| < 0.1.

---

## 7. Transient Thermal Constraints

### 7.1 Thermal Penetration Depth

The unsteady heat conduction equation is:

```
∂T/∂t = α · ∂²T/∂x²
```

where `α = k/(ρ·c_p)` is the thermal diffusivity (m²/s).

By dimensional analysis, a thermal disturbance propagates a distance:

```
δ ~ √(α · dt)
```

in one time step `dt`. In an implicit scheme, whose stability places no
limit on `dt`, an element much larger than δ cannot resolve the thermal
wavefront within a step and smears it. The resolution guideline is:

```
Δx_max ≤ C · √(α · dt),    C ≈ 1.5–2.0   (implicit schemes)
```

In Fourier-number form (§7.2) this reads `Fo = α·dt/Δx² ≥ 1/C²`: a **lower**
bound on Fo. It is not a stability limit, and it is not applied to explicit
schemes (§7.4 shows why it cannot be).

**Implementation:** `transient.py`, `penetration_depth()`:
```python
alpha = k / (rho * cp)
pen = math.sqrt(alpha * dt)
max_dx_m = safety_factor * pen
```

### 7.2 Fourier Number and Explicit Stability

The Fourier number is the dimensionless ratio of the diffusion time scale
to the element time scale:

```
Fo = α · dt / (Δx)²
```

Rearranging for minimum element size at a given Fo limit:

```
Δx_min = √(α · dt / Fo_max)
```

**Stability limits:**

```
Explicit (forward Euler), 1D:   Fo ≤ 0.5
Explicit (forward Euler), 3D:   Fo ≤ 1/6    [≈ 0.167]
Implicit (backward Euler):      Unconditionally stable,
                                 but Fo ≤ 5 for accuracy
```

The 3D limit (1/6) arises because heat diffuses in three independent
directions simultaneously; with uniform mesh spacing, each direction
contributes Fo/3, and the effective 1D Fourier number in any direction
must satisfy Fo_1D ≤ 0.5.

**Implementation:** `transient.py`, `fourier_number_limit()`:
```python
alpha = k / (rho * cp)
min_dx_m = math.sqrt(alpha * dt / fo_max)
```

### 7.3 Drive-Cycle Constraint

In automotive thermal analysis, boundary conditions change on a time scale
`τ_bc` (e.g., 10–30 s between operating point transitions in a drive cycle).
The mesh must resolve the thermal response within each segment:

```
Δx_max ≤ √(α · τ_bc)
```

This is independent of the solver time-step — it constrains the mesh to
ensure the spatial field can represent the thermal state at each transition.

**Example:** Steel (α ≈ 1.2 × 10⁻⁵ m²/s), τ_bc = 1 s:
```
Δx_max ≈ √(1.2 × 10⁻⁵ × 1) = 3.46 × 10⁻³ m ≈ 3.5 mm
```

### 7.4 Combined Constraints

The constraints bound the element from both sides:

```
Δx ≥ Δx_min,Fourier = √(α·dt / Fo_max)                 [lower bound]
Δx ≤ Δx_drive_cycle  = √(α·τ_bc)                         [upper bound]
Δx ≤ Δx_penetration  = C·√(α·dt)      (implicit only)   [upper bound]
```

The binding constraint is the tightest applicable upper bound, provided it
is not below the Fourier minimum.

**Why the penetration bound is not applied to explicit schemes.** The
explicit (forward-Euler, lumped) update in 1D is

```
T_i^(n+1) = Fo·T_(i−1)^n + (1 − 2·Fo)·T_i^n + Fo·T_(i+1)^n
```

For Fo ≤ 1/2 its three weights are non-negative and sum to one: each new
value is a weighted average of old ones, so the update creates no new
extremum and cannot oscillate, however small Fo is. An element larger than
the per-step penetration depth (Fo < 1) is the normal explicit regime: it
costs steps, not accuracy. The upper bounds on an explicit scheme's
element are the drive-cycle limit and the steady constraints (§2–§5).

**The ratio that does not depend on dt.** Both √dt bounds scale alike:

```
Δx_min,Fourier / Δx_penetration = 1 / (C·√Fo_max)
```

so their window is non-empty only if `C²·Fo_max ≥ 1`, at every dt or at
none. With C = 1, the explicit limits Fo_max = 1/2 and 1/6 give ratios of
√2 ≈ 1.414 and √6 ≈ 2.449: applied to an explicit scheme, the penetration
bound leaves no feasible element at any time step, and the advice "reduce
dt" cannot change that. The implicit guideline pair C = 2, Fo_max = 5 gives
`C²·Fo_max = 20`, a window from 0.447 to 2 × √(α·dt).

**Remedies that change the result.**

| Conflict | Remedy | Why it closes the conflict |
|---|---|---|
| Fourier minimum > drive-cycle bound | `dt ≤ Fo_max · τ_bc` | √(α·dt/Fo_max) ≤ √(α·τ_bc) ⇔ dt ≤ Fo_max·τ_bc, and the drive-cycle bound does not depend on dt |
| Fourier minimum > governing size Δx_gov (batch) | `dt ≤ Fo_max · Δx_gov² / α` | the minimum scales with √dt and Δx_gov does not shrink with dt (an implicit penetration bound at the new dt is C·√Fo_max·Δx_gov ≥ Δx_gov when `C²·Fo_max ≥ 1`) |
| Implicit, `C²·Fo_max < 1` | `C ≥ 1/√Fo_max`, or `Fo_max ≥ 1/C²` | the ratio above does not depend on dt, so no dt closes it |

**Worked example.** `steel_mild` (k = 54, ρ = 7833, c_p = 465, so
α = 1.4826 × 10⁻⁵ m²/s), explicit, Fo_max = 0.5:

```
dt = 0.01 s:  Δx ≥ 0.5445 mm
dt = 1 s:     Δx ≥ 5.445 mm
dt = 100 s:   Δx ≥ 54.45 mm
```

With no τ_bc no transient upper bound applies, and the steady constraints
set the size. With τ_bc = 20 s the drive-cycle bound is √(α·20) = 17.22 mm,
so dt = 100 s conflicts; the remedy dt ≤ 0.5 × 20 = 10 s brings the minimum
down to exactly 17.22 mm.

**Implementation:** `transient.py`, `combined_transient_limits()` (the
window, `feasible`, and `conflict` with its `remedies`); `batch.py`,
`process_part()`, which uses the transient upper bound as a size candidate
and reports a Fourier minimum above the governing size as a
`TRANSIENT_CONFLICT` warning carrying the dt remedy.

---

## 8. Convective Heat Transfer Coefficient Correlations

### 8.1 Fundamental Relationships

Heat transfer from/to a fluid is characterized by the Nusselt number:

```
Nu = h · L / k_fluid
```

where:
- `h` = convection coefficient (W/m²·K)
- `L` = characteristic length (m)
- `k_fluid` = thermal conductivity of the fluid (W/m·K)

Properties are evaluated at the **film temperature**:

```
T_film = (T_surface + T_fluid) / 2
```

### 8.2 Air Properties

For air at atmospheric pressure (polynomial fits valid 250–700 K):

```
k_air = 0.0241 + 7.0×10⁻⁵ · (T − 300)     [W/m·K, linear fit]
ν     = 1.5×10⁻⁵ · (T/300)^1.7              [m²/s, power-law fit]
Pr    ≈ 0.71                                  [nearly constant]
β     = 1/T                                   [ideal gas, 1/K]
ρ     = 101325 / (287.05 · T)                [ideal gas, kg/m³]
```

**Implementation:** `h_estimator.py` lines 61–76.

### 8.3 Forced Convection — Flat Plate

**Reynolds number:**
```
Re = V · L / ν
```

**Laminar regime** (Re < 5 × 10⁵):
```
Nu = 0.664 · Re^0.5 · Pr^(1/3)
```
(Blasius boundary layer solution)

**Turbulent regime** (Re ≥ 5 × 10⁵):
```
Nu = 0.037 · Re^0.8 · Pr^(1/3)
```
(Colburn analogy / modified Reynolds analogy)

**Geometry note:** These are flat-plate correlations. The package uses them
for automotive components as a first approximation. For cylindrical geometry
(exhaust pipes), Churchill-Bernstein is more appropriate but would add
complexity for marginal improvement in a mesh-sizing context.

**Implementation:** `h_estimator.py` lines 271–282:
```python
if Re < 5e5:
    Nu = 0.664 * (Re ** 0.5) * Pr_term    # laminar
else:
    Nu = 0.037 * (Re ** 0.8) * Pr_term    # turbulent
h_raw = Nu * props["k_air"] / char_length
```

### 8.4 Natural Convection — Buoyancy-Driven

**Grashof number:**
```
Gr = g · β · |ΔT| · L³ / ν²
```

**Rayleigh number:**
```
Ra = Gr · Pr
```

**Vertical surface (characteristic length = height):**
```
Laminar  (Ra < 10⁹):   Nu = 0.59 · Ra^(1/4)
Turbulent (Ra ≥ 10⁹):  Nu = 0.10 · Ra^(1/3)
```
(Simplified Churchill-Chu; coefficients from Incropera & DeWitt)

**Horizontal surface, hot side up (L = Area/Perimeter):**
```
Laminar  (Ra < 10⁷):   Nu = 0.54 · Ra^(1/4)
Turbulent (Ra ≥ 10⁷):  Nu = 0.15 · Ra^(1/3)
```
(McAdams correlations)

**Horizontal surface, hot side down:**
```
Nu = 0.27 · Ra^(1/4)   [always laminar — gravitationally stable]
```

The hot-side-down configuration suppresses buoyant plumes, producing
weak convection. h is typically very low (3–10 W/m²K).

**Implementation:** `h_estimator.py` lines 349–377.

### 8.5 Richardson Number — Regime Classification

When both forced flow and buoyancy are present, the **Richardson number**
determines which dominates:

```
Ri = Gr / Re²
```

**Physical interpretation:** Ri compares buoyancy forces (Gr) to inertia
forces (Re²). When Ri is small, forced flow overwhelms buoyancy.

**Regime classification:**

```
Ri < 0.1     →  Forced convection dominates
                 Use forced correlation only

0.1 ≤ Ri ≤ 10 → Mixed convection
                 Both effects are significant
                 Use asymptotic blending

Ri > 10       →  Natural convection dominates
                 Use natural correlation only
```

**Implementation:** `h_estimator.py` lines 440–453.

### 8.6 Mixed Convection — Churchill & Usagi Asymptotic Blending

For assisting mixed convection (buoyancy aids forced flow), the standard
blending formula is:

```
h_mixed = (h_forced³ + h_natural³)^(1/3)
```

The cubic exponent (n = 3) is the standard choice for assisting flows
(Churchill & Usagi, 1972). This formula:
- Tends to the larger of the two for disparate magnitudes
- Smoothly blends when both are of similar order
- Avoids the discontinuity of simply switching between correlations

**Note:** For opposing flows (buoyancy opposes forced flow), some references
use n = 4 or subtract the natural component. The code uses the assisting
formulation for all cases, which is conservative for mesh sizing.

**Implementation:** `h_estimator.py` line 537:
```python
h_raw = (h_frc**3 + h_nat**3) ** (1.0 / 3.0)
```

### 8.7 External h Cap

For external surfaces at typical off-road analysis speeds (~10 mph),
the practical maximum h is bounded:

```
H_EXTERNAL_TYPICAL = 60 W/m²K    [typical max at ~10 mph]
H_EXTERNAL_MAX     = 100 W/m²K   [absolute cap for external surfaces]
```

Internal flows (exhaust gas, clutch outlet) are **not** capped — internal
forced convection can legitimately exceed 100 W/m²K.

---

## 9. Solver Selection Advisory Logic

### 9.1 Steady-State vs. Transient Analysis

The decision depends on whether the dominant convection mode produces
inherently unsteady heat transfer:

```
Laminar natural convection  →  Steady plumes       →  SS solve OK
Turbulent forced convection →  Time-averaged Nu OK  →  SS solve OK
Turbulent natural convection → Oscillating plumes   →  Transient needed
Mixed convection             → Depends on Ri        →  Check case-by-case
```

### 9.2 Turbulent Natural Convection — Why It Breaks Steady-State

At high Rayleigh numbers (Ra >> Ra_crit), buoyant plumes:
- Form at the surface, rise, cool, and detach periodically
- Create oscillating temperature and velocity fields
- Have time-varying local h that deviates significantly from time-averaged Nu

A steady-state solver using time-averaged Nu produces a single temperature
field that may not converge (residuals oscillate) or may converge to a
solution that misrepresents peak temperatures and thermal stresses.

### 9.3 Ra Thresholds for Turbulence

The critical Rayleigh number depends on surface orientation:

```
Vertical:          Ra_crit = 10⁹
Horizontal hot-up: Ra_crit = 10⁷    [lower — plumes form more easily]
Horizontal hot-down: Ra_crit = ∞    [always laminar — stable stratification]
```

**Key finding from parametric studies:** For horizontal hot-up dead zones,
the Ra = 10⁷ threshold corresponds to characteristic lengths of only
~130–190 mm. Nearly ALL horizontal underbody parts exceed this — meaning
turbulent natural convection is the norm, not the exception, in dead zones.

### 9.4 Forced Flow Suppression of Turbulent Buoyancy

When Ri < 0.1, forced convection dominates and buoyancy effects are
negligible. Even a modest vehicle speed (~10 mph / 4.5 m/s) typically
produces Ri << 0.1 for most automotive scenarios.

**Key finding:** Critical velocity for buoyancy suppression is typically
3–9 mph. At the 10 mph analysis baseline, all tested off-road scenarios
showed Ri < 0.1 (forced-dominant, SS-safe).

### 9.5 Decision Tree Implementation

The `solver_advisory()` function evaluates:

1. Is natural convection turbulent? (Ra > Ra_crit for the orientation)
2. Is forced convection present and dominant? (Ri < 0.1)
3. What is the combined regime?

Severity levels:

```
"info"     — Laminar natural or turbulent forced only; SS is fine
"caution"  — Mixed regime, or forced-dominant with turbulent natural
             background; monitor convergence
"warning"  — Turbulent natural dominates; transient solve recommended
```

### 9.6 Summary Table: Solver Selection

| Flow Type | Regime | SS OK? | Severity | Rationale |
|-----------|--------|--------|----------|-----------|
| Forced only | Laminar | Yes | info | Stable boundary layer |
| Forced only | Turbulent | Yes | info | Time-averaged Nu valid |
| Natural only | Laminar | Yes | info | Stable plumes |
| Natural only | Turbulent | No | warning | Oscillating plumes |
| Mixed (Ri<0.1) | Forced dom. | Yes | info/caution | Forced overwhelms buoyancy |
| Mixed (0.1<Ri<10) | Comparable | Monitor | caution | Both contribute |
| Mixed (Ri>10) | Natural dom. | Check Ra | warning if turbulent | Buoyancy controls regime |

---

## 10. Fin Theory Lateral Gradient Constraint (v0.4)

### 10.1 Problem Statement

The boundary-driven conduction constraint (Section 1) sizes elements to resolve
through-thickness gradients at the boundary. However, a localised heat source
on a thin conductive plate also creates lateral (in-plane) temperature gradients
that decay exponentially away from the source. If the surface mesh is too coarse,
it smears these hot spots and under-predicts peak temperatures.

### 10.2 Extended Surface (Fin) Theory

From classical fin theory, a localised heat input on a thin plate of
conductivity `k` and thickness `t` loses heat to the surrounding fluid
with total coefficient `h_total`. The temperature perturbation decays
as:

```
θ(x) = θ₀ · exp(−m · x)
```

where `m` is the fin parameter:

```
m = √(h_total / (k · t))
```

The characteristic decay length is `1/m`. The linearised radiation contribution
to `h_total` is:

```
h_total = h_conv + h_rad = h + 4·ε·σ·T³
```

where `4·ε·σ·T³` is the derivative of the Stefan-Boltzmann law evaluated at
the local surface temperature (same quantity appearing in Section 3).

### 10.3 Mesh Constraint

To capture the exponential decay within a fraction `f` of the characteristic
length (so the mesh resolves ~95% of the gradient when `f = 1/3`):

```
Δx_lateral ≤ f / m = f · √(k · t / h_total)
```

The default `f = 1/3` is conservative; `f = 1/4` is tighter for critical areas.

**Code:** `BoundaryDrivenConductionCalculator.lateral_gradient_limit()` in `conduction.py`

### 10.4 Example

Steel bracket: `k = 45 W/m·K`, `h = 20 W/m²·K`, `t = 3 mm`, `ε = 0.92`,
`T_s = 473 K`:

```
h_rad = 4 × 0.92 × 5.67e-8 × 473³ = 22.1 W/m²·K
h_total = 20 + 22.1 = 42.1 W/m²·K
m = √(42.1 / (45 × 0.003)) = 17.66 1/m
decay_length = 1/17.66 = 56.6 mm
Δx_lateral = 56.6 / 3 = 18.9 mm
```

This is wired into `batch.py` as the `lateral_gradient` constraint.

---

## 11. Cell Biot Number Heuristic (v0.4)

### 11.1 Motivation

The classical Biot number `Bi = h·L_c/k` determines whether a shell or
3D solid mesh is needed (Section 2). When 3D elements are required, the
number of elements across the thickness must be sufficient to prevent
numerical oscillation.

### 11.2 Cell Biot Number

Define the cell Biot number as the Biot number of a single element:

```
Bi_cell = h · Δx / k
```

where `Δx = thickness / n` is the element size across the thickness.
For numerical stability (particularly with explicit time integration),
we require:

```
Bi_cell ≤ 0.1
```

Substituting:

```
h · (thickness / n) / k ≤ 0.1
n ≥ h · thickness / (0.1 · k) = 10 · h · (thickness/2) / (0.5 · k)
```

Using `L_c = thickness/2` and `Bi = h·L_c/k`:

```
n ≥ 10 · 2 · Bi = 20 · Bi
```

The factor of 20 (rather than the simpler 10) accounts for the half-thickness
definition of the Biot number and provides margin for non-uniform h distributions
across the element face.

### 11.3 Implementation

```
n = max(2, min(10, ceil(20 · Bi)))
```

Clamped to [2, 10] for practical meshing limits.

**Code:** `ConvectionMeshCalculator.biot_number()` in `convection.py`

---

## 12. Gap View Factor F₁₂ in Multilayer Shields (v0.4)

### 12.1 Classical Parallel-Plate Emissivity

Section 6 derived the effective gap emissivity for infinite parallel plates:

```
ε_eff = 1 / (1/ε_g1 + 1/ε_g2 − 1)
```

This assumes a view factor `F₁₂ = 1` between the gap faces — valid for
infinite parallel plates but not for real automotive shields that have
finite extent, offset geometry, or non-parallel surfaces.

### 12.2 General Two-Surface Enclosure

For a two-surface enclosure with gray, diffuse surfaces, the net radiative
exchange is:

```
q_gap = σ·(T₁⁴ − T₂⁴) / (1/ε_g1 + 1/ε_g2 − 2 + 1/F₁₂)
```

This reduces to the classical formula when `F₁₂ = 1` (note: `−2 + 1/1 = −1`).

The generalised effective emissivity is therefore:

```
ε_eff = 1 / (1/ε_g1 + 1/ε_g2 − 2 + 1/F₁₂)
```

### 12.3 Physical Interpretation

When `F₁₂ < 1`, the gap faces "see" less of each other, reducing radiative
coupling across the gap. This means:
- Lower `ε_eff` → less radiation crosses the gap
- Inner layer retains more heat → higher T₁
- Outer layer receives less radiation → lower T₂
- Larger temperature difference `ΔT = T₁ − T₂`

### 12.4 Typical Values

| Geometry | F₁₂ |
|---|---|
| Infinite parallel plates | 1.0 |
| Typical automotive offset shield | 0.80–0.90 |
| Small or highly non-parallel gap | 0.50–0.75 |
| Deeply concave surface into flat | 0.30–0.60 |

Default: `F₁₂ = 1.0` (backward compatible; use ~0.85 for typical automotive).

**Code:** `MultilayerShieldCalculator.solve_temperatures(f12=...)` in `shields.py`

---

## 13. Opposing Mixed Convection (v0.4)

### 13.1 Assisting vs. Opposing Flow

Section 8 introduced the asymptotic combination for mixed convection:

```
h_mixed = (h_forced³ + h_natural³)^(1/3)
```

This is valid for **assisting** flow, where forced and buoyant flows act
in the same direction (e.g., upward forced flow over a hot vertical surface).

In **opposing** flow, the forced flow acts against the buoyant plume.
This occurs, for example, when:
- Forced flow pushes downward over a hot horizontal surface (hot-side-down)
- The buoyancy-driven plume rises against the forced flow direction

### 13.2 Opposing Formula

For opposing flow, the asymptotic combination subtracts the modes:

```
h_opposing = |h_forced³ − h_natural³|^(1/3)
```

This can produce very small values when the two modes nearly cancel.
A physical lower bound prevents singularity:

```
h_opposing = max(|h_f³ − h_n³|^(1/3), k_air / L)
```

The `k_air / L` term represents pure conduction across the characteristic
length — the minimum possible heat transfer when all convective motion
is suppressed by the opposing flows.

### 13.3 Detection

Opposing flow is detected by orientation:
- `horizontal_down` with forced flow → opposing (buoyancy lifts off the surface,
  forced flow pushes fluid toward it)
- `vertical` and `horizontal_up` → assisting (default)

### 13.4 Implementation

In `estimate_h()`:
```python
if is_opposing:
    h_raw = max(abs(h_f**3 - h_n**3)**(1/3), k_air / L)
else:
    h_raw = (h_f**3 + h_n**3)**(1/3)
```

**Code:** `estimate_h()` in `h_estimator.py`

---

## Summary of Key Design Equations

**Mesh sizing (conduction):**
```
Δx_max = k · ΔT_max / |h(T_s − T_∞) + εσ(T_s⁴ − T_surr⁴)|
```

**Lateral gradient (fin theory):**
```
Δx_lateral ≤ (1/3) · √(k · t / h_total)
h_total = h + 4·ε·σ·T³
```

**Biot number and cell Biot element count:**
```
Bi = h · L_c / k     [L_c = thickness/2]
n = max(2, min(10, ceil(20 · Bi)))    [ensures Bi_cell ≤ 0.1]
```

**Radiation sensitivity:**
```
dq/dT = 4εσT³
Δx_max = δq_allow / (4εσT³ · |dT/dx|)
```

**Fourier number (explicit stability):**
```
Fo = α · dt / (Δx)² ≤ 0.5 (1D) or ≤ 1/6 (3D)
Δx_min = √(α · dt / Fo_max)
```

**Single shield Newton-Raphson:**
```
F(T_s) = ε_in·σ·(T_exh⁴ − T_s⁴) − ε_out·σ·(T_s⁴ − T_surr⁴)
         − (h_in + h_out)·(T_s − T_fluid)
F'(T_s) = −4σ·T_s³·(ε_in + ε_out) − (h_in + h_out)
T_s^{n+1} = T_s^n − F/F'
```

**Gap view factor (multilayer shield):**
```
ε_eff = 1 / (1/ε_g1 + 1/ε_g2 − 2 + 1/F₁₂)
```

**Natural convection (vertical):**
```
Nu = 0.59 · Ra^(1/4)  [Ra < 10⁹, laminar]
Nu = 0.10 · Ra^(1/3)  [Ra ≥ 10⁹, turbulent]
```

**Mixed convection (Churchill & Usagi):**
```
Assisting:  h_mixed = (h_forced³ + h_natural³)^(1/3)
Opposing:   h_mixed = max(|h_forced³ − h_natural³|^(1/3), k_air/L)
```

**Richardson number:**
```
Ri = Gr / Re²
Ri < 0.1: forced;  0.1 ≤ Ri ≤ 10: mixed;  Ri > 10: natural
```

---

## 14. Aerodynamic Boundary Layer Mesh Constraints

### 14.1 Motivation

Thermal mesh sizing (Sections 1–13) determines element sizes needed to resolve
temperature gradients within the solid and at the solid–fluid boundary.  However,
when the surface mesh feeds a coupled CFD–thermal solver, the volume mesh near
walls must also resolve the velocity and temperature boundary layers in the fluid.
The surface element size constrains the inflation (prism) layer footprint, and
excessively coarse surface elements produce degenerate volume cells that smear the
flow field — and therefore the convective heat transfer coefficient that the
thermal solver depends on.

### 14.2 First Cell Height from Target y⁺

The non-dimensional wall distance y⁺ relates the first cell height y₁ to the
friction velocity u_τ:

```
y⁺ = y₁ · u_τ / ν
→  y₁ = y⁺ · ν / u_τ
```

where the friction velocity is obtained from the wall shear stress:

```
u_τ = √(τ_w / ρ)
τ_w = ½ · C_f · ρ · U²
```

**Skin friction coefficient (flat plate):**

```
Laminar  (Re_x < 5×10⁵):  C_f = 0.664 · Re_x^(−0.5)    [Blasius]
Turbulent (Re_x ≥ 5×10⁵):  C_f = 0.058 · Re_x^(−0.2)    [Schlichting]
```

The local Reynolds number Re_x = U·x/ν uses the distance from the leading edge
(or an equivalent characteristic length for non-plate geometries).

### 14.3 Leading-Edge Singularity Protection

As x → 0, Re_x → 0, causing C_f → ∞ and y₁ → 0.  This is handled by:

1. **Reynolds number floor:** Re_x = max(U·x/ν, Re_x_min) with Re_x_min = 1000
2. **Absolute y₁ floor:** y₁ = max(calculated, 1×10⁻⁶ m)

Re_x_min = 1000 corresponds to ~1 mm at 22 m/s in air at 350 K, which is
effectively the smallest length scale relevant to vehicle-level meshing.

### 14.4 Boundary Layer Thickness

```
Laminar:   δ = 5.0  · x · Re_x^(−0.5)
Turbulent: δ = 0.37 · x · Re_x^(−0.2)
```

### 14.5 Inflation Layer Count

Given first cell height y₁, boundary layer thickness δ, and geometric growth
ratio r (typically 1.2), the number of prism layers to span δ:

```
n = ⌈ ln(1 + (δ/y₁)·(r − 1)) / ln(r) ⌉
```

with a minimum of 3 layers.  The outermost prism height is:

```
y_last = y₁ · r^(n−1)
```

and the total inflation stack height (geometric series):

```
H_prism = y₁ · (r^n − 1) / (r − 1)
```

### 14.6 Surface Mesh Constraint — Two-Regime Model

The surface element size dx must be compatible with the inflation stack to avoid
degenerate cell quality.  The constraint differs by flow regime:

**External forced (attached, unidirectional flow):**

The critical quality metric is the volume transition from the outermost prism
cell into the isotropic tet/poly core.  The constraint is:

```
dx ≤ min(AR_transition · y_last,  C_BL · δ)
```

where AR_transition ≈ 3 controls the aspect ratio at the prism-to-tet transition,
and C_BL ≈ 0.3 ensures enough lateral points to resolve streamwise BL variations.
Prism layers are always recommended in this regime.

**Mixed / unknown (recirculating, impinging, buoyancy-dominated):**

The surface constraint uses only the boundary layer fraction:

```
dx ≤ C_BL · δ
```

The prism aspect ratio constraint (`dx ≤ AR_max · y₁`) is **not applied** in this
regime.  This is based on the analysis that prism layers provide no meaningful
accuracy improvement over isotropic tets in low-Re mixed/buoyancy zones (see
Section 14.9 for the full justification).

### 14.7 Buoyancy Velocity Scale

For the mixed/unknown regime, the freestream velocity U may be zero or irrelevant.
The characteristic velocity is the maximum of U and the buoyancy velocity:

```
V_buoy = √(g · β · ΔT · L)
U_eff  = max(U, V_buoy, 0.1 m/s)
```

where β = 1/T_film (ideal gas), ΔT = |T_surf − T_fluid|, and L is the
characteristic length (e.g. gap height between exhaust and cargo bed).

For the exhaust scenario (ΔT ≈ 200 K, T_film ≈ 450 K, L = 0.1 m):

```
V_buoy = √(9.81 · (1/450) · 200 · 0.1) ≈ 0.66 m/s
```

This is an order of magnitude below a 22 m/s vehicle speed, producing much
thicker y₁ and correspondingly relaxed surface mesh — matching the physical
reality of gentle, low-shear buoyancy-driven flows.

### 14.8 Practical Guidelines

| Parameter | External Forced | Mixed / Unknown |
|-----------|----------------|-----------------|
| y⁺ target | 30–300 (wall functions) | 30–100 |
| Growth ratio | 1.2 | 1.2 |
| Surface dx | min(AR_trans · y_last, C·δ) | C·δ (BL fraction only) |
| Velocity scale | U (freestream) | max(U, V_buoy) |
| Prism layers | Yes (always) | No (isotropic tets to wall) |

Wall-resolved (y⁺ = 1) meshing is not recommended for full-vehicle underbody
thermal analysis.  It produces sub-millimetre surface elements that explode
cell counts without proportionate accuracy gains in the thermal solution.

### 14.9 Prism Layer Recommendation — Low-Re Justification

The calculator outputs a `prisms_recommended` flag that evaluates whether
structured inflation layers provide meaningful accuracy improvement over
isotropic tets at the wall.  For the mixed/unknown regime, prisms are generally
**not recommended**.  The justification rests on two independent arguments:

#### 14.9.1 First Cell Height vs. Surface Mesh Size

At low velocities (0.1–1 m/s typical of buoyancy-dominated zones), the y⁺
formalism produces first cell heights of 10–60 mm — vastly exceeding the
surface mesh size of 1–5 mm set by the BL fraction constraint.  This means
isotropic tets at the thermal mesh size already contain 5–10+ cells within the
y⁺-derived "first cell" zone, over-resolving the wall-normal gradient.  Prism
layers would be coarser than the existing tet resolution at the wall.

#### 14.9.2 Wall Temperature Gradient Magnitude

The finite-volume cross-diffusion error from tet non-orthogonality at the wall
is:

```
Error ∝ k_geo · ∇φ
```

where `k_geo` is the non-orthogonality factor (~0.1–0.2 for reasonable quality
tets) and `∇φ = dT/dy` is the wall temperature gradient.

For the two regimes:

```
External forced (exhaust, 22 m/s):  dT/dy ~ 10⁶–10⁷ K/m
Mixed/buoyancy (dead zone, 0.5 m/s): dT/dy ~ 10³–10⁴ K/m
```

The gradient is ~500× smaller in the mixed/buoyancy zone.  Even with the
geometric penalty of tet non-orthogonality (~10–20% error in h), the absolute
error is within the ±30% physics uncertainty of the h correlation itself in
these mixed/buoyancy regions.

#### 14.9.3 The "Kicking the Can" Argument

Prism layers don't eliminate non-orthogonality — they move it from the wall to
the prism-to-tet transition (typically layer 5–10).  In high-Re attached flow,
the temperature gradient decays exponentially through the boundary layer, so the
error at the transition is much smaller than at the wall.  Prisms provide a real
benefit.

In low-Re buoyancy-dominated flow, the temperature profile across a gap is
nearly linear.  The gradient at layer 5 is essentially the same as at the wall.
The non-orthogonality error is simply relocated, not reduced.  The volume
expansion ratio (ER_v) at the transition introduces additional truncation error
that offsets any benefit from wall-normal alignment.

#### 14.9.4 Volume Expansion Ratio (ER_v)

For the external forced regime where prisms are recommended, the calculator
reports the volume expansion ratio at the prism-to-tet transition:

```
ER_v = (dx_surface / y_last)³
```

This estimates the volume jump from the outermost prism cell to the first
isotropic tet.  The truncation error of 2nd-order schemes degrades as ER_v
increases:

```
ER_v ≤ 5  : stable, 2nd-order accuracy maintained
5 < ER_v ≤ 10 : marginal, gradient truncation at interface
ER_v > 10 : unstable, diagonal dominance loss likely
```

For the mixed/unknown regime where prisms are not recommended, ER_v is not
computed (reported as "n/a").

---

## Summary of Key Design Equations (Updated)

All constraints from Sections 1–14:

**Conduction (boundary-driven):**
```
dx_max = k · ΔT_max / |q''_conv + q''_rad|
```

**Fin theory (lateral gradient):**
```
dx_lateral ≤ (1/3) · √(k·t / h_total)
h_total = h_conv + 4·ε·σ·T³
```

**Biot number (through-thickness):**
```
Bi = h·L/k;  Bi > 0.1 → 3D Solid
n ≥ ceil(20·Bi) elements through thickness
```

**Radiation (T⁴ sensitivity):**
```
dx_rad = ΔF_allow / (4·ε·σ·T³ · ∇T_spatial)
```

**Transient (Fourier number):**
```
dx_fo = √(α · dt / Fo_max)
```

**Shield solver (Newton–Raphson):**
```
F(T_s) = ε_in·σ·(T_exh⁴ − T_s⁴) − ε_out·σ·(T_s⁴ − T_surr⁴)
         − (h_in + h_out)·(T_s − T_fluid)
F'(T_s) = −4σ·T_s³·(ε_in + ε_out) − (h_in + h_out)
T_s^{n+1} = T_s^n − F/F'
```

**Gap view factor (multilayer shield):**
```
ε_eff = 1 / (1/ε_g1 + 1/ε_g2 − 2 + 1/F₁₂)
```

**Natural convection (vertical):**
```
Nu = 0.59 · Ra^(1/4)  [Ra < 10⁹, laminar]
Nu = 0.10 · Ra^(1/3)  [Ra ≥ 10⁹, turbulent]
```

**Mixed convection (Churchill & Usagi):**
```
Assisting:  h_mixed = (h_forced³ + h_natural³)^(1/3)
Opposing:   h_mixed = max(|h_forced³ − h_natural³|^(1/3), k_air/L)
```

**Richardson number:**
```
Ri = Gr / Re²
Ri < 0.1: forced;  0.1 ≤ Ri ≤ 10: mixed;  Ri > 10: natural
```

**Aerodynamic boundary layer (first cell):**
```
y₁ = y⁺ · ν / u_τ
u_τ = U · √(C_f / 2)
C_f = 0.058 · Re_x^(−0.2)  [turbulent]
```

**Aerodynamic boundary layer (surface mesh):**
```
External forced:  dx ≤ min(AR_trans · y_last, C_BL · δ)   [prisms recommended]
Mixed/unknown:    dx ≤ C_BL · δ                            [isotropic tets to wall]
```

**Volume expansion ratio (prism-to-tet transition):**
```
ER_v = (dx_surface / y_last)³
ER_v ≤ 5: stable;  5 < ER_v ≤ 10: marginal;  ER_v > 10: unstable
```

---

*End of Mathematical Derivations Document*
