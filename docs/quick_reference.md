# Thermal Mesh Calculators — Quick Reference

## Getting Started

- **No install needed** — pure Python 3.6+
- **Import pattern:** `from thermal_mesh_calculators.batch import process_batch, summary_table`
- **Temperature units:** Kelvin internally, output includes `_C` keys for Celsius
- **Length units:** meters internally, output in `_mm` keys for millimeters

---

## Common Workflows

### 1. "I need mesh size for a single exhaust component"

Use `BoundaryDrivenConductionCalculator.max_mesh_size()` with known T_surf, h, emissivity.

```python
from thermal_mesh_calculators.conduction import BoundaryDrivenConductionCalculator

result = BoundaryDrivenConductionCalculator.max_mesh_size(
    k=52.0,              # W/m K — cast iron
    h=200.0,             # W/m^2 K — exhaust internal
    t_surf=1073.15,      # K (800°C)
    t_fluid=373.15,      # K (100°C underhood air)
    epsilon=0.85,        # oxidised iron
    t_surr=373.15,       # K radiation sink
    max_dt=50.0,         # K — accuracy target per element
)
print(f"Max element size: {result['max_dx_mm']:.2f} mm")
```

### 2. "I need mesh size for a structural part (non-exhaust)"

Same calculator, but T_surf is estimated (±30 K typical).

```python
result = BoundaryDrivenConductionCalculator.max_mesh_size(
    k=54.0,              # W/m K — mild steel
    h=12.0,              # W/m^2 K — engine bay dead zone
    t_surf=473.15,       # K (200°C estimate)
    t_fluid=353.15,      # K (80°C underhood)
    epsilon=0.92,        # painted surface
    t_surr=353.15,       # K
    max_dt=20.0,         # K — stricter accuracy for structural
)
```

### 3. "Should I use shell or solid elements?"

Use `ConvectionMeshCalculator.biot_number()` — Biot number governs the decision.

```python
from thermal_mesh_calculators.convection import ConvectionMeshCalculator

result = ConvectionMeshCalculator.biot_number(
    h=15.0,              # W/m^2 K
    k=54.0,              # W/m K
    thickness=0.003,     # m (3 mm)
)
# If biot < 0.1 → shell mesh (lumped capacitance OK)
# If biot ≥ 0.1 → 3D solid mesh needed
```

For a wall at known temperatures, `BoundaryDrivenConductionCalculator.size_wall()`
gives the number of cells through it, with radiation included:

```python
from thermal_mesh_calculators.conduction import BoundaryDrivenConductionCalculator

sized = BoundaryDrivenConductionCalculator.size_wall(
    k=0.25, h=40.0,                  # W/m K, W/m^2 K — glass-filled nylon
    t_surf=393.15, t_fluid=353.15,   # K (120°C surface, 80°C air)
    epsilon=0.92, t_surr=353.15,     # K
    max_dt=5.0,                      # K — accuracy target per element
    thickness_m=0.003,               # m (3 mm wall)
)
# sized["n_cells"] == 5, sized["regime"] == "resolve", round(sized["biot"], 2) == 0.61
```

`regime` is `"thermally_thin"` (Bi < 0.1: shell), `"resolve"` (0.1 ≤ Bi < 1),
`"steep_gradient"` (Bi ≥ 1) or `"near_equilibrium"`. This Biot number uses the whole
wall thickness and includes radiation, so it differs from `biot_number()`'s.

**Rule of thumb table:**

| Material Type | k (W/m K) | Typical h (W/m² K) | Biot | Recommendation |
|---|---|---|---|---|
| Steel/Iron | 50–60 | 5–15 | << 0.1 | 2D Shell |
| Aluminium | 150–170 | 5–15 | << 0.1 | 2D Shell |
| Plastic (unfilled) | 0.2–0.3 | 8–15 | >> 0.1 | 3D Solid if h > 25 |
| Thin gauge metal | — | any | << 0.1 | Always Shell |
| < 1 mm thick | — | any | << 0.1 | Always Shell |

### 4. "What's the mesh size for a heat shield?"

Use `SingleLayerShieldCalculator.mesh_size()` — temperature is **solved**, not input.

```python
from thermal_mesh_calculators.shields import SingleLayerShieldCalculator

result = SingleLayerShieldCalculator.mesh_size(
    t_exh=1073.15,       # K (800°C exhaust surface, known)
    t_fluid=353.15,      # K (80°C underhood)
    t_surr=353.15,       # K (radiation sink)
    eps_in=0.40,         # aluminised facing exhaust
    eps_out=0.40,        # aluminised facing ambient
    h_in=12.0,           # W/m^2 K — natural convection near exhaust
    h_out=8.0,           # W/m^2 K — stagnant air side
    k=52.0,              # W/m K — mild steel
    thickness=0.0008,    # m (0.8 mm)
    max_dt=30.0,         # K
)
# Returns equilibrium T_shield and mesh size for each surface
```

### 5. "I have a dual-wall shield with an air gap"

Use `MultilayerShieldCalculator.mesh_sizes()` for dimpled/offset shields.

```python
from thermal_mesh_calculators.shields import MultilayerShieldCalculator

result = MultilayerShieldCalculator.mesh_sizes(
    t_exh=1073.15,       # K (exhaust side, known)
    t_fluid=353.15,      # K (ambient air)
    t_surr=353.15,       # K
    eps_in=0.40,         # inner surface emissivity
    eps_out=0.40,        # outer surface emissivity
    eps_gap=0.40,        # gap surface (both layers)
    h_in=12.0,           # W/m^2 K — near exhaust
    h_gap=8.0,           # W/m^2 K — stagnant air gap
    h_out=5.0,           # W/m^2 K — cavity/cargo side
    k1=52.0, k2=52.0,    # W/m K — layer materials
    thickness1=0.0008, thickness2=0.0008,  # m — layer thickness
    gap_thickness=0.015, # m — air gap spacing
    max_dt=30.0,         # K
)
# Returns T1, T2, T_gap and individual mesh sizes
```

### 6. "I need transient mesh constraints"

Use `TransientMeshCalculator.combined_transient_limits()` — find the binding constraint.

```python
from thermal_mesh_calculators.transient import TransientMeshCalculator

result = TransientMeshCalculator.combined_transient_limits(
    k=52.0, rho=7200.0, cp=460.0,  # cast iron diffusivity
    dt=0.5,              # solver time-step (s)
    fo_max=0.5,          # explicit stability limit
    safety_factor=1.0,
    tau_bc=20.0,         # drive-cycle segment duration (s)
)
# Returns recommended_dx_mm and binding_constraint name
# binding_constraint tells you which limit is the tightest
```

### 7. "I want to estimate h without CFD"

Use `estimate_h()` with zone physics, or `h_from_velocity()` with known velocity.

```python
from thermal_mesh_calculators.h_estimator import estimate_h, h_from_velocity

# Method A: zone-based (proper physics)
h_result = estimate_h(
    t_surf=473.15,           # K (200°C)
    t_fluid=353.15,          # K (80°C)
    velocity_ms=0.0,         # m/s (dead zone → natural convection)
    orientation="vertical",  # or "horizontal_up", "horizontal_down"
)
print(f"h = {h_result['h']:.1f} W/m^2 K")

# Method B: direct velocity (flat-plate forced convection)
h_direct = h_from_velocity(
    velocity_ms=3.0,         # m/s (cooling pack flow)
    t_film=363.15,           # film temperature (K)
    char_length_m=0.1,       # characteristic length (m)
)
```

**When to use each:**
- **estimate_h()** — when you have surface temps and zone info (preferred for automation)
- **h_from_velocity()** — when you have measured/CFD velocity and want quick HTC

### 8. "I need to process hundreds of parts from a BOM"

Use `process_batch()` with a list of part dicts.

```python
from thermal_mesh_calculators.batch import process_batch, summary_table

project = {
    "t_fluid_K": 353.15,      # 80°C underhood
    "t_surr_K": 353.15,       # radiation sink
    "t_exh_K": 1073.15,       # 800°C exhaust (for shields)
    "max_dt": None,           # use class defaults
    "dt": 0.5,                # transient time-step (s)
    "fo_max": 0.5,            # explicit solver
    "tau_bc": 20.0,           # drive-cycle segment (s)
}

parts = [
    {
        "part_id": "STR-001",
        "material": "steel_mild",
        "surface": "painted",
        "component_class": "structural",
        "convection_zone": "engine_beside",
        "thickness_mm": 3.0,
        "t_surf_K": 473.15,    # 200°C
        "char_length_mm": 200.0,
    },
    {
        "part_id": "SH-001",
        "material": "steel_mild",
        "component_class": "shield",
        "convection_zone_in": "exhaust_beside",
        "convection_zone_out": "engine_beside",
        "thickness_mm": 0.8,
        "t_surf_K": 0,         # ignored for shields
        "surface_in": "aluminised",
        "surface_out": "aluminised",
    },
]

results = process_batch(parts, project)
print(summary_table(results))
```

**Part dict structure:**
- `part_id`, `material`, `component_class` — required
- `convection_zone` — zone name (see [Convection Zones](#convection-zones))
- `thickness_mm` — wall thickness
- `t_surf_K` — surface temperature (ignored for shields; solved instead)
- `surface` or `surface_in`/`surface_out` — optional emissivity specification
- `h_override` — optional specific h instead of zone lookup
- `char_length_mm` — optional for h estimation

### 9. "How do I know if steady-state is appropriate?"

Check the `solver_advisory` output in batch results.

```python
advisory = result["solver_advisory"]
print(f"Severity: {advisory['severity']}")  # "info", "caution", "warning"
print(f"Steady-state OK: {advisory['steady_state_ok']}")
print(f"Reasoning: {advisory['message']}")
```

**Quick decision tree:**

| Severity | Recommendation |
|---|---|
| `info` | Proceed with steady-state |
| `caution` | Steady-state acceptable, but monitor solver convergence |
| `warning` | **Strongly consider transient** — turbulent natural convection flagged |

---

## Material & Surface Lookup

### Available Materials (43 total)

List by category with material key names for part dicts:

**Steels (5):**
`steel_mild`, `steel_high_strength`, `steel_galvanised`, `steel_stainless_304`, `steel_stainless_409`

**Aluminium (4):**
`aluminium_6061`, `aluminium_5052`, `cast_aluminium`, `cast_aluminium_a380`

**Cast Iron (2):**
`cast_iron`, `cast_iron_ductile`

**Other Metals (6):**
`copper`, `brass`, `magnesium_az91`, `titanium_6al4v`, `zinc`, `nickel_alloy`

**Plastics — Glass-Filled (5):**
`plastic_pa66_gf30`, `plastic_pa6_gf30`, `plastic_pp_gf30`, `plastic_pbt_gf30`, `plastic_pps_gf40`

**Plastics — Unfilled (8):**
`plastic_pa66`, `plastic_pp`, `plastic_hdpe`, `plastic_abs`, `plastic_pc`, `plastic_pc_abs`, `plastic_acetal_pom`, `plastic_pet`

**Rubbers (6):**
`rubber_epdm`, `rubber_silicone`, `rubber_nbr`, `rubber_natural`, `rubber_cr`, `rubber_fkm`

**Composites & Specialty (7):**
`composite_smc`, `composite_cfrp`, `glass_soda_lime`, `ceramic_alumina`, `ceramic_cordierite`, `insulation_fibreglass`, `insulation_ceramic_blanket`

### Available Surface Treatments (30 total)

**Bare/Polished Metals:**
`bare_metal`, `polished_aluminium` (ε=0.05), `polished_steel` (ε=0.07), `polished_copper` (ε=0.04), `cast_aluminium_bare` (ε=0.25)

**Oxide Layers:**
`lightly_oxidised` (ε=0.50), `heavily_oxidised` (ε=0.85), `cast_iron_oxidised` (ε=0.80), `aluminium_oxidised` (ε=0.25), `copper_oxidised` (ε=0.78), `stainless_weathered` (ε=0.85), `rusted_steel` (ε=0.70)

**Coatings (Industrial):**
`aluminised` (ε=0.40), `painted` (ε=0.92), `painted_gloss_white` (ε=0.90), `painted_matte_black` (ε=0.95), `anodised` (ε=0.80), `galvanised_new` (ε=0.23), `galvanised_weathered` (ε=0.88), `chrome_plated` (ε=0.06), `nickel_plated` (ε=0.12), `zinc_plated` (ε=0.20), `e_coat` (ε=0.92), `powder_coat` (ε=0.92), `ceramic_coating` (ε=0.60)

**Non-Metals:**
`plastic` (ε=0.92), `rubber` (ε=0.90), `glass` (ε=0.90), `composite_smc` (ε=0.90), `fabric_woven` (ε=0.90)

### Emissivity Fallback Logic

If no surface treatment is specified in the part dict:

| Material Class | Default ε |
|---|---|
| Non-metal (plastic, rubber, ceramic, glass) | 0.90 |
| Aluminium (any alloy) | 0.30 |
| All other metals (steel, cast iron, copper, etc.) | 0.73 |

---

## Convection Zones

All zones are in `zones.CONVECTION_ZONES`. Use the key name in `convection_zone` field.

**16 zones with regime, velocity, and typical air temperature range:**

| Zone | Regime | Velocity (m/s) | Typical T_air (°C) | Orientation |
|---|---|---|---|---|
| `cooling_pack_downstream` | forced | 3.0 | 80–100 | vertical |
| `front_end_edges` | forced | 4.5 | 40–50 | vertical |
| `front_end_dead_zone` | natural | 0.0 | 40–60 | vertical |
| `cabin_dead_zone` | natural | 0.0 | 40–50 | vertical |
| `cabin_tunnel` | natural | 0.0 | 40–50 | horizontal_up |
| `engine_bay_edges` | mixed | 2.0 | 40–50 | vertical |
| `engine_beside` | natural | 0.0 | 60–70 | vertical |
| `engine_below` | natural | 0.0 | 50–60 | horizontal_down |
| `engine_above` | natural | 0.0 | 80–100 | horizontal_up |
| `exhaust_beside` | natural | 0.0 | 90–150 | vertical |
| `exhaust_below` | natural | 0.0 | 90–150 | horizontal_down |
| `exhaust_above` | natural | 0.0 | 150–350 | horizontal_up |
| `clutch_outlet_downstream` | forced | 5.0 | 90–110 | vertical |
| `exhaust_internal` | forced | 30.0 | 500–900 | vertical |
| `shield_gap_confined` | natural | 0.0 | 80–200 | vertical |
| `cargo_bed_exterior` | mixed | 2.0 | 40–80 | horizontal_up |

**Legacy aliases:** `engine_bay_dead_zone` → `engine_beside`, `near_exhaust_natural` → `exhaust_beside`

---

## Warning Codes

Common warnings from the batch processor:

| Code | Severity | Meaning | Action |
|---|---|---|---|
| `TURB_NAT_HORIZ_UP` | caution | Horizontal hot-side-up dead zone > 150 mm | May trigger turbulent natural convection; monitor solver |
| `TURB_NAT_VERTICAL` | caution | Vertical dead zone > 600 mm | Onset of turbulent buoyancy; consider transient |
| `PLASTIC_BIOT_MARGINAL` | info | Non-metal Biot ≥ 0.1 | 3D solid mesh required (not critical) |
| `SOLVER_TRANSIENT_RECOMMENDED` | warning | Turbulent natural convection flagged | Strongly consider transient analysis |
| `LOW_VELOCITY_HIGH_DT` | info | Low forced velocity with high temp rise | Verify zone velocity is realistic |

---

## Key Thresholds (from Parametric Studies)

Critical thresholds discovered in the package's regime studies:

| Phenomenon | Threshold | Implication |
|---|---|---|
| Horizontal hot-up turbulent onset | L > 130–190 mm (Ra > 1e7) | Natural convection becomes turbulent |
| Vertical turbulent onset | L > 600 mm (Ra > 1e9) | Vertical dead zones turn turbulent at large scale |
| Horizontal hot-down | Always laminar | Weakest natural convection mode (h very low) |
| Plastic Biot crossover | h ≈ 25 W/m² K | For typical thickness; 3D effects appear above this h |
| Wind suppression of buoyancy | v ≈ 3–9 mph (1.3–4 m/s) | Natural convection suppressed by forced flow |
| Metal Biot | Bi << 0.1 | Conduction always governs mesh (shell mesh valid) |

---

## Units Reference

| Quantity | Input Unit | Output Key |
|---|---|---|
| Temperature | K (Kelvin) | `_K` (Kelvin), `_C` (Celsius) |
| Length | m (meters) | `_mm` (millimeters) |
| Thermal conductivity | W/m K | — |
| Heat transfer coefficient | W/m² K | — |
| Heat flux | W/m² | — |
| Density | kg/m³ | — |
| Specific heat | J/kg K | — |
| Emissivity | 0–1 (dimensionless) | — |

**Quick conversions:**
- 273.15 K = 0°C, 373.15 K = 100°C
- 1 mm = 0.001 m

---

## Examples Quick Start

See `/examples/batch_example.py` for a complete side-by-side vehicle analysis with 12 parts (exhaust, structural, shields, plastics, rubbers).

Run: `python -m examples.batch_example`

---

## Common Errors & Fixes

| Error | Cause | Fix |
|---|---|---|
| `KeyError: "Unknown material 'X'"` | Typo in material name | Use `list_materials()` to see available names |
| `KeyError: "Unknown surface treatment 'X'"` | Typo in surface name | Use `list_surface_treatments()` |
| `KeyError: "Unknown convection zone 'X'"` | Typo in zone name | Use `list_zones()` or check legacy aliases |
| `ValueError: "Provide either (h_in, h_out) or h_total"` | Shield solver missing HTC | Supply both `h_in` and `h_out` or single `h_total` |
| `AttributeError: 'NoneType'` in transient | Missing `k`, `rho`, `cp` for material | Verify material is in `MATERIALS` dict; check spelling |
| Mesh size = `inf` or NaN | q_total = 0 W/m² (no heat transfer) | Check temperature difference is non-zero |

