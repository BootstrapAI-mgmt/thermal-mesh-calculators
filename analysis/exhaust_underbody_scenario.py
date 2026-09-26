"""
Exhaust Underbody Scenario Analysis
====================================

Off-road vehicle underbody exhaust system with shields, cargo shield, and cargo bed.

Layout (temperatures in °C, from rear forward):
  825 → 800 → [bend] → 775 → 727 (UNSHIELDED) → [bend] → 750 → 675 → silencer (675/650)

Components:
  - Exhaust pipes + silencer: 2mm SS409, heavily oxidised (eps=0.85)
  - Exhaust shields: 1mm SS409, aluminised inner (eps=0.40) / lightly oxidised outer (eps=0.50)
    On all upper halves EXCEPT the 727°C region
  - Cargo shield: 3.5mm aluminium Al 5052, bare both sides (eps=0.25)
  - Cargo bed: 6.6mm PP GF30 plastic (eps=0.92)

Geometry:
  - Ground to exhaust centreline: 1m (0.5m to plastic skid plate)
  - Exhaust pipes to cargo shield: 190mm
  - Silencer to cargo shield: 77mm
  - Cargo shield to cargo bed: 67mm
  - Exhaust lengths: 400mm (to 1st bend), 300mm (to 2nd bend), 350mm (to silencer mid)
  - Silencer: 300mm square base × 200mm height
  - Pipe OD: ~60mm (assumed)

Shield temperature assumptions (engineering estimates):
  - Exhaust pipe shields: ~300°C cooler than adjacent pipe
  - Silencer shield: ~425°C
  - Cargo shield: ~100-125°C (slightly hotter near silencer)

Key findings:
  1. Solver confirms ~307°C average drop across exhaust shields (validates estimate)
  2. View factor F12 is CRITICAL for cargo shield — parallel-plate assumption (F12=1)
     gives absurd 280-580°C; realistic F12 (0.15 for pipes, 0.55 for silencer) gives
     ~100-160°C matching engineering expectation
  3. Unshielded 727°C section is the thermal design driver
  4. Cargo bed (plastic) needs 3D solid mesh (Bi > 0.1) with very fine elements
  5. Exhaust shields: lateral gradient (fin theory) is the binding constraint, not conduction
  6. All metallic components are shell-meshable (Bi << 0.1)

Run from the repository root, as a script or as a module:
    python analysis/exhaust_underbody_scenario.py
    python -m analysis.exhaust_underbody_scenario
"""

import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from thermal_mesh_calculators.shields import SingleLayerShieldCalculator
from thermal_mesh_calculators.conduction import BoundaryDrivenConductionCalculator
from thermal_mesh_calculators.convection import ConvectionMeshCalculator
from thermal_mesh_calculators.h_estimator import estimate_h
from thermal_mesh_calculators.boundary_layer import BoundaryLayerCalculator
from thermal_mesh_calculators.constants import STEFAN_BOLTZMANN

# ============================================================
# Global boundary conditions
# ============================================================
T_FLUID = 353.15    # 80°C underhood air (K)
T_SURR = 333.15     # 60°C radiation sink (K)

# ============================================================
# Material properties
# ============================================================
K_SS409 = 25.0      # SS409 thermal conductivity (W/m·K)
K_AL5052 = 138.0    # Al 5052 thermal conductivity (W/m·K)
K_PP_GF30 = 0.22    # PP GF30 thermal conductivity (W/m·K)

# ============================================================
# Exhaust pipe temperatures (°C)
# ============================================================
PIPE_SECTIONS = {
    "pipe_A_rear":      {"t_C": 825, "length_mm": 400, "note": "rear, to 1st bend"},
    "pipe_B":           {"t_C": 800, "length_mm": 200, "note": "mid"},
    "pipe_C_postbend":  {"t_C": 775, "length_mm": 100, "note": "after 1st bend"},
    "pipe_D_unshielded":{"t_C": 727, "length_mm": 200, "note": "NO SHIELD"},
    "pipe_E_postbend2": {"t_C": 750, "length_mm": 150, "note": "after 2nd bend"},
    "pipe_F_presilencer":{"t_C": 675, "length_mm": 200, "note": "pre-silencer"},
    "silencer_hot":     {"t_C": 675, "length_mm": 300, "note": "silencer inlet side"},
    "silencer_cool":    {"t_C": 650, "length_mm": 300, "note": "silencer outlet side"},
}

SHIELDED_SECTIONS = {k: v for k, v in PIPE_SECTIONS.items()
                     if k != "pipe_D_unshielded"}


def run_analysis():
    """Run the complete exhaust underbody scenario analysis."""
    results = {"exhaust_shields": {}, "cargo_shield": {}, "cargo_bed": {}, "mesh_sizes": {}}

    # ----------------------------------------------------------------
    # 1. Exhaust shield temperatures
    # ----------------------------------------------------------------
    print("=" * 80)
    print("1. EXHAUST SHIELD TEMPERATURES")
    print("=" * 80)
    print(f"\n  {'Section':<25s} {'T_pipe(°C)':>10s} {'T_shield(°C)':>12s} {'ΔT(°C)':>8s}")
    print("  " + "-" * 60)

    for name, sec in SHIELDED_SECTIONS.items():
        t_exh = sec["t_C"] + 273.15
        r = SingleLayerShieldCalculator.solve_temperature(
            t_exh=t_exh, t_fluid=T_FLUID, t_surr=T_SURR,
            eps_in=0.40, eps_out=0.50,  # aluminised / lightly oxidised
            h_in=12.0, h_out=20.0,
        )
        results["exhaust_shields"][name] = r
        delta = sec["t_C"] - r["t_shield_C"]
        print(f"  {name:<25s} {sec['t_C']:>10d} {r['t_shield_C']:>12.1f} {delta:>8.1f}")

    # ----------------------------------------------------------------
    # 2. Cargo shield temperatures (with view factor correction)
    # ----------------------------------------------------------------
    print("\n" + "=" * 80)
    print("2. CARGO SHIELD TEMPERATURES (view-factor corrected)")
    print("=" * 80)

    cargo_zones = [
        ("shielded_worst",  510, 0.50, 0.15, 190, 10.0, 8.0, "Worst shielded pipe"),
        ("shielded_mid",    445, 0.50, 0.15, 190, 10.0, 8.0, "Mid shielded pipe"),
        ("silencer",        380, 0.50, 0.55,  77, 12.0, 8.0, "Silencer zone"),
        ("unshielded_727",  727, 0.85, 0.15, 190, 12.0, 8.0, "UNSHIELDED 727 pipe"),
    ]

    print(f"\n  {'Zone':<25s} {'T_src(°C)':>10s} {'F12':>6s} {'eps_eff':>8s} {'T_cargo(°C)':>12s}")
    print("  " + "-" * 65)

    for name, t_src_C, eps_src, F12, gap, h_in, h_out, desc in cargo_zones:
        t_src_K = t_src_C + 273.15
        eps_cargo_in = 0.25  # bare aluminium
        eps_eff = 1.0 / (1.0/eps_src + 1.0/eps_cargo_in - 2.0 + 1.0/F12)

        r = SingleLayerShieldCalculator.solve_temperature(
            t_exh=t_src_K, t_fluid=T_FLUID, t_surr=T_SURR,
            eps_in=eps_eff, eps_out=0.25,
            h_in=h_in, h_out=h_out,
        )
        results["cargo_shield"][name] = {"result": r, "desc": desc, "eps_eff": eps_eff}
        print(f"  {desc:<25s} {t_src_C:>10d} {F12:>6.2f} {eps_eff:>8.4f} {r['t_shield_C']:>12.1f}")

    # ----------------------------------------------------------------
    # 3. Cargo bed temperatures
    # ----------------------------------------------------------------
    print("\n" + "=" * 80)
    print("3. CARGO BED TEMPERATURES (67mm above cargo shield)")
    print("=" * 80)

    F_BED = 0.80  # parallel plates at 67mm
    print(f"\n  {'Zone':<25s} {'T_cargo(°C)':>12s} {'T_bed(°C)':>10s} {'OK?':>6s}")
    print("  " + "-" * 58)

    for name, data in results["cargo_shield"].items():
        t_cargo_K = data["result"]["t_shield_K"]
        eps_eff_bed = 1.0 / (1.0/0.25 + 1.0/0.92 - 2.0 + 1.0/F_BED)

        r = SingleLayerShieldCalculator.solve_temperature(
            t_exh=t_cargo_K, t_fluid=T_FLUID, t_surr=T_SURR,
            eps_in=eps_eff_bed, eps_out=0.92,
            h_in=6.0, h_out=8.0,
        )
        results["cargo_bed"][name] = r
        ok = "✓" if r["t_shield_C"] < 130 else "⚠"
        print(f"  {data['desc']:<25s} {data['result']['t_shield_C']:>12.1f} "
              f"{r['t_shield_C']:>10.1f} {ok:>6s}")

    # ----------------------------------------------------------------
    # 4. Boundary layer analysis
    # ----------------------------------------------------------------
    print("\n" + "=" * 80)
    print("4. BOUNDARY LAYER ANALYSIS")
    print("=" * 80)
    results["boundary_layer"] = {}

    # Define BL scenarios for each component group
    # Vehicle speed: 80 km/h = 22.2 m/s for external-forced surfaces
    V_VEHICLE = 22.2  # m/s (~80 km/h)

    bl_scenarios = [
        # (label, U, x, regime, t_surf_C, delta_t_buoy, ar_max, y+, description)
        ("exhaust_pipe_external", V_VEHICLE, 0.4, "external_forced",
         750, 0, 5.0, 30.0, "Exhaust pipes (freestream underbody)"),
        ("exhaust_shield_external", V_VEHICLE, 0.3, "external_forced",
         450, 0, 5.0, 30.0, "Exhaust shields (outer face, freestream)"),
        ("cargo_shield_top_mixed", 0.0, 0.067, "mixed_unknown",
         150, 0, 3.0, 50.0, "Cargo shield top (gap to bed, buoyancy)"),
        ("cargo_shield_bottom_mixed", 0.0, 0.190, "mixed_unknown",
         150, 0, 3.0, 50.0, "Cargo shield bottom (gap to exhaust, buoyancy)"),
        ("cargo_bed_underside", 0.0, 0.067, "mixed_unknown",
         90, 0, 3.0, 50.0, "Cargo bed underside (gap to cargo shield)"),
        ("cargo_bed_topside", V_VEHICLE, 0.5, "external_forced",
         90, 0, 5.0, 30.0, "Cargo bed topside (if exposed to underbody flow)"),
    ]

    print(f"\n  {'Component':<35s} {'Regime':<18s} {'U_eff':>6s} {'y1(mm)':>8s} "
          f"{'n_lay':>6s} {'δ(mm)':>7s} {'dx_BL':>7s}")
    print("  " + "-" * 92)

    for (label, U, x, regime, t_surf_C, dt_buoy, ar_max,
         yp, desc) in bl_scenarios:
        t_surf_K = t_surf_C + 273.15
        bl = BoundaryLayerCalculator.estimate_mesh(
            U=U, x=x, t_fluid=T_FLUID,
            y_plus_target=yp, growth_ratio=1.2,
            regime=regime,
            delta_t_buoyancy=dt_buoy,
            ar_max_prism=ar_max,
            bl_fraction=0.3,
            t_surf=t_surf_K,
        )
        results["boundary_layer"][label] = bl
        print(f"  {desc:<35s} {regime:<18s} {bl['U_effective']:>6.2f} "
              f"{bl['y1_mm']:>8.4f} {bl['n_layers']:>6d} "
              f"{bl['delta_mm']:>7.2f} {bl['max_dx_surface_mm']:>7.2f}")

    # ----------------------------------------------------------------
    # 5. Combined mesh sizing (thermal + aero BL)
    # ----------------------------------------------------------------
    print("\n" + "=" * 80)
    print("5. COMBINED MESH SIZE SUMMARY (thermal + aero boundary layer)")
    print("=" * 80)
    print(f"\n  {'Component':<30s} {'Mat':>7s} {'t(mm)':>6s} {'T(°C)':>7s} "
          f"{'dx_therm':>8s} {'dx_BL':>7s} {'dx_gov':>7s} {'Constraint':>14s} {'Mesh':>10s}")
    print("  " + "-" * 102)

    # Exhaust pipes — external forced BL
    bl_pipe = results["boundary_layer"]["exhaust_pipe_external"]
    for name, sec in PIPE_SECTIONS.items():
        t_K = sec["t_C"] + 273.15
        h_est = estimate_h(velocity=0.0, t_surf=t_K, t_fluid=T_FLUID,
                           char_length=0.05, orientation="vertical")
        cond = BoundaryDrivenConductionCalculator.max_mesh_size(
            k=K_SS409, h=h_est["h"], t_surf=t_K, t_fluid=T_FLUID,
            epsilon=0.85, t_surr=T_SURR, max_dt=10.0)
        lat = BoundaryDrivenConductionCalculator.lateral_gradient_limit(
            k=K_SS409, h=h_est["h"], thickness_m=0.002,
            epsilon=0.85, t_surf=t_K)
        dx_therm = min(cond["max_dx_mm"], lat["max_dx_mm"])
        therm_gov = "conduction" if cond["max_dx_mm"] <= lat["max_dx_mm"] else "lateral"
        dx_bl = bl_pipe["max_dx_surface_mm"]
        dx = min(dx_therm, dx_bl)
        gov = therm_gov if dx_therm <= dx_bl else "aero_BL"
        results["mesh_sizes"][name] = {"dx_mm": dx, "dx_thermal_mm": dx_therm,
                                       "dx_bl_mm": dx_bl, "constraint": gov}
        tag = " *" if name == "pipe_D_unshielded" else ""
        print(f"  {name:<30s} {'SS409':>7s} {'2.0':>6s} {sec['t_C']:>7d} "
              f"{dx_therm:>8.1f} {dx_bl:>7.1f} {dx:>7.1f} {gov:>14s} {'Shell':>10s}{tag}")

    # Exhaust shields — external forced BL on outer face
    bl_shield = results["boundary_layer"]["exhaust_shield_external"]
    for name, r in results["exhaust_shields"].items():
        t_K = r["t_shield_K"]
        lat = BoundaryDrivenConductionCalculator.lateral_gradient_limit(
            k=K_SS409, h=32.0, thickness_m=0.001, epsilon=0.50, t_surf=t_K)
        dx_therm = lat["max_dx_mm"]
        dx_bl = bl_shield["max_dx_surface_mm"]
        dx = min(dx_therm, dx_bl)
        gov = "lateral" if dx_therm <= dx_bl else "aero_BL"
        label = f"shield_{name}"
        results["mesh_sizes"][label] = {"dx_mm": dx, "dx_thermal_mm": dx_therm,
                                        "dx_bl_mm": dx_bl, "constraint": gov}
        print(f"  {label:<30s} {'SS409':>7s} {'1.0':>6s} {r['t_shield_C']:>7.0f} "
              f"{dx_therm:>8.1f} {dx_bl:>7.1f} {dx:>7.1f} {gov:>14s} {'Shell':>10s}")

    # Cargo shield — mixed/unknown BL (bottom face toward exhaust)
    bl_cargo_bot = results["boundary_layer"]["cargo_shield_bottom_mixed"]
    bl_cargo_top = results["boundary_layer"]["cargo_shield_top_mixed"]
    for name, data in results["cargo_shield"].items():
        t_K = data["result"]["t_shield_K"]
        lat = BoundaryDrivenConductionCalculator.lateral_gradient_limit(
            k=K_AL5052, h=18.0, thickness_m=0.0035, epsilon=0.25, t_surf=t_K)
        dx_therm = lat["max_dx_mm"]
        # Use the tighter of the two faces' BL constraints
        dx_bl = min(bl_cargo_bot["max_dx_surface_mm"],
                    bl_cargo_top["max_dx_surface_mm"])
        dx = min(dx_therm, dx_bl)
        gov = "lateral" if dx_therm <= dx_bl else "aero_BL"
        label = f"cargo_shield_{name}"
        results["mesh_sizes"][label] = {"dx_mm": dx, "dx_thermal_mm": dx_therm,
                                        "dx_bl_mm": dx_bl, "constraint": gov}
        print(f"  {label:<30s} {'Al5052':>7s} {'3.5':>6s} {data['result']['t_shield_C']:>7.0f} "
              f"{dx_therm:>8.1f} {dx_bl:>7.1f} {dx:>7.1f} {gov:>14s} {'Shell':>10s}")

    # Cargo bed — mixed/unknown BL on underside
    bl_bed = results["boundary_layer"]["cargo_bed_underside"]
    for name, r in results["cargo_bed"].items():
        t_K = r["t_shield_K"]
        cond = BoundaryDrivenConductionCalculator.max_mesh_size(
            k=K_PP_GF30, h=8.0, t_surf=t_K, t_fluid=T_FLUID,
            epsilon=0.92, t_surr=T_SURR, max_dt=15.0)
        lat = BoundaryDrivenConductionCalculator.lateral_gradient_limit(
            k=K_PP_GF30, h=14.0, thickness_m=0.0066,
            epsilon=0.92, t_surf=t_K)
        biot = ConvectionMeshCalculator.biot_number(h=8.0, k=K_PP_GF30, thickness=0.0066)
        dx_therm = min(cond["max_dx_mm"], lat["max_dx_mm"])
        therm_gov = "conduction" if cond["max_dx_mm"] <= lat["max_dx_mm"] else "lateral"
        dx_bl = bl_bed["max_dx_surface_mm"]
        dx = min(dx_therm, dx_bl)
        gov = therm_gov if dx_therm <= dx_bl else "aero_BL"
        label = f"cargo_bed_{name}"
        results["mesh_sizes"][label] = {"dx_mm": dx, "dx_thermal_mm": dx_therm,
                                        "dx_bl_mm": dx_bl, "constraint": gov,
                                        "mesh_type": biot["mesh_type"],
                                        "biot": biot["biot"]}
        tag = " ***" if r["t_shield_C"] > 130 else ""
        print(f"  {label:<30s} {'PP-GF30':>7s} {'6.6':>6s} {r['t_shield_C']:>7.0f} "
              f"{dx_therm:>8.1f} {dx_bl:>7.1f} {dx:>7.1f} {gov:>14s} {biot['mesh_type']:>10s}{tag}")

    # ----------------------------------------------------------------
    # 6. Boundary layer recommendations summary
    # ----------------------------------------------------------------
    print("\n" + "=" * 80)
    print("6. INFLATION LAYER RECOMMENDATIONS")
    print("=" * 80)
    print(f"\n  {'Surface':<35s} {'y1(mm)':>8s} {'Layers':>7s} {'Growth':>7s} "
          f"{'δ(mm)':>7s} {'H_prism':>8s} {'y+(tgt)':>8s} {'Prisms?':>8s} {'ER_v':>6s}")
    print("  " + "-" * 105)

    for label, bl in results["boundary_layer"].items():
        desc_map = {
            "exhaust_pipe_external": "Exhaust pipes (outer)",
            "exhaust_shield_external": "Exhaust shields (outer)",
            "cargo_shield_top_mixed": "Cargo shield (top, to bed)",
            "cargo_shield_bottom_mixed": "Cargo shield (bottom, to exh)",
            "cargo_bed_underside": "Cargo bed (underside)",
            "cargo_bed_topside": "Cargo bed (topside)",
        }
        desc = desc_map.get(label, label)
        prism_flag = "YES" if bl["prisms_recommended"] else "NO"
        er_v_str = f"{bl['volume_expansion_ratio']:.1f}" if bl["er_v_status"] != "n/a" else "n/a"
        print(f"  {desc:<35s} {bl['y1_mm']:>8.4f} {bl['n_layers']:>7d} "
              f"{bl['growth_ratio']:>7.1f} {bl['delta_mm']:>7.2f} "
              f"{bl['total_prism_mm']:>8.2f} {bl['y_plus_target']:>8.0f} "
              f"{prism_flag:>8s} {er_v_str:>6s}")

    # Prism recommendation rationale
    print("\n  Prism Recommendation Rationale:")
    print("  " + "-" * 80)
    for label, bl in results["boundary_layer"].items():
        desc_map = {
            "exhaust_pipe_external": "Exhaust pipes (outer)",
            "exhaust_shield_external": "Exhaust shields (outer)",
            "cargo_shield_top_mixed": "Cargo shield (top, to bed)",
            "cargo_shield_bottom_mixed": "Cargo shield (bottom, to exh)",
            "cargo_bed_underside": "Cargo bed (underside)",
            "cargo_bed_topside": "Cargo bed (topside)",
        }
        desc = desc_map.get(label, label)
        print(f"  {desc}: {bl['prism_reason']}")

    return results


if __name__ == "__main__":
    results = run_analysis()
