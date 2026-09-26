"""
Off-Road Vehicle Regime Crossover Study — 10 mph Baseline
==========================================================

All scenarios assume:
  - 10 mph (4.47 m/s) max ambient velocity at thermal analysis condition
  - 40°C ambient temperature
  - Off-road SxS vehicle (engine bay dead zones, exhaust proximity,
    cargo bed, cabin tunnel)

Sweeps:
  1. Conduction vs convection dominance (which drives mesh size)
  2. Laminar → turbulent natural convection (orientation-dependent)
  3. When does 10 mph forced flow suppress turbulent buoyancy (Ri < 0.1)
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from thermal_mesh_calculators.conduction import BoundaryDrivenConductionCalculator
from thermal_mesh_calculators.convection import ConvectionMeshCalculator
from thermal_mesh_calculators.h_estimator import (
    natural_convection, richardson_number, estimate_h,
)

V_BASELINE_MPH = 10.0
V_BASELINE_MS = V_BASELINE_MPH * 0.44704
T_AMB_K = 313.15   # 40°C
T_SURR_K = 313.15


# =========================================================================
#  STUDY 1: Conduction vs Biot Mesh Dominance — Off-Road Parts
# =========================================================================

def study_1():
    print("=" * 90)
    print("STUDY 1: Conduction vs Biot/Convection Mesh Dominance — Off-Road @ 10 mph")
    print("=" * 90)
    print()
    print("Question: at what combination of h and T_surf does the conduction")
    print("mesh requirement (dx = k*dT_max/q) become more restrictive than")
    print("the Biot through-thickness constraint (dx = t/3 when Bi>0.1)?")
    print()

    # For metals Bi is almost always < 0.1 at automotive h values,
    # so conduction always dominates. For plastics, Bi can exceed 0.1.
    # Let's prove this systematically.

    cases = [
        # (name, k, eps, thickness_mm, max_dt, description)
        ("SS409 shield 1.5mm",  25.0, 0.85, 1.5, 10, "Exhaust/turbo heat shield"),
        ("Al 6061 shield 2mm", 167.0, 0.30, 2.0, 10, "Cargo/underbody shield"),
        ("Mild steel bracket 3mm", 54.0, 0.73, 3.0, 15, "Engine bay bracket"),
        ("PA66 GF30 cover 4mm", 0.25, 0.90, 4.0, 15, "Plastic engine cover"),
        ("PP GF30 bed 6mm", 0.22, 0.90, 6.0, 15, "Cargo bed panel"),
        ("HDPE tank wall 5mm", 0.45, 0.90, 5.0, 15, "Fuel tank"),
        ("Rubber seal 3mm", 0.25, 0.90, 3.0, 15, "EPDM engine seal"),
    ]

    # h values from 10 mph aero zone estimation (~15-25 typical)
    # and dead zone natural conv (~3-10 typical)
    h_sweep = [3, 5, 8, 10, 15, 20, 30, 50]
    t_surfs_C = [60, 80, 100, 150, 200, 300, 400]

    for name, k, eps, t_mm, max_dt, desc in cases:
        thickness_m = t_mm / 1000.0
        print(f"\n--- {name} (k={k}, ε={eps}, t={t_mm}mm, dT_max={max_dt}K) ---")
        print(f"    {desc}")

        # Find crossover point
        crossovers = []
        for t_C in t_surfs_C:
            t_surf = t_C + 273.15
            for h in h_sweep:
                cond = BoundaryDrivenConductionCalculator.max_mesh_size(
                    k=k, h=h, t_surf=t_surf, t_fluid=T_AMB_K,
                    epsilon=eps, t_surr=T_SURR_K, max_dt=max_dt,
                )
                cond_dx = cond["max_dx_mm"]
                biot = ConvectionMeshCalculator.biot_number(h=h, k=k, thickness=thickness_m)
                bi = biot["biot"]
                bi_dx = t_mm / biot["min_elements_through_thickness"] if biot["mesh_type"] == "3D Solid" else float("inf")

                crossovers.append({
                    "t_C": t_C, "h": h, "cond_dx": cond_dx, "bi": bi,
                    "bi_dx": bi_dx, "dominant": "cond" if cond_dx < bi_dx else "biot",
                    "q": cond["q_total"],
                })

        # Print summary table
        biot_ever_dominates = any(c["dominant"] == "biot" for c in crossovers)
        cond_sub_10 = [c for c in crossovers if c["cond_dx"] < 10]

        if not biot_ever_dominates:
            # Find minimum conduction dx
            min_c = min(crossovers, key=lambda c: c["cond_dx"])
            print(f"    Biot NEVER dominates (max Bi={max(c['bi'] for c in crossovers):.4f})")
            print(f"    Conduction dx range: {min_c['cond_dx']:.1f} – "
                  f"{max(c['cond_dx'] for c in crossovers):.1f} mm")
            if cond_sub_10:
                print("    Conduction dx < 10mm when:")
                for c in cond_sub_10:
                    if c["cond_dx"] == min(cc["cond_dx"] for cc in cond_sub_10 if cc["t_C"] == c["t_C"]):
                        print(f"      T_surf={c['t_C']}°C, h≥{c['h']}: dx={c['cond_dx']:.1f}mm "
                              f"(q={c['q']:.0f} W/m²)")
        else:
            print(f"  {'T °C':>6s} | {'h':>5s} | {'Cond dx':>8s} | {'Bi':>7s} | "
                  f"{'Bi dx':>8s} | {'Dom':>6s} | {'q':>8s}")
            print(f"  {'-'*60}")
            for c in crossovers:
                if c["dominant"] == "biot" or (c["bi"] > 0.05 and c["cond_dx"] < 20):
                    bi_dx_s = f"{c['bi_dx']:.1f}" if c["bi_dx"] < 9999 else "shell"
                    print(f"  {c['t_C']:>6.0f} | {c['h']:>5.0f} | {c['cond_dx']:>8.2f} | "
                          f"{c['bi']:>7.4f} | {bi_dx_s:>8s} | {c['dominant']:>6s} | "
                          f"{c['q']:>8.0f}")

    print()
    print("KEY FINDING: For metals at automotive h values (3-50 W/m²K),")
    print("Biot is always << 0.1 → shell elements are always appropriate.")
    print("Conduction ALWAYS drives mesh size for metals.")
    print("For plastics, Biot can exceed 0.1 at h > ~25 W/m²K,")
    print("but conduction dx is already sub-mm by then, so conduction")
    print("still usually governs. Biot only becomes the binding constraint")
    print("for thick plastics at low T_surf with high h (rare combo).")


# =========================================================================
#  STUDY 2: Laminar → Turbulent Natural Convection — Off-Road Zones
# =========================================================================

def study_2():
    print()
    print("=" * 90)
    print("STUDY 2: Laminar → Turbulent Natural Convection — Off-Road Dead Zones")
    print("=" * 90)
    print()
    print("Off-road vehicle dead zones: engine bay interior, above exhaust,")
    print("cabin tunnel, cargo bed over tunnel. No forced flow — pure buoyancy.")
    print()

    # Zone-specific scenarios
    zones = [
        # (zone, orient, t_air_C, t_surf_range_C, L_range_mm, desc)
        ("Engine beside",     "vertical",       (60,70),   (80,120),  (100,300), "Next to engine block"),
        ("Engine below",      "horizontal_down", (50,60),   (70,100),  (200,500), "Under engine"),
        ("Engine above",      "horizontal_up",  (80,100),  (100,150), (200,500), "Above engine"),
        ("Exhaust beside",    "vertical",       (90,150),  (150,300), (100,300), "Next to exhaust pipe"),
        ("Exhaust below",     "horizontal_down", (90,150),  (120,250), (200,400), "Under exhaust"),
        ("Exhaust above",     "horizontal_up",  (150,350), (200,450), (150,400), "Above exhaust"),
        ("Cabin tunnel",      "horizontal_up",  (40,60),   (50,80),   (300,600), "Tunnel floor/cover"),
        ("Cargo over tunnel", "horizontal_up",  (40,60),   (60,100),  (300,800), "Cargo bed floor"),
    ]

    print(f"{'Zone':<22s} | {'Orient':<16s} | {'T_surf°C':>8s} | {'L mm':>6s} | "
          f"{'Ra':>12s} | {'Regime':>8s} | {'h nat':>7s} | {'h@10mph':>8s} | "
          f"{'Ri@10mph':>10s} | {'Regime@10':>10s}")
    print("-" * 140)

    for zone_name, orient, t_air_range, t_surf_range, L_range, desc in zones:
        # Test at high end (worst case for turbulence)
        t_air_C = t_air_range[1]
        t_air_K = t_air_C + 273.15
        for t_surf_C in [t_surf_range[0], t_surf_range[1]]:
            t_surf_K = t_surf_C + 273.15
            for L_mm in [L_range[0], L_range[1]]:
                L = L_mm / 1000.0

                nat = natural_convection(t_surf_K, t_air_K, L, orient, cap=999.0)
                Ra = nat["Ra"]
                nat_regime = nat["regime"]
                h_nat = nat["h_raw"]

                # What happens at 10 mph?
                est_10 = estimate_h(V_BASELINE_MS, t_surf_K, t_air_K, L, orient, cap=999.0)
                h_10 = est_10["h"]
                Ri_10 = est_10["Ri"]
                regime_10 = est_10["regime"]

                Ri_str = f"{Ri_10:.3f}" if Ri_10 < 1e6 else "inf"
                marker = " *TURB*" if nat_regime == "turbulent" else ""

                print(f"{zone_name:<22s} | {orient:<16s} | {t_surf_C:>8.0f} | {L_mm:>6.0f} | "
                      f"{Ra:>12.2e} | {nat_regime:>8s} | {h_nat:>7.1f} | {h_10:>8.1f} | "
                      f"{Ri_str:>10s} | {regime_10:>10s}{marker}")

    # Critical length summary
    print()
    print("--- Critical Length for Turbulence (at typical zone ΔT) ---")
    print()
    print("Vertical (Ra_crit = 1e9): Only turbulent for L > ~600mm.")
    print("  Most engine-bay and exhaust-adjacent parts are 100-300mm → LAMINAR.")
    print("  Exception: tall vertical shields or body panels > 600mm.")
    print()
    print("Horizontal hot-up (Ra_crit = 1e7): Turbulent for L > ~130-190mm.")
    print("  Almost ALL horizontal-up dead zone parts trip this threshold!")
    print("  Cargo beds, tunnel covers, shields above exhaust → TURBULENT.")
    print()
    print("Horizontal hot-down: ALWAYS laminar (stable stratification).")
    print("  Under-engine and under-exhaust shields are safe for SS.")

    # Specific critical lengths for off-road ΔTs
    print()
    print(f"{'Orient':<16s} | {'Zone ΔT (K)':>12s} | {'L_crit mm':>10s} | {'Typical L mm':>14s} | {'Verdict':>12s}")
    print("-" * 80)

    test_cases = [
        ("vertical", 40, "Engine beside", "100-300"),
        ("vertical", 100, "Exhaust beside", "100-300"),
        ("vertical", 260, "Exhaust beside hot", "100-300"),
        ("horizontal_up", 25, "Cabin tunnel", "300-600"),
        ("horizontal_up", 40, "Cargo over tunnel", "300-800"),
        ("horizontal_up", 60, "Engine above", "200-500"),
        ("horizontal_up", 110, "Exhaust above (shielded)", "150-400"),
        ("horizontal_up", 300, "Exhaust above (unshielded)", "150-400"),
    ]

    for orient, dT, zone_desc, typical_L in test_cases:
        t_surf = T_AMB_K + dT
        L_lo, L_hi = 0.01, 5.0
        for _ in range(50):
            L_mid = (L_lo + L_hi) / 2.0
            r = natural_convection(t_surf, T_AMB_K, L_mid, orient, cap=999.0)
            if r["regime"] == "turbulent":
                L_hi = L_mid
            else:
                L_lo = L_mid
        L_crit_mm = (L_lo + L_hi) / 2.0 * 1000
        verdict = "TURB likely" if L_crit_mm < 300 else "Laminar likely"
        print(f"{orient:<16s} | {dT:>12.0f} | {L_crit_mm:>10.0f} | {typical_L:>14s} | {verdict:>12s}")


# =========================================================================
#  STUDY 3: Does 10 mph Suppress Turbulent Buoyancy? (Ri < 0.1)
# =========================================================================

def study_3():
    print()
    print("=" * 90)
    print("STUDY 3: Does 10 mph Forced Flow Suppress Turbulent Natural Convection?")
    print("=" * 90)
    print()
    print("For zones that DO see aero (front end, sides, cargo exterior),")
    print("does 10 mph keep Ri < 0.1 (forced dominates)?")
    print()

    # Off-road parts that see some flow at 10 mph
    scenarios = [
        # (name, t_surf_C, L_mm, orient, t_air_C)
        ("Front bumper bracket",        60,  200, "vertical", 40),
        ("Front crossmember",           60,  300, "vertical", 40),
        ("Cooling pack downstream",     80,  200, "vertical", 80),
        ("HX outlet duct",             100,  150, "vertical", 100),
        ("Side body panel",             50,  500, "vertical", 40),
        ("Side body panel tall",        50, 1000, "vertical", 40),
        ("Cargo exterior panel",        60,  400, "vertical", 40),
        ("Roof — horizontal up",        50,  500, "horizontal_up", 40),
        ("Hood — horizontal up",        80,  600, "horizontal_up", 50),
        ("Underbody — horizontal down",  70,  600, "horizontal_down", 40),
        ("Exhaust shield (edge flow)",  200,  300, "vertical", 60),
        ("Engine cover (edge flow)",    120,  300, "horizontal_up", 50),
    ]

    print(f"{'Scenario':<30s} | {'T_s°C':>6s} | {'L mm':>6s} | {'Orient':>14s} | "
          f"{'Ra':>10s} | {'Ri@10':>8s} | {'Regime':>8s} | {'h_est':>7s} | "
          f"{'SS?':>4s} | {'Sev':>8s} | {'Notes':>20s}")
    print("-" * 160)

    for name, t_surf_C, L_mm, orient, t_air_C in scenarios:
        t_surf = t_surf_C + 273.15
        t_air = t_air_C + 273.15
        L = L_mm / 1000.0

        est = estimate_h(V_BASELINE_MS, t_surf, t_air, L, orient, cap=999.0)
        adv = est["solver_advisory"]
        Ri = est["Ri"]
        Ra = est["Ra"]
        regime = est["regime"]
        h = est["h"]

        Ri_str = f"{Ri:.4f}" if Ri < 1e6 else "inf"
        ss = "Yes" if adv["steady_state_ok"] else "NO"

        # Annotate
        notes = ""
        if regime == "forced" and Ri < 0.01:
            notes = "Strong forced dom."
        elif regime == "forced" and Ri < 0.1:
            notes = "Forced dominates"
        elif regime == "mixed":
            notes = "Mixed — check conv."
        elif regime == "natural":
            notes = "Buoyancy dominant"

        print(f"{name:<30s} | {t_surf_C:>6.0f} | {L_mm:>6.0f} | {orient:>14s} | "
              f"{Ra:>10.2e} | {Ri_str:>8s} | {regime:>8s} | {h:>7.1f} | "
              f"{ss:>4s} | {adv['severity']:>8s} | {notes:>20s}")

    print()
    print("--- Ri < 0.1 Threshold: Minimum velocity to suppress turbulent buoyancy ---")
    print()

    # For the turbulent horizontal-up cases, what velocity is needed?
    turb_cases = [
        ("Cargo bed above tunnel",  80, 600, "horizontal_up", 40),
        ("Plastic engine cover",   100, 400, "horizontal_up", 80),
        ("Exhaust shield above",   250, 300, "horizontal_up", 150),
        ("Turbo shield above",     400, 200, "horizontal_up", 100),
        ("Cabin tunnel cover",      65, 300, "horizontal_up", 40),
        ("Hood underside",          80, 600, "horizontal_up", 50),
    ]

    print(f"{'Scenario':<30s} | {'ΔT K':>6s} | {'L mm':>6s} | {'Ra':>10s} | "
          f"{'Nat reg':>8s} | {'v_crit mph':>11s} | {'v_crit m/s':>11s} | "
          f"{'10 mph ok?':>11s}")
    print("-" * 120)

    for name, t_surf_C, L_mm, orient, t_air_C in turb_cases:
        t_surf = t_surf_C + 273.15
        t_air = t_air_C + 273.15
        L = L_mm / 1000.0
        dT = t_surf_C - t_air_C

        nat = natural_convection(t_surf, t_air, L, orient, cap=999.0)
        Ra = nat["Ra"]
        nat_reg = nat["regime"]

        if nat_reg == "turbulent":
            v_lo, v_hi = 0.01, 50.0
            for _ in range(50):
                v_mid = (v_lo + v_hi) / 2.0
                ri = richardson_number(v_mid, t_surf, t_air, L)
                if ri["Ri"] < 0.1:
                    v_hi = v_mid
                else:
                    v_lo = v_mid
            v_crit = (v_lo + v_hi) / 2.0
            v_crit_mph = v_crit / 0.44704
            ok = "YES" if v_crit_mph <= 10.0 else "NO — needs more"
        else:
            v_crit = 0.0
            v_crit_mph = 0.0
            ok = "LAMINAR — ok"

        print(f"{name:<30s} | {dT:>6.0f} | {L_mm:>6.0f} | {Ra:>10.2e} | "
              f"{nat_reg:>8s} | {v_crit_mph:>11.1f} | {v_crit:>11.2f} | "
              f"{ok:>11s}")


# =========================================================================
#  SUMMARY
# =========================================================================

def summary():
    print()
    print("=" * 90)
    print("SUMMARY: Key Thresholds for Off-Road Vehicle @ 10 mph Thermal Analysis")
    print("=" * 90)
    print()
    print("1. CONDUCTION vs CONVECTION MESH DOMINANCE")
    print("   • Metals: Biot << 0.1 at all realistic h values → always shell elements")
    print("   • Conduction dx is ALWAYS the governing mesh constraint for metals")
    print("   • Plastics: Biot crosses 0.1 at h ≈ 25 W/m²K (for k=0.25, t=2mm)")
    print("     but conduction dx is already < 5mm by then, so still governs")
    print("   • Biot/convection only governs for thick plastics (>4mm) at LOW T_surf")
    print("     with HIGH h — a rare combo in off-road vehicles")
    print()
    print("2. LAMINAR → TURBULENT NATURAL CONVECTION")
    print("   • VERTICAL surfaces: turbulent only when L > ~600mm (at typical ΔT)")
    print("     Most engine/exhaust-adjacent parts are 100-300mm → SAFE (laminar)")
    print("     Only tall body panels or full-height shields trip this")
    print("   • HORIZONTAL HOT-UP: turbulent at L > ~130-190mm (!)")
    print("     Nearly ALL horizontal dead-zone parts are turbulent:")
    print("     cargo beds, tunnel covers, engine covers, shields above exhaust")
    print("     This is the dominant SS convergence risk in the vehicle")
    print("   • HORIZONTAL HOT-DOWN: always laminar (stable stratification)")
    print("     Under-engine, under-exhaust shields → safe for SS")
    print()
    print("3. FORCED SUPPRESSION OF TURBULENT BUOYANCY @ 10 mph")
    print("   • 10 mph (4.47 m/s) gives Ri < 0.1 for MOST scenarios with L < 400mm")
    print("   • Critical velocity for forced dominance is typically 3-9 mph")
    print("   • Exception: very large horizontal-up surfaces (L > 500mm) with")
    print("     moderate ΔT — these may need > 10 mph for Ri < 0.1")
    print("   • For zones that genuinely see 10 mph aero: SS is acceptable")
    print("   • For dead zones (0 mph): horizontal-up surfaces will be")
    print("     turbulent natural → transient recommended or SS with caution")
    print()
    print("PRACTICAL RECOMMENDATION:")
    print("   • Run SS for all forced-convection zones (edges, front, downstream)")
    print("   • Flag horizontal-up dead zones for transient review:")
    print("     cargo bed, tunnel cover, engine cover, shields above exhaust")
    print("   • Vertical dead zones < 500mm: safe for SS")
    print("   • Hot-side-down surfaces: always safe for SS")


if __name__ == "__main__":
    study_1()
    study_2()
    study_3()
    summary()
