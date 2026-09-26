"""
Regime Crossover Parametric Study
==================================

Sweeps realistic automotive parameter ranges to identify:
  1. Conduction vs convection mesh dominance crossover
  2. Laminar → turbulent natural convection boundaries
  3. SS-acceptable turbulent natural conv (forced dominates)

Uses the actual thermal_mesh_calculators functions — no approximations.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from thermal_mesh_calculators.conduction import BoundaryDrivenConductionCalculator
from thermal_mesh_calculators.convection import ConvectionMeshCalculator
from thermal_mesh_calculators.h_estimator import (
    natural_convection,
    richardson_number,
    estimate_h,
    solver_advisory,
)
from thermal_mesh_calculators.batch import MATERIALS


# =========================================================================
#  STUDY 1: Conduction vs Convection Mesh Dominance
# =========================================================================

def study_1_conduction_vs_convection():
    """
    When does the conduction mesh requirement (driven by boundary flux)
    become more restrictive than the convection mesh requirement (driven
    by Biot number / through-thickness resolution)?

    Conduction dx_max = k * dT_max / q_boundary
    Convection says: if Bi < 0.1, shell is fine (no through-thickness constraint).
                     if Bi > 0.1, need 3+ elements through thickness.

    So the question is: for a given material + thickness, at what h and
    T_surf does conduction dx become smaller than thickness/3 (Biot-driven)?
    """
    print("=" * 80)
    print("STUDY 1: Conduction vs Convection Mesh Dominance Crossover")
    print("=" * 80)
    print()

    # Representative materials
    materials = {
        "steel_mild":          MATERIALS["steel_mild"],
        "steel_stainless_409": MATERIALS["steel_stainless_409"],
        "aluminium_6061":      MATERIALS["aluminium_6061"],
        "plastic_pa66_gf30":   MATERIALS["plastic_pa66_gf30"],
    }

    # Sweep parameters
    h_values = [2, 5, 10, 15, 20, 30, 50, 75, 100]  # W/m²K
    t_surfs_C = [60, 100, 150, 200, 300, 500]  # surface temp, °C
    t_fluid_K = 313.15   # 40°C ambient
    t_surr_K = 313.15    # 40°C radiation sink
    max_dt = 10.0        # K per element
    thickness_mm = 2.0   # typical shield/bracket

    print(f"Fixed: t_fluid={t_fluid_K-273.15:.0f}°C, t_surr={t_surr_K-273.15:.0f}°C, "
          f"max_dT={max_dt}K, thickness={thickness_mm}mm")
    print()

    for mat_name, mat in materials.items():
        k = mat["k"]
        eps = 0.73 if "steel" in mat_name else (0.30 if "aluminium" in mat_name else 0.90)
        thickness_m = thickness_mm / 1000.0

        print(f"--- {mat_name} (k={k} W/mK, ε={eps}) ---")
        print(f"{'T_surf °C':>10s} | {'h':>5s} | {'Cond dx':>8s} | {'Biot':>7s} | "
              f"{'Bi mesh':>8s} | {'Dominant':>12s} | {'q_total':>10s}")
        print("-" * 80)

        for t_surf_C in t_surfs_C:
            t_surf = t_surf_C + 273.15
            for h in h_values:
                # Conduction mesh limit
                cond = BoundaryDrivenConductionCalculator.max_mesh_size(
                    k=k, h=h, t_surf=t_surf, t_fluid=t_fluid_K,
                    epsilon=eps, t_surr=t_surr_K, max_dt=max_dt,
                )
                cond_dx = cond["max_dx_mm"]

                # Biot number → through-thickness constraint
                biot = ConvectionMeshCalculator.biot_number(
                    h=h, k=k, thickness=thickness_m,
                )
                bi_val = biot["biot"]
                if biot["mesh_type"] == "3D Solid":
                    bi_dx = thickness_mm / biot["min_elements_through_thickness"]
                else:
                    bi_dx = float("inf")  # shell — no through-thickness limit

                dominant = "conduction" if cond_dx < bi_dx else "biot/conv"
                if bi_dx == float("inf"):
                    dominant = "conduction" if cond_dx < 999 else "neither"

                q_total = cond["q_total"]

                # Only print interesting crossover region
                if cond_dx < 50 or bi_dx < 50:
                    print(f"{t_surf_C:>10.0f} | {h:>5.0f} | {cond_dx:>8.2f} | "
                          f"{bi_val:>7.4f} | {bi_dx:>8.2f} | {dominant:>12s} | "
                          f"{q_total:>10.1f}")

        print()


# =========================================================================
#  STUDY 2: Laminar → Turbulent Natural Convection Boundaries
# =========================================================================

def study_2_natural_convection_transition():
    """
    Map the Ra transition boundaries for each orientation.

    Vertical:        Ra_crit = 1e9
    Horizontal up:   Ra_crit = 1e7
    Horizontal down: always laminar

    For a given ΔT, what characteristic length triggers turbulence?
    For a given length, what ΔT triggers turbulence?
    """
    print("=" * 80)
    print("STUDY 2: Laminar → Turbulent Natural Convection Boundaries")
    print("=" * 80)
    print()

    t_fluid_K = 313.15  # 40°C

    # Sweep: ΔT × char_length × orientation
    delta_Ts = [20, 40, 60, 80, 100, 150, 200, 300, 400, 500]  # K
    char_lengths_mm = [50, 100, 150, 200, 300, 500, 750, 1000]  # mm
    orientations = ["vertical", "horizontal_up", "horizontal_down"]

    for orient in orientations:
        print(f"\n--- Orientation: {orient} ---")
        print(f"{'ΔT (K)':>8s} | {'L (mm)':>8s} | {'Ra':>12s} | {'Regime':>10s} | "
              f"{'h (W/m²K)':>10s} | {'Nu':>8s}")
        print("-" * 72)

        for dT in delta_Ts:
            t_surf = t_fluid_K + dT
            for L_mm in char_lengths_mm:
                L = L_mm / 1000.0
                result = natural_convection(
                    t_surf=t_surf, t_fluid=t_fluid_K,
                    char_length=L, orientation=orient, cap=999.0,
                )
                Ra = result["Ra"]
                regime = result["regime"]
                h = result["h_raw"]
                Nu = result["Nu"]

                # Only print near transitions and key data points
                if orient == "vertical" and (1e8 < Ra < 1e10):
                    marker = " <<<" if regime == "turbulent" else ""
                    print(f"{dT:>8.0f} | {L_mm:>8.0f} | {Ra:>12.2e} | "
                          f"{regime:>10s} | {h:>10.2f} | {Nu:>8.1f}{marker}")
                elif orient == "horizontal_up" and (1e6 < Ra < 1e8):
                    marker = " <<<" if regime == "turbulent" else ""
                    print(f"{dT:>8.0f} | {L_mm:>8.0f} | {Ra:>12.2e} | "
                          f"{regime:>10s} | {h:>10.2f} | {Nu:>8.1f}{marker}")
                elif orient == "horizontal_down":
                    if L_mm in [100, 300, 500] and dT in [40, 100, 300, 500]:
                        print(f"{dT:>8.0f} | {L_mm:>8.0f} | {Ra:>12.2e} | "
                              f"{regime:>10s} | {h:>10.2f} | {Nu:>8.1f}")

    # Now find critical lengths for typical automotive ΔT values
    print("\n\n--- Critical Characteristic Length for Turbulence Onset ---")
    print(f"{'Orient':>15s} | {'ΔT (K)':>8s} | {'T_surf °C':>10s} | "
          f"{'L_crit (mm)':>12s} | {'L_crit (m)':>10s}")
    print("-" * 70)

    for orient in ["vertical", "horizontal_up"]:
        for dT in [20, 40, 60, 80, 100, 150, 200, 300, 500]:
            t_surf = t_fluid_K + dT
            # Binary search for critical length
            L_lo, L_hi = 0.01, 5.0  # meters
            for _ in range(50):
                L_mid = (L_lo + L_hi) / 2.0
                r = natural_convection(t_surf, t_fluid_K, L_mid, orient, cap=999.0)
                if r["regime"] == "turbulent":
                    L_hi = L_mid
                else:
                    L_lo = L_mid
            L_crit = (L_lo + L_hi) / 2.0
            t_surf_C = t_surf - 273.15
            print(f"{orient:>15s} | {dT:>8.0f} | {t_surf_C:>10.0f} | "
                  f"{L_crit*1000:>12.0f} | {L_crit:>10.3f}")


# =========================================================================
#  STUDY 3: SS-Acceptable Turbulent Natural Conv (Forced Dominates)
# =========================================================================

def study_3_forced_suppresses_turbulent_natural():
    """
    When turbulent natural convection is present but forced convection
    dominates (Ri < 0.1), steady-state solve is acceptable.

    Question: for a given buoyancy scenario (ΔT, L, orientation), what
    minimum forced velocity keeps Ri < 0.1? And what's Ri at typical
    automotive velocities?
    """
    print("\n" + "=" * 80)
    print("STUDY 3: Forced Convection Suppression of Turbulent Natural Conv")
    print("=" * 80)
    print()

    t_fluid_K = 313.15  # 40°C

    # Typical automotive scenarios
    scenarios = [
        # (name, t_surf_C, L_mm, orient, zone_desc)
        ("Exhaust shield above", 250, 300, "horizontal_up",
         "Aluminium shield 300mm above exhaust"),
        ("Exhaust shield beside", 150, 200, "vertical",
         "SS409 shield 200mm beside exhaust"),
        ("Engine cover above", 100, 400, "horizontal_up",
         "Plastic cover 400mm above engine"),
        ("Cargo bed floor", 80, 600, "horizontal_up",
         "Plastic cargo bed 600mm over tunnel"),
        ("A-pillar bracket", 60, 150, "vertical",
         "Steel bracket near engine bay edge"),
        ("Upper cargo shield", 120, 500, "horizontal_down",
         "Aluminium shield, hot side down"),
        ("Exhaust manifold zone", 400, 300, "vertical",
         "SS304 bracket near exhaust manifold"),
        ("Turbo heat shield", 350, 200, "horizontal_up",
         "SS409 shield above turbo"),
    ]

    print("--- Scenario Analysis ---")
    print()

    for name, t_surf_C, L_mm, orient, desc in scenarios:
        t_surf = t_surf_C + 273.15
        L = L_mm / 1000.0

        # Natural convection (no forced flow)
        nat = natural_convection(t_surf, t_fluid_K, L, orient, cap=999.0)
        Ra = nat["Ra"]
        nat_regime = nat["regime"]
        h_nat = nat["h_raw"]

        print(f"Scenario: {name}")
        print(f"  {desc}")
        print(f"  T_surf={t_surf_C}°C, L={L_mm}mm, orient={orient}")
        print(f"  Natural conv: Ra={Ra:.2e}, regime={nat_regime}, h_nat={h_nat:.1f} W/m²K")

        # Sweep velocities to find Ri < 0.1 threshold
        velocities_mph = [0, 2, 5, 10, 15, 20, 30]
        print(f"  {'v (mph)':>8s} | {'v (m/s)':>8s} | {'Ri':>10s} | {'Regime':>8s} | "
              f"{'h_est':>8s} | {'SS ok?':>7s} | {'Severity':>10s}")
        print(f"  {'-'*72}")

        for v_mph in velocities_mph:
            v_ms = v_mph * 0.44704
            if v_ms < 0.01:
                # Pure natural
                adv = solver_advisory(natural_result=nat, regime="natural",
                                      orientation=orient)
                print(f"  {v_mph:>8.0f} | {v_ms:>8.2f} | {'inf':>10s} | "
                      f"{'natural':>8s} | {h_nat:>8.1f} | "
                      f"{'Yes' if adv['steady_state_ok'] else 'NO':>7s} | "
                      f"{adv['severity']:>10s}")
            else:
                est = estimate_h(v_ms, t_surf, t_fluid_K, L, orient, cap=999.0)
                adv = est["solver_advisory"]
                Ri = est["Ri"]
                Ri_str = f"{Ri:.3f}" if Ri < 1e6 else "inf"
                print(f"  {v_mph:>8.0f} | {v_ms:>8.2f} | {Ri_str:>10s} | "
                      f"{est['regime']:>8s} | {est['h']:>8.1f} | "
                      f"{'Yes' if adv['steady_state_ok'] else 'NO':>7s} | "
                      f"{adv['severity']:>10s}")

        # Find minimum velocity for Ri < 0.1 (forced dominance)
        if nat_regime == "turbulent":
            # Binary search for critical velocity
            v_lo, v_hi = 0.01, 50.0
            for _ in range(50):
                v_mid = (v_lo + v_hi) / 2.0
                ri = richardson_number(v_mid, t_surf, t_fluid_K, L)
                if ri["Ri"] < 0.1:
                    v_hi = v_mid
                else:
                    v_lo = v_mid
            v_crit_ms = (v_lo + v_hi) / 2.0
            v_crit_mph = v_crit_ms / 0.44704
            print(f"  → Min velocity for forced dominance (Ri<0.1): "
                  f"{v_crit_ms:.2f} m/s = {v_crit_mph:.1f} mph")
        else:
            print("  → Natural conv is laminar — SS ok without forced flow")

        print()


# =========================================================================
#  STUDY 4: Real Vehicle Worked Examples
# =========================================================================

def study_4_vehicle_examples():
    """
    Concrete automotive scenarios tying all three studies together.
    """
    print("=" * 80)
    print("STUDY 4: Integrated Vehicle Examples")
    print("=" * 80)
    print()

    t_amb = 313.15  # 40°C ambient
    t_surr = 313.15

    examples = [
        {
            "name": "Cargo shield — aluminium 2mm, above exhaust tunnel",
            "material": "aluminium_6061",
            "k": 167.0, "eps": 0.30,
            "t_surf_C": 150, "thickness_mm": 2.0,
            "L_mm": 400, "orient": "horizontal_up",
            "zone_velocity_mph": 0,  # dead zone
            "max_dt": 10.0,
        },
        {
            "name": "Cargo bed — plastic PP 6mm, above tunnel",
            "material": "plastic_pp_gf30",
            "k": 0.22, "eps": 0.90,
            "t_surf_C": 80, "thickness_mm": 6.0,
            "L_mm": 600, "orient": "horizontal_up",
            "zone_velocity_mph": 0,
            "max_dt": 15.0,
        },
        {
            "name": "Exhaust heat shield — SS409 1.5mm, beside exhaust",
            "material": "steel_stainless_409",
            "k": 25.0, "eps": 0.85,
            "t_surf_C": 300, "thickness_mm": 1.5,
            "L_mm": 250, "orient": "vertical",
            "zone_velocity_mph": 0,
            "max_dt": 10.0,
        },
        {
            "name": "Engine bay bracket — mild steel 3mm, bay edge",
            "material": "steel_mild",
            "k": 54.0, "eps": 0.73,
            "t_surf_C": 80, "thickness_mm": 3.0,
            "L_mm": 150, "orient": "vertical",
            "zone_velocity_mph": 5,  # some edge flow
            "max_dt": 15.0,
        },
        {
            "name": "Front crossmember — galv steel 2mm, aero zone",
            "material": "steel_galvanised",
            "k": 52.0, "eps": 0.23,
            "t_surf_C": 60, "thickness_mm": 2.0,
            "L_mm": 200, "orient": "vertical",
            "zone_velocity_mph": 10,  # full aero
            "max_dt": 15.0,
        },
        {
            "name": "Turbo heat shield — SS409 1mm, above turbo",
            "material": "steel_stainless_409",
            "k": 25.0, "eps": 0.85,
            "t_surf_C": 400, "thickness_mm": 1.0,
            "L_mm": 200, "orient": "horizontal_up",
            "zone_velocity_mph": 0,
            "max_dt": 10.0,
        },
        {
            "name": "Cabin tunnel cover — plastic PA66 GF30 4mm",
            "material": "plastic_pa66_gf30",
            "k": 0.25, "eps": 0.90,
            "t_surf_C": 65, "thickness_mm": 4.0,
            "L_mm": 300, "orient": "horizontal_up",
            "zone_velocity_mph": 0,
            "max_dt": 15.0,
        },
        {
            "name": "Exhaust clamp — SS409 2mm, near manifold",
            "material": "steel_stainless_409",
            "k": 25.0, "eps": 0.85,
            "t_surf_C": 500, "thickness_mm": 2.0,
            "L_mm": 100, "orient": "vertical",
            "zone_velocity_mph": 0,
            "max_dt": 10.0,
        },
    ]

    print(f"{'Example':<55s} | {'Cond dx':>8s} | {'Bi dx':>8s} | {'Dom':>10s} | "
          f"{'h':>6s} | {'Ra':>10s} | {'Nat':>6s} | {'SS?':>4s} | {'Sev':>8s}")
    print("-" * 140)

    for ex in examples:
        t_surf = ex["t_surf_C"] + 273.15
        k = ex["k"]
        eps = ex["eps"]
        thickness_m = ex["thickness_mm"] / 1000.0
        L = ex["L_mm"] / 1000.0
        v_ms = ex["zone_velocity_mph"] * 0.44704

        # h estimation
        est = estimate_h(v_ms, t_surf, t_amb, L, ex["orient"], cap=999.0)
        h = est["h"]
        adv = est["solver_advisory"]

        # Conduction mesh
        cond = BoundaryDrivenConductionCalculator.max_mesh_size(
            k=k, h=h, t_surf=t_surf, t_fluid=t_amb,
            epsilon=eps, t_surr=t_surr, max_dt=ex["max_dt"],
        )
        cond_dx = cond["max_dx_mm"]

        # Biot
        biot = ConvectionMeshCalculator.biot_number(h=h, k=k, thickness=thickness_m)
        if biot["mesh_type"] == "3D Solid":
            bi_dx = ex["thickness_mm"] / biot["min_elements_through_thickness"]
        else:
            bi_dx = float("inf")

        dominant = "conduction" if cond_dx < bi_dx else "biot"
        if bi_dx == float("inf"):
            dominant = "conduction"

        Ra = est.get("Ra", 0)
        nat_reg = est.get("natural_regime", "n/a")
        ss_ok = "Yes" if adv["steady_state_ok"] else "NO"

        bi_dx_str = f"{bi_dx:.1f}" if bi_dx < 9999 else "shell"

        print(f"{ex['name']:<55s} | {cond_dx:>8.1f} | {bi_dx_str:>8s} | "
              f"{dominant:>10s} | {h:>6.1f} | {Ra:>10.2e} | "
              f"{nat_reg:>6s} | {ss_ok:>4s} | {adv['severity']:>8s}")

    print()
    print("--- Detailed Breakdown ---")
    print()

    for ex in examples:
        t_surf = ex["t_surf_C"] + 273.15
        k = ex["k"]
        eps = ex["eps"]
        thickness_m = ex["thickness_mm"] / 1000.0
        L = ex["L_mm"] / 1000.0
        v_ms = ex["zone_velocity_mph"] * 0.44704

        est = estimate_h(v_ms, t_surf, t_amb, L, ex["orient"], cap=999.0)
        h = est["h"]
        adv = est["solver_advisory"]

        cond = BoundaryDrivenConductionCalculator.max_mesh_size(
            k=k, h=h, t_surf=t_surf, t_fluid=t_amb,
            epsilon=eps, t_surr=t_surr, max_dt=ex["max_dt"],
        )

        biot = ConvectionMeshCalculator.biot_number(h=h, k=k, thickness=thickness_m)

        print(f"  {ex['name']}")
        print(f"    Material: {ex['material']}, k={k}, ε={eps}")
        print(f"    T_surf={ex['t_surf_C']}°C, t={ex['thickness_mm']}mm, "
              f"L={ex['L_mm']}mm, orient={ex['orient']}")
        print(f"    Flow: {ex['zone_velocity_mph']} mph → h={h:.1f} W/m²K "
              f"(regime={est['regime']}, {est['dominant_mode']})")
        print(f"    Ra={est.get('Ra',0):.2e}, natural_regime={est.get('natural_regime','n/a')}")
        print(f"    Conduction dx={cond['max_dx_mm']:.1f}mm "
              f"(q_total={cond['q_total']:.0f} W/m², "
              f"conv={cond['q_conv']:.0f}, rad={cond['q_rad']:.0f})")
        print(f"    Biot={biot['biot']:.4f} → {biot['mesh_type']}")
        print(f"    Solver: SS_ok={adv['steady_state_ok']}, "
              f"severity={adv['severity']}")
        if not adv["steady_state_ok"]:
            print(f"    ⚠ {adv['reason']}")
        print()


# =========================================================================
#  Main
# =========================================================================

if __name__ == "__main__":
    study_1_conduction_vs_convection()
    study_2_natural_convection_transition()
    study_3_forced_suppresses_turbulent_natural()
    study_4_vehicle_examples()
