#!/usr/bin/env python3
"""
Batch Mesh Sizing — Off-Road Side-by-Side Vehicle
===================================================

Demonstrates the batch processor on a representative set of underhood
and cargo area components for a side-by-side off-road vehicle.

Per-part input is minimal:
    part_id, material, surface treatment, component_class,
    convection_zone, thickness_mm, t_surf_K

Project-level settings are defined once and applied to all parts.

Run:  python -m examples.batch_example
"""

from thermal_mesh_calculators.batch import (
    process_batch,
    summary_table,
    list_materials,
    list_surface_treatments,
)
from thermal_mesh_calculators.zones import list_zones


def main():
    # ------------------------------------------------------------------
    #  Project-level settings (set once per analysis)
    # ------------------------------------------------------------------
    project = {
        "t_fluid_K": 353.15,        # 80 C underhood air
        "t_surr_K": 353.15,         # radiation sink ~ underhood air
        "t_exh_K": 1073.15,         # 800 C exhaust surface (shields)
        "max_dt": None,             # use class defaults
        "dt": 0.5,                  # transient time-step (s)
        "fo_max": 0.5,              # explicit solver
        "safety_factor": 1.0,
        "tau_bc": 20.0,             # drive-cycle segment (s)
    }

    # ------------------------------------------------------------------
    #  Component list — this is what a BOM spreadsheet would produce
    # ------------------------------------------------------------------
    parts = [
        # --- Exhaust system ---
        {
            "part_id": "EXH-001 Exhaust Manifold",
            "material": "cast_iron",
            "surface": "cast_iron_oxidised",
            "component_class": "exhaust",
            "convection_zone": "exhaust_internal",
            "thickness_mm": 6.0,
            "t_surf_K": 1073.15,    # 800 C — fixed BC
        },
        {
            "part_id": "EXH-002 Downpipe",
            "material": "steel_stainless_304",
            "surface": "heavily_oxidised",
            "component_class": "exhaust",
            "convection_zone": "exhaust_internal",
            "thickness_mm": 2.0,
            "t_surf_K": 973.15,     # 700 C
        },

        # --- Engine bay structural ---
        {
            "part_id": "STR-001 Engine Mount Bracket",
            "material": "steel_mild",
            "surface": "painted",
            "component_class": "structural",
            "convection_zone": "engine_bay_dead_zone",
            "thickness_mm": 3.0,
            "t_surf_K": 473.15,     # 200 C estimate
            "char_length_mm": 200.0,
        },
        {
            "part_id": "STR-002 Oil Pan",
            "material": "cast_aluminium",
            "surface": "cast_aluminium_bare",
            "component_class": "structural",
            "convection_zone": "engine_bay_dead_zone",
            "thickness_mm": 4.0,
            "t_surf_K": 423.15,     # 150 C
            "char_length_mm": 300.0,
        },

        # --- Front end ---
        {
            "part_id": "FE-001 Radiator Support",
            "material": "plastic_pa66_gf30",
            "component_class": "structural",
            "convection_zone": "cooling_pack_downstream",
            "thickness_mm": 3.5,
            "t_surf_K": 383.15,     # 110 C
            "char_length_mm": 250.0,
            # No "surface" key → falls back to plastic default (eps=0.92)
        },
        {
            "part_id": "FE-002 Charge Air Cooler Bracket",
            "material": "aluminium_6061",
            "surface": "anodised",
            "component_class": "structural",
            "convection_zone": "cooling_pack_downstream",
            "thickness_mm": 2.5,
            "t_surf_K": 393.15,     # 120 C
        },

        # --- Cabin / tunnel ---
        {
            "part_id": "CAB-001 Tunnel Crossmember",
            "material": "steel_mild",
            "surface": "painted",
            "component_class": "structural",
            "convection_zone": "cabin_tunnel",
            "thickness_mm": 2.0,
            "t_surf_K": 373.15,     # 100 C (mild, mostly dead air)
            "char_length_mm": 400.0,
        },

        # --- Heat shields ---
        {
            "part_id": "SH-001 Manifold Shield",
            "material": "steel_mild",
            "component_class": "shield",
            "convection_zone_in": "near_exhaust_natural",
            "convection_zone_out": "engine_bay_dead_zone",
            "thickness_mm": 0.8,
            "t_surf_K": 0,           # ignored — solved by shield calculator
            "surface_in": "aluminised",
            "surface_out": "aluminised",
        },
        {
            "part_id": "SH-002 Cargo Bed Shield",
            "material": "steel_mild",
            "component_class": "shield",
            "convection_zone_in": "near_exhaust_natural",
            "convection_zone_out": "shield_gap_confined",
            "thickness_mm": 0.6,
            "t_surf_K": 0,
            "surface_in": "aluminised",
            "surface_out": "lightly_oxidised",
        },

        # --- Multilayer shield ---
        {
            "part_id": "SH-003 Dual-Wall Cargo Shield",
            "material": "steel_mild",
            "component_class": "multilayer_shield",
            "convection_zone_in": "near_exhaust_natural",
            "convection_zone_out": "shield_gap_confined",
            "thickness_mm": 0.6,
            "t_surf_K": 0,
            "surface_in": "aluminised",
            "surface_out": "aluminised",
            "surface_g1": "aluminised",
            "surface_g2": "aluminised",
            "h_gap": 12.0,
        },

        # --- Plastic cargo bed ---
        {
            "part_id": "CARGO-001 Bed Liner",
            "material": "plastic_hdpe",
            "component_class": "structural",
            "convection_zone": "cargo_bed_exterior",
            "thickness_mm": 5.0,
            "t_surf_K": 373.15,     # 100 C estimate
            "char_length_mm": 350.0,
        },

        # --- Rubber / soft parts ---
        {
            "part_id": "SOFT-001 Exhaust Grommet",
            "material": "rubber_epdm",
            "component_class": "exhaust_adjacent",
            "convection_zone": "near_exhaust_natural",
            "thickness_mm": 8.0,
            "t_surf_K": 423.15,     # 150 C — near exhaust
            "char_length_mm": 150.0,
        },
    ]

    # ------------------------------------------------------------------
    #  Run batch
    # ------------------------------------------------------------------
    print("=" * 80)
    print("  AVAILABLE LOOKUPS")
    print("=" * 80)
    print(f"  Materials ({len(list_materials())}):  {', '.join(list_materials())}")
    print(f"  Surfaces ({len(list_surface_treatments())}):   {', '.join(list_surface_treatments())}")
    print(f"  Zones ({len(list_zones())}):      {', '.join(list_zones())}")
    print()

    results = process_batch(parts, project)

    print("=" * 80)
    print("  BATCH MESH SIZING RESULTS — Off-Road Side-by-Side")
    print("=" * 80)
    print()
    print(summary_table(results))
    print()

    # ------------------------------------------------------------------
    #  Detailed breakdown for a few interesting parts
    # ------------------------------------------------------------------
    for r in results:
        if r.get("error"):
            continue
        if r["part_id"] in (
            "STR-001 Engine Mount Bracket",
            "EXH-001 Exhaust Manifold",
            "FE-001 Radiator Support",
            "CARGO-001 Bed Liner",
            "SH-001 Manifold Shield",
            "SH-003 Dual-Wall Cargo Shield",
            "SOFT-001 Exhaust Grommet",
        ):
            print("-" * 60)
            print(f"  {r['part_id']}")
            print("-" * 60)
            print(f"    Material:    {r['material']}")
            print(f"    Class:       {r['component_class']}")
            print(f"    h used:      {r['h_used']}")
            print(f"    eps used:    {r.get('eps_used', '—')}")

            if r.get("conduction"):
                c = r["conduction"]
                print(f"    Conduction:  dx={c['max_dx_mm']:.2f} mm  "
                      f"(q={c['q_total']:.0f} W/m^2, rad_frac={c['rad_fraction']:.1%})")

            if r.get("biot"):
                b = r["biot"]
                print(f"    Biot:        {b['biot']:.4f} → {b['mesh_type']}")

            if r.get("radiation"):
                rd = r["radiation"]
                print(f"    Radiation:   dx={rd['max_dx_mm']:.2f} mm  "
                      f"(dq/dT={rd['dq_dt']:.1f} W/m^2K)")

            if r.get("shield"):
                s = r["shield"]
                if "t_shield_C" in s:
                    print(f"    Shield:      T={s['t_shield_C']:.1f} C, "
                          f"dx={s['max_dx_mm']:.2f} mm")
                elif "t1_C" in s:
                    print(f"    Shield L1:   T={s['t1_C']:.1f} C, "
                          f"dx={s['layer1_max_dx_mm']:.2f} mm")
                    print(f"    Shield L2:   T={s['t2_C']:.1f} C, "
                          f"dx={s['layer2_max_dx_mm']:.2f} mm")

            if r.get("transient"):
                t = r["transient"]
                print(f"    Transient:   dx={t['recommended_dx_mm']:.2f} mm  "
                      f"({t['binding_constraint']})")

            if r.get("all_constraints"):
                print(f"    Constraints: ", end="")
                for i, c in enumerate(r["all_constraints"]):
                    sep = "  |  " if i > 0 else ""
                    print(f"{sep}{c['source']}={c['dx_mm']:.2f}", end="")
                print()

            # h estimation details (non-shield parts)
            if r.get("h_estimation"):
                he = r["h_estimation"]
                print(f"    h method:    {he.get('method', '?')} "
                      f"({he.get('regime', '?')})")
                if he.get("details"):
                    d = he["details"]
                    if "Ri" in d:
                        print(f"    Richardson:  Ri={d['Ri']:.3f}  "
                              f"(h_forced={d.get('h_forced', 0):.1f}, "
                              f"h_natural={d.get('h_natural', 0):.1f})")

            # Solver advisory
            if r.get("solver_advisory"):
                sa = r["solver_advisory"]
                ss_str = "YES" if sa["steady_state_ok"] else "NO"
                print(f"    Solver:      SS={ss_str}  "
                      f"severity={sa['severity']}  "
                      f"({sa.get('natural_regime', '?')} nat / "
                      f"{sa.get('forced_regime', '?')} forced)")

            # Warnings
            if r.get("warnings"):
                for w in r["warnings"]:
                    icon = "!!" if w["severity"] == "warning" else " ?"
                    print(f"    [{icon}] {w['code']}: {w['message']}")

            print(f"    >> GOVERNING: {r['governing_dx_mm']:.2f} mm "
                  f"({r['governing_constraint']})")
            print()

    # ------------------------------------------------------------------
    #  Warnings summary
    # ------------------------------------------------------------------
    all_warnings = []
    for r in results:
        for w in r.get("warnings", []):
            all_warnings.append((r.get("part_id", "?"), w))

    if all_warnings:
        print("=" * 80)
        print("  WARNINGS SUMMARY")
        print("=" * 80)
        for pid, w in all_warnings:
            icon = "!!" if w["severity"] == "warning" else " ?"
            print(f"  [{icon}] {pid}")
            print(f"       {w['message']}")
        print()


if __name__ == "__main__":
    main()
