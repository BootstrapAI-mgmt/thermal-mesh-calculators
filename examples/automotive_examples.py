#!/usr/bin/env python3
"""
Automotive Thermal Mesh Sizing — Worked Examples
=================================================

Demonstrates all calculators against realistic underhood scenarios.
Run directly:  python -m examples.automotive_examples
"""

from thermal_mesh_calculators import (
    BoundaryDrivenConductionCalculator,
    ConvectionMeshCalculator,
    RadiationMeshCalculator,
    SingleLayerShieldCalculator,
    MultilayerShieldCalculator,
    TransientMeshCalculator,
)


def separator(title: str) -> None:
    print(f"\n{'='*60}")
    print(f"  {title}")
    print(f"{'='*60}\n")


def main():
    # ------------------------------------------------------------------
    # 1. CONDUCTION — Exhaust manifold (fixed surface temp)
    # ------------------------------------------------------------------
    separator("1. Boundary-Driven Conduction — Exhaust Manifold")

    calc = BoundaryDrivenConductionCalculator()
    result = calc.max_mesh_size(
        k=45.0,            # Cast iron / steel (W/m K)
        h=150.0,           # Underhood forced convection
        t_surf=1073.15,    # 800 C manifold surface (fixed BC)
        t_fluid=353.15,    # 80 C underhood air
        epsilon=0.85,      # Oxidised steel
        t_surr=353.15,     # 80 C surroundings
        max_dt=10.0,       # 10 K accuracy target per element
    )
    print(f"  Convection flux:   {result['q_conv']:.0f} W/m^2")
    print(f"  Radiation flux:    {result['q_rad']:.0f} W/m^2")
    print(f"  Radiation fraction: {result['rad_fraction']*100:.1f}%")
    print(f"  Total flux:        {result['q_total']:.0f} W/m^2")
    print(f"  >> Max element size: {result['max_dx_mm']:.2f} mm")

    # ------------------------------------------------------------------
    # 2. CONDUCTION — Plastic intake manifold (low k, moderate temp)
    # ------------------------------------------------------------------
    separator("2. Boundary-Driven Conduction — Plastic Intake Manifold")

    result2 = calc.max_mesh_size(
        k=0.25,            # Glass-filled nylon (W/m K)
        h=40.0,            # Moderate underhood convection
        t_surf=393.15,     # 120 C surface
        t_fluid=353.15,    # 80 C ambient
        epsilon=0.92,      # Plastic emissivity
        t_surr=353.15,
        max_dt=5.0,        # Tighter accuracy for low-k material
    )
    print(f"  Convection flux:   {result2['q_conv']:.0f} W/m^2")
    print(f"  Radiation flux:    {result2['q_rad']:.0f} W/m^2")
    print(f"  Total flux:        {result2['q_total']:.0f} W/m^2")
    print(f"  >> Max element size: {result2['max_dx_mm']:.2f} mm")
    print("  NOTE: Low conductivity drives very fine mesh even at")
    print("  moderate temperatures.")

    # ------------------------------------------------------------------
    # 3. CONVECTION — Biot number evaluation
    # ------------------------------------------------------------------
    separator("3. Convection — Biot Number (Shell vs. Solid Decision)")

    conv_calc = ConvectionMeshCalculator()

    # Thin stamped bracket
    bi_thin = conv_calc.biot_number(h=100.0, k=45.0, thickness=0.002)
    print(f"  2mm steel bracket:  Bi = {bi_thin['biot']:.4f}")
    print(f"    -> {bi_thin['mesh_type']}")
    print(f"    -> {bi_thin['rationale']}")
    print()

    # Thick cast housing
    bi_thick = conv_calc.biot_number(h=250.0, k=45.0, thickness=0.015)
    print(f"  15mm cast housing:  Bi = {bi_thick['biot']:.4f}")
    print(f"    -> {bi_thick['mesh_type']}")
    print(f"    -> {bi_thick['rationale']}")

    # ------------------------------------------------------------------
    # 4. CONVECTION — h-gradient and CFD mapping
    # ------------------------------------------------------------------
    separator("4. Convection — h-Gradient & CFD Mapping Limits")

    h_grad = conv_calc.h_gradient_mesh_limit(
        h_max=800.0,       # Impingement zone
        h_min=50.0,        # Far-field
        gradient_length_mm=15.0,
    )
    print("  h varies 50 → 800 W/m^2K over 15 mm:")
    print(f"    h ratio: {h_grad['h_ratio']:.1f}x")
    print(f"    Elements needed: {h_grad['elements_needed']}")
    print(f"    >> Max element size: {h_grad['max_dx_mm']:.2f} mm")
    print()

    cfd_limit = conv_calc.cfd_mapping_mesh_limit(
        fluid_wall_face_mm=2.0,
        mapping_ratio=4.0,
    )
    print("  CFD wall face = 2 mm, 4:1 mapping ratio:")
    print(f"    >> Max solid surface element: {cfd_limit:.1f} mm")

    # ------------------------------------------------------------------
    # 5. RADIATION — High-temp sensitivity
    # ------------------------------------------------------------------
    separator("5. Radiation — T^4 Sensitivity at Different Temperatures")

    rad_calc = RadiationMeshCalculator()

    for t_c in [200, 500, 800, 1000]:
        t_k = t_c + 273.15
        sens = rad_calc.flux_sensitivity(t_k, emissivity=0.85)
        r = rad_calc.max_mesh_size(
            t_local=t_k,
            emissivity=0.85,
            allowable_flux_error=500.0,
            spatial_gradient=2000.0,   # 2 K/mm
        )
        print(f"  {t_c:>5d} C:  dq/dT = {sens:>8.1f} W/m^2K"
              f"  |  max dT/elem = {r['max_dt_element']:.1f} K"
              f"  |  max dx = {r['max_dx_mm']:.1f} mm")

    # ------------------------------------------------------------------
    # 6. RADIATION — View factor curvature limit
    # ------------------------------------------------------------------
    separator("6. Radiation — Curvature / View Factor Limit")

    for r_mm in [25, 50, 100, 200]:
        vf_limit = rad_calc.view_factor_curvature_limit(r_mm, max_facet_angle_deg=15.0)
        print(f"  R = {r_mm:>3d} mm, 15 deg facet limit"
              f"  ->  max edge = {vf_limit:.1f} mm")

    # ------------------------------------------------------------------
    # 7. SINGLE-LAYER SHIELD — Aluminised steel
    # ------------------------------------------------------------------
    separator("7. Single-Layer Shield — Aluminised Steel")

    shield_calc = SingleLayerShieldCalculator()

    # v0.2 API: separate h_in / h_out for asymmetric convection
    s1 = shield_calc.mesh_size(
        k=45.0,
        max_dt=15.0,
        t_exh=1073.15,     # 800 C exhaust
        t_fluid=353.15,    # 80 C air
        t_surr=353.15,
        h_in=30.0,         # Exhaust-side convection (sheltered)
        h_out=30.0,        # Ambient-side convection (open)
        eps_in=0.4,        # Aluminised
        eps_out=0.4,
    )
    print(f"  Shield temp:  {s1['t_shield_C']:.1f} C")
    print(f"  Rad absorbed: {s1['q_rad_in']:.0f} W/m^2")
    print(f"  Rad emitted:  {s1['q_rad_out']:.0f} W/m^2")
    print(f"  Conv (in):    {s1['q_conv_in']:.0f} W/m^2")
    print(f"  Conv (out):   {s1['q_conv_out']:.0f} W/m^2")
    print(f"  >> Max element size: {s1['max_dx_mm']:.2f} mm")
    print()

    # Demonstrate asymmetric convection — sheltered inner vs exposed outer
    s1b = shield_calc.mesh_size(
        k=45.0,
        max_dt=15.0,
        t_exh=1073.15,
        t_fluid=353.15,
        t_surr=353.15,
        h_in=15.0,         # Sheltered exhaust gap — low velocity
        h_out=60.0,        # Exposed to underhood airflow
        eps_in=0.4,
        eps_out=0.4,
    )
    print("  Asymmetric case (h_in=15, h_out=60):")
    print(f"    Shield temp:  {s1b['t_shield_C']:.1f} C")
    print(f"    Conv (in):    {s1b['q_conv_in']:.0f} W/m^2")
    print(f"    Conv (out):   {s1b['q_conv_out']:.0f} W/m^2")
    print(f"    >> Max element size: {s1b['max_dx_mm']:.2f} mm")

    # ------------------------------------------------------------------
    # 8. MULTILAYER SHIELD — Dual-wall with air gap
    # ------------------------------------------------------------------
    separator("8. Multilayer Shield — Dual-Wall with Air Gap")

    multi_calc = MultilayerShieldCalculator()
    m = multi_calc.mesh_sizes(
        k_metal=45.0,
        max_dt=15.0,
        t_exh=1073.15,
        t_fluid=353.15,
        t_surr=353.15,
        h_in=30.0,
        h_out=30.0,
        h_gap=15.0,         # Air gap + contact points
        eps_in=0.4,
        eps_out=0.4,
        eps_g1=0.4,
        eps_g2=0.4,
    )
    print(f"  Layer 1 (exhaust side): {m['t1_C']:.1f} C"
          f"  ->  max dx = {m['layer1_max_dx_mm']:.2f} mm")
    print(f"  Layer 2 (ambient side): {m['t2_C']:.1f} C"
          f"  ->  max dx = {m['layer2_max_dx_mm']:.2f} mm")
    print(f"  Delta T across gap:     {m['delta_T_C']:.1f} C")
    print(f"  Converged: {m['converged']} in {m['iterations']} iterations")
    print()
    print(f"  KEY INSIGHT: Layer 2 allows {m['layer2_max_dx_mm']/m['layer1_max_dx_mm']:.1f}x")
    print("  coarser mesh than Layer 1 — significant node count savings.")

    # ------------------------------------------------------------------
    # 9. MULTILAYER — High-emissivity (oxidised) comparison
    # ------------------------------------------------------------------
    separator("9. Multilayer Shield — Oxidised Steel (eps=0.85)")

    m2 = multi_calc.mesh_sizes(
        k_metal=45.0,
        max_dt=15.0,
        t_exh=1073.15,
        t_fluid=353.15,
        t_surr=353.15,
        h_in=30.0,
        h_out=30.0,
        h_gap=15.0,
        eps_in=0.85,
        eps_out=0.85,
        eps_g1=0.85,
        eps_g2=0.85,
    )
    print(f"  Layer 1: {m2['t1_C']:.1f} C  ->  max dx = {m2['layer1_max_dx_mm']:.2f} mm")
    print(f"  Layer 2: {m2['t2_C']:.1f} C  ->  max dx = {m2['layer2_max_dx_mm']:.2f} mm")
    print(f"  Delta T: {m2['delta_T_C']:.1f} C")
    print()
    print("  COMPARISON: Higher emissivity increases radiation coupling,")
    print("  raising Layer 1 heat load and requiring finer mesh.")

    # ------------------------------------------------------------------
    # 10. TRANSIENT — Penetration depth and Fourier number
    # ------------------------------------------------------------------
    separator("10. Transient Mesh Constraints — Steel Exhaust Bracket")

    # Steel bracket: k=45, rho=7800, cp=500, dt=0.5 s, explicit solver
    # (Fo <= 0.5).  Stability bounds the element from BELOW; the per-step
    # penetration depth does not bound an explicit scheme (it is an
    # implicit resolution guideline), so the drive cycle sets the upper
    # bound.
    tc = TransientMeshCalculator()

    pen = tc.penetration_depth(k=45.0, rho=7800.0, cp=500.0, dt=0.5)
    print(f"  Thermal diffusivity: {pen['alpha']:.3e} m^2/s")
    print(f"  Penetration depth (dt=0.5s): {pen['max_dx_mm']:.2f} mm"
          f"  (implicit guideline; not applied here)")

    fo = tc.fourier_number_limit(k=45.0, rho=7800.0, cp=500.0, dt=0.5, fo_max=0.5)
    print(f"  Fourier stability min dx:    {fo['min_dx_mm']:.2f} mm  (Fo <= {fo['fo_max']})")

    combo = tc.combined_transient_limits(
        k=45.0, rho=7800.0, cp=500.0, dt=0.5,
        fo_max=0.5, safety_factor=1.0, tau_bc=20.0,
    )
    print(f"  Drive-cycle limit (tau=20s): {combo['drive_cycle_max_dx_mm']:.2f} mm")
    print(f"  >> Feasible window ({combo['scheme']}): "
          f"{combo['fourier_min_dx_mm']:.2f} to {combo['max_dx_mm']:.2f} mm")
    print(f"  >> Recommended dx: {combo['recommended_dx_mm']:.2f} mm")
    print(f"     Binding constraint: {combo['binding_constraint']}")

    # A time step too long for the drive cycle: the conflict carries a
    # remedy that changes the result (dt <= Fo_max * tau_bc).
    slow = tc.combined_transient_limits(
        k=45.0, rho=7800.0, cp=500.0, dt=50.0,
        fo_max=0.5, safety_factor=1.0, tau_bc=20.0,
    )
    print()
    print(f"  dt = 50 s: {slow['advice']}")
    fixed = tc.combined_transient_limits(
        k=45.0, rho=7800.0, cp=500.0,
        dt=slow["conflict"]["remedies"][0]["max_value"],
        fo_max=0.5, safety_factor=1.0, tau_bc=20.0,
    )
    print(f"  With that dt: feasible = {fixed['feasible']}, "
          f"recommended dx = {fixed['recommended_dx_mm']:.2f} mm")

    # ------------------------------------------------------------------
    # 11. TRANSIENT — Implicit solver comparison (plastic part)
    # ------------------------------------------------------------------
    separator("11. Transient — Implicit Solver (Plastic Intake, dt=2s)")

    # Plastic: k=0.25, rho=1400, cp=1600
    combo_plastic = tc.combined_transient_limits(
        k=0.25, rho=1400.0, cp=1600.0, dt=2.0,
        fo_max=5.0,             # implicit solver — accuracy limit
        safety_factor=2.0,      # looser penetration bound for implicit
        tau_bc=30.0,
    )
    print(f"  Diffusivity: {combo_plastic['alpha']:.3e} m^2/s")
    print(f"  Penetration depth limit: {combo_plastic['penetration_max_dx_mm']:.3f} mm")
    print(f"  Fourier accuracy min dx: {combo_plastic['fourier_min_dx_mm']:.3f} mm")
    print(f"  Drive-cycle limit:       {combo_plastic['drive_cycle_max_dx_mm']:.3f} mm")
    print(f"  >> Recommended dx: {combo_plastic['recommended_dx_mm']:.3f} mm")
    print(f"     Binding constraint: {combo_plastic['binding_constraint']}")
    print()
    print("  NOTE: Low thermal diffusivity in plastics makes transient")
    print("  constraints much tighter than for metals.")


if __name__ == "__main__":
    main()
