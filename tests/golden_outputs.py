"""
Outputs captured before the internal-flow correlations were added.

test_examples_regression.py holds the package to them: the two example scripts'
printed output byte for byte, and estimate_h() over a grid of inputs with
its floats written to 12 significant digits.  Regenerate them only for a
change that is meant to move a result, and say so in CHANGELOG.md.
"""

AUTOMOTIVE_EXAMPLES_STDOUT = (
    '\n'
    '============================================================\n'
    '  1. Boundary-Driven Conduction — Exhaust Manifold\n'
    '============================================================\n'
    '\n'
    '  Convection flux:   108000 W/m^2\n'
    '  Radiation flux:    63171 W/m^2\n'
    '  Radiation fraction: 36.9%\n'
    '  Total flux:        171171 W/m^2\n'
    '  >> Max element size: 2.63 mm\n'
    '\n'
    '============================================================\n'
    '  2. Boundary-Driven Conduction — Plastic Intake Manifold\n'
    '============================================================\n'
    '\n'
    '  Convection flux:   1600 W/m^2\n'
    '  Radiation flux:    435 W/m^2\n'
    '  Total flux:        2035 W/m^2\n'
    '  >> Max element size: 0.61 mm\n'
    '  NOTE: Low conductivity drives very fine mesh even at\n'
    '  moderate temperatures.\n'
    '\n'
    '============================================================\n'
    '  3. Convection — Biot Number (Shell vs. Solid Decision)\n'
    '============================================================\n'
    '\n'
    '  2mm steel bracket:  Bi = 0.0022\n'
    '    -> 2D Shell\n'
    '    -> Bi = 0.0022 < 0.1 — internal conduction resistance is negligible relative to surface convection.  A single shell element across the thickness is physically valid.\n'
    '\n'
    '  15mm cast housing:  Bi = 0.0417\n'
    '    -> 2D Shell\n'
    '    -> Bi = 0.0417 < 0.1 — internal conduction resistance is negligible relative to surface convection.  A single shell element across the thickness is physically valid.\n'
    '\n'
    '============================================================\n'
    '  4. Convection — h-Gradient & CFD Mapping Limits\n'
    '============================================================\n'
    '\n'
    '  h varies 50 → 800 W/m^2K over 15 mm:\n'
    '    h ratio: 16.0x\n'
    '    Elements needed: 6\n'
    '    >> Max element size: 2.50 mm\n'
    '\n'
    '  CFD wall face = 2 mm, 4:1 mapping ratio:\n'
    '    >> Max solid surface element: 8.0 mm\n'
    '\n'
    '============================================================\n'
    '  5. Radiation — T^4 Sensitivity at Different Temperatures\n'
    '============================================================\n'
    '\n'
    '    200 C:  dq/dT =     20.4 W/m^2K  |  max dT/elem = 24.5 K  |  max dx = 12.2 mm\n'
    '    500 C:  dq/dT =     89.1 W/m^2K  |  max dT/elem = 5.6 K  |  max dx = 2.8 mm\n'
    '    800 C:  dq/dT =    238.3 W/m^2K  |  max dT/elem = 2.1 K  |  max dx = 1.0 mm\n'
    '   1000 C:  dq/dT =    397.8 W/m^2K  |  max dT/elem = 1.3 K  |  max dx = 0.6 mm\n'
    '\n'
    '============================================================\n'
    '  6. Radiation — Curvature / View Factor Limit\n'
    '============================================================\n'
    '\n'
    '  R =  25 mm, 15 deg facet limit  ->  max edge = 6.5 mm\n'
    '  R =  50 mm, 15 deg facet limit  ->  max edge = 13.1 mm\n'
    '  R = 100 mm, 15 deg facet limit  ->  max edge = 26.1 mm\n'
    '  R = 200 mm, 15 deg facet limit  ->  max edge = 52.2 mm\n'
    '\n'
    '============================================================\n'
    '  7. Single-Layer Shield — Aluminised Steel\n'
    '============================================================\n'
    '\n'
    '  Shield temp:  416.4 C\n'
    '  Rad absorbed: 24954 W/m^2\n'
    '  Rad emitted:  4773 W/m^2\n'
    '  Conv (in):    10091 W/m^2\n'
    '  Conv (out):   10091 W/m^2\n'
    '  >> Max element size: 19.26 mm\n'
    '\n'
    '  Asymmetric case (h_in=15, h_out=60):\n'
    '    Shield temp:  377.4 C\n'
    '    Conv (in):    4461 W/m^2\n'
    '    Conv (out):   17846 W/m^2\n'
    '    >> Max element size: 22.15 mm\n'
    '\n'
    '============================================================\n'
    '  8. Multilayer Shield — Dual-Wall with Air Gap\n'
    '============================================================\n'
    '\n'
    '  Layer 1 (exhaust side): 520.6 C  ->  max dx = 19.68 mm\n'
    '  Layer 2 (ambient side): 282.0 C  ->  max dx = 85.88 mm\n'
    '  Delta T across gap:     238.6 C\n'
    '  Converged: True in 5 iterations\n'
    '\n'
    '  KEY INSIGHT: Layer 2 allows 4.4x\n'
    '  coarser mesh than Layer 1 — significant node count savings.\n'
    '\n'
    '============================================================\n'
    '  9. Multilayer Shield — Oxidised Steel (eps=0.85)\n'
    '============================================================\n'
    '\n'
    '  Layer 1: 606.7 C  ->  max dx = 13.28 mm\n'
    '  Layer 2: 405.5 C  ->  max dx = 35.08 mm\n'
    '  Delta T: 201.2 C\n'
    '\n'
    '  COMPARISON: Higher emissivity increases radiation coupling,\n'
    '  raising Layer 1 heat load and requiring finer mesh.\n'
    '\n'
    '============================================================\n'
    '  10. Transient Mesh Constraints — Steel Exhaust Bracket\n'
    '============================================================\n'
    '\n'
    '  Thermal diffusivity: 1.154e-05 m^2/s\n'
    '  Penetration depth (dt=0.5s): 2.40 mm  (implicit guideline; not applied here)\n'
    '  Fourier stability min dx:    3.40 mm  (Fo <= 0.5)\n'
    '  Drive-cycle limit (tau=20s): 15.19 mm\n'
    '  >> Feasible window (explicit): 3.40 to 15.19 mm\n'
    '  >> Recommended dx: 15.19 mm\n'
    '     Binding constraint: drive_cycle\n'
    '\n'
    '  dt = 50 s: The Fourier stability minimum (33.97 mm, Fo <= 0.5) exceeds the drive-cycle bound (15.19 mm): no element satisfies both at dt = 50 s. Remedy: reduce dt to <= 10 s.\n'
    '  With that dt: feasible = True, recommended dx = 15.19 mm\n'
    '\n'
    '============================================================\n'
    '  11. Transient — Implicit Solver (Plastic Intake, dt=2s)\n'
    '============================================================\n'
    '\n'
    '  Diffusivity: 1.116e-07 m^2/s\n'
    '  Penetration depth limit: 0.945 mm\n'
    '  Fourier accuracy min dx: 0.211 mm\n'
    '  Drive-cycle limit:       1.830 mm\n'
    '  >> Recommended dx: 0.945 mm\n'
    '     Binding constraint: penetration_depth\n'
    '\n'
    '  NOTE: Low thermal diffusivity in plastics makes transient\n'
    '  constraints much tighter than for metals.\n'
)

BATCH_EXAMPLE_STDOUT = (
    '================================================================================\n'
    '  AVAILABLE LOOKUPS\n'
    '================================================================================\n'
    '  Materials (43):  aluminium_5052, aluminium_6061, brass, cast_aluminium, cast_aluminium_a380, cast_iron, cast_iron_ductile, ceramic_alumina, ceramic_cordierite, composite_cfrp, composite_smc, copper, glass_soda_lime, insulation_ceramic_blanket, insulation_fibreglass, magnesium_az91, nickel_alloy, plastic_abs, plastic_acetal_pom, plastic_hdpe, plastic_pa66, plastic_pa66_gf30, plastic_pa6_gf30, plastic_pbt_gf30, plastic_pc, plastic_pc_abs, plastic_pet, plastic_pp, plastic_pp_gf30, plastic_pps_gf40, rubber_cr, rubber_epdm, rubber_fkm, rubber_natural, rubber_nbr, rubber_silicone, steel_galvanised, steel_high_strength, steel_mild, steel_stainless_304, steel_stainless_409, titanium_6al4v, zinc\n'
    '  Surfaces (30):   aluminised, aluminium_oxidised, anodised, bare_metal, cast_aluminium_bare, cast_iron_oxidised, ceramic_coating, chrome_plated, composite_smc, copper_oxidised, e_coat, fabric_woven, galvanised_new, galvanised_weathered, glass, heavily_oxidised, lightly_oxidised, nickel_plated, painted, painted_gloss_white, painted_matte_black, plastic, polished_aluminium, polished_copper, polished_steel, powder_coat, rubber, rusted_steel, stainless_weathered, zinc_plated\n'
    '  Zones (18):      cabin_dead_zone, cabin_tunnel, cargo_bed_exterior, clutch_outlet_downstream, cooling_pack_downstream, engine_above, engine_bay_dead_zone, engine_bay_edges, engine_below, engine_beside, exhaust_above, exhaust_below, exhaust_beside, exhaust_internal, front_end_dead_zone, front_end_edges, near_exhaust_natural, shield_gap_confined\n'
    '\n'
    '================================================================================\n'
    '  BATCH MESH SIZING RESULTS — Off-Road Side-by-Side\n'
    '================================================================================\n'
    '\n'
    'Part ID                        Material               Class               dx (mm)               Governing      Biot\n'
    '-------------------------------------------------------------------------------------------------------------------\n'
    'EXH-001 Exhaust Manifold       cast_iron              exhaust                4.82              conduction    0.0039\n'
    'EXH-002 Downpipe               steel_stainless_304    exhaust                1.91              conduction    0.0042\n'
    'STR-001 Engine Mount Bracket   steel_mild             structural            17.22               transient    0.0002\n'
    'STR-002 Oil Pan                cast_aluminium         structural            35.14               transient    0.0001\n'
    'FE-001 Radiator Support        plastic_pa66_gf30      structural             1.49               transient    0.0900\n'
    'FE-002 Charge Air Cooler Bracket aluminium_6061         structural            37.09               transient    0.0002\n'
    'CAB-001 Tunnel Crossmember     steel_mild             structural            17.22               transient    0.0001\n'
    'SH-001 Manifold Shield         steel_mild             shield                17.22               transient    0.0001\n'
    'SH-002 Cargo Bed Shield        steel_mild             shield                17.22               transient    0.0001\n'
    'SH-003 Dual-Wall Cargo Shield  steel_mild             multilayer_shield     16.95               radiation    0.0001\n'
    'CARGO-001 Bed Liner            plastic_hdpe           structural             2.22               transient    0.0492\n'
    'SOFT-001 Exhaust Grommet       rubber_epdm            exhaust_adjacent       1.51               transient    0.1008\n'
    '\n'
    '===================================================================================================================\n'
    '  WARNINGS\n'
    '===================================================================================================================\n'
    '  [ ?] EXH-001 Exhaust Manifold: The air-property fits cover film temperatures of 250-700 K; the h estimate (film temperature 713 K) used air properties evaluated at 700 K. Treat the result as extrapolated.\n'
    '  [!!] EXH-002 Downpipe: Transient (explicit, dt = 0.5 s): the smallest element this time step allows, 1.993 mm (Fo <= 0.5, stability limit), exceeds the governing size, 1.907 mm. Remedy: reduce dt to <= 0.4577 s.\n'
    "  [ ?] FE-001 Radiator Support: No surface treatment ('surface') and no emissivity ('epsilon') given; assumed eps = 0.90, the non-metal fallback for material 'plastic_pa66_gf30'. Set 'surface' to one of list_surface_treatments(), or give 'epsilon'.\n"
    '  [ ?] FE-001 Radiator Support: Non-metal part (Bi=0.090) approaching or exceeding Bi=0.1 threshold. Through-thickness gradients may require 3D solid mesh instead of shell. h=12.9 W/m²K.\n'
    '  [!!] CAB-001 Tunnel Crossmember: Horizontal hot-up dead zone with L=400 mm > 150 mm threshold. Natural convection is likely turbulent (Ra > 1e7). Recommend transient solve or verify with solver_advisory.\n'
    '  [!!] CAB-001 Tunnel Crossmember: Turbulent natural convection (Ra=5.70e+07, threshold=1e+07 for horizontal_up). Buoyant plumes are inherently unsteady — steady-state solver may not converge or may give misleading time-averaged results. Recommend transient solve with time-averaged output, or use relaxation with caution in steady-state.\n'
    "  [ ?] SH-003 Dual-Wall Cargo Shield: No 'f12' given (gap view factor; 1.0 is parallel plates, offset or curved gaps are lower); assumed f12 = 1. Give 'f12' on the part.\n"
    "  [ ?] CARGO-001 Bed Liner: No surface treatment ('surface') and no emissivity ('epsilon') given; assumed eps = 0.90, the non-metal fallback for material 'plastic_hdpe'. Set 'surface' to one of list_surface_treatments(), or give 'epsilon'.\n"
    "  [ ?] SOFT-001 Exhaust Grommet: No surface treatment ('surface') and no emissivity ('epsilon') given; assumed eps = 0.90, the non-metal fallback for material 'rubber_epdm'. Set 'surface' to one of list_surface_treatments(), or give 'epsilon'.\n"
    '  [!!] SOFT-001 Exhaust Grommet: Non-metal part (Bi=0.101) approaching or exceeding Bi=0.1 threshold. Through-thickness gradients may require 3D solid mesh instead of shell. h=6.3 W/m²K.\n'
    '\n'
    '------------------------------------------------------------\n'
    '  EXH-001 Exhaust Manifold\n'
    '------------------------------------------------------------\n'
    '    Material:    cast_iron\n'
    '    Class:       exhaust\n'
    '    h used:      67.16762350746704\n'
    '    eps used:    0.8\n'
    '    Conduction:  dx=4.82 mm  (q=107816 W/m^2, rad_frac=55.1%)\n'
    '    Biot:        0.0039 → 2D Shell\n'
    '    Radiation:   dx=13.35 mm  (dq/dT=224.2 W/m^2K)\n'
    '    Transient:   dx=17.72 mm  (drive_cycle)\n'
    '    Constraints: conduction=4.82  |  lateral_gradient=10.91  |  radiation=13.35  |  transient=17.72\n'
    '    h method:    correlation (forced)\n'
    '    Richardson:  Ri=0.001  (h_forced=67.2, h_natural=11.2)\n'
    '    Solver:      SS=YES  severity=info  (laminar nat / laminar forced)\n'
    '    [ ?] FILM_TEMP_OUT_OF_RANGE: The air-property fits cover film temperatures of 250-700 K; the h estimate (film temperature 713 K) used air properties evaluated at 700 K. Treat the result as extrapolated.\n'
    '    >> GOVERNING: 4.82 mm (conduction)\n'
    '\n'
    '------------------------------------------------------------\n'
    '  STR-001 Engine Mount Bracket\n'
    '------------------------------------------------------------\n'
    '    Material:    steel_mild\n'
    '    Class:       structural\n'
    '    h used:      6.627331225108186\n'
    '    eps used:    0.92\n'
    '    Conduction:  dx=311.74 mm  (q=2598 W/m^2, rad_frac=69.4%)\n'
    '    Biot:        0.0002 → 2D Shell\n'
    '    Radiation:   dx=188.06 mm  (dq/dT=22.1 W/m^2K)\n'
    '    Transient:   dx=17.22 mm  (drive_cycle)\n'
    '    Constraints: transient=17.22  |  lateral_gradient=25.03  |  radiation=188.06  |  conduction=311.74\n'
    '    h method:    correlation (natural)\n'
    '    Richardson:  Ri=inf  (h_forced=0.0, h_natural=6.6)\n'
    '    Solver:      SS=YES  severity=info  (laminar nat / n/a forced)\n'
    '    >> GOVERNING: 17.22 mm (transient)\n'
    '\n'
    '------------------------------------------------------------\n'
    '  FE-001 Radiator Support\n'
    '------------------------------------------------------------\n'
    '    Material:    plastic_pa66_gf30\n'
    '    Class:       structural\n'
    '    h used:      12.853423296747291\n'
    '    eps used:    0.9\n'
    '    Conduction:  dx=5.42 mm  (q=692 W/m^2, rad_frac=44.2%)\n'
    '    Biot:        0.0900 → 2D Shell\n'
    '    Radiation:   dx=6.30 mm  (dq/dT=11.5 W/m^2K)\n'
    '    Transient:   dx=1.49 mm  (drive_cycle)\n'
    '    Constraints: transient=1.49  |  lateral_gradient=2.00  |  conduction=5.42  |  radiation=6.30\n'
    '    h method:    correlation (forced)\n'
    '    Richardson:  Ri=0.022  (h_forced=12.9, h_natural=4.5)\n'
    '    Solver:      SS=YES  severity=info  (laminar nat / laminar forced)\n'
    "    [ ?] SURFACE_DEFAULTED: No surface treatment ('surface') and no emissivity ('epsilon') given; assumed eps = 0.90, the non-metal fallback for material 'plastic_pa66_gf30'. Set 'surface' to one of list_surface_treatments(), or give 'epsilon'.\n"
    '    [ ?] PLASTIC_BIOT_MARGINAL: Non-metal part (Bi=0.090) approaching or exceeding Bi=0.1 threshold. Through-thickness gradients may require 3D solid mesh instead of shell. h=12.9 W/m²K.\n'
    '    >> GOVERNING: 1.49 mm (transient)\n'
    '\n'
    '------------------------------------------------------------\n'
    '  SH-001 Manifold Shield\n'
    '------------------------------------------------------------\n'
    '    Material:    steel_mild\n'
    '    Class:       shield\n'
    "    h used:      {'h_in': 20.0, 'h_out': 15.0}\n"
    "    eps used:    {'eps_in': 0.4, 'eps_out': 0.4}\n"
    '    Biot:        0.0001 → 2D Shell\n'
    '    Radiation:   dx=21.47 mm  (dq/dT=41.3 W/m^2K)\n'
    '    Shield:      T=496.0 C, dx=26.59 mm\n'
    '    Transient:   dx=17.22 mm  (drive_cycle)\n'
    '    Constraints: transient=17.22  |  radiation=21.47  |  shield=26.59\n'
    '    >> GOVERNING: 17.22 mm (transient)\n'
    '\n'
    '------------------------------------------------------------\n'
    '  SH-003 Dual-Wall Cargo Shield\n'
    '------------------------------------------------------------\n'
    '    Material:    steel_mild\n'
    '    Class:       multilayer_shield\n'
    "    h used:      {'h_in': 20.0, 'h_out': 15.0, 'h_gap': 12.0}\n"
    "    eps used:    {'eps_in': 0.4, 'eps_out': 0.4, 'eps_g1': 0.4, 'eps_g2': 0.4}\n"
    '    Biot:        0.0001 → 2D Shell\n'
    '    Radiation:   dx=16.95 mm  (dq/dT=57.0 W/m^2K)\n'
    '    Shield L1:   T=583.3 C, dx=28.99 mm\n'
    '    Shield L2:   T=368.3 C, dx=103.72 mm\n'
    '    Transient:   dx=17.22 mm  (drive_cycle)\n'
    '    Constraints: radiation=16.95  |  transient=17.22  |  shield_layer1=28.99  |  shield_layer2=103.72\n'
    "    [ ?] SHIELD_INPUT_DEFAULTED: No 'f12' given (gap view factor; 1.0 is parallel plates, offset or curved gaps are lower); assumed f12 = 1. Give 'f12' on the part.\n"
    '    >> GOVERNING: 16.95 mm (radiation)\n'
    '\n'
    '------------------------------------------------------------\n'
    '  CARGO-001 Bed Liner\n'
    '------------------------------------------------------------\n'
    '    Material:    plastic_hdpe\n'
    '    Class:       structural\n'
    '    h used:      8.864616160843356\n'
    '    eps used:    0.9\n'
    '    Conduction:  dx=18.10 mm  (q=373 W/m^2, rad_frac=52.5%)\n'
    '    Biot:        0.0492 → 2D Shell\n'
    '    Radiation:   dx=22.75 mm  (dq/dT=10.6 W/m^2K)\n'
    '    Transient:   dx=2.22 mm  (drive_cycle)\n'
    '    Constraints: transient=2.22  |  lateral_gradient=3.58  |  conduction=18.10  |  radiation=22.75\n'
    '    h method:    correlation (forced)\n'
    '    Richardson:  Ri=0.047  (h_forced=8.9, h_natural=4.1)\n'
    '    Solver:      SS=YES  severity=caution  (turbulent nat / laminar forced)\n'
    "    [ ?] SURFACE_DEFAULTED: No surface treatment ('surface') and no emissivity ('epsilon') given; assumed eps = 0.90, the non-metal fallback for material 'plastic_hdpe'. Set 'surface' to one of list_surface_treatments(), or give 'epsilon'.\n"
    '    >> GOVERNING: 2.22 mm (transient)\n'
    '\n'
    '------------------------------------------------------------\n'
    '  SOFT-001 Exhaust Grommet\n'
    '------------------------------------------------------------\n'
    '    Material:    rubber_epdm\n'
    '    Class:       exhaust_adjacent\n'
    '    h used:      6.301776958178754\n'
    '    eps used:    0.9\n'
    '    Conduction:  dx=1.95 mm  (q=1283 W/m^2, rad_frac=65.6%)\n'
    '    Biot:        0.1008 → 3D Solid\n'
    '    Radiation:   dx=389.51 mm  (dq/dT=15.5 W/m^2K)\n'
    '    Transient:   dx=1.51 mm  (drive_cycle)\n'
    '    Constraints: transient=1.51  |  conduction=1.95  |  biot_through_thickness=2.67  |  lateral_gradient=3.20  |  radiation=389.51\n'
    '    h method:    correlation (natural)\n'
    '    Richardson:  Ri=inf  (h_forced=0.0, h_natural=6.3)\n'
    '    Solver:      SS=YES  severity=info  (laminar nat / n/a forced)\n'
    "    [ ?] SURFACE_DEFAULTED: No surface treatment ('surface') and no emissivity ('epsilon') given; assumed eps = 0.90, the non-metal fallback for material 'rubber_epdm'. Set 'surface' to one of list_surface_treatments(), or give 'epsilon'.\n"
    '    [!!] PLASTIC_BIOT_MARGINAL: Non-metal part (Bi=0.101) approaching or exceeding Bi=0.1 threshold. Through-thickness gradients may require 3D solid mesh instead of shell. h=6.3 W/m²K.\n'
    '    >> GOVERNING: 1.51 mm (transient)\n'
    '\n'
    '================================================================================\n'
    '  WARNINGS SUMMARY\n'
    '================================================================================\n'
    '  [ ?] EXH-001 Exhaust Manifold\n'
    '       The air-property fits cover film temperatures of 250-700 K; the h estimate (film temperature 713 K) used air properties evaluated at 700 K. Treat the result as extrapolated.\n'
    '  [!!] EXH-002 Downpipe\n'
    '       Transient (explicit, dt = 0.5 s): the smallest element this time step allows, 1.993 mm (Fo <= 0.5, stability limit), exceeds the governing size, 1.907 mm. Remedy: reduce dt to <= 0.4577 s.\n'
    '  [ ?] FE-001 Radiator Support\n'
    "       No surface treatment ('surface') and no emissivity ('epsilon') given; assumed eps = 0.90, the non-metal fallback for material 'plastic_pa66_gf30'. Set 'surface' to one of list_surface_treatments(), or give 'epsilon'.\n"
    '  [ ?] FE-001 Radiator Support\n'
    '       Non-metal part (Bi=0.090) approaching or exceeding Bi=0.1 threshold. Through-thickness gradients may require 3D solid mesh instead of shell. h=12.9 W/m²K.\n'
    '  [!!] CAB-001 Tunnel Crossmember\n'
    '       Horizontal hot-up dead zone with L=400 mm > 150 mm threshold. Natural convection is likely turbulent (Ra > 1e7). Recommend transient solve or verify with solver_advisory.\n'
    '  [!!] CAB-001 Tunnel Crossmember\n'
    '       Turbulent natural convection (Ra=5.70e+07, threshold=1e+07 for horizontal_up). Buoyant plumes are inherently unsteady — steady-state solver may not converge or may give misleading time-averaged results. Recommend transient solve with time-averaged output, or use relaxation with caution in steady-state.\n'
    '  [ ?] SH-003 Dual-Wall Cargo Shield\n'
    "       No 'f12' given (gap view factor; 1.0 is parallel plates, offset or curved gaps are lower); assumed f12 = 1. Give 'f12' on the part.\n"
    '  [ ?] CARGO-001 Bed Liner\n'
    "       No surface treatment ('surface') and no emissivity ('epsilon') given; assumed eps = 0.90, the non-metal fallback for material 'plastic_hdpe'. Set 'surface' to one of list_surface_treatments(), or give 'epsilon'.\n"
    '  [ ?] SOFT-001 Exhaust Grommet\n'
    "       No surface treatment ('surface') and no emissivity ('epsilon') given; assumed eps = 0.90, the non-metal fallback for material 'rubber_epdm'. Set 'surface' to one of list_surface_treatments(), or give 'epsilon'.\n"
    '  [!!] SOFT-001 Exhaust Grommet\n'
    '       Non-metal part (Bi=0.101) approaching or exceeding Bi=0.1 threshold. Through-thickness gradients may require 3D solid mesh instead of shell. h=6.3 W/m²K.\n'
    '\n'
)

ESTIMATE_H_GRID_TEXT = (
    "0.0 400.0 300.0 0.05 vertical: {Ra: 654587.521306, Re: 0, Ri: inf, dominant_mode: 'natural convection (vertical)', film_in_range: True, forced_regime: 'n/a', h: 9.26366880706, h_forced: 0, h_natural: 9.26366880706, natural_regime: 'laminar', regime: 'natural', solver_advisory: {Ra: 654587.521306, Re: 0, forced_regime: 'n/a', natural_regime: 'laminar', reason: 'Laminar convection — steady-state solve is appropriate.', severity: 'info', steady_state_ok: True, transient_advisory: False}, t_film_K: 350}\n"
    "0.0 400.0 300.0 0.05 horizontal_up: {Ra: 654587.521306, Re: 0, Ri: inf, dominant_mode: 'natural convection (horizontal_up)', film_in_range: True, forced_regime: 'n/a', h: 8.4786121285, h_forced: 0, h_natural: 8.4786121285, natural_regime: 'laminar', regime: 'natural', solver_advisory: {Ra: 654587.521306, Re: 0, forced_regime: 'n/a', natural_regime: 'laminar', reason: 'Laminar convection — steady-state solve is appropriate.', severity: 'info', steady_state_ok: True, transient_advisory: False}, t_film_K: 350}\n"
    "0.0 400.0 300.0 0.05 horizontal_down: {Ra: 654587.521306, Re: 0, Ri: inf, dominant_mode: 'natural convection (horizontal_down)', film_in_range: True, forced_regime: 'n/a', h: 4.23930606425, h_forced: 0, h_natural: 4.23930606425, natural_regime: 'laminar', regime: 'natural', solver_advisory: {Ra: 654587.521306, Re: 0, forced_regime: 'n/a', natural_regime: 'laminar', reason: 'Laminar convection — steady-state solve is appropriate.', severity: 'info', steady_state_ok: True, transient_advisory: False}, t_film_K: 350}\n"
    "0.0 400.0 300.0 0.4 vertical: {Ra: 335148810.909, Re: 0, Ri: inf, dominant_mode: 'natural convection (vertical)', film_in_range: True, forced_regime: 'n/a', h: 5.50821042819, h_forced: 0, h_natural: 5.50821042819, natural_regime: 'laminar', regime: 'natural', solver_advisory: {Ra: 335148810.909, Re: 0, forced_regime: 'n/a', natural_regime: 'laminar', reason: 'Laminar convection — steady-state solve is appropriate.', severity: 'info', steady_state_ok: True, transient_advisory: False}, t_film_K: 350}\n"
    "0.0 400.0 300.0 0.4 horizontal_up: {Ra: 335148810.909, Re: 0, Ri: inf, dominant_mode: 'natural convection (horizontal_up)', film_in_range: True, forced_regime: 'n/a', h: 7.18929400018, h_forced: 0, h_natural: 7.18929400018, natural_regime: 'turbulent', regime: 'natural', solver_advisory: {Ra: 335148810.909, Re: 0, forced_regime: 'n/a', natural_regime: 'turbulent', reason: 'Turbulent natural convection (Ra=3.35e+08, threshold=1e+07 for horizontal_up). Buoyant plumes are inherently unsteady — steady-state solver may not converge or may give misleading time-averaged results. Recommend transient solve with time-averaged output, or use relaxation with caution in steady-state.', severity: 'warning', steady_state_ok: False, transient_advisory: True}, t_film_K: 350}\n"
    "0.0 400.0 300.0 0.4 horizontal_down: {Ra: 335148810.909, Re: 0, Ri: inf, dominant_mode: 'natural convection (horizontal_down)', film_in_range: True, forced_regime: 'n/a', h: 2.52070646714, h_forced: 0, h_natural: 2.52070646714, natural_regime: 'laminar', regime: 'natural', solver_advisory: {Ra: 335148810.909, Re: 0, forced_regime: 'n/a', natural_regime: 'laminar', reason: 'Laminar convection — steady-state solve is appropriate.', severity: 'info', steady_state_ok: True, transient_advisory: False}, t_film_K: 350}\n"
    "0.0 900.0 330.0 0.05 vertical: {Ra: 312385.928901, Re: 0, Ri: inf, dominant_mode: 'natural convection (vertical)', film_in_range: True, forced_regime: 'n/a', h: 12.8743875558, h_forced: 0, h_natural: 12.8743875558, natural_regime: 'laminar', regime: 'natural', solver_advisory: {Ra: 312385.928901, Re: 0, forced_regime: 'n/a', natural_regime: 'laminar', reason: 'Laminar convection — steady-state solve is appropriate.', severity: 'info', steady_state_ok: True, transient_advisory: False}, t_film_K: 615}\n"
    "0.0 900.0 330.0 0.05 horizontal_up: {Ra: 312385.928901, Re: 0, Ri: inf, dominant_mode: 'natural convection (horizontal_up)', film_in_range: True, forced_regime: 'n/a', h: 11.783337763, h_forced: 0, h_natural: 11.783337763, natural_regime: 'laminar', regime: 'natural', solver_advisory: {Ra: 312385.928901, Re: 0, forced_regime: 'n/a', natural_regime: 'laminar', reason: 'Laminar convection — steady-state solve is appropriate.', severity: 'info', steady_state_ok: True, transient_advisory: False}, t_film_K: 615}\n"
    "0.0 900.0 330.0 0.05 horizontal_down: {Ra: 312385.928901, Re: 0, Ri: inf, dominant_mode: 'natural convection (horizontal_down)', film_in_range: True, forced_regime: 'n/a', h: 5.89166888148, h_forced: 0, h_natural: 5.89166888148, natural_regime: 'laminar', regime: 'natural', solver_advisory: {Ra: 312385.928901, Re: 0, forced_regime: 'n/a', natural_regime: 'laminar', reason: 'Laminar convection — steady-state solve is appropriate.', severity: 'info', steady_state_ok: True, transient_advisory: False}, t_film_K: 615}\n"
    "0.0 900.0 330.0 0.4 vertical: {Ra: 159941595.597, Re: 0, Ri: inf, dominant_mode: 'natural convection (vertical)', film_in_range: True, forced_regime: 'n/a', h: 7.65515664135, h_forced: 0, h_natural: 7.65515664135, natural_regime: 'laminar', regime: 'natural', solver_advisory: {Ra: 159941595.597, Re: 0, forced_regime: 'n/a', natural_regime: 'laminar', reason: 'Laminar convection — steady-state solve is appropriate.', severity: 'info', steady_state_ok: True, transient_advisory: False}, t_film_K: 615}\n"
    "0.0 900.0 330.0 0.4 horizontal_up: {Ra: 159941595.597, Re: 0, Ri: inf, dominant_mode: 'natural convection (horizontal_up)', film_in_range: True, forced_regime: 'n/a', h: 9.39413465799, h_forced: 0, h_natural: 9.39413465799, natural_regime: 'turbulent', regime: 'natural', solver_advisory: {Ra: 159941595.597, Re: 0, forced_regime: 'n/a', natural_regime: 'turbulent', reason: 'Turbulent natural convection (Ra=1.60e+08, threshold=1e+07 for horizontal_up). Buoyant plumes are inherently unsteady — steady-state solver may not converge or may give misleading time-averaged results. Recommend transient solve with time-averaged output, or use relaxation with caution in steady-state.', severity: 'warning', steady_state_ok: False, transient_advisory: True}, t_film_K: 615}\n"
    "0.0 900.0 330.0 0.4 horizontal_down: {Ra: 159941595.597, Re: 0, Ri: inf, dominant_mode: 'natural convection (horizontal_down)', film_in_range: True, forced_regime: 'n/a', h: 3.50320727655, h_forced: 0, h_natural: 3.50320727655, natural_regime: 'laminar', regime: 'natural', solver_advisory: {Ra: 159941595.597, Re: 0, forced_regime: 'n/a', natural_regime: 'laminar', reason: 'Laminar convection — steady-state solve is appropriate.', severity: 'info', steady_state_ok: True, transient_advisory: False}, t_film_K: 615}\n"
    "0.5 400.0 300.0 0.05 vertical: {Ra: 654587.521306, Re: 1282.44635414, Ri: 0.560571428571, dominant_mode: 'mixed-assisting (Ri=0.56)', film_in_range: True, forced_regime: 'laminar', h: 13.3897218737, h_forced: 11.7097119069, h_natural: 9.26366880706, natural_regime: 'laminar', regime: 'mixed', solver_advisory: {Ra: 654587.521306, Re: 1282.44635414, forced_regime: 'laminar', natural_regime: 'laminar', reason: 'Mixed convection regime (Ri between 0.1 and 10). Both forced and buoyant effects are significant. Verify convergence carefully in steady-state, or consider transient solve for more reliable results.', severity: 'caution', steady_state_ok: True, transient_advisory: False}, t_film_K: 350}\n"
    "0.5 400.0 300.0 0.05 horizontal_up: {Ra: 654587.521306, Re: 1282.44635414, Ri: 0.560571428571, dominant_mode: 'mixed-assisting (Ri=0.56)', film_in_range: True, forced_regime: 'laminar', h: 13.0356134608, h_forced: 11.7097119069, h_natural: 8.4786121285, natural_regime: 'laminar', regime: 'mixed', solver_advisory: {Ra: 654587.521306, Re: 1282.44635414, forced_regime: 'laminar', natural_regime: 'laminar', reason: 'Mixed convection regime (Ri between 0.1 and 10). Both forced and buoyant effects are significant. Verify convergence carefully in steady-state, or consider transient solve for more reliable results.', severity: 'caution', steady_state_ok: True, transient_advisory: False}, t_film_K: 350}\n"
    "0.5 400.0 300.0 0.05 horizontal_down: {Ra: 654587.521306, Re: 1282.44635414, Ri: 0.560571428571, dominant_mode: 'mixed-opposing (Ri=0.56)', film_in_range: True, forced_regime: 'laminar', h: 11.5214899899, h_forced: 11.7097119069, h_natural: 4.23930606425, natural_regime: 'laminar', regime: 'mixed', solver_advisory: {Ra: 654587.521306, Re: 1282.44635414, forced_regime: 'laminar', natural_regime: 'laminar', reason: 'Mixed convection regime (Ri between 0.1 and 10). Both forced and buoyant effects are significant. Verify convergence carefully in steady-state, or consider transient solve for more reliable results.', severity: 'caution', steady_state_ok: True, transient_advisory: False}, t_film_K: 350}\n"
    "0.5 400.0 300.0 0.4 vertical: {Ra: 335148810.909, Re: 10259.5708331, Ri: 4.48457142857, dominant_mode: 'mixed-assisting (Ri=4.48)', film_in_range: True, forced_regime: 'laminar', h: 6.19784509241, h_forced: 4.14000834756, h_natural: 5.50821042819, natural_regime: 'laminar', regime: 'mixed', solver_advisory: {Ra: 335148810.909, Re: 10259.5708331, forced_regime: 'laminar', natural_regime: 'laminar', reason: 'Mixed convection regime (Ri between 0.1 and 10). Both forced and buoyant effects are significant. Verify convergence carefully in steady-state, or consider transient solve for more reliable results.', severity: 'caution', steady_state_ok: True, transient_advisory: False}, t_film_K: 350}\n"
    "0.5 400.0 300.0 0.4 horizontal_up: {Ra: 335148810.909, Re: 10259.5708331, Ri: 4.48457142857, dominant_mode: 'mixed-assisting (Ri=4.48)', film_in_range: True, forced_regime: 'laminar', h: 7.62053455663, h_forced: 4.14000834756, h_natural: 7.18929400018, natural_regime: 'turbulent', regime: 'mixed', solver_advisory: {Ra: 335148810.909, Re: 10259.5708331, forced_regime: 'laminar', natural_regime: 'turbulent', reason: 'Turbulent natural convection (Ra=3.35e+08, threshold=1e+07 for horizontal_up). Buoyant plumes are inherently unsteady — steady-state solver may not converge or may give misleading time-averaged results. Recommend transient solve with time-averaged output, or use relaxation with caution in steady-state.', severity: 'warning', steady_state_ok: False, transient_advisory: True}, t_film_K: 350}\n"
    "0.5 400.0 300.0 0.4 horizontal_down: {Ra: 335148810.909, Re: 10259.5708331, Ri: 4.48457142857, dominant_mode: 'mixed-opposing (Ri=4.48)', film_in_range: True, forced_regime: 'laminar', h: 3.80161294413, h_forced: 4.14000834756, h_natural: 2.52070646714, natural_regime: 'laminar', regime: 'mixed', solver_advisory: {Ra: 335148810.909, Re: 10259.5708331, forced_regime: 'laminar', natural_regime: 'laminar', reason: 'Mixed convection regime (Ri between 0.1 and 10). Both forced and buoyant effects are significant. Verify convergence carefully in steady-state, or consider transient solve for more reliable results.', severity: 'caution', steady_state_ok: True, transient_advisory: False}, t_film_K: 350}\n"
    "0.5 900.0 330.0 0.05 vertical: {Ra: 312385.928901, Re: 491.889085151, Ri: 1.81843902439, dominant_mode: 'mixed-assisting (Ri=1.82)', film_in_range: True, forced_regime: 'laminar', h: 15.7634445442, h_forced: 12.1261468376, h_natural: 12.8743875558, natural_regime: 'laminar', regime: 'mixed', solver_advisory: {Ra: 312385.928901, Re: 491.889085151, forced_regime: 'laminar', natural_regime: 'laminar', reason: 'Mixed convection regime (Ri between 0.1 and 10). Both forced and buoyant effects are significant. Verify convergence carefully in steady-state, or consider transient solve for more reliable results.', severity: 'caution', steady_state_ok: True, transient_advisory: False}, t_film_K: 615}\n"
    "0.5 900.0 330.0 0.05 horizontal_up: {Ra: 312385.928901, Re: 491.889085151, Ri: 1.81843902439, dominant_mode: 'mixed-assisting (Ri=1.82)', film_in_range: True, forced_regime: 'laminar', h: 15.0651271675, h_forced: 12.1261468376, h_natural: 11.783337763, natural_regime: 'laminar', regime: 'mixed', solver_advisory: {Ra: 312385.928901, Re: 491.889085151, forced_regime: 'laminar', natural_regime: 'laminar', reason: 'Mixed convection regime (Ri between 0.1 and 10). Both forced and buoyant effects are significant. Verify convergence carefully in steady-state, or consider transient solve for more reliable results.', severity: 'caution', steady_state_ok: True, transient_advisory: False}, t_film_K: 615}\n"
    "0.5 900.0 330.0 0.05 horizontal_down: {Ra: 312385.928901, Re: 491.889085151, Ri: 1.81843902439, dominant_mode: 'mixed-opposing (Ri=1.82)', film_in_range: True, forced_regime: 'laminar', h: 11.6435936482, h_forced: 12.1261468376, h_natural: 5.89166888148, natural_regime: 'laminar', regime: 'mixed', solver_advisory: {Ra: 312385.928901, Re: 491.889085151, forced_regime: 'laminar', natural_regime: 'laminar', reason: 'Mixed convection regime (Ri between 0.1 and 10). Both forced and buoyant effects are significant. Verify convergence carefully in steady-state, or consider transient solve for more reliable results.', severity: 'caution', steady_state_ok: True, transient_advisory: False}, t_film_K: 615}\n"
    "0.5 900.0 330.0 0.4 vertical: {Ra: 159941595.597, Re: 3935.11268121, Ri: 14.5475121951, dominant_mode: 'natural convection (vertical, Ra=1.6e+08)', film_in_range: True, forced_regime: 'laminar', h: 7.65515664135, h_forced: 4.28724032927, h_natural: 7.65515664135, natural_regime: 'laminar', regime: 'natural', solver_advisory: {Ra: 159941595.597, Re: 3935.11268121, forced_regime: 'laminar', natural_regime: 'laminar', reason: 'Laminar convection — steady-state solve is appropriate.', severity: 'info', steady_state_ok: True, transient_advisory: False}, t_film_K: 615}\n"
    "0.5 900.0 330.0 0.4 horizontal_up: {Ra: 159941595.597, Re: 3935.11268121, Ri: 14.5475121951, dominant_mode: 'natural convection (horizontal_up, Ra=1.6e+08)', film_in_range: True, forced_regime: 'laminar', h: 9.39413465799, h_forced: 4.28724032927, h_natural: 9.39413465799, natural_regime: 'turbulent', regime: 'natural', solver_advisory: {Ra: 159941595.597, Re: 3935.11268121, forced_regime: 'laminar', natural_regime: 'turbulent', reason: 'Turbulent natural convection (Ra=1.60e+08, threshold=1e+07 for horizontal_up). Buoyant plumes are inherently unsteady — steady-state solver may not converge or may give misleading time-averaged results. Recommend transient solve with time-averaged output, or use relaxation with caution in steady-state.', severity: 'warning', steady_state_ok: False, transient_advisory: True}, t_film_K: 615}\n"
    "0.5 900.0 330.0 0.4 horizontal_down: {Ra: 159941595.597, Re: 3935.11268121, Ri: 14.5475121951, dominant_mode: 'natural convection (horizontal_down, Ra=1.6e+08)', film_in_range: True, forced_regime: 'laminar', h: 3.50320727655, h_forced: 4.28724032927, h_natural: 3.50320727655, natural_regime: 'laminar', regime: 'natural', solver_advisory: {Ra: 159941595.597, Re: 3935.11268121, forced_regime: 'laminar', natural_regime: 'laminar', reason: 'Laminar convection — steady-state solve is appropriate.', severity: 'info', steady_state_ok: True, transient_advisory: False}, t_film_K: 615}\n"
    "4.47 400.0 300.0 0.05 vertical: {Ra: 654587.521306, Re: 11465.070406, Ri: 0.00701384107537, dominant_mode: 'forced convection (Re=11465)', film_in_range: True, forced_regime: 'laminar', h: 35.0118427866, h_forced: 35.0118427866, h_natural: 9.26366880706, natural_regime: 'laminar', regime: 'forced', solver_advisory: {Ra: 654587.521306, Re: 11465.070406, forced_regime: 'laminar', natural_regime: 'laminar', reason: 'Laminar convection — steady-state solve is appropriate.', severity: 'info', steady_state_ok: True, transient_advisory: False}, t_film_K: 350}\n"
    "4.47 400.0 300.0 0.05 horizontal_up: {Ra: 654587.521306, Re: 11465.070406, Ri: 0.00701384107537, dominant_mode: 'forced convection (Re=11465)', film_in_range: True, forced_regime: 'laminar', h: 35.0118427866, h_forced: 35.0118427866, h_natural: 8.4786121285, natural_regime: 'laminar', regime: 'forced', solver_advisory: {Ra: 654587.521306, Re: 11465.070406, forced_regime: 'laminar', natural_regime: 'laminar', reason: 'Laminar convection — steady-state solve is appropriate.', severity: 'info', steady_state_ok: True, transient_advisory: False}, t_film_K: 350}\n"
    "4.47 400.0 300.0 0.05 horizontal_down: {Ra: 654587.521306, Re: 11465.070406, Ri: 0.00701384107537, dominant_mode: 'forced convection (Re=11465)', film_in_range: True, forced_regime: 'laminar', h: 35.0118427866, h_forced: 35.0118427866, h_natural: 4.23930606425, natural_regime: 'laminar', regime: 'forced', solver_advisory: {Ra: 654587.521306, Re: 11465.070406, forced_regime: 'laminar', natural_regime: 'laminar', reason: 'Laminar convection — steady-state solve is appropriate.', severity: 'info', steady_state_ok: True, transient_advisory: False}, t_film_K: 350}\n"
    "4.47 400.0 300.0 0.4 vertical: {Ra: 335148810.909, Re: 91720.5632483, Ri: 0.056110728603, dominant_mode: 'forced convection (Re=91721)', film_in_range: True, forced_regime: 'laminar', h: 12.3785557281, h_forced: 12.3785557281, h_natural: 5.50821042819, natural_regime: 'laminar', regime: 'forced', solver_advisory: {Ra: 335148810.909, Re: 91720.5632483, forced_regime: 'laminar', natural_regime: 'laminar', reason: 'Laminar convection — steady-state solve is appropriate.', severity: 'info', steady_state_ok: True, transient_advisory: False}, t_film_K: 350}\n"
    "4.47 400.0 300.0 0.4 horizontal_up: {Ra: 335148810.909, Re: 91720.5632483, Ri: 0.056110728603, dominant_mode: 'forced convection (Re=91721)', film_in_range: True, forced_regime: 'laminar', h: 12.3785557281, h_forced: 12.3785557281, h_natural: 7.18929400018, natural_regime: 'turbulent', regime: 'forced', solver_advisory: {Ra: 335148810.909, Re: 91720.5632483, forced_regime: 'laminar', natural_regime: 'turbulent', reason: 'Forced convection dominates (Ri<0.1), but natural convection is turbulent (Ra=3.35e+08). If forced flow is intermittent or uncertain, buoyancy oscillations could affect results. Monitor convergence.', severity: 'caution', steady_state_ok: True, transient_advisory: False}, t_film_K: 350}\n"
    "4.47 400.0 300.0 0.4 horizontal_down: {Ra: 335148810.909, Re: 91720.5632483, Ri: 0.056110728603, dominant_mode: 'forced convection (Re=91721)', film_in_range: True, forced_regime: 'laminar', h: 12.3785557281, h_forced: 12.3785557281, h_natural: 2.52070646714, natural_regime: 'laminar', regime: 'forced', solver_advisory: {Ra: 335148810.909, Re: 91720.5632483, forced_regime: 'laminar', natural_regime: 'laminar', reason: 'Laminar convection — steady-state solve is appropriate.', severity: 'info', steady_state_ok: True, transient_advisory: False}, t_film_K: 350}\n"
    "4.47 900.0 330.0 0.05 vertical: {Ra: 312385.928901, Re: 4397.48842125, Ri: 0.0227522161713, dominant_mode: 'forced convection (Re=4397)', film_in_range: True, forced_regime: 'laminar', h: 36.2569762655, h_forced: 36.2569762655, h_natural: 12.8743875558, natural_regime: 'laminar', regime: 'forced', solver_advisory: {Ra: 312385.928901, Re: 4397.48842125, forced_regime: 'laminar', natural_regime: 'laminar', reason: 'Laminar convection — steady-state solve is appropriate.', severity: 'info', steady_state_ok: True, transient_advisory: False}, t_film_K: 615}\n"
    "4.47 900.0 330.0 0.05 horizontal_up: {Ra: 312385.928901, Re: 4397.48842125, Ri: 0.0227522161713, dominant_mode: 'forced convection (Re=4397)', film_in_range: True, forced_regime: 'laminar', h: 36.2569762655, h_forced: 36.2569762655, h_natural: 11.783337763, natural_regime: 'laminar', regime: 'forced', solver_advisory: {Ra: 312385.928901, Re: 4397.48842125, forced_regime: 'laminar', natural_regime: 'laminar', reason: 'Laminar convection — steady-state solve is appropriate.', severity: 'info', steady_state_ok: True, transient_advisory: False}, t_film_K: 615}\n"
    "4.47 900.0 330.0 0.05 horizontal_down: {Ra: 312385.928901, Re: 4397.48842125, Ri: 0.0227522161713, dominant_mode: 'forced convection (Re=4397)', film_in_range: True, forced_regime: 'laminar', h: 36.2569762655, h_forced: 36.2569762655, h_natural: 5.89166888148, natural_regime: 'laminar', regime: 'forced', solver_advisory: {Ra: 312385.928901, Re: 4397.48842125, forced_regime: 'laminar', natural_regime: 'laminar', reason: 'Laminar convection — steady-state solve is appropriate.', severity: 'info', steady_state_ok: True, transient_advisory: False}, t_film_K: 615}\n"
    "4.47 900.0 330.0 0.4 vertical: {Ra: 159941595.597, Re: 35179.90737, Ri: 0.182017729371, dominant_mode: 'mixed-assisting (Ri=0.18)', film_in_range: True, forced_regime: 'laminar', h: 13.6708917864, h_forced: 12.8187768913, h_natural: 7.65515664135, natural_regime: 'laminar', regime: 'mixed', solver_advisory: {Ra: 159941595.597, Re: 35179.90737, forced_regime: 'laminar', natural_regime: 'laminar', reason: 'Mixed convection regime (Ri between 0.1 and 10). Both forced and buoyant effects are significant. Verify convergence carefully in steady-state, or consider transient solve for more reliable results.', severity: 'caution', steady_state_ok: True, transient_advisory: False}, t_film_K: 615}\n"
    "4.47 900.0 330.0 0.4 horizontal_up: {Ra: 159941595.597, Re: 35179.90737, Ri: 0.182017729371, dominant_mode: 'mixed-assisting (Ri=0.18)', film_in_range: True, forced_regime: 'laminar', h: 14.3182627246, h_forced: 12.8187768913, h_natural: 9.39413465799, natural_regime: 'turbulent', regime: 'mixed', solver_advisory: {Ra: 159941595.597, Re: 35179.90737, forced_regime: 'laminar', natural_regime: 'turbulent', reason: 'Turbulent natural convection (Ra=1.60e+08, threshold=1e+07 for horizontal_up). Buoyant plumes are inherently unsteady — steady-state solver may not converge or may give misleading time-averaged results. Recommend transient solve with time-averaged output, or use relaxation with caution in steady-state.', severity: 'warning', steady_state_ok: False, transient_advisory: True}, t_film_K: 615}\n"
    "4.47 900.0 330.0 0.4 horizontal_down: {Ra: 159941595.597, Re: 35179.90737, Ri: 0.182017729371, dominant_mode: 'mixed-opposing (Ri=0.18)', film_in_range: True, forced_regime: 'laminar', h: 12.7309633379, h_forced: 12.8187768913, h_natural: 3.50320727655, natural_regime: 'laminar', regime: 'mixed', solver_advisory: {Ra: 159941595.597, Re: 35179.90737, forced_regime: 'laminar', natural_regime: 'laminar', reason: 'Mixed convection regime (Ri between 0.1 and 10). Both forced and buoyant effects are significant. Verify convergence carefully in steady-state, or consider transient solve for more reliable results.', severity: 'caution', steady_state_ok: True, transient_advisory: False}, t_film_K: 615}\n"
    "15.0 400.0 300.0 0.05 vertical: {Ra: 654587.521306, Re: 38473.3906243, Ri: 0.000622857142857, dominant_mode: 'forced convection (Re=38473)', film_in_range: True, forced_regime: 'laminar', h: 64.1367335331, h_forced: 64.1367335331, h_natural: 9.26366880706, natural_regime: 'laminar', regime: 'forced', solver_advisory: {Ra: 654587.521306, Re: 38473.3906243, forced_regime: 'laminar', natural_regime: 'laminar', reason: 'Laminar convection — steady-state solve is appropriate.', severity: 'info', steady_state_ok: True, transient_advisory: False}, t_film_K: 350}\n"
    "15.0 400.0 300.0 0.05 horizontal_up: {Ra: 654587.521306, Re: 38473.3906243, Ri: 0.000622857142857, dominant_mode: 'forced convection (Re=38473)', film_in_range: True, forced_regime: 'laminar', h: 64.1367335331, h_forced: 64.1367335331, h_natural: 8.4786121285, natural_regime: 'laminar', regime: 'forced', solver_advisory: {Ra: 654587.521306, Re: 38473.3906243, forced_regime: 'laminar', natural_regime: 'laminar', reason: 'Laminar convection — steady-state solve is appropriate.', severity: 'info', steady_state_ok: True, transient_advisory: False}, t_film_K: 350}\n"
    "15.0 400.0 300.0 0.05 horizontal_down: {Ra: 654587.521306, Re: 38473.3906243, Ri: 0.000622857142857, dominant_mode: 'forced convection (Re=38473)', film_in_range: True, forced_regime: 'laminar', h: 64.1367335331, h_forced: 64.1367335331, h_natural: 4.23930606425, natural_regime: 'laminar', regime: 'forced', solver_advisory: {Ra: 654587.521306, Re: 38473.3906243, forced_regime: 'laminar', natural_regime: 'laminar', reason: 'Laminar convection — steady-state solve is appropriate.', severity: 'info', steady_state_ok: True, transient_advisory: False}, t_film_K: 350}\n"
    "15.0 400.0 300.0 0.4 vertical: {Ra: 335148810.909, Re: 307787.124994, Ri: 0.00498285714286, dominant_mode: 'forced convection (Re=307787)', film_in_range: True, forced_regime: 'laminar', h: 22.6757596022, h_forced: 22.6757596022, h_natural: 5.50821042819, natural_regime: 'laminar', regime: 'forced', solver_advisory: {Ra: 335148810.909, Re: 307787.124994, forced_regime: 'laminar', natural_regime: 'laminar', reason: 'Laminar convection — steady-state solve is appropriate.', severity: 'info', steady_state_ok: True, transient_advisory: False}, t_film_K: 350}\n"
    "15.0 400.0 300.0 0.4 horizontal_up: {Ra: 335148810.909, Re: 307787.124994, Ri: 0.00498285714286, dominant_mode: 'forced convection (Re=307787)', film_in_range: True, forced_regime: 'laminar', h: 22.6757596022, h_forced: 22.6757596022, h_natural: 7.18929400018, natural_regime: 'turbulent', regime: 'forced', solver_advisory: {Ra: 335148810.909, Re: 307787.124994, forced_regime: 'laminar', natural_regime: 'turbulent', reason: 'Forced convection dominates (Ri<0.1), but natural convection is turbulent (Ra=3.35e+08). If forced flow is intermittent or uncertain, buoyancy oscillations could affect results. Monitor convergence.', severity: 'caution', steady_state_ok: True, transient_advisory: False}, t_film_K: 350}\n"
    "15.0 400.0 300.0 0.4 horizontal_down: {Ra: 335148810.909, Re: 307787.124994, Ri: 0.00498285714286, dominant_mode: 'forced convection (Re=307787)', film_in_range: True, forced_regime: 'laminar', h: 22.6757596022, h_forced: 22.6757596022, h_natural: 2.52070646714, natural_regime: 'laminar', regime: 'forced', solver_advisory: {Ra: 335148810.909, Re: 307787.124994, forced_regime: 'laminar', natural_regime: 'laminar', reason: 'Laminar convection — steady-state solve is appropriate.', severity: 'info', steady_state_ok: True, transient_advisory: False}, t_film_K: 350}\n"
    "15.0 900.0 330.0 0.05 vertical: {Ra: 312385.928901, Re: 14756.6725545, Ri: 0.00202048780488, dominant_mode: 'forced convection (Re=14757)', film_in_range: True, forced_regime: 'laminar', h: 66.4176415858, h_forced: 66.4176415858, h_natural: 12.8743875558, natural_regime: 'laminar', regime: 'forced', solver_advisory: {Ra: 312385.928901, Re: 14756.6725545, forced_regime: 'laminar', natural_regime: 'laminar', reason: 'Laminar convection — steady-state solve is appropriate.', severity: 'info', steady_state_ok: True, transient_advisory: False}, t_film_K: 615}\n"
    "15.0 900.0 330.0 0.05 horizontal_up: {Ra: 312385.928901, Re: 14756.6725545, Ri: 0.00202048780488, dominant_mode: 'forced convection (Re=14757)', film_in_range: True, forced_regime: 'laminar', h: 66.4176415858, h_forced: 66.4176415858, h_natural: 11.783337763, natural_regime: 'laminar', regime: 'forced', solver_advisory: {Ra: 312385.928901, Re: 14756.6725545, forced_regime: 'laminar', natural_regime: 'laminar', reason: 'Laminar convection — steady-state solve is appropriate.', severity: 'info', steady_state_ok: True, transient_advisory: False}, t_film_K: 615}\n"
    "15.0 900.0 330.0 0.05 horizontal_down: {Ra: 312385.928901, Re: 14756.6725545, Ri: 0.00202048780488, dominant_mode: 'forced convection (Re=14757)', film_in_range: True, forced_regime: 'laminar', h: 66.4176415858, h_forced: 66.4176415858, h_natural: 5.89166888148, natural_regime: 'laminar', regime: 'forced', solver_advisory: {Ra: 312385.928901, Re: 14756.6725545, forced_regime: 'laminar', natural_regime: 'laminar', reason: 'Laminar convection — steady-state solve is appropriate.', severity: 'info', steady_state_ok: True, transient_advisory: False}, t_film_K: 615}\n"
    "15.0 900.0 330.0 0.4 vertical: {Ra: 159941595.597, Re: 118053.380436, Ri: 0.016163902439, dominant_mode: 'forced convection (Re=118053)', film_in_range: True, forced_regime: 'laminar', h: 23.4821823779, h_forced: 23.4821823779, h_natural: 7.65515664135, natural_regime: 'laminar', regime: 'forced', solver_advisory: {Ra: 159941595.597, Re: 118053.380436, forced_regime: 'laminar', natural_regime: 'laminar', reason: 'Laminar convection — steady-state solve is appropriate.', severity: 'info', steady_state_ok: True, transient_advisory: False}, t_film_K: 615}\n"
    "15.0 900.0 330.0 0.4 horizontal_up: {Ra: 159941595.597, Re: 118053.380436, Ri: 0.016163902439, dominant_mode: 'forced convection (Re=118053)', film_in_range: True, forced_regime: 'laminar', h: 23.4821823779, h_forced: 23.4821823779, h_natural: 9.39413465799, natural_regime: 'turbulent', regime: 'forced', solver_advisory: {Ra: 159941595.597, Re: 118053.380436, forced_regime: 'laminar', natural_regime: 'turbulent', reason: 'Forced convection dominates (Ri<0.1), but natural convection is turbulent (Ra=1.60e+08). If forced flow is intermittent or uncertain, buoyancy oscillations could affect results. Monitor convergence.', severity: 'caution', steady_state_ok: True, transient_advisory: False}, t_film_K: 615}\n"
    "15.0 900.0 330.0 0.4 horizontal_down: {Ra: 159941595.597, Re: 118053.380436, Ri: 0.016163902439, dominant_mode: 'forced convection (Re=118053)', film_in_range: True, forced_regime: 'laminar', h: 23.4821823779, h_forced: 23.4821823779, h_natural: 3.50320727655, natural_regime: 'laminar', regime: 'forced', solver_advisory: {Ra: 159941595.597, Re: 118053.380436, forced_regime: 'laminar', natural_regime: 'laminar', reason: 'Laminar convection — steady-state solve is appropriate.', severity: 'info', steady_state_ok: True, transient_advisory: False}, t_film_K: 615}\n"
)
