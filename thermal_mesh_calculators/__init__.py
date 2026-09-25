"""
Thermal Mesh Calculators for Automotive CAE

Physics-driven mesh sizing tools that determine maximum element sizes
based on material properties, boundary conditions, and heat transfer physics
rather than arbitrary rules of thumb.

Modules:
    conduction  - Boundary-driven conduction mesh sizing (Fourier's Law)
    convection  - Biot number evaluation, y+ mapping, h-gradient resolution
    radiation   - T^4 sensitivity and view factor mesh constraints
    shields     - Single-layer and multilayer heat shield solvers (Newton-Raphson)
    transient   - Penetration depth, Fourier number, and drive-cycle constraints
    zones       - Convection zone lookup and spatial gradient estimation
    h_estimator - Convective HTC from correlations (forced, natural, mixed)
    boundary_layer - Aerodynamic boundary layer mesh constraints
    batch       - Bulk component processing (BOM → mesh sizes)
"""

from thermal_mesh_calculators.conduction import BoundaryDrivenConductionCalculator
from thermal_mesh_calculators.convection import ConvectionMeshCalculator
from thermal_mesh_calculators.radiation import RadiationMeshCalculator
from thermal_mesh_calculators.shields import (
    SingleLayerShieldCalculator,
    MultilayerShieldCalculator,
)
from thermal_mesh_calculators.transient import TransientMeshCalculator
from thermal_mesh_calculators.batch import (
    process_batch,
    process_part,
    summary_table,
    MATERIALS,
    SURFACE_TREATMENTS,
)
from thermal_mesh_calculators.zones import (
    CONVECTION_ZONES,
    get_zone,
    get_conservative_h,
    get_zone_air_temp,
    estimate_spatial_gradient,
)
from thermal_mesh_calculators.boundary_layer import BoundaryLayerCalculator
from thermal_mesh_calculators.h_estimator import (
    air_properties,
    forced_convection_flat_plate,
    natural_convection,
    richardson_number,
    estimate_h,
    solver_advisory,
    h_from_velocity,
    H_EXTERNAL_MAX,
    H_EXTERNAL_TYPICAL,
)

__version__ = "0.6.2.post1"
