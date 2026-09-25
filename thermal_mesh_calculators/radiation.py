"""
Radiation Mesh Calculator
=========================

Radiation is often the least sensitive heat transfer mode for typical
automotive structural meshes at moderate temperatures.  However, it
becomes the *dominant* mesh driver in two regimes:

1.  High-temperature components (exhaust manifolds, turbo housings,
    catalytic converters) where T^4 dependence amplifies small spatial
    temperature differences into large flux differences.

2.  Complex view-factor geometries (nested shields, concave surfaces)
    where geometric fidelity of the mesh facets directly controls the
    accuracy of the calculated view factors.

Physics:
    Stefan-Boltzmann:  q'' = eps * sigma * (T_s^4 - T_surr^4)

    The derivative dq/dT = 4 * eps * sigma * T^3 gives the local
    sensitivity of radiant flux to temperature.  At eps=0.8, this is
    ~4.9 W/m^2 K at 300 K and ~181 W/m^2 K at 1000 K — a 37x increase
    (T^3 scaling).  A 10 K uncertainty across an element face causes a
    ~49 W/m^2 flux error at room temperature but ~1814 W/m^2 error
    at 1000 K.

    The mesh must therefore be refined in regions where:
    (a) temperatures are high AND spatial gradients exist, and
    (b) surface curvature is high enough that coarse faceting would
        distort view factors.

Emissivity boundaries:
    Where surface treatment changes abruptly (bare metal → coated,
    aluminised → oxidised), the mesh edges must align with the
    material boundary.  This is a meshing topology constraint, not
    a density constraint, but is flagged here for completeness.
"""

import math
from thermal_mesh_calculators.constants import STEFAN_BOLTZMANN
from thermal_mesh_calculators._guards import (
    require_fraction,
    require_non_negative,
    require_positive,
)


class RadiationMeshCalculator:
    """
    Evaluates radiation-driven mesh size constraints.
    """

    @staticmethod
    def flux_sensitivity(t_local: float, emissivity: float) -> float:
        """
        Local sensitivity of radiant flux to temperature (W/m^2 per K).

        This is dq/dT = 4 * eps * sigma * T^3.

        Parameters
        ----------
        t_local : float     — local surface temperature (K)
        emissivity : float  — surface emissivity (0-1)

        Returns
        -------
        float — dq/dT (W/m^2 K)

        Raises ValueError for t_local not > 0 or emissivity outside [0, 1].
        """
        require_positive("t_local", t_local, "K")
        require_fraction("emissivity", emissivity)
        return 4.0 * emissivity * STEFAN_BOLTZMANN * (t_local ** 3)

    @staticmethod
    def max_mesh_size(
        t_local: float,
        emissivity: float,
        allowable_flux_error: float,
        spatial_gradient: float,
    ) -> dict:
        """
        Maximum element size to keep the radiation linearisation error
        across an element face within an acceptable bound.

        The element spans a temperature range of (dT/dx * dx).  The
        resulting flux error is approximately (dq/dT) * (dT/dx * dx).
        Bounding this by the allowable error and solving for dx:

            dx_max = allowable_flux_error / (dq/dT * dT/dx)

        Parameters
        ----------
        t_local : float
            Local surface temperature (K).
        emissivity : float
            Surface emissivity (0-1).
        allowable_flux_error : float
            Acceptable flux error per element face (W/m^2).
            Typical values: 100-500 W/m^2 for automotive exhaust,
            50-100 W/m^2 for precision thermal management.
        spatial_gradient : float
            Expected spatial temperature gradient (K/m).

        Returns
        -------
        dict with keys:
            max_dx_mm       : float — maximum element size (mm)
            dq_dt           : float — flux sensitivity (W/m^2 K)
            max_dt_element  : float — max delta-T across element (K)

        Raises
        ------
        ValueError
            t_local or allowable_flux_error not > 0, emissivity outside
            [0, 1], or spatial_gradient < 0.
        """
        require_positive("t_local", t_local, "K")
        require_fraction("emissivity", emissivity)
        require_positive("allowable_flux_error", allowable_flux_error, "W/m^2")
        require_non_negative("spatial_gradient", spatial_gradient, "K/m")
        if spatial_gradient == 0:
            return {
                "max_dx_mm": float("inf"),
                "dq_dt": 0.0,
                "max_dt_element": float("inf"),
            }

        dq_dt = 4.0 * emissivity * STEFAN_BOLTZMANN * (t_local ** 3)
        max_dt_element = allowable_flux_error / dq_dt if dq_dt > 0 else float("inf")
        max_dx_m = max_dt_element / spatial_gradient

        return {
            "max_dx_mm": max_dx_m * 1000.0,
            "dq_dt": dq_dt,
            "max_dt_element": max_dt_element,
        }

    @staticmethod
    def view_factor_curvature_limit(
        radius_mm: float,
        max_facet_angle_deg: float = 15.0,
    ) -> float:
        """
        Maximum element edge length on a curved surface so that the
        flat-facet approximation stays within a given angular tolerance.

        A chord of length L on a circle of radius R subtends an angle
        theta = 2 * arcsin(L / 2R).  Inverting:

            L_max = 2 R sin(theta_max / 2)

        Parameters
        ----------
        radius_mm : float
            Local radius of curvature (mm).
        max_facet_angle_deg : float
            Maximum acceptable angle subtended by a single facet (deg).
            15 deg is a reasonable default; tighter for concave surfaces
            facing a heat source.

        Returns
        -------
        float — maximum element edge length (mm)
        """
        theta_rad = math.radians(max_facet_angle_deg)
        return 2.0 * radius_mm * math.sin(theta_rad / 2.0)
