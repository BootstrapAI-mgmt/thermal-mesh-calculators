"""
Convection Mesh Calculator
==========================

The solid-side mesh requirement for convection is dictated by the fluid
behaviour at the boundary.  There are two independent constraints:

1.  Biot number (Bi = h Lc / k)
    Determines whether the part needs through-thickness 3D elements or
    whether a 2D shell is physically valid.

    - Bi << 0.1  →  lumped capacitance applies, shell mesh OK
    - Bi >  0.1  →  internal temperature gradients exist, need 3D solid
                     mesh with multiple elements across thickness

2.  Surface h-gradient resolution
    Where the convective heat transfer coefficient varies sharply over
    a short distance (flow impingement, separation, stagnation lines),
    the solid surface mesh must be fine enough to capture those spatial
    variations.  If the solid mesh is coarser than the spatial wavelength
    of h, the solver smears the peak flux and under-predicts local temps.

    The practical constraint is:

        dx_solid  <=  L_h_gradient

    where L_h_gradient is the characteristic length over which h changes
    significantly (e.g. halves or doubles).

y+ considerations for CFD-coupled workflows:
    - Wall-function CFD (y+ ~ 30-300):  fluid wall faces are relatively
      large.  Solid mesh can be similar or slightly finer.
    - Resolved boundary layer (y+ ~ 1):  fluid mesh is extremely fine
      near the wall.  Mapping those gradients onto a coarse solid mesh
      introduces interpolation errors and artificially diffuses heat flux.
      Rule of thumb: solid surface element size should not exceed ~3-5x
      the fluid wall face size to avoid significant mapping loss.
"""

import math
from thermal_mesh_calculators._guards import require_non_negative, require_positive


class ConvectionMeshCalculator:
    """
    Evaluates convection-driven mesh constraints on the solid side.
    """

    @staticmethod
    def biot_number(
        h: float,
        k: float,
        thickness: float,
    ) -> dict:
        """
        Evaluate the Biot number to determine shell vs. solid meshing.

        Parameters
        ----------
        h : float
            Convective HTC (W/m^2 K).
        k : float
            Solid thermal conductivity (W/m K).
        thickness : float
            Material / wall thickness (m).

        Returns
        -------
        dict with keys:
            biot          : float  — Biot number
            mesh_type     : str    — "2D Shell" or "3D Solid"
            min_elements  : int    — minimum elements through thickness
                                     (1 for shell, ceil(20*Bi) for solid,
                                      clamped to [2, 10])
            rationale     : str

        Raises
        ------
        ValueError
            k not > 0, or h < 0.
        """
        require_positive("k", k, "W/m K")
        require_non_negative("h", h, "W/m^2 K")
        lc = thickness / 2.0
        bi = (h * lc) / k

        if bi < 0.1:
            return {
                "biot": bi,
                "mesh_type": "2D Shell",
                "min_elements_through_thickness": 1,
                "rationale": (
                    f"Bi = {bi:.4f} < 0.1 — internal conduction resistance "
                    f"is negligible relative to surface convection.  A single "
                    f"shell element across the thickness is physically valid."
                ),
            }
        else:
            # Cell Biot heuristic: n >= 20*Bi ensures Bi_cell = h*dx/k <= 0.1
            # per element, preventing numerical oscillation in explicit solvers.
            # The factor of 20 (rather than 10) accounts for the half-thickness
            # definition of Bi and provides margin for non-uniform h.
            n_elem = max(2, min(10, math.ceil(20 * bi)))
            return {
                "biot": bi,
                "mesh_type": "3D Solid",
                "min_elements_through_thickness": n_elem,
                "rationale": (
                    f"Bi = {bi:.4f} >= 0.1 — significant through-thickness "
                    f"temperature gradient exists.  Use at least {n_elem} "
                    f"elements across the {thickness*1000:.1f} mm thickness "
                    f"(element size ~ {thickness/n_elem*1000:.2f} mm)."
                ),
            }

    @staticmethod
    def h_gradient_mesh_limit(
        h_max: float,
        h_min: float,
        gradient_length_mm: float,
    ) -> dict:
        """
        Estimate the surface mesh size needed to resolve a spatial
        variation in the heat transfer coefficient.

        Parameters
        ----------
        h_max : float
            Peak HTC in the gradient region (W/m^2 K).
        h_min : float
            Trough HTC in the gradient region (W/m^2 K).
        gradient_length_mm : float
            Physical distance over which h transitions from h_min to
            h_max (mm).

        Returns
        -------
        dict with keys:
            max_dx_mm     : float — recommended max element size (mm)
            h_ratio       : float — h_max / h_min
            elements_needed : int — minimum elements across gradient zone
        """
        h_ratio = h_max / h_min if h_min > 0 else float("inf")

        # We want at least 4 elements across the gradient zone to
        # capture the shape, more if the ratio is extreme.
        n_min = max(4, math.ceil(2 * math.log(h_ratio + 1)))
        max_dx = gradient_length_mm / n_min

        return {
            "max_dx_mm": max_dx,
            "h_ratio": h_ratio,
            "elements_needed": n_min,
        }

    @staticmethod
    def cfd_mapping_mesh_limit(
        fluid_wall_face_mm: float,
        mapping_ratio: float = 4.0,
    ) -> float:
        """
        Maximum solid surface element size to avoid significant
        interpolation loss when mapping CFD wall data.

        Parameters
        ----------
        fluid_wall_face_mm : float
            Characteristic size of the CFD wall face (mm).
        mapping_ratio : float
            Max allowable ratio of solid-to-fluid face size.
            Default 4.0 is a practical upper bound; 2.0 is ideal.

        Returns
        -------
        float — max solid surface element size (mm)
        """
        return fluid_wall_face_mm * mapping_ratio
