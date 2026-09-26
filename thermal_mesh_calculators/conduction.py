"""
Boundary-Driven Conduction Mesh Calculator
==========================================

Physics basis:
    At steady state the conduction flux through the solid must equal the
    surface energy balance (convection + radiation out).  Rather than
    requiring q'' as an input (which is a solver *output*), we substitute
    the boundary condition directly:

        q''_cond = h (T_s - T_inf) + eps sigma (T_s^4 - T_surr^4)

    Discretising Fourier's law across the first element at the boundary:

        q''_cond  ~  k * dT_max / dx_max

    Solving for dx_max:

        dx_max = k * dT_max / |q''_conv + q''_rad|

    This yields the maximum element size in terms of quantities that are
    either known inputs (material props, BCs) or can be reasonably
    estimated (surface temperature).

Notes for automotive workflows:
    - Exhaust component surface temperatures are typically fixed in the
      solver, so T_s is known exactly for those parts.
    - Non-exhaust structural components can usually be estimated within
      +/- 30 K from experience or a coarse preliminary run.
    - Heat shields are the main floating unknown; use the shield solvers
      in shields.py for those.

Transient problems:
    The thermal penetration depth bounds the element size one time step
    dt can resolve:

        dx <= sqrt(alpha * dt)

    where alpha = k / (rho * cp) is the thermal diffusivity.
    BoundaryDrivenConductionCalculator.transient_penetration_depth()
    returns this bound.  The full set of transient constraints (the
    penetration depth, a resolution guideline for implicit schemes; the
    Fourier-number limit; the drive-cycle limit) is in transient.py,
    TransientMeshCalculator.
"""

import math
from thermal_mesh_calculators.constants import STEFAN_BOLTZMANN
from thermal_mesh_calculators._guards import (
    require_fraction,
    require_non_negative,
    require_positive,
    require_temperatures,
)

# Biot-number regime thresholds for size_wall() (classical lumped-capacitance
# screen: Bi < 0.1 is thermally thin; Bi > 1 means the wall resistance
# dominates and the 1-D steady derivation itself is suspect).
THIN_BIOT = 0.1
STEEP_BIOT = 1.0

# Below this net flux — relative to the larger individual transfer term — the
# sizing relation has a pole that the physics does not (convection and
# radiation nearly cancel; the part is simply near equilibrium).
NEAR_EQUILIBRIUM_FRACTION = 0.02


class BoundaryDrivenConductionCalculator:
    """
    Computes the maximum conduction element size at a boundary using the
    surface energy balance to eliminate the need for a-priori knowledge of
    the conduction heat flux.
    """

    @staticmethod
    def max_mesh_size(
        k: float,
        h: float,
        t_surf: float,
        t_fluid: float,
        epsilon: float,
        t_surr: float,
        max_dt: float,
    ) -> dict:
        """
        Parameters
        ----------
        k : float
            Thermal conductivity of the solid (W/m K).
        h : float
            Convective heat transfer coefficient at the boundary (W/m^2 K).
        t_surf : float
            Surface temperature (K).  For exhaust components this is the
            fixed boundary condition.  For other parts it is an estimate.
        t_fluid : float
            Free-stream / bulk fluid temperature (K).
        epsilon : float
            Surface emissivity (0-1).
        t_surr : float
            Surrounding radiation sink temperature (K).
        max_dt : float
            Maximum allowable temperature drop across a single element (K).
            This is the user's accuracy target — smaller values yield
            finer meshes and better gradient resolution.

        Returns
        -------
        dict with keys:
            max_dx_mm : float   — maximum element size (mm)
            q_conv    : float   — convection flux component (W/m^2)
            q_rad     : float   — radiation flux component (W/m^2)
            q_total   : float   — total boundary flux magnitude (W/m^2)
            rad_fraction : float — fraction of total flux from radiation

        Raises
        ------
        ValueError
            k, max_dt or a temperature not > 0, h < 0, or epsilon outside
            [0, 1].
        """
        require_positive("k", k, "W/m K")
        require_non_negative("h", h, "W/m^2 K")
        require_temperatures(t_surf=t_surf, t_fluid=t_fluid, t_surr=t_surr)
        require_fraction("epsilon", epsilon)
        require_positive("max_dt", max_dt, "K")

        q_conv = h * (t_surf - t_fluid)
        q_rad = epsilon * STEFAN_BOLTZMANN * (t_surf**4 - t_surr**4)
        q_total = abs(q_conv + q_rad)

        if q_total == 0:
            return {
                "max_dx_mm": float("inf"),
                "q_conv": 0.0,
                "q_rad": 0.0,
                "q_total": 0.0,
                "rad_fraction": 0.0,
            }

        max_dx_m = (k * max_dt) / q_total
        rad_frac = abs(q_rad) / (abs(q_conv) + abs(q_rad)) if (abs(q_conv) + abs(q_rad)) > 0 else 0.0

        return {
            "max_dx_mm": max_dx_m * 1000.0,
            "q_conv": q_conv,
            "q_rad": q_rad,
            "q_total": q_total,
            "rad_fraction": rad_frac,
        }

    @staticmethod
    def size_wall(
        k: float,
        h: float,
        t_surf: float,
        t_fluid: float,
        epsilon: float,
        t_surr: float,
        max_dt: float,
        thickness_m: float,
        max_cells: int = 100,
    ) -> dict:
        """
        Thickness-aware through-thickness sizing: cell count + Biot number.

        ``max_mesh_size`` returns a raw length with no knowledge of the part
        being sized, which produces cell sizes tens of times larger than the
        object whenever the wall is thermally thin (e.g. 85.8 mm for a 1 mm
        aluminium heat-shield layer).  That raw dx is the relation reporting
        "thermally thin", not a mesh size.  This method closes the reporting
        gap: it clamps against the actual wall thickness and leads with the
        two quantities that support a modelling decision — ``n_cells`` and
        the Biot number — retaining the raw dx only for traceability.

        Rearranged, the sizing relation is the Biot number with extra steps:

            dT_wall = q'' L / k = Bi (Ts - Tf)
            N_cells = dT_wall / dT_max = Bi (Ts - Tf) / dT_max

        so this is a thermal-thinness screen (which solids need
        through-thickness resolution vs. shell conduction), not an error
        estimator: under the derivation's own assumptions the profile is
        linear and one cell captures it exactly.

        Parameters
        ----------
        k : float
            Thermal conductivity of the solid (W/m K) — at the relevant
            temperature (carbon steel falls ~50 -> ~30 W/m K from 20 C to
            800 C, which moves dx by the same factor).
        h : float
            Convective heat transfer coefficient at the boundary (W/m^2 K).
        t_surf : float
            Surface temperature (K).
        t_fluid : float
            Free-stream / bulk fluid temperature (K).
        epsilon : float
            Surface emissivity (0-1).
        t_surr : float
            Surrounding radiation sink temperature (K).
        max_dt : float
            Maximum allowable temperature drop across a single element (K).
            A user tolerance, not physics — it sets resolution, not accuracy.
        thickness_m : float
            Wall thickness (m).  Required — the entire defect class this
            method fixes came from sizing without it.
        max_cells : int
            Cap on the reported through-thickness cell count (default 100).

        Returns
        -------
        dict with keys:
            n_cells          : int   — clamped through-thickness cell count
            biot             : float — h_eff * L / k (inf when Ts == Tf)
            regime           : str   — "thermally_thin" | "resolve" |
                                       "steep_gradient" | "near_equilibrium"
            dx_used_mm       : float — thickness / n_cells (the size implied)
            dx_raw_mm        : float — unclamped relation (traceability only)
            dx_exceeds_part  : bool  — raw dx larger than the wall itself
            dt_wall_K        : float — implied through-thickness drop
            h_effective      : float — q'' / (Ts - Tf), incl. radiation
            q_conv, q_rad, q_total, rad_fraction — as in ``max_mesh_size``
            recommendation   : str   — the modelling decision, in words
            warnings         : list[str]
        """
        if thickness_m <= 0:
            raise ValueError("thickness_m must be positive")
        if k <= 0:
            raise ValueError("k must be positive")
        if max_dt <= 0:
            raise ValueError("max_dt must be positive")

        base = BoundaryDrivenConductionCalculator.max_mesh_size(
            k=k, h=h, t_surf=t_surf, t_fluid=t_fluid,
            epsilon=epsilon, t_surr=t_surr, max_dt=max_dt,
        )
        q_conv = base["q_conv"]
        q_rad = base["q_rad"]
        thickness_mm = thickness_m * 1000.0
        warnings: list = []

        # --- Pole guard ------------------------------------------------
        # When convection and radiation oppose (a warm part inside a hotter
        # enclosure, or a cooldown) the net flux crosses zero and the raw dx
        # diverges.  There is no physical singularity there — the part is
        # near equilibrium.  Guard on the net relative to the larger
        # individual term, not on an absolute threshold, so the guard is
        # narrow.
        net = q_conv + q_rad
        scale = max(abs(q_conv), abs(q_rad))
        if scale == 0.0 or abs(net) < NEAR_EQUILIBRIUM_FRACTION * scale:
            warnings.append(
                "Net surface flux is near zero relative to the individual "
                "transfer terms (convection and radiation nearly cancel). "
                "The sizing relation has a pole here that the physics does "
                "not. Falling back to a single cell; size this region "
                "geometrically or from the transient instead."
            )
            return {
                "n_cells": 1,
                "biot": 0.0,
                "regime": "near_equilibrium",
                "dx_used_mm": thickness_mm,
                "dx_raw_mm": float("inf"),
                "dx_exceeds_part": True,
                "dt_wall_K": 0.0,
                "h_effective": 0.0,
                "q_conv": q_conv,
                "q_rad": q_rad,
                "q_total": base["q_total"],
                "rad_fraction": base["rad_fraction"],
                "recommendation": (
                    "NEAR EQUILIBRIUM — convection and radiation cancel. Do "
                    "not size from the steady balance. Use one cell, or size "
                    "from the transient."
                ),
                "warnings": warnings,
            }

        q = abs(net)
        dx_raw_m = k * max_dt / q
        dt_wall = q * thickness_m / k

        # Effective coefficient including radiation. A bookkeeping device for
        # the Biot number, not a film coefficient: it lumps radiation driven
        # by (Ts^4 - Tsurr^4) into a coefficient referenced to Tf, and it
        # scales as Ts^3, so Bi is not a constant of the problem.
        driving = t_surf - t_fluid
        h_eff = q / driving if abs(driving) > 1e-9 else float("inf")
        biot = (
            h_eff * thickness_m / k if math.isfinite(h_eff) else float("inf")
        )

        # --- The clamp --------------------------------------------------
        n_cells = max(1, math.ceil(thickness_m / dx_raw_m))
        if n_cells > max_cells:
            warnings.append(
                f"Requested resolution implies {n_cells} cells through a "
                f"{thickness_mm:.2f} mm wall; capped at {max_cells}. Either "
                f"max_dt={max_dt} K is tighter than the problem needs, or "
                f"the 1-D steady assumption has stopped being appropriate."
            )
            n_cells = max_cells

        if abs(driving) < 1e-9:
            warnings.append(
                "Ts equals Tf, so the effective coefficient and Biot number "
                "are undefined. Reported Biot is infinite by convention."
            )

        if math.isfinite(biot) and biot < THIN_BIOT:
            regime = "thermally_thin"
            recommendation = (
                f"THERMALLY THIN (Bi = {biot:.4f} < {THIN_BIOT}). "
                f"Through-thickness drop is {dt_wall:.2f} K. Use shell "
                f"conduction or a lumped region; do not spend cells "
                f"resolving this."
            )
        elif math.isfinite(biot) and biot < STEEP_BIOT:
            regime = "resolve"
            recommendation = (
                f"RESOLVE (Bi = {biot:.3f}). Through-thickness drop is "
                f"{dt_wall:.2f} K. Use {n_cells} cells."
            )
        else:
            regime = "steep_gradient"
            recommendation = (
                f"STEEP GRADIENT (Bi = {biot:.3f}). Resolve with {n_cells} "
                f"cells, and check that the 1-D steady assumption still "
                f"holds — at this Biot number lateral spreading and "
                f"transient response usually matter."
            )

        if dx_raw_m > thickness_m:
            warnings.append(
                f"Unclamped relation gives dx = {dx_raw_m * 1000.0:.1f} mm "
                f"for a {thickness_mm:.2f} mm wall "
                f"({dx_raw_m / thickness_m:.0f}x the part). That is the "
                f"relation reporting 'thermally thin', not a cell size. "
                f"Clamped to {n_cells} cell(s)."
            )

        # --- Physical consistency check ---------------------------------
        # If the implied through-wall drop exceeds the total surface-to-fluid
        # driving difference, the assumed flux cannot be sustained across
        # this wall. The wall's own resistance dominates, so Ts is NOT
        # independent of the wall and the prescribed-Ts assumption has broken
        # down.
        if abs(driving) > 1e-9 and dt_wall > abs(driving):
            warnings.append(
                f"PHYSICALLY INCONSISTENT: implied through-wall drop "
                f"({dt_wall:.0f} K) exceeds the total surface-to-fluid "
                f"difference ({abs(driving):.0f} K). The wall resistance "
                f"dominates the surface resistance (Bi = {biot:.2f}), so the "
                f"prescribed surface temperature cannot be independent of "
                f"the wall. Solve the wall and its surface condition "
                f"together rather than sizing from an assumed Ts."
            )

        if base["rad_fraction"] > 0.6:
            warnings.append(
                f"Radiation carries {base['rad_fraction'] * 100:.0f}% of the "
                f"surface flux, so the result is dominated by an assumed "
                f"emissivity. Aluminized surfaces run eps ~0.1 fresh and "
                f">0.5 oxidised — a 5x swing in the dominant term. Treat dx "
                f"as a band, not a number."
            )

        return {
            "n_cells": n_cells,
            "biot": biot,
            "regime": regime,
            "dx_used_mm": thickness_mm / n_cells,
            "dx_raw_mm": dx_raw_m * 1000.0,
            "dx_exceeds_part": dx_raw_m > thickness_m,
            "dt_wall_K": dt_wall,
            "h_effective": h_eff,
            "q_conv": q_conv,
            "q_rad": q_rad,
            "q_total": base["q_total"],
            "rad_fraction": base["rad_fraction"],
            "recommendation": recommendation,
            "warnings": warnings,
        }

    @staticmethod
    def lateral_gradient_limit(
        k: float,
        h: float,
        thickness_m: float,
        epsilon: float = 0.0,
        t_surf: float = 300.0,
        fraction: float = 1.0 / 3.0,
    ) -> dict:
        """
        Maximum surface element size to resolve lateral temperature gradients.

        Based on fin theory: a localised heat source on a thin conductive
        plate produces a temperature field that decays exponentially with
        characteristic length 1/m, where:

            m = sqrt(h_total / (k * t))

        and h_total includes the linearised radiation contribution:

            h_total = h_conv + h_rad = h + 4*eps*sigma*T^3

        The surface mesh should capture this decay within ``fraction``
        of the decay length to avoid smearing hot spots.

            dx_lateral <= fraction / m = fraction * sqrt(k * t / h_total)

        Physics from classical extended surface (fin) theory.  Requires
        no new inputs beyond what the conduction calculator already has.

        Parameters
        ----------
        k : float
            Thermal conductivity of the solid (W/m K).
        h : float
            Convective HTC at the boundary (W/m^2 K).
        thickness_m : float
            Wall thickness (m).
        epsilon : float
            Surface emissivity (0–1).  Default 0 ignores radiation.
        t_surf : float
            Surface temperature (K).  Used for linearised radiation h.
        fraction : float
            Fraction of the decay length to use as mesh limit.
            1/3 is conservative (captures ~95% of gradient); 1/4 is
            tighter.  Default 1/3.

        Returns
        -------
        dict with keys:
            max_dx_mm          : float — maximum lateral element size (mm)
            h_total            : float — combined convective + radiative h (W/m^2 K)
            h_rad_linearised   : float — 4*eps*sigma*T^3 (W/m^2 K)
            fin_parameter_m    : float — m = sqrt(h_total / (k*t)) (1/m)
            decay_length_mm    : float — 1/m in mm

        Raises
        ------
        ValueError
            k or t_surf not > 0, h < 0, or epsilon outside [0, 1].
        """
        require_positive("k", k, "W/m K")
        require_non_negative("h", h, "W/m^2 K")
        require_fraction("epsilon", epsilon)
        require_positive("t_surf", t_surf, "K")

        h_rad = 4.0 * epsilon * STEFAN_BOLTZMANN * (t_surf ** 3)
        h_total = h + h_rad

        if h_total <= 0 or thickness_m <= 0 or k <= 0:
            return {
                "max_dx_mm": float("inf"),
                "h_total": h_total,
                "h_rad_linearised": h_rad,
                "fin_parameter_m": 0.0,
                "decay_length_mm": float("inf"),
            }

        m = math.sqrt(h_total / (k * thickness_m))
        decay_length_m = 1.0 / m
        max_dx_m = fraction * decay_length_m

        return {
            "max_dx_mm": max_dx_m * 1000.0,
            "h_total": h_total,
            "h_rad_linearised": h_rad,
            "fin_parameter_m": m,
            "decay_length_mm": decay_length_m * 1000.0,
        }

    @staticmethod
    def transient_penetration_depth(k: float, rho: float, cp: float, dt: float) -> float:
        """
        Maximum element size to resolve a transient thermal wave.

        Parameters
        ----------
        k : float    — conductivity (W/m K)
        rho : float  — density (kg/m^3)
        cp : float   — specific heat (J/kg K)
        dt : float   — time-step size (s)

        Returns
        -------
        float — maximum element size (mm) to resolve the thermal front
        """
        require_positive("k", k, "W/m K")
        require_positive("rho", rho, "kg/m^3")
        require_positive("cp", cp, "J/kg K")
        require_positive("dt", dt, "s")
        alpha = k / (rho * cp)
        return math.sqrt(alpha * dt) * 1000.0
