"""
1-D Wall Mesh-Convergence Check Against CalculiX
================================================

Checks the wall sizing rule against a finite-element solver: a plane wall,
its two faces held at fixed temperatures, meshed at the cell count the
package's sizing rule gives and at two refinements of it, solved by
CalculiX (``ccx``), and compared with the exact steady answer.

The problem
    A plane wall of thickness L, the face x = 0 held at T_hot and the face
    x = L at T_cold, no heat source, conductivity falling linearly with
    temperature:

        k(T) = k_1 + (k_2 - k_1) (T - T_1) / (T_2 - T_1)

    through the two points (T_1, k_1) and (T_2, k_2).

Why the conductivity varies
    With a constant k the steady profile is a straight line, and linear
    elements reproduce a straight line exactly at any mesh (size_wall()'s
    docstring says the same of its own relation): a coarse mesh would pass
    and the check could never fail.  A conductivity that varies with
    temperature bends the profile.  The default case takes carbon steel's
    fall from about 50 W/m K at 20 C to about 30 W/m K at 800 C (the
    figures size_wall()'s docstring gives) across a wall held at 800 C and
    100 C.

The exact answer
    The Kirchhoff transform theta(T) = integral of k(T') dT' from T_1 makes
    the problem linear: theta varies linearly with x between theta(T_hot)
    and theta(T_cold), the flux is q'' = (theta(T_hot) - theta(T_cold)) / L,
    and T(x) is theta's inverse, the root of a quadratic.

The rule's cells
    The sizing rule needs the flux the boundary delivers; here it is q''
    exactly.  The package's sizing hook (size_for_fea()) is given the wall,
    the smallest conductivity across it (the most demanding) and the
    coefficient h = q'' / (T_hot - T_cold) that delivers q'', and returns
    the rule's cell count n.  The check meshes the wall with n, 2n and 4n
    cells.

The comparison
    Each mesh is a row of 8-node bricks (C3D8, one through the width), its
    end faces held at T_hot and T_cold, solved steady with the
    temperature-dependent conductivity, the Newton iteration driven to
    1e-10 of the field (ccx's default criteria stop it about 1 K short).  The mesh's temperature field is the
    piecewise-linear interpolation of the solver's nodal temperatures
    through the thickness; it is compared with the exact profile at sample
    points across the wall, and the largest difference is divided by
    T_hot - T_cold.  Every level must lie within the band (2 %) and the
    error must fall with refinement.  A planted coarse mesh — one cell
    through the wall — misses the band, and the check fails it by name.

Run from the repository root:

    python analysis/mesh_convergence_1d_wall.py              # the study
    python analysis/mesh_convergence_1d_wall.py --cells 1    # one mesh

The study exits 0 when the rule's three levels pass and the planted coarse
mesh fails; ``--cells N`` exits 0 when that mesh lies within the band and 1
when it does not.  Without ``ccx`` on PATH it says so and exits 2.

Standard library only; ``ccx`` runs in a temporary directory that is
removed afterwards.
"""

import argparse
import math
import os
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from thermal_mesh_calculators.fea_sizing import size_for_fea  # noqa: E402

# The default wall: 10 mm of carbon steel, 800 C on one face, 100 C on the
# other; k 50 W/m K at 20 C falling linearly to 30 W/m K at 800 C.
DEFAULT_CASE = {
    "thickness_m": 0.010,
    "t_hot_K": 1073.15,
    "t_cold_K": 373.15,
    "t_1_K": 293.15,
    "k_1": 50.0,
    "t_2_K": 1073.15,
    "k_2": 30.0,
    "rho": 7850.0,       # not used by a steady solve; ccx reads a material whole
    "cp": 460.0,
    "max_dt_K": 20.0,    # the rule's temperature step per cell
    "width_m": 0.001,    # the bar's cross-section (any; the field is 1-D)
}
BAND = 0.02              # the largest profile error allowed, of T_hot - T_cold
SAMPLES = 400            # sample intervals across the wall
CCX_TIMEOUT_S = 60


# ---------------------------------------------------------------------------
#  The exact answer
# ---------------------------------------------------------------------------

def conductivity(t: float, case: dict = DEFAULT_CASE) -> float:
    """k(T), linear through (t_1_K, k_1) and (t_2_K, k_2) (W/m K)."""
    slope = (case["k_2"] - case["k_1"]) / (case["t_2_K"] - case["t_1_K"])
    return case["k_1"] + slope * (t - case["t_1_K"])


def kirchhoff(t: float, case: dict = DEFAULT_CASE) -> float:
    """theta(T) = integral of k from t_1_K to T (W/m)."""
    slope = (case["k_2"] - case["k_1"]) / (case["t_2_K"] - case["t_1_K"])
    u = t - case["t_1_K"]
    return case["k_1"] * u + 0.5 * slope * u * u


def heat_flux(case: dict = DEFAULT_CASE) -> float:
    """The steady flux from the hot face to the cold one (W/m^2)."""
    return ((kirchhoff(case["t_hot_K"], case) - kirchhoff(case["t_cold_K"], case))
            / case["thickness_m"])


def exact_temperature(x: float, case: dict = DEFAULT_CASE) -> float:
    """T(x): theta is linear in x; T is theta's inverse (K)."""
    length = case["thickness_m"]
    th_hot = kirchhoff(case["t_hot_K"], case)
    th_cold = kirchhoff(case["t_cold_K"], case)
    theta = th_hot + (th_cold - th_hot) * x / length
    slope = (case["k_2"] - case["k_1"]) / (case["t_2_K"] - case["t_1_K"])
    k1 = case["k_1"]
    # 0.5 slope u^2 + k1 u - theta = 0, in the form that stays accurate as
    # slope -> 0: u = 2 theta / (k1 + sqrt(k1^2 + 2 slope theta)).
    u = 2.0 * theta / (k1 + math.sqrt(k1 * k1 + 2.0 * slope * theta))
    return case["t_1_K"] + u


# ---------------------------------------------------------------------------
#  The rule's cell count
# ---------------------------------------------------------------------------

def rule_sizing(case: dict = DEFAULT_CASE) -> dict:
    """The sizing hook's answer for the wall: the smallest conductivity
    across it, and the coefficient that delivers the exact flux."""
    k_min = min(conductivity(case["t_hot_K"], case),
                conductivity(case["t_cold_K"], case))
    q = heat_flux(case)
    drop = case["t_hot_K"] - case["t_cold_K"]
    return size_for_fea({
        "k": k_min, "rho": case["rho"], "cp": case["cp"],
        "thickness_mm": case["thickness_m"] * 1000.0,
        "h": q / drop, "t_surf_K": case["t_hot_K"],
        "t_fluid_K": case["t_cold_K"], "t_surr_K": case["t_cold_K"],
        "epsilon": 0.0, "max_dt": case["max_dt_K"], "max_cells": 10000,
    })


# ---------------------------------------------------------------------------
#  The CalculiX deck and its output
# ---------------------------------------------------------------------------

def _real(value: float) -> str:
    """A real for the deck.  ccx reads at most 20 characters of a number
    (a 21-character 1.000000000000000e-03 is read as 1.0), so 13
    significant digits, 18 characters."""
    return f"{value:.12e}"


def write_deck(path: str, n_cells: int, case: dict = DEFAULT_CASE) -> None:
    """A row of n_cells C3D8 bricks through the thickness, the end faces
    held at t_hot_K and t_cold_K, a steady heat-transfer step, and the
    nodal temperatures printed to the .dat file."""
    length, w = case["thickness_m"], case["width_m"]
    lines = ["*HEADING", f"1-D steady wall of {n_cells} cells", "*NODE, NSET=NALL"]
    for i in range(n_cells + 1):
        x = length * i / n_cells
        for j, (y, z) in enumerate(((0.0, 0.0), (w, 0.0), (w, w), (0.0, w))):
            lines.append(f"{4 * i + j + 1}, {_real(x)}, {_real(y)}, {_real(z)}")
    lines.append("*ELEMENT, TYPE=C3D8, ELSET=EALL")
    for i in range(n_cells):
        a, b = 4 * i, 4 * (i + 1)
        # bottom face (z = 0) counter-clockwise from (x_i, 0), then the top
        lines.append(f"{i + 1}, {a + 1}, {b + 1}, {b + 2}, {a + 2}, "
                     f"{a + 4}, {b + 4}, {b + 3}, {a + 3}")
    last = 4 * n_cells
    lines += [
        "*NSET, NSET=HOT", "1, 2, 3, 4",
        "*NSET, NSET=COLD", f"{last + 1}, {last + 2}, {last + 3}, {last + 4}",
        "*MATERIAL, NAME=WALL",
        "*CONDUCTIVITY",
        f"{_real(case['k_1'])}, {_real(case['t_1_K'])}",
        f"{_real(case['k_2'])}, {_real(case['t_2_K'])}",
        "*DENSITY", _real(case["rho"]),
        "*SPECIFIC HEAT", _real(case["cp"]),
        "*SOLID SECTION, ELSET=EALL, MATERIAL=WALL",
        "*INITIAL CONDITIONS, TYPE=TEMPERATURE",
        f"NALL, {_real(0.5 * (case['t_hot_K'] + case['t_cold_K']))}",
        "*STEP",
        "*HEAT TRANSFER, STEADY STATE",
        "1., 1.",
        # The conductivity varies, so the solve is a Newton iteration, and
        # ccx's default criteria stop it with a temperature correction of
        # about 1 K still to come (0.4 % of the field), an error that no
        # refinement removes.  These ask for corrections and residuals of
        # 1e-10 of the field, which Newton reaches in a few more iterations.
        "*CONTROLS, PARAMETERS=FIELD",
        "1.e-10, 1.e-10, 1.e-2, 1.e-15, 0.02, 1.e-5, 1.e-3, 1.e-8",
        "*BOUNDARY",
        f"HOT, 11, 11, {_real(case['t_hot_K'])}",
        f"COLD, 11, 11, {_real(case['t_cold_K'])}",
        "*NODE PRINT, NSET=NALL",
        "NT",
        "*END STEP",
    ]
    with open(path, "w", encoding="ascii", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")


def read_dat_temperatures(text: str) -> dict:
    """Node number -> temperature from a ccx .dat file's NT block."""
    temps = {}
    in_block = False
    for line in text.splitlines():
        if "temperatures" in line:
            in_block = True
            continue
        parts = line.split()
        if in_block and len(parts) == 2:
            try:
                temps[int(parts[0])] = float(parts[1])
            except ValueError:
                in_block = False
    if not temps:
        raise ValueError("no nodal temperatures found in the .dat output")
    return temps


def station_temperatures(temps: dict, n_cells: int) -> list:
    """The temperature at each of the n_cells + 1 stations through the
    thickness: the mean of its four nodes, which a 1-D field makes equal."""
    out = []
    for i in range(n_cells + 1):
        four = [temps[4 * i + j + 1] for j in range(4)]
        if max(four) - min(four) > 1e-6 * max(abs(t) for t in four):
            raise ValueError(f"station {i}: the four nodes disagree: {four}")
        out.append(sum(four) / 4.0)
    return out


def solve(n_cells: int, case: dict = DEFAULT_CASE) -> list:
    """Run ccx on the deck for n_cells; the station temperatures (K)."""
    ccx = shutil.which("ccx")
    if ccx is None:
        raise FileNotFoundError("ccx (CalculiX) is not on PATH")
    work = tempfile.mkdtemp(prefix="wall1d-")
    try:
        write_deck(os.path.join(work, "wall.inp"), n_cells, case)
        env = dict(os.environ, OMP_NUM_THREADS="1")
        run = subprocess.run([ccx, "-i", "wall"], cwd=work, env=env,
                             capture_output=True, text=True,
                             timeout=CCX_TIMEOUT_S)
        dat = os.path.join(work, "wall.dat")
        if run.returncode != 0 or not os.path.isfile(dat):
            raise RuntimeError(f"ccx failed on {n_cells} cells (exit "
                               f"{run.returncode}): {run.stdout[-400:]}")
        with open(dat, encoding="ascii", errors="replace") as fh:
            return station_temperatures(read_dat_temperatures(fh.read()),
                                        n_cells)
    finally:
        shutil.rmtree(work, ignore_errors=True)


# ---------------------------------------------------------------------------
#  The comparison
# ---------------------------------------------------------------------------

def profile_error(stations: list, case: dict = DEFAULT_CASE) -> dict:
    """The mesh's piecewise-linear field against the exact profile.

    Returns max_error (the largest difference at SAMPLES + 1 points across
    the wall and at every station, over T_hot - T_cold) and nodal_error
    (the largest difference at the stations alone, likewise)."""
    n = len(stations) - 1
    length = case["thickness_m"]
    drop = case["t_hot_K"] - case["t_cold_K"]

    def field(x):
        s = min(n - 1, int(x / length * n))
        x0 = length * s / n
        frac = (x - x0) / (length / n)
        return stations[s] + (stations[s + 1] - stations[s]) * frac

    points = {length * i / SAMPLES for i in range(SAMPLES + 1)}
    points |= {length * i / n for i in range(n + 1)}
    worst = max(abs(field(x) - exact_temperature(x, case)) for x in points)
    nodal = max(abs(stations[i] - exact_temperature(length * i / n, case))
                for i in range(n + 1))
    return {"max_error": worst / drop, "nodal_error": nodal / drop}


def check_mesh(name: str, n_cells: int, stations: list,
               case: dict = DEFAULT_CASE, band: float = BAND) -> dict:
    """One mesh against the band: its errors, and a verdict naming it."""
    err = profile_error(stations, case)
    within = err["max_error"] <= band
    verdict = (f"{'PASS' if within else 'FAIL'}: {name} "
               f"({n_cells} cell{'s' if n_cells != 1 else ''}, "
               f"{case['thickness_m'] * 1000.0 / n_cells:.4g} mm): profile "
               f"error {err['max_error'] * 100:.4g} % of T_hot - T_cold "
               f"{'within' if within else 'misses'} the {band * 100:g} % band")
    return {"name": name, "n_cells": n_cells, "within_band": within,
            "verdict": verdict, **err}


def convergence_study(case: dict = DEFAULT_CASE, band: float = BAND,
                      solver=solve) -> dict:
    """The rule's three levels and the planted coarse mesh.

    Returns the rule's sizing, the four checks, the observed orders of
    convergence between the levels, and passed: True when the three levels
    lie within the band with the error falling at each refinement and the
    planted coarse mesh misses it."""
    sized = rule_sizing(case)
    n = sized["n_cells"]
    levels = [("the rule's mesh", n), ("refined x2", 2 * n),
              ("refined x4", 4 * n)]
    checks = [check_mesh(name, cells, solver(cells, case), case, band)
              for name, cells in levels]
    planted = check_mesh("planted coarse mesh", 1, solver(1, case), case, band)
    errors = [c["max_error"] for c in checks]
    orders = [math.log(errors[i] / errors[i + 1], 2)
              for i in range(len(errors) - 1)
              if errors[i] > 0 and errors[i + 1] > 0]
    falling = all(errors[i + 1] < errors[i] for i in range(len(errors) - 1))
    passed = (all(c["within_band"] for c in checks) and falling
              and not planted["within_band"])
    return {"sizing": sized, "checks": checks, "planted": planted,
            "orders": orders, "falling": falling, "passed": passed}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="1-D wall mesh-convergence check against CalculiX")
    parser.add_argument("--cells", type=int, default=None,
                        help="check one mesh of this many cells instead")
    args = parser.parse_args(argv)
    if shutil.which("ccx") is None:
        print("ccx (CalculiX) is not on PATH; nothing was checked.")
        return 2
    case = DEFAULT_CASE
    print(f"Wall: {case['thickness_m'] * 1000:g} mm, {case['t_hot_K']} K -> "
          f"{case['t_cold_K']} K, k {case['k_1']:g} W/m K at {case['t_1_K']} K "
          f"-> {case['k_2']:g} at {case['t_2_K']} K; q'' = "
          f"{heat_flux(case):.6g} W/m^2; band {BAND * 100:g} %")
    if args.cells is not None:
        result = check_mesh(f"{args.cells}-cell mesh", args.cells,
                            solve(args.cells, case), case)
        print(result["verdict"])
        print(f"  nodal error {result['nodal_error'] * 100:.3g} %")
        return 0 if result["within_band"] else 1
    study = convergence_study(case)
    sized = study["sizing"]
    print(f"Rule: {sized['n_cells']} cells of {sized['char_size_mm']:.4g} mm "
          f"({sized['binding_constraint']}; max_dt = {case['max_dt_K']:g} K)")
    for c in study["checks"] + [study["planted"]]:
        print(f"  {c['verdict']}; nodal error {c['nodal_error'] * 100:.3g} %")
    print("  observed order between levels: "
          + ", ".join(f"{p:.2f}" for p in study["orders"]))
    print("CHECK " + ("PASSED: the rule's meshes converge within the band and "
                      "the planted coarse mesh fails it" if study["passed"]
                      else "FAILED"))
    return 0 if study["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
