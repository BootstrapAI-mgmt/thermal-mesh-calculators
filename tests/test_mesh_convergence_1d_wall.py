"""
The 1-D wall mesh-convergence check (analysis/mesh_convergence_1d_wall.py).

The exact answer, the rule's cell count, the deck and the comparison are
tested without a solver, so they run everywhere.  The tests that run
CalculiX skip, by name, where ``ccx`` is not on PATH.
"""

import contextlib
import importlib.util
import io
import math
import shutil
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT = REPO_ROOT / "analysis" / "mesh_convergence_1d_wall.py"


def _load():
    if not SCRIPT.is_file():
        if (REPO_ROOT / ".git").exists():
            pytest.fail("analysis/mesh_convergence_1d_wall.py is missing")
        pytest.skip("analysis/ is not in this tree (the sdist ships none)")
    spec = importlib.util.spec_from_file_location("mesh_convergence_1d_wall",
                                                  SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def wall():
    return _load()


needs_ccx = pytest.mark.skipif(
    shutil.which("ccx") is None,
    reason="CalculiX (ccx) is not on PATH: the solver half of the 1-D wall "
           "mesh-convergence check runs where ccx is installed")


def _exact_stations(wall, n):
    length = wall.DEFAULT_CASE["thickness_m"]
    return [wall.exact_temperature(length * i / n) for i in range(n + 1)]


# --- The exact answer ---------------------------------------------------------

def test_the_exact_profile_meets_both_faces(wall):
    case = wall.DEFAULT_CASE
    assert wall.exact_temperature(0.0) == pytest.approx(case["t_hot_K"],
                                                       rel=1e-14)
    assert wall.exact_temperature(case["thickness_m"]) == pytest.approx(
        case["t_cold_K"], rel=1e-14)


def test_the_flux_by_hand(wall):
    # k(T) = 50 - (20/780)(T - 293.15); theta(T) = 50 u - (10/780) u^2,
    # u = T - 293.15:
    #   theta(1073.15) = 50 * 780 - (10/780) * 780^2 = 39000 - 7800 = 31200
    #   theta(373.15)  = 50 * 80 - (10/780) * 80^2   = 4000 - 82.0512821
    #                  = 3917.9487179
    #   q'' = (31200 - 3917.9487179) / 0.010 m       = 2728205.128 W/m^2
    assert wall.heat_flux() == pytest.approx(2728205.128205128, rel=1e-12)


def test_the_exact_profile_carries_one_flux(wall):
    # Independently of the Kirchhoff inversion: k(T(x)) dT/dx, by central
    # differences of the profile, is -q'' everywhere across the wall.
    q = wall.heat_flux()
    length = wall.DEFAULT_CASE["thickness_m"]
    step = length * 1e-5
    for i in range(1, 10):
        x = length * i / 10
        slope = (wall.exact_temperature(x + step)
                 - wall.exact_temperature(x - step)) / (2 * step)
        assert wall.conductivity(wall.exact_temperature(x)) * slope == \
            pytest.approx(-q, rel=1e-6)


def test_the_profile_bends(wall):
    # Mid-wall the exact temperature lies below the straight line between
    # the faces (723.15 K): k falls with T, so the hot side carries the
    # flux on a steeper gradient.  By hand: theta(L/2) = (31200 +
    # 3917.9487) / 2 = 17558.9744; (10/780) u^2 - 50 u + 17558.9744 = 0,
    # smaller root (50 - sqrt(1599.53978)) / (20/780) = 390.22437511, so
    # T = 683.37437511 K.
    assert wall.exact_temperature(0.005) == pytest.approx(683.3743751,
                                                         rel=1e-9)


# --- The rule's cell count ----------------------------------------------------

def test_the_rule_cell_count_by_hand(wall):
    # The smallest k across the wall is k(1073.15 K) = 30 W/m K; with
    # max_dt = 20 K the rule's cell is 30 * 20 / 2728205.128 m = 0.21993 mm,
    # and 10 mm / 0.21993 mm = 45.47 -> 46 cells of 10/46 = 0.21739 mm.
    sized = wall.rule_sizing()
    assert sized["n_cells"] == 46
    assert sized["char_size_mm"] == pytest.approx(10.0 / 46, rel=1e-12)
    assert sized["binding_constraint"] == "conduction"


# --- The deck and the solver's output -----------------------------------------

def test_the_deck(wall, tmp_path):
    path = tmp_path / "wall.inp"
    wall.write_deck(str(path), 3)
    text = path.read_text(encoding="ascii")
    node_lines = text.split("*NODE, NSET=NALL\n")[1].split("*ELEMENT")[0]
    assert len(node_lines.strip().splitlines()) == 4 * (3 + 1)
    elements = text.split("ELSET=EALL\n")[1].split("*NSET")[0]
    assert elements.strip().splitlines() == [
        "1, 1, 5, 6, 2, 4, 8, 7, 3",
        "2, 5, 9, 10, 6, 8, 12, 11, 7",
        "3, 9, 13, 14, 10, 12, 16, 15, 11",
    ]
    assert "*HEAT TRANSFER, STEADY STATE" in text
    assert "*CONTROLS, PARAMETERS=FIELD" in text
    # ccx reads at most 20 characters of a number (the title under
    # *HEADING is text).
    data = text.split("*NODE, NSET=NALL\n", 1)[1]
    for line in data.splitlines():
        if not line.startswith("*"):
            for field in line.split(","):
                assert len(field.strip()) <= 20, line


def test_the_dat_reader(wall):
    sample = ("\n                        S T E P       1\n\n"
              " temperatures for set NALL and time  0.1000000E+01\n\n"
              "         1  1.073150E+03\n         2  1.073150E+03\n"
              "         3  1.073150E+03\n         4  1.073150E+03\n"
              "         5  6.833743E+02\n         6  6.833743E+02\n"
              "         7  6.833743E+02\n         8  6.833743E+02\n"
              "         9  3.731500E+02\n        10  3.731500E+02\n"
              "        11  3.731500E+02\n        12  3.731500E+02\n")
    temps = wall.read_dat_temperatures(sample)
    assert wall.station_temperatures(temps, 2) == [1073.15, 683.3743, 373.15]
    with pytest.raises(ValueError, match="no nodal temperatures"):
        wall.read_dat_temperatures("nothing here\n")


# --- The comparison, without a solver ------------------------------------------

def test_one_cell_misses_the_band_and_the_rule_mesh_does_not(wall):
    # With the exact nodal values, one cell is the straight line between
    # the faces: 723.15 K mid-wall against 683.37 K, and more where the
    # bend is largest — over 5 % of the 700 K drop.
    coarse = wall.profile_error(_exact_stations(wall, 1))
    assert coarse["max_error"] > 0.05 > wall.BAND
    assert coarse["nodal_error"] == 0.0
    fine = wall.profile_error(_exact_stations(wall, 46))
    assert fine["max_error"] < 1e-4


def test_the_check_names_the_failing_mesh(wall):
    result = wall.check_mesh("planted coarse mesh", 1,
                             _exact_stations(wall, 1))
    assert result["within_band"] is False
    assert result["verdict"].startswith("FAIL: planted coarse mesh (1 cell")
    assert "misses the 2 % band" in result["verdict"]


def test_the_study_logic_with_exact_values(wall):
    study = wall.convergence_study(
        solver=lambda n, case: _exact_stations(wall, n))
    assert study["passed"] is True
    assert [c["n_cells"] for c in study["checks"]] == [46, 92, 184]
    assert study["planted"]["within_band"] is False
    assert all(1.8 < p < 2.2 for p in study["orders"])


def test_the_study_fails_when_the_error_stops_falling(wall):
    # A solver whose answer is off by a fixed amount at every mesh (what an
    # iteration stopped short gives) does not converge: the check fails.
    def stalled(n, case):
        return [t - 0.2 if 0 < i < n else t
                for i, t in enumerate(_exact_stations(wall, n))]
    assert wall.convergence_study(solver=stalled)["passed"] is False


# --- The solver: CalculiX ------------------------------------------------------

@needs_ccx
def test_ccx_the_rule_meshes_converge_and_the_planted_mesh_fails(wall):
    study = wall.convergence_study()
    for check in study["checks"]:
        assert check["within_band"], check["verdict"]
        assert check["nodal_error"] < 1e-5, check["verdict"]
    assert study["falling"]
    assert all(1.8 < p < 2.2 for p in study["orders"]), study["orders"]
    assert study["planted"]["within_band"] is False
    assert study["planted"]["verdict"].startswith("FAIL: planted coarse mesh")
    assert study["passed"] is True


@needs_ccx
def test_ccx_one_mesh_from_the_command_line(wall):
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        assert wall.main(["--cells", "1"]) == 1
        assert wall.main(["--cells", "46"]) == 0
    lines = out.getvalue().splitlines()
    assert any(line.startswith("FAIL: 1-cell mesh") for line in lines)
    assert any(line.startswith("PASS: 46-cell mesh") for line in lines)


@needs_ccx
def test_ccx_matches_the_exact_answer_at_the_nodes(wall):
    stations = wall.solve(4)
    exact = _exact_stations(wall, 4)
    for got, want in zip(stations, exact):
        assert math.isclose(got, want, rel_tol=2e-6), (got, want)
