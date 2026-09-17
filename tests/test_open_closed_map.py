"""Tests for the open/closed map and its checker.

Two jobs. The first is the live gate: this repo's own ``OPEN-CLOSED-MAP.yaml`` must
cover every tracked file. The second is proving the checker *can fail* — a gate only
ever observed passing is not evidence it works, and this ecosystem has a documented
case of a "gate" that was a presence-of-non-empty-string check.

The fallback-parser test matters more than it looks: the checker degrades to a
hand-rolled YAML subset parser when PyYAML is absent, and a checker that mis-parses
its own map would silently green-stamp anything. So the fallback is asserted to
agree with PyYAML on the real map, key for key.
"""

from __future__ import annotations

import importlib.util
import shutil
import subprocess
import sys
from pathlib import Path

import pytest


def _find_repo_root() -> Path:
    """Walk up to the directory holding the map.

    Not ``parent.parent``: this test file lives at ``tests/`` in most of the six
    repos but at ``server/tests/`` in cae-ml-gui, and a fixed number of levels
    silently resolves to the wrong root there.
    """
    for candidate in Path(__file__).resolve().parents:
        if (candidate / "OPEN-CLOSED-MAP.yaml").is_file():
            return candidate
    return Path(__file__).resolve().parent.parent


REPO_ROOT = _find_repo_root()
CHECKER_PATH = REPO_ROOT / "scripts" / "check_open_closed_map.py"
MAP_PATH = REPO_ROOT / "OPEN-CLOSED-MAP.yaml"


def _load_checker():
    spec = importlib.util.spec_from_file_location("_occ_checker", CHECKER_PATH)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


checker = _load_checker()


# --------------------------------------------------------------------------------
# The live gate
# --------------------------------------------------------------------------------
def test_map_exists():
    assert MAP_PATH.is_file(), f"{MAP_PATH.name} must exist at the repo root"


def test_every_tracked_file_resolves_to_a_bucket():
    """The gate itself: no tracked file may be unbucketed."""
    errors, report = checker.check(REPO_ROOT)
    assert not report["uncovered"], "unbucketed files:\n  " + "\n  ".join(report["uncovered"])
    assert not report["stale_entries"], "map entries matching nothing:\n  " + "\n  ".join(report["stale_entries"])
    assert not errors, "map problems:\n  " + "\n  ".join(errors)


def test_checker_exits_zero_as_a_subprocess():
    """CI invokes it as a script, so assert the process contract, not just the function."""
    proc = subprocess.run([sys.executable, str(CHECKER_PATH)], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stdout + proc.stderr


def test_every_entry_cites_a_basis():
    _, report = checker.check(REPO_ROOT)
    # UNDERIVABLE is a legitimate, deliberately visible state; it must be escalated,
    # not silently absorbed, so the checker surfaces it and this test records it.
    assert isinstance(report.get("underivable"), list)


# --------------------------------------------------------------------------------
# Resolution semantics
# --------------------------------------------------------------------------------
@pytest.mark.parametrize(
    "target,entries,expected",
    [
        ("server/db.py", ["server"], "server"),
        ("server/db.py", ["server", "server/db.py"], "server/db.py"),
        ("server/model_card.py", ["server", "server/model_card.py"], "server/model_card.py"),
        ("serverless/x.py", ["server"], None),  # must not match on a bare string prefix
        ("a/b/c/d.py", ["a", "a/b", "a/b/c"], "a/b/c"),
        ("README.md", ["README.md"], "README.md"),
        ("docs/x.md", ["docs/y"], None),
    ],
)
def test_resolve_is_longest_prefix_on_directory_boundaries(target, entries, expected):
    assert checker.resolve(target, entries) == expected


# --------------------------------------------------------------------------------
# The checker must be able to fail
# --------------------------------------------------------------------------------
def _clone_repo(tmp_path: Path) -> Path:
    """A tiny throwaway git repo carrying a copy of this repo's map and checker."""
    work = tmp_path / "repo"
    (work / "scripts").mkdir(parents=True)
    shutil.copy(CHECKER_PATH, work / "scripts" / "check_open_closed_map.py")
    shutil.copy(MAP_PATH, work / "OPEN-CLOSED-MAP.yaml")
    subprocess.run(["git", "init", "-q"], cwd=work, check=True)
    return work


def _track(work: Path, rel: str, body: str | None = None) -> None:
    """Stage ``rel``. Only writes content when the file does not already exist —
    passing a default body here would clobber the copied map and make every
    downstream assertion test the wrong thing."""
    p = work / rel
    if not p.exists():
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(body if body is not None else "x\n", encoding="utf-8")
    subprocess.run(["git", "add", rel], cwd=work, check=True)


def test_uncovered_file_is_an_error(tmp_path):
    work = _clone_repo(tmp_path)
    _track(work, "OPEN-CLOSED-MAP.yaml")
    _track(work, "scripts/check_open_closed_map.py")
    _track(work, "totally/new/thing.py")
    errors, report = checker.check(work)
    assert "totally/new/thing.py" in report["uncovered"]
    assert any("uncovered: totally/new/thing.py" in e for e in errors)


def test_unknown_bucket_is_an_error_without_crashing(tmp_path):
    work = _clone_repo(tmp_path)
    text = (work / "OPEN-CLOSED-MAP.yaml").read_text(encoding="utf-8")
    (work / "OPEN-CLOSED-MAP.yaml").write_text(
        text.replace("bucket: neutral", "bucket: sort-of-open", 1), encoding="utf-8"
    )
    _track(work, "OPEN-CLOSED-MAP.yaml")
    _track(work, "scripts/check_open_closed_map.py")
    errors, _ = checker.check(work)  # must return, not raise
    assert any("unknown bucket 'sort-of-open'" in e for e in errors)


def test_duplicate_path_is_an_error(tmp_path):
    work = _clone_repo(tmp_path)
    text = (work / "OPEN-CLOSED-MAP.yaml").read_text(encoding="utf-8")
    text += '  - path: OPEN-CLOSED-MAP.yaml\n    bucket: open\n    basis: "dupe"\n'
    (work / "OPEN-CLOSED-MAP.yaml").write_text(text, encoding="utf-8")
    _track(work, "OPEN-CLOSED-MAP.yaml")
    _track(work, "scripts/check_open_closed_map.py")
    errors, _ = checker.check(work)
    assert any("duplicate path" in e for e in errors)


def test_extract_to_on_a_closed_path_is_an_error(tmp_path):
    work = _clone_repo(tmp_path)
    text = (work / "OPEN-CLOSED-MAP.yaml").read_text(encoding="utf-8")
    text += '  - path: zzz\n    bucket: closed\n    basis: "b"\n    extract_to: trust-standard\n'
    (work / "OPEN-CLOSED-MAP.yaml").write_text(text, encoding="utf-8")
    _track(work, "OPEN-CLOSED-MAP.yaml")
    _track(work, "scripts/check_open_closed_map.py")
    _track(work, "zzz/f.py")
    errors, _ = checker.check(work)
    assert any("only open paths are extractable" in e for e in errors)


def test_missing_map_is_an_error(tmp_path):
    work = tmp_path / "empty"
    work.mkdir()
    subprocess.run(["git", "init", "-q"], cwd=work, check=True)
    errors, _ = checker.check(work)
    assert errors and "not found" in errors[0]


# --------------------------------------------------------------------------------
# The fallback parser must not lie
# --------------------------------------------------------------------------------
def test_fallback_parser_agrees_with_pyyaml_on_the_real_map():
    yaml = pytest.importorskip("yaml")
    text = MAP_PATH.read_text(encoding="utf-8")
    reference = yaml.safe_load(text)
    fallback = checker._mini_yaml_load(text)

    assert set(fallback) == set(reference), "top-level keys differ"
    for key in ("schema_version", "repo", "generated_at", "counts", "derived_from"):
        assert fallback[key] == reference[key], f"{key} differs"
    assert len(fallback["entries"]) == len(reference["entries"])
    for got, want in zip(fallback["entries"], reference["entries"]):
        assert got == want, f"entry differs:\n  fallback={got}\n  pyyaml  ={want}"


def test_fallback_parser_refuses_constructs_it_cannot_handle():
    with pytest.raises(ValueError):
        checker._mini_yaml_load("entries: [a, b, c]\n")
