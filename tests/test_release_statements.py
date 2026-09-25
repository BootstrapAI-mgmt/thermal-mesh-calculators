"""Release statements in the docs must name releases that exist.

README.md, CHANGELOG.md and .github/copilot-instructions.md tell a reader which
version to install and which tag to pin. Those statements drifted from what was
published: the README's git pin named a tag this repository does not carry, so
its documented install command failed, and CHANGELOG.md linked a release page
that answers 404.

The checks hold each statement to two sources of truth kept in the tree, so they
run offline, on the standard library alone:

* ``__version__`` in ``thermal_mesh_calculators/__init__.py`` -- no statement may
  name a version newer than the one the package declares;
* CHANGELOG.md -- its ``## [x.y.z]`` headings are the released versions, and its
  Keep a Changelog link definitions, ``[x.y.z]: <repository>/releases/tag/vx.y.z``,
  are the release tags this repository carries. A test cannot see the remote, so
  those definitions stand in for it: add one when a tag is pushed, not before.

The test counts that CLAUDE.md and .github/copilot-instructions.md state (in
total, per part and, in CLAUDE.md's tree, per file) are held to this tree's own
suite, collected in a subprocess by ``pytest --collect-only -qq``, which prints
one "path: count" line per test file and runs no test. One of those documents
was left stating an earlier size of the suite when tests were added.

The second half plants each defect in synthetic text and asserts that it is
reported: a check only ever observed passing is not evidence that it works.
"""

from __future__ import annotations

import ast
import functools
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest


def _find_repo_root() -> Path:
    """Walk up to the directory that holds CHANGELOG.md beside the package."""
    here = Path(__file__).resolve()
    for candidate in here.parents:
        if (candidate / "CHANGELOG.md").is_file() and (candidate / "thermal_mesh_calculators").is_dir():
            return candidate
    return here.parent.parent


REPO_ROOT = _find_repo_root()

# The documents whose git pins are checked. The sdist prunes .github/ (MANIFEST.in), so a
# packager running the sdist's own suite sees that one reported as a skip, not an error.
DOCUMENTS = ("README.md", "CHANGELOG.md", ".github/copilot-instructions.md")

_PLAIN_VERSION = re.compile(r"(\d+)\.(\d+)\.(\d+)")
_HEADING = re.compile(r"^## \[([^\]]+)\]", re.MULTILINE)
_LINK_DEFINITION = re.compile(r"^\[([^\]]+)\]:[ \t]*(\S+)", re.MULTILINE)
# A version token runs to the first character a version cannot contain, so a pin such as
# ``@v0.7.0rc1`` is read whole and rejected, never silently read as ``0.7.0``.
_TOKEN = r"(\d[0-9A-Za-z.+-]*)"
_GIT_REF = re.compile("@v" + _TOKEN)
_TAG_URL = re.compile("/releases/tag/v" + _TOKEN)
_ON_PYPI_SINCE = re.compile("On PyPI since v?" + _TOKEN)
_REPOSITORY_URL = re.compile(r'^Repository\s*=\s*"([^"]+)"', re.MULTILINE)


# --------------------------------------------------------------------------------
# The checks
# --------------------------------------------------------------------------------
def parse_version(token: str) -> tuple[int, ...] | None:
    """``"0.6.2"`` -> ``(0, 6, 2)``; anything but a plain x.y.z -> ``None``."""
    match = _PLAIN_VERSION.fullmatch(token)
    return tuple(int(part) for part in match.groups()) if match else None


def _dotted(version: tuple[int, ...]) -> str:
    return ".".join(str(part) for part in version)


def _find_tokens(pattern: re.Pattern, text: str) -> list[tuple[int, str]]:
    """``(line number, token)`` per match; a sentence-ending full stop is not part of it."""
    return [
        (lineno, match.group(1).rstrip("."))
        for lineno, line in enumerate(text.splitlines(), start=1)
        for match in pattern.finditer(line)
    ]


def released_versions(changelog: str) -> set[tuple[int, ...]]:
    """Versions with a ``## [x.y.z]`` heading; ``## [Unreleased]`` is not a release."""
    parsed = (parse_version(label) for label in _HEADING.findall(changelog))
    return {version for version in parsed if version is not None}


def version_problems(
    where: str, token: str, declared: tuple[int, ...], released: set[tuple[int, ...]]
) -> list[str]:
    """The rule for a named version: released, and no newer than ``__version__``."""
    version = parse_version(token)
    if version is None:
        return [f"{where} -- {token!r} is not a plain x.y.z release version"]
    problems = []
    if version > declared:
        problems.append(f"{where} -- newer than __version__ {_dotted(declared)}")
    if version not in released:
        problems.append(f"{where} -- CHANGELOG.md has no '## [{token}]' heading")
    return problems


def read_tag_links(
    changelog: str, repository: str, declared: tuple[int, ...]
) -> tuple[set[tuple[int, ...]], list[str]]:
    """CHANGELOG.md's release-tag link definitions: the versions they declare tagged, and
    every problem found. A malformed link declares nothing."""
    released = released_versions(changelog)
    prefix = repository.rstrip("/") + "/releases/tag/v"
    tagged: set[tuple[int, ...]] = set()
    problems: list[str] = []
    for label, url in _LINK_DEFINITION.findall(changelog):
        if not _TAG_URL.search(url):
            continue  # not a tag link (a compare link, say)
        where = f"CHANGELOG.md link [{label}]: {url}"
        found = []
        if url != prefix + label:
            found.append(f"{where} -- expected {prefix}{label}")
        found += version_problems(where, label, declared, released)
        if found:
            problems += found
        else:
            tagged.add(parse_version(label))
    return tagged, problems


def git_ref_problems(
    name: str, text: str, changelog: str, repository: str, declared: tuple[int, ...]
) -> list[str]:
    """Every ``@vX.Y.Z`` git ref names a released version no newer than ``__version__``,
    and one whose release tag CHANGELOG.md links -- a tag this repository carries."""
    released = released_versions(changelog)
    tagged, _ = read_tag_links(changelog, repository, declared)
    problems = []
    for lineno, token in _find_tokens(_GIT_REF, text):
        where = f"{name}:{lineno}: git ref v{token}"
        found = version_problems(where, token, declared, released)
        if not found and parse_version(token) not in tagged:
            found.append(
                f"{where} -- CHANGELOG.md links no v{token} release tag; "
                "pin a tag this repository carries"
            )
        problems += found
    return problems


def on_pypi_since_problems(readme: str, changelog: str, declared: tuple[int, ...]) -> list[str]:
    """README's "On PyPI since" names a released version no newer than ``__version__``."""
    found = _find_tokens(_ON_PYPI_SINCE, readme)
    if not found:
        return [
            "README.md has no 'On PyPI since <version>' statement; if the wording changed, "
            "change _ON_PYPI_SINCE with it rather than let this check pass on nothing"
        ]
    released = released_versions(changelog)
    problems = []
    for lineno, token in found:
        where = f"README.md:{lineno}: On PyPI since {token}"
        problems += version_problems(where, token, declared, released)
    return problems


# The documents that state test counts, and the two test files that are governance
# rather than physics.
COUNT_DOCUMENTS = ("CLAUDE.md", ".github/copilot-instructions.md")
_MAP_FILE = "test_open_closed_map.py"
_RELEASE_FILE = "test_release_statements.py"

_ACROSS = re.compile(r"(\d+) tests across (\d+) (?:test )?files")
_TOTAL = re.compile(r"(\d+) tests total|pytest tests/ -q\s+# (\d+) passed")
_PHYSICS = re.compile(r"(\d+)(?: physics| pytest unit tests)")
_MAP_PART = re.compile(r"(\d+)\** for the open/closed map checker")
_RELEASE_PART = re.compile(r"(\d+)\** for the release statements")
_PER_FILE = re.compile(r"^\s*(test_\w+\.py)\s+# (\d+) tests\b")


def count_problems(name: str, text: str, counts: dict[str, int]) -> list[str]:
    """Every test count ``text`` states that differs from ``counts``, the number of
    tests collected per test file (keyed by file name). A document that states no
    count at all is reported too, so the check cannot pass on nothing."""
    total = sum(counts.values())
    physics = total - counts.get(_MAP_FILE, 0) - counts.get(_RELEASE_FILE, 0)
    problems: list[str] = []
    stated = 0
    for lineno, line in enumerate(text.splitlines(), start=1):
        found = []
        for match in _ACROSS.finditer(line):
            found.append(("tests", match.group(1), total))
            found.append(("test files", match.group(2), len(counts)))
        found += [("tests", a or b, total) for a, b in _TOTAL.findall(line)]
        found += [("physics and input-validation tests", n, physics)
                  for n in _PHYSICS.findall(line)]
        found += [(f"tests in {_MAP_FILE}", n, counts.get(_MAP_FILE, 0))
                  for n in _MAP_PART.findall(line)]
        found += [(f"tests in {_RELEASE_FILE}", n, counts.get(_RELEASE_FILE, 0))
                  for n in _RELEASE_PART.findall(line)]
        per_file = _PER_FILE.match(line)
        if per_file:
            found.append((f"tests in {per_file.group(1)}", per_file.group(2),
                          counts.get(per_file.group(1), 0)))
        for what, said, actual in found:
            stated += 1
            if int(said) != actual:
                problems.append(f"{name}:{lineno}: says {said} {what}; the suite has {actual}")
    if not stated:
        problems.append(
            f"{name} states no test count; if the wording changed, change the patterns "
            "with it rather than let this check pass on nothing"
        )
    return problems


# --------------------------------------------------------------------------------
# The live gate
# --------------------------------------------------------------------------------
def _read(relpath: str) -> str:
    return (REPO_ROOT / relpath).read_text(encoding="utf-8")


def _declared_version() -> tuple[int, ...]:
    """``__version__`` read statically from this tree, as pyproject.toml has setuptools do."""
    module = ast.parse(_read("thermal_mesh_calculators/__init__.py"))
    for node in module.body:
        if isinstance(node, ast.Assign) and any(
            isinstance(target, ast.Name) and target.id == "__version__" for target in node.targets
        ):
            value = ast.literal_eval(node.value)
            version = parse_version(value) if isinstance(value, str) else None
            assert version is not None, f"__version__ = {value!r} is not a plain x.y.z version"
            return version
    raise AssertionError("thermal_mesh_calculators/__init__.py assigns no __version__")


def _repository() -> str:
    """This repository's URL, from pyproject.toml's ``[project.urls] Repository``."""
    match = _REPOSITORY_URL.search(_read("pyproject.toml"))
    assert match, "pyproject.toml declares no [project.urls] Repository"
    return match.group(1)


@pytest.mark.parametrize("relpath", DOCUMENTS)
def test_every_git_ref_names_a_released_and_tagged_version(relpath):
    if not (REPO_ROOT / relpath).is_file():
        if (REPO_ROOT / ".git").exists():
            pytest.fail(f"{relpath} is missing from this checkout; update DOCUMENTS if it moved")
        pytest.skip(f"{relpath} is not in this tree (the sdist prunes .github/, see MANIFEST.in)")
    problems = git_ref_problems(
        relpath, _read(relpath), _read("CHANGELOG.md"), _repository(), _declared_version()
    )
    assert not problems, "\n".join(problems)


def test_changelog_tag_links_are_well_formed():
    _, problems = read_tag_links(_read("CHANGELOG.md"), _repository(), _declared_version())
    assert not problems, "\n".join(problems)


def test_readme_on_pypi_since_names_a_released_version():
    problems = on_pypi_since_problems(_read("README.md"), _read("CHANGELOG.md"), _declared_version())
    assert not problems, "\n".join(problems)


_COLLECTED = re.compile(r"^(.+\.py): (\d+)$", re.MULTILINE)


@functools.lru_cache(maxsize=None)
def _collected_counts() -> dict[str, int]:
    """Tests per test file in this tree's suite, from ``pytest --collect-only -qq`` run
    in a subprocess (collecting runs no test)."""
    env = {key: value for key, value in os.environ.items() if key != "PYTEST_ADDOPTS"}
    run = subprocess.run(
        [sys.executable, "-m", "pytest", "--collect-only", "-qq", "-p", "no:cacheprovider", "tests"],
        cwd=REPO_ROOT, env=env, capture_output=True, text=True, timeout=300,
    )
    counts = {
        path.replace("\\", "/").rsplit("/", 1)[-1]: int(count)
        for path, count in _COLLECTED.findall(run.stdout)
    }
    assert run.returncode == 0 and counts, run.stdout + run.stderr
    return counts


@pytest.mark.parametrize("relpath", COUNT_DOCUMENTS)
def test_documented_test_counts_match_the_suite(relpath):
    if not (REPO_ROOT / relpath).is_file():
        if (REPO_ROOT / ".git").exists():
            pytest.fail(f"{relpath} is missing from this checkout; update COUNT_DOCUMENTS if it moved")
        pytest.skip(f"{relpath} is not in this tree (the sdist ships neither CLAUDE.md nor .github/)")
    problems = count_problems(relpath, _read(relpath), _collected_counts())
    assert not problems, "\n".join(problems)


# --------------------------------------------------------------------------------
# Proving each check can fail
# --------------------------------------------------------------------------------
_REPO = "https://github.com/example-org/example-repo"
_TAG_LINK_062 = "[0.6.2]: " + _REPO + "/releases/tag/v0.6.2"
_CHANGELOG = (
    "# Changelog\n\n## [Unreleased]\n\n## [0.6.2] - 2026-09-17\n\n"
    "## [0.6.1] - 2026-09-16\n\n## [0.6.0] - 2026-08-12\n\n" + _TAG_LINK_062 + "\n"
)
_DECLARED = (0, 6, 2)


def _pin(token: str) -> str:
    return 'pip install "example @ git+' + _REPO + ".git@v" + token + '"\n'


def test_consistent_statements_report_nothing():
    """Guards the tests below: each must fail for its planted reason, not for every input."""
    assert git_ref_problems("README.md", _pin("0.6.2"), _CHANGELOG, _REPO, _DECLARED) == []
    assert read_tag_links(_CHANGELOG, _REPO, _DECLARED) == ({(0, 6, 2)}, [])
    assert on_pypi_since_problems("On PyPI since 0.6.1.\n", _CHANGELOG, _DECLARED) == []


@pytest.mark.parametrize(
    "token, expected",
    [
        ("0.6.1", "links no v0.6.1 release tag"),  # released but never tagged here
        ("0.5.9", "no '## [0.5.9]' heading"),
        ("9.9.9", "newer than __version__ 0.6.2"),
        ("0.7.0rc1", "is not a plain x.y.z release version"),
    ],
)
def test_a_bad_git_ref_is_reported(token, expected):
    problems = git_ref_problems("README.md", _pin(token), _CHANGELOG, _REPO, _DECLARED)
    assert any(expected in problem for problem in problems), problems


@pytest.mark.parametrize(
    "link, expected",
    [
        (
            "[0.6.2]: https://github.com/example-org/elsewhere/releases/tag/v0.6.2",
            "expected " + _REPO + "/releases/tag/v0.6.2",
        ),
        ("[0.6.2]: " + _REPO + "/releases/tag/v0.6.1", "expected " + _REPO + "/releases/tag/v0.6.2"),
        ("[0.5.9]: " + _REPO + "/releases/tag/v0.5.9", "no '## [0.5.9]' heading"),
    ],
)
def test_a_bad_changelog_tag_link_is_reported_and_declares_nothing(link, expected):
    changelog = _CHANGELOG.replace(_TAG_LINK_062, link)
    assert link in changelog
    tagged, problems = read_tag_links(changelog, _REPO, _DECLARED)
    assert any(expected in problem for problem in problems), problems
    assert tagged == set()


@pytest.mark.parametrize(
    "readme, expected",
    [
        ("On PyPI since 0.5.0.\n", "no '## [0.5.0]' heading"),
        ("Install it from PyPI.\n", "no 'On PyPI since <version>' statement"),
    ],
)
def test_a_bad_on_pypi_since_statement_is_reported(readme, expected):
    problems = on_pypi_since_problems(readme, _CHANGELOG, _DECLARED)
    assert any(expected in problem for problem in problems), problems


_SUITE = {"test_a.py": 3, _MAP_FILE: 2, _RELEASE_FILE: 1}  # 6 tests, 3 files, 3 physics
_COUNTS_TEXT = (
    "    test_a.py      # 3 tests\n"
    "                   # 6 tests total\n"
    "6 tests across 3 test files (3 physics + 2 for the open/closed map checker"
    " + 1 for the release statements).\n"
    "python -m pytest tests/ -q   # 6 passed\n"
)


def test_consistent_counts_report_nothing():
    """Guards the test below: each case must fail for its planted count alone."""
    assert count_problems("DOC.md", _COUNTS_TEXT, _SUITE) == []


@pytest.mark.parametrize(
    "fresh, stale, expected",
    [
        ("# 6 tests total", "# 7 tests total", "says 7 tests; the suite has 6"),
        ("across 3 test files", "across 4 test files", "says 4 test files; the suite has 3"),
        ("(3 physics", "(5 physics", "says 5 physics and input-validation tests; the suite has 3"),
        ("2 for the open", "4 for the open", f"says 4 tests in {_MAP_FILE}; the suite has 2"),
        ("1 for the release", "4 for the release", f"says 4 tests in {_RELEASE_FILE}; the suite has 1"),
        ("test_a.py      # 3 tests", "test_a.py      # 9 tests", "says 9 tests in test_a.py; the suite has 3"),
        ("# 6 passed", "# 254 passed", "says 254 tests; the suite has 6"),
    ],
)
def test_a_stale_count_is_reported(fresh, stale, expected):
    text = _COUNTS_TEXT.replace(fresh, stale)
    assert stale in text
    problems = count_problems("DOC.md", text, _SUITE)
    assert any(expected in problem for problem in problems), problems


def test_a_document_that_states_no_count_is_reported():
    problems = count_problems("DOC.md", "The suite is large.\n", _SUITE)
    assert any("states no test count" in problem for problem in problems), problems
