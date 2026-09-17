#!/usr/bin/env python3
"""Falsifier for ``OPEN-CLOSED-MAP.yaml``: every tracked file must resolve to a bucket.

The map records, per path, which side of the open/closed line the code sits on and
which phrase of ``MASTER-PRODUCTIZATION-STRATEGY.md`` §4.3 put it there. That map is
only worth anything if it cannot silently go stale, so this script is the gate:

* **Every tracked file must resolve.** A file that no map entry covers is an error.
  This is the point of the whole exercise — a new file cannot land unbucketed.
* **Every map entry must cover something.** An entry whose path matches no tracked
  file is an error too (the reverse falsifier: it catches a map that still describes
  files somebody deleted).
* **Structure is checked, not assumed.** Unknown bucket values, missing required
  keys, duplicate paths, bad ``extract_to`` targets, and a header ``counts`` block
  that disagrees with the entries are all errors.

Resolution is **longest-prefix wins**, directory-boundary aware. ``server/`` covers
``server/db.py``; an explicit ``server/model_card.py`` entry beats it. That rule is
mechanical: two people resolving the same file always get the same answer.

``counts`` in the header counts *entries* per bucket, not files. It is a property of
the map file itself, so it only changes when the map changes — adding a source file
under an already-covered directory does not force a bookkeeping edit. The per-bucket
*file* tallies are computed fresh and printed in the summary.

Usage::

    python scripts/check_open_closed_map.py            # check this repo
    python scripts/check_open_closed_map.py --json     # machine-readable
    python scripts/check_open_closed_map.py --repo ../other

Exit code 0 when the map is complete and well-formed, 1 otherwise.

This file is duplicated verbatim into each of the six ecosystem repos on purpose:
each repo has to be independently checkable, including after the open core is split
out into its own repository.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any

MAP_FILENAME = "OPEN-CLOSED-MAP.yaml"

VALID_BUCKETS = frozenset({"open", "closed", "neutral"})
VALID_EXTRACT_TARGETS = frozenset({"trust-standard"})

REQUIRED_HEADER_KEYS = ("schema_version", "repo", "derived_from", "generated_at", "counts")
REQUIRED_ENTRY_KEYS = ("path", "bucket", "basis")
OPTIONAL_ENTRY_KEYS = ("extract_to", "notes")

SUPPORTED_SCHEMA_VERSION = 1


# --------------------------------------------------------------------------------------
# YAML loading
# --------------------------------------------------------------------------------------
def _load_yaml(text: str) -> Any:
    """Parse the map.

    Prefers PyYAML. Falls back to :func:`_mini_yaml_load`, a deliberately narrow
    parser for the exact grammar the map is written in, so the gate still runs in a
    repo that does not carry a PyYAML dependency (``thermal-mesh-calculators`` is
    dependency-free on purpose). The fallback raises on anything it does not
    recognise rather than guessing — a checker that mis-parses is worse than one
    that refuses to run.
    """
    try:
        import yaml  # type: ignore
    except ImportError:
        return _mini_yaml_load(text)
    return yaml.safe_load(text)


def _scalar(raw: str, lineno: int) -> Any:
    raw = raw.strip()
    if raw in ("null", "~", ""):
        return None
    if raw in ("true", "false"):
        return raw == "true"
    if len(raw) >= 2 and raw[0] == raw[-1] and raw[0] in "\"'":
        body = raw[1:-1]
        if raw[0] == '"':
            body = body.replace('\\"', '"').replace("\\\\", "\\")
        else:
            body = body.replace("''", "'")
        return body
    try:
        return int(raw)
    except ValueError:
        pass
    if any(ch in raw for ch in "{}[]#&*!|>%@`"):
        raise ValueError(f"line {lineno}: unsupported YAML construct in fallback parser: {raw!r}")
    return raw


def _split_kv(line: str, lineno: int) -> tuple[str, str]:
    """Split ``key: value`` outside of quotes."""
    in_q: str | None = None
    i = 0
    while i < len(line):
        ch = line[i]
        if in_q:
            if ch == "\\" and in_q == '"':
                i += 2
                continue
            if ch == in_q:
                in_q = None
        elif ch in "\"'":
            in_q = ch
        elif ch == ":" and (i + 1 == len(line) or line[i + 1] in " \t"):
            return line[:i].strip(), line[i + 1 :].strip()
        i += 1
    raise ValueError(f"line {lineno}: expected 'key: value', got {line!r}")


def _mini_yaml_load(text: str) -> Any:
    """Strict parser for the restricted grammar the map uses.

    Supports exactly: top-level ``key: scalar``, one level of nested mapping, and
    lists of flat mappings (``- key: scalar`` then aligned ``key: scalar``).
    Anything else raises.
    """
    root: dict[str, Any] = {}
    stack: list[tuple[int, Any]] = [(-1, root)]
    cur_list: list[dict[str, Any]] | None = None
    cur_item: dict[str, Any] | None = None
    list_indent = -1

    for lineno, raw_line in enumerate(text.splitlines(), start=1):
        if not raw_line.strip() or raw_line.lstrip().startswith("#"):
            continue
        line = raw_line.split(" #")[0].rstrip() if " #" in raw_line else raw_line.rstrip()
        indent = len(line) - len(line.lstrip())
        body = line.strip()

        if body.startswith("- "):
            if cur_list is None or indent != list_indent:
                raise ValueError(f"line {lineno}: list item outside a known list context")
            cur_item = {}
            cur_list.append(cur_item)
            k, v = _split_kv(body[2:], lineno)
            cur_item[k] = _scalar(v, lineno)
            continue

        if cur_item is not None and indent == list_indent + 2:
            k, v = _split_kv(body, lineno)
            cur_item[k] = _scalar(v, lineno)
            continue

        cur_list, cur_item, list_indent = None, None, -1
        while stack and indent <= stack[-1][0]:
            stack.pop()
        if not stack:
            raise ValueError(f"line {lineno}: bad indentation")
        parent = stack[-1][1]
        k, v = _split_kv(body, lineno)
        if v == "":
            nxt = _peek_next(text, lineno)
            if nxt is not None and nxt.strip().startswith("- "):
                new_list: list[dict[str, Any]] = []
                parent[k] = new_list
                cur_list = new_list
                list_indent = len(nxt) - len(nxt.lstrip())
            else:
                child: dict[str, Any] = {}
                parent[k] = child
                stack.append((indent, child))
        else:
            parent[k] = _scalar(v, lineno)
    return root


def _peek_next(text: str, after_lineno: int) -> str | None:
    for line in text.splitlines()[after_lineno:]:
        if line.strip() and not line.lstrip().startswith("#"):
            return line
    return None


# --------------------------------------------------------------------------------------
# Resolution
# --------------------------------------------------------------------------------------
def _norm(path: str) -> str:
    """Normalise to a repo-relative POSIX path with no leading ``./`` and no trailing ``/``.

    Strips the ``./`` *prefix*, not the characters ``.`` and ``/`` — ``str.lstrip("./")``
    would turn ``.gitignore`` into ``gitignore`` and collide with a real ``gitignore/``
    directory.
    """
    p = path.replace("\\", "/").strip()
    while p.startswith("./"):
        p = p[2:]
    return p.rstrip("/")


def resolve(file_path: str, entry_paths: Sequence[str]) -> str | None:
    """Longest-prefix, directory-boundary-aware match. Returns the winning entry path."""
    target = _norm(file_path)
    best: str | None = None
    for ep in entry_paths:
        if target == ep or target.startswith(ep + "/"):
            if best is None or len(ep) > len(best):
                best = ep
    return best


def tracked_files(repo: Path) -> list[str]:
    out = subprocess.run(
        ["git", "-C", str(repo), "ls-files"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    return [_norm(line) for line in out.splitlines() if line.strip()]


# --------------------------------------------------------------------------------------
# Checking
# --------------------------------------------------------------------------------------
def check(repo: Path) -> tuple[list[str], dict[str, Any]]:
    """Return ``(errors, report)``. Empty ``errors`` means the map is sound."""
    errors: list[str] = []
    report: dict[str, Any] = {"repo": str(repo), "map": MAP_FILENAME}

    map_path = repo / MAP_FILENAME
    if not map_path.is_file():
        return [f"{MAP_FILENAME} not found at {map_path}"], report

    try:
        data = _load_yaml(map_path.read_text(encoding="utf-8"))
    except Exception as exc:  # noqa: BLE001 - any parse failure is a gate failure
        return [f"{MAP_FILENAME} failed to parse: {exc}"], report

    if not isinstance(data, dict):
        return [f"{MAP_FILENAME} top level must be a mapping, got {type(data).__name__}"], report

    for key in REQUIRED_HEADER_KEYS:
        if key not in data:
            errors.append(f"header: missing required key {key!r}")
    if data.get("schema_version") != SUPPORTED_SCHEMA_VERSION:
        errors.append(f"header: schema_version must be {SUPPORTED_SCHEMA_VERSION}, got {data.get('schema_version')!r}")
    report["repo_declared"] = data.get("repo")

    entries = data.get("entries")
    if not isinstance(entries, list) or not entries:
        errors.append("entries: must be a non-empty list")
        return errors, report

    # --- per-entry structural validation -------------------------------------------
    seen: dict[str, int] = {}
    clean: list[dict[str, Any]] = []
    for idx, entry in enumerate(entries):
        where = f"entries[{idx}]"
        if not isinstance(entry, dict):
            errors.append(f"{where}: must be a mapping, got {type(entry).__name__}")
            continue
        for key in REQUIRED_ENTRY_KEYS:
            if key not in entry or entry[key] in (None, ""):
                errors.append(f"{where}: missing or empty required key {key!r}")
        unknown = set(entry) - set(REQUIRED_ENTRY_KEYS) - set(OPTIONAL_ENTRY_KEYS)
        if unknown:
            errors.append(f"{where}: unknown key(s) {sorted(unknown)}")

        bucket = entry.get("bucket")
        if bucket not in VALID_BUCKETS:
            errors.append(
                f"{where} ({entry.get('path')!r}): unknown bucket {bucket!r}; allowed {sorted(VALID_BUCKETS)}"
            )

        target = entry.get("extract_to")
        if target is not None and target not in VALID_EXTRACT_TARGETS:
            errors.append(
                f"{where} ({entry.get('path')!r}): unknown extract_to {target!r}; "
                f"allowed {sorted(VALID_EXTRACT_TARGETS)} or null"
            )
        if target is not None and bucket != "open":
            errors.append(
                f"{where} ({entry.get('path')!r}): extract_to={target!r} on a {bucket!r} path; "
                "only open paths are extractable"
            )

        path = entry.get("path")
        if isinstance(path, str) and path:
            npath = _norm(path)
            if npath in seen:
                errors.append(f"{where}: duplicate path {npath!r} (first seen at entries[{seen[npath]}])")
            else:
                seen[npath] = idx
            clean.append({**entry, "path": npath})

    if not clean:
        return errors, report

    # --- header counts describe the map itself ---------------------------------------
    declared = data.get("counts")
    if isinstance(declared, dict):
        actual = {b: sum(1 for e in clean if e.get("bucket") == b) for b in sorted(VALID_BUCKETS)}
        actual["total"] = len(clean)
        for key, want in declared.items():
            got = actual.get(key)
            if got is None:
                errors.append(f"header counts: unknown key {key!r}; allowed {sorted(actual)}")
            elif want != got:
                errors.append(f"header counts.{key}: declares {want}, entries contain {got}")
        report["entry_counts"] = actual
    else:
        errors.append("header: 'counts' must be a mapping")

    # --- coverage: the falsifier -------------------------------------------------------
    files = tracked_files(repo)
    entry_paths = [e["path"] for e in clean]
    by_path = {e["path"]: e for e in clean}

    uncovered: list[str] = []
    # An invalid bucket is already an error above; tally it under "<invalid>" rather
    # than crashing, so the operator sees the full report instead of a traceback.
    file_counts = {b: 0 for b in sorted(VALID_BUCKETS)}
    used: dict[str, int] = {p: 0 for p in entry_paths}
    for f in files:
        winner = resolve(f, entry_paths)
        if winner is None:
            uncovered.append(f)
            continue
        used[winner] += 1
        bucket = by_path[winner].get("bucket")
        key = bucket if bucket in VALID_BUCKETS else "<invalid>"
        file_counts[key] = file_counts.get(key, 0) + 1

    for f in uncovered:
        errors.append(f"uncovered: {f}")

    stale = sorted(p for p, n in used.items() if n == 0)
    for p in stale:
        errors.append(f"stale entry: {p!r} matches no tracked file")

    report["tracked_files"] = len(files)
    report["file_counts"] = file_counts
    report["uncovered"] = uncovered
    report["stale_entries"] = stale
    report["extractable"] = sorted(e["path"] for e in clean if e.get("extract_to"))
    report["underivable"] = sorted(e["path"] for e in clean if str(e.get("basis", "")).strip() == "UNDERIVABLE")
    return errors, report


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--repo", default=None, help="repo root to check (default: the repo this script lives in)")
    ap.add_argument("--json", action="store_true", dest="as_json", help="emit the report as JSON")
    args = ap.parse_args(argv)

    repo = Path(args.repo).resolve() if args.repo else Path(__file__).resolve().parent.parent
    errors, report = check(repo)
    report["ok"] = not errors
    report["errors"] = errors

    if args.as_json:
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0 if not errors else 1

    name = report.get("repo_declared") or repo.name
    print(f"open/closed map check — {name}")
    print(f"  map:            {MAP_FILENAME}")
    if "tracked_files" in report:
        print(f"  tracked files:  {report['tracked_files']}")
        fc = report["file_counts"]
        for bucket in sorted(fc):
            print(f"    {bucket:<8} {fc[bucket]:>5} files")
        print(f"  map entries:    {report['entry_counts']['total']}")
        if report["extractable"]:
            print(f"  extractable:    {len(report['extractable'])} path(s) -> trust-standard")
        if report["underivable"]:
            print(f"  UNDERIVABLE:    {len(report['underivable'])} path(s) need an owner ruling")
            for p in report["underivable"]:
                print(f"      {p}")

    if errors:
        print(f"\nFAIL — {len(errors)} problem(s):")
        for e in errors:
            print(f"  - {e}")
        return 1

    print("\nOK — every tracked file resolves to exactly one bucket.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
