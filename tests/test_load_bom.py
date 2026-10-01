"""
load_bom(): a BOM file in, part dicts out, every problem of every row kept
===========================================================================

The gate: a BOM with good rows and bad rows goes through load_bom() and
process_batch(); every bad row comes back as an error result naming its row
and each problem (key, code, reason), and every good row is sized.  No bad
row may be dropped, and the first one may not stop the load.

The files are written by the tests themselves, so the sdist's suite needs
nothing beside it.
"""

import io
import json

import pytest

from thermal_mesh_calculators.batch import process_batch, process_part
from thermal_mesh_calculators.intake import load_bom


PROJECT = {"t_fluid_K": 353.15, "t_surr_K": 353.15, "t_exh_K": 1073.15}

# Rows 2, 4 and 6 are good; row 3 has one problem, row 5 has two.  Row 2's
# material and row 4's zone are spellings the alias tables resolve.
GATE_CSV = (
    "part_id,material,component_class,convection_zone,convection_zone_in,"
    "convection_zone_out,thickness_mm,t_surf_K,surface,surface_in,"
    "surface_out,x_note\n"
    "BRK-1,Mild steel,structural,engine_beside,,,3.0,473.15,painted,,,"
    "a bracket\n"
    "BRK-2,steel_mild,structural,engine_beside,,,-2,473.15,painted,,,"
    "negative thickness\n"
    "SH-1,steel_mild,shield,,beside exhaust,engine_beside,0.8,,,"
    "aluminised,aluminised,\n"
    "PNL-1,unobtanium,structural,cabin_tunnel,,,2.0,hot,painted,,,"
    "two problems\n"
    "FL-1,,fluid,front_end_edges,,,,,,,,the air\n"
)


def _write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def _problems(part):
    return [(p["key"], p["code"]) for p in part.get("bom_errors", [])]


# --- The gate ---

class TestBadRowsAreReportedAndGoodRowsSize:

    def test_load_keeps_every_row_and_marks_the_bad_ones(self, tmp_path):
        parts = load_bom(_write(tmp_path, "bom.csv", GATE_CSV))
        assert [p["bom_row"] for p in parts] == [2, 3, 4, 5, 6]
        assert [p["part_id"] for p in parts] == [
            "BRK-1", "BRK-2", "SH-1", "PNL-1", "FL-1"]
        assert _problems(parts[1]) == [("thickness_mm", "OUT_OF_RANGE")]
        assert _problems(parts[3]) == [("material", "UNKNOWN_NAME"),
                                       ("t_surf_K", "WRONG_TYPE")]
        assert all("bom_errors" not in parts[i] for i in (0, 2, 4))

    def test_batch_reports_each_bad_row_and_sizes_the_good_ones(self, tmp_path):
        results = process_batch(load_bom(_write(tmp_path, "bom.csv", GATE_CSV)),
                                PROJECT)
        assert len(results) == 5
        good = [r for r in results if r["error"] is None]
        bad = [r for r in results if r["error"] is not None]
        assert [r["part_id"] for r in good] == ["BRK-1", "SH-1", "FL-1"]
        assert all(0 < r["governing_dx_mm"] < float("inf") for r in good)
        assert [(r["bom_row"], r["part_id"]) for r in bad] == [
            (3, "BRK-2"), (5, "PNL-1")]
        assert [(p["key"], p["code"]) for p in bad[0]["problems"]] == [
            ("thickness_mm", "OUT_OF_RANGE")]
        assert [(p["key"], p["code"]) for p in bad[1]["problems"]] == [
            ("material", "UNKNOWN_NAME"), ("t_surf_K", "WRONG_TYPE")]
        assert bad[0]["error"] == (
            "BOM row 3 did not load: thickness_mm must be > 0 mm, got -2.0")
        assert bad[1]["error"].startswith(
            "BOM row 5 did not load: 2 problems with part 'PNL-1':")
        assert "t_surf_K must be a number, got 'hot'" in bad[1]["error"]

    def test_a_good_row_sizes_as_the_same_part_written_by_hand(self, tmp_path):
        parts = load_bom(_write(tmp_path, "bom.csv", GATE_CSV))
        by_hand = {"part_id": "BRK-1", "material": "steel_mild",
                   "component_class": "structural",
                   "convection_zone": "engine_beside", "thickness_mm": 3.0,
                   "t_surf_K": 473.15, "surface": "painted",
                   "x_note": "a bracket", "bom_row": 2}
        assert parts[0] == by_hand
        assert process_part(parts[0], PROJECT) == process_part(by_hand, PROJECT)

    def test_the_names_were_resolved(self, tmp_path):
        parts = load_bom(_write(tmp_path, "bom.csv", GATE_CSV))
        assert parts[0]["material"] == "steel_mild"
        assert parts[2]["convection_zone_in"] == "exhaust_beside"


# --- CSV ---

class TestCsv:

    def test_a_clean_bom_has_no_problems_and_echoes_its_rows(self, tmp_path):
        text = ("part_id,material,component_class,convection_zone,"
                "thickness_mm,t_surf_K,surface\n"
                "A,steel_mild,structural,engine_beside,3,473.15,painted\n"
                "B,cast_iron,exhaust,exhaust_internal,6,1173.15,"
                "cast_iron_oxidised\n")
        parts = load_bom(_write(tmp_path, "clean.csv", text))
        assert all("bom_errors" not in p for p in parts)
        results = process_batch(parts, PROJECT)
        assert [(r["bom_row"], r["error"]) for r in results] == [
            (2, None), (3, None)]

    def test_blank_rows_are_skipped_and_rows_keep_their_line(self, tmp_path):
        text = ("part_id,component_class,convection_zone\n"
                "\n"
                "F-1,fluid,front_end_edges\n"
                ",,\n"
                "F-2,fluid,cabin_tunnel\n")
        parts = load_bom(_write(tmp_path, "gaps.csv", text))
        assert [(p["part_id"], p["bom_row"]) for p in parts] == [
            ("F-1", 3), ("F-2", 5)]

    def test_an_unknown_column_is_a_problem_of_every_row(self, tmp_path):
        text = ("part_id,component_class,convection_zone,velocity\n"
                "F-1,fluid,front_end_edges,3\n"
                "F-2,fluid,cabin_tunnel,\n"
                "F-3,fluid,cabin_tunnel,4\n")
        parts = load_bom(_write(tmp_path, "typo.csv", text))
        assert [_problems(p) for p in parts] == [
            [("velocity", "UNKNOWN_KEY")], [], [("velocity", "UNKNOWN_KEY")]]
        assert "Keys starting with 'x_' are carried through unread" in (
            parts[0]["bom_errors"][0]["reason"])

    def test_extra_values_are_a_malformed_row(self, tmp_path):
        text = ("part_id,component_class,convection_zone\n"
                "F-1,fluid,front_end_edges,surplus\n")
        [part] = load_bom(_write(tmp_path, "wide.csv", text))
        assert _problems(part) == [(None, "MALFORMED_ROW")]

    def test_a_reserved_column_is_refused(self, tmp_path):
        text = ("part_id,component_class,convection_zone,bom_row\n"
                "F-1,fluid,front_end_edges,7\n")
        [part] = load_bom(_write(tmp_path, "reserved.csv", text))
        assert _problems(part) == [("bom_row", "RESERVED_KEY")]
        assert part["bom_row"] == 2

    def test_a_repeated_part_id_names_the_first_row(self, tmp_path):
        text = ("part_id,component_class,convection_zone\n"
                "F-1,fluid,front_end_edges\n"
                "F-1,fluid,cabin_tunnel\n")
        parts = load_bom(_write(tmp_path, "dup.csv", text))
        assert _problems(parts[0]) == []
        assert _problems(parts[1]) == [("part_id", "DUPLICATE_PART_ID")]
        assert "already used by row 2" in parts[1]["bom_errors"][0]["reason"]

    def test_cells_are_read_as_the_schema_types(self, tmp_path):
        text = ("part_id,material,component_class,convection_zone_in,"
                "convection_zone_out,thickness_mm,shield_max_iter,h_gap\n"
                "ML-1,steel_mild,Multilayer Shield,exhaust_beside,"
                "shield_gap_confined,6e-1,80,12\n"
                "ML-2,steel_mild,multilayer_shield,exhaust_beside,"
                "shield_gap_confined,0.6,80.0,12\n")
        parts = load_bom(_write(tmp_path, "types.csv", text))
        assert parts[0]["component_class"] == "multilayer_shield"
        assert parts[0]["thickness_mm"] == 0.6
        assert parts[0]["shield_max_iter"] == 80
        assert isinstance(parts[0]["shield_max_iter"], int)
        assert _problems(parts[1]) == [("shield_max_iter", "WRONG_TYPE")]

    def test_names_as_written_when_not_resolving(self, tmp_path):
        parts = load_bom(_write(tmp_path, "bom.csv", GATE_CSV),
                         resolve_names=False)
        assert parts[0]["material"] == "Mild steel"
        assert _problems(parts[0]) == [("material", "UNKNOWN_NAME")]
        assert _problems(parts[2]) == [("convection_zone_in", "UNKNOWN_NAME")]

    def test_a_family_name_gets_its_candidates(self, tmp_path):
        text = ("part_id,material,component_class,convection_zone,"
                "thickness_mm,t_surf_K\n"
                "P-1,aluminium,structural,engine_beside,2,400\n")
        [part] = load_bom(_write(tmp_path, "family.csv", text))
        [problem] = part["bom_errors"]
        assert problem["key"] == "material"
        assert ("'aluminium' could be any of: aluminium_5052, aluminium_6061, "
                "cast_aluminium, cast_aluminium_a380") in problem["reason"]

    def test_a_header_naming_a_column_twice_is_not_a_bom(self, tmp_path):
        with pytest.raises(ValueError, match="duplicate column"):
            load_bom(_write(tmp_path, "twice.csv",
                            "part_id,material,material\nA,b,c\n"))

    def test_an_empty_file_and_a_bare_header_load_nothing(self, tmp_path):
        assert load_bom(_write(tmp_path, "empty.csv", "")) == []
        assert load_bom(_write(tmp_path, "header.csv", "part_id,material\n")) == []

    def test_an_open_file_a_byte_order_mark_and_bytes(self):
        text = "part_id,component_class,convection_zone\nF-1,fluid,front_end_edges\n"
        from_text = load_bom(io.StringIO("﻿" + text), fmt="csv")
        from_bytes = load_bom(io.BytesIO(("﻿" + text).encode("utf-8")),
                              fmt="csv")
        assert from_text == from_bytes == [
            {"part_id": "F-1", "component_class": "fluid",
             "convection_zone": "front_end_edges", "bom_row": 2}]

    def test_text_that_is_not_utf8_is_refused(self, tmp_path):
        path = tmp_path / "latin1.csv"
        path.write_bytes("part_id\nGr\xfc\xdfe\n".encode("latin-1"))
        with pytest.raises(ValueError, match="not UTF-8"):
            load_bom(path)


# --- JSON ---

class TestJson:

    def test_an_array_of_parts(self, tmp_path):
        rows = [
            {"part_id": "A", "material": "Aluminum 6061",
             "component_class": "structural", "convection_zone": "Beside Engine",
             "thickness_mm": 2.0, "t_surf_K": 400.0, "surface": "anodised"},
            "not an object",
            {"part_id": "C", "component_class": "fluid",
             "convection_zone": "front_end_edges", "thickness_mm": "3"},
            {"part_id": "A", "component_class": "fluid",
             "convection_zone": "cabin_tunnel"},
        ]
        parts = load_bom(_write(tmp_path, "bom.json", json.dumps(rows)))
        assert [p["bom_row"] for p in parts] == [1, 2, 3, 4]
        assert parts[0]["material"] == "aluminium_6061"
        assert parts[0]["convection_zone"] == "engine_beside"
        assert _problems(parts[0]) == []
        assert _problems(parts[1]) == [(None, "NOT_AN_OBJECT")]
        # A fluid region does not read thickness_mm, so its type is not
        # checked; JSON values are not coerced, so the string stays.
        assert parts[2]["thickness_mm"] == "3" and _problems(parts[2]) == []
        assert _problems(parts[3]) == [("part_id", "DUPLICATE_PART_ID")]
        results = process_batch(parts, PROJECT)
        assert [r["error"] is None for r in results] == [True, False, True, False]
        assert results[1]["part_id"] == "UNKNOWN" and results[1]["bom_row"] == 2

    def test_a_json_string_is_not_a_number(self, tmp_path):
        rows = [{"part_id": "A", "material": "steel_mild",
                 "component_class": "structural",
                 "convection_zone": "engine_beside", "thickness_mm": "3",
                 "t_surf_K": 400.0}]
        [part] = load_bom(_write(tmp_path, "s.json", json.dumps(rows)))
        assert _problems(part) == [("thickness_mm", "WRONG_TYPE")]

    def test_an_object_with_parts_and_a_version(self, tmp_path):
        doc = {"schema_version": 1, "parts": [
            {"part_id": "F", "component_class": "fluid",
             "convection_zone": "front_end_edges"}]}
        [part] = load_bom(_write(tmp_path, "v.json", json.dumps(doc)))
        assert _problems(part) == []

    @pytest.mark.parametrize("doc, match", [
        ({"schema_version": 2, "parts": []}, "schema_version 2 is not supported"),
        ({"parts": [], "project": {}}, "found project"),
        ({"rows": []}, "a JSON BOM is an array"),
        ("parts", "a JSON BOM is an array of part objects, got str"),
    ])
    def test_a_document_that_is_not_a_bom(self, tmp_path, doc, match):
        with pytest.raises(ValueError, match=match):
            load_bom(_write(tmp_path, "x.json", json.dumps(doc)))

    def test_invalid_json(self, tmp_path):
        with pytest.raises(ValueError, match="not valid JSON"):
            load_bom(_write(tmp_path, "bad.json", "[{"))


# --- Format ---

class TestFormat:

    def test_the_extension_decides_and_fmt_overrides(self, tmp_path):
        path = _write(tmp_path, "parts.txt",
                      "part_id,component_class,convection_zone\n"
                      "F,fluid,front_end_edges\n")
        with pytest.raises(ValueError, match="pass fmt='csv' or fmt='json'"):
            load_bom(path)
        assert load_bom(path, fmt="csv")[0]["part_id"] == "F"

    def test_an_unknown_format(self, tmp_path):
        with pytest.raises(ValueError, match="fmt must be 'csv' or 'json'"):
            load_bom(io.StringIO(""), fmt="xlsx")

    def test_upper_case_extension(self, tmp_path):
        path = _write(tmp_path, "PARTS.CSV",
                      "part_id,component_class,convection_zone\n"
                      "F,fluid,front_end_edges\n")
        assert load_bom(path)[0]["bom_row"] == 2
