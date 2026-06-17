import pandas as pd
import pytest

from mges_to_gff import (
    junction_relative_coords,
    to_gff_entries,
    load_type_map,
    load_assignments,
)

# columns emitted by scripts/assign_junctions.py (results/mges_to_junctions/*.csv)
ASSIGN_COLS = "id,junction,iso,in_core,in_acc,ib,ie,jcb,jce,jab,jae,j_strand"


# --- junction_relative_coords: the coordinate transform -----------------------


def test_forward_non_wrapping_is_plain_offset():
    # junction [20, 80] on a length-100 circle (len 60); element [30, 50] inside it
    assert junction_relative_coords(30, 50, 20, 80, True, 100) == (10, 30, 60, False)


def test_reverse_strand_flips_within_interval():
    # same geometry, reverse-oriented edge: p -> interval_len - p, ends swapped
    assert junction_relative_coords(30, 50, 20, 80, False, 100) == (30, 50, 60, False)


def test_forward_wrapping_junction():
    # junction [80, 20] wraps the origin (len 40); element [90, 15] also wraps
    assert junction_relative_coords(90, 15, 80, 20, True, 100) == (10, 35, 40, False)


def test_reverse_wrapping_junction_flips():
    # same as above, reverse edge: 40-35, 40-10 -> (5, 30)
    assert junction_relative_coords(90, 15, 80, 20, False, 100) == (5, 30, 40, False)


def test_partial_when_element_overhangs_boundary():
    # element [10, 50] starts left of junction [20, 80] -> partial, start clipped to 0
    assert junction_relative_coords(10, 50, 20, 80, True, 100) == (0, 30, 60, True)


def test_origin_wrapping_real_junction_regression():
    # the wrapping edge verified by hand against the pipeline output:
    # jcb=4841508 > jce=39466, L=4898796, reverse strand
    assert junction_relative_coords(
        4855474, 4898795, 4841508, 39466, False, 4898796
    ) == (39467, 82788, 96754, False)


# --- to_gff_entries: GFF3 record shape ----------------------------------------


def test_to_gff_entries_shape_and_attributes():
    ann = [
        {
            "iso": "isoA",
            "ann_id": "x1",
            "start": 10,
            "end": 30,
            "strand": True,
            "is_partial": False,
            "type": "IS",
        }
    ]
    (e,) = to_gff_entries(ann)
    assert e["seqid"] == "isoA"
    # source and type both carry the tool category; strand is unused (".")
    assert e["source"] == "IS" and e["type"] == "IS"
    # 0-based half-open [10, 30) -> GFF3 1-based end-inclusive [11, 30]
    assert e["start"] == 11 and e["end"] == 30
    assert e["score"] == "." and e["strand"] == "." and e["phase"] == "."
    assert e["attributes"] == {"ID": "x1", "is_partial": False, "tool": "IS"}


# --- load_type_map: id -> type across tools -----------------------------------


def test_load_type_map_merges_across_files(tmp_path):
    f1 = tmp_path / "genomad.csv"
    f1.write_text("id,iso,beg,end,type\na,iso1,1,2,prophage\n")
    f2 = tmp_path / "ISEScan.csv"
    f2.write_text("id,iso,beg,end,type\nb,iso1,3,4,IS\n")
    assert load_type_map([str(f1), str(f2)]) == {"a": "prophage", "b": "IS"}


# --- load_assignments: concat + filter, tolerate empty tool tables ------------


def test_load_assignments_filters_to_junction(tmp_path):
    f = tmp_path / "genomad.csv"
    f.write_text(
        ASSIGN_COLS + "\n"
        "x,J1,iso1,subset,subset,10,20,5,30,8,25,True\n"
        "y,J2,iso1,subset,subset,1,2,0,5,1,4,True\n"
    )
    df = load_assignments([str(f)], "J1")
    assert list(df["id"]) == ["x"]


def test_load_assignments_skips_empty_file(tmp_path):
    # a tool with zero assignments writes an empty (zero-byte) file -> EmptyDataError
    good = tmp_path / "genomad.csv"
    good.write_text(ASSIGN_COLS + "\nx,J1,iso1,subset,subset,10,20,5,30,8,25,True\n")
    empty = tmp_path / "ISEScan.csv"
    empty.write_text("")
    df = load_assignments([str(good), str(empty)], "J1")
    assert list(df["id"]) == ["x"]


def test_load_assignments_all_empty_returns_empty(tmp_path):
    empty = tmp_path / "ISEScan.csv"
    empty.write_text("")
    assert load_assignments([str(empty)], "J1").empty
