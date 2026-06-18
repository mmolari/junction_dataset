import pytest

from utils import load_junction_positions, merged_relative_coords

# Two edges, three rows; covers strand bool parsing and the nested {edge: {iso}} shape.
CSV = """edge,iso,left_start,left_end,right_start,right_end,strand
e1,isoA,100,200,300,400,True
e1,isoB,4841508,4855335,29951,39466,False
e2,isoA,10,20,30,40,True
"""


def _write(tmp_path):
    p = tmp_path / "jp.csv"
    p.write_text(CSV)
    return str(p)


def test_nested_by_edge_then_iso(tmp_path):
    # {edge: {iso: ...}} preserving every (edge, iso) pair
    d = load_junction_positions(_write(tmp_path))
    assert set(d) == {"e1", "e2"}
    assert set(d["e1"]) == {"isoA", "isoB"}


def test_tuple_order_and_strand_type(tmp_path):
    # (left_start, left_end, right_start, right_end, strand) with a real bool
    d = load_junction_positions(_write(tmp_path))
    assert d["e1"]["isoA"] == (100, 200, 300, 400, True)
    assert d["e1"]["isoB"][4] is False  # not the truthy string "False"


def test_index_single_junction(tmp_path):
    # the pattern extract_junctions uses: index the nested dict by edge id
    one = load_junction_positions(_write(tmp_path))["e1"]
    assert set(one) == {"isoA", "isoB"}


def test_edge_filter_materialises_only_that_junction(tmp_path):
    # passing edge= keeps the nested shape but only the requested junction's rows
    d = load_junction_positions(_write(tmp_path), edge="e1")
    assert set(d) == {"e1"}
    assert d["e1"]["isoA"] == (100, 200, 300, 400, True)


def test_non_bool_strand_rejected(tmp_path):
    # a stray strand value makes pandas read the column as object, not bool;
    # bool("anything-nonempty") is True, so guard against the silent flip.
    bad = "edge,iso,left_start,left_end,right_start,right_end,strand\ne1,isoA,1,2,3,4,fwd\n"
    p = tmp_path / "bad.csv"
    p.write_text(bad)
    with pytest.raises(AssertionError):
        load_junction_positions(str(p))


# --- merged_relative_coords: per-part mapping into the extracted interval ---
# A single-part feature has parts == [(start, end)], so these also pin down the
# behaviour that the per-part loop must keep identical to a plain transform.


def test_single_part_fully_inside():
    # contained feature: exact relative coords, not partial
    assert merged_relative_coords([(1200, 1500)], 1000, 2000, 10000) == (200, 500, False)


def test_single_part_overhangs_edge():
    # extends past the interval end -> clipped to interval_len, flagged partial
    assert merged_relative_coords([(1800, 2200)], 1000, 2000, 10000) == (800, 1000, True)


def test_single_part_outside_returns_none():
    # no overlap at all -> feature dropped
    assert merged_relative_coords([(3000, 4000)], 1000, 2000, 10000) is None


def test_origin_wrapping_gene_outside_interval_dropped():
    # The bug: ggt wraps the origin (join 4992021..4992088 + 0..1667); its envelope
    # is the whole genome and used to fill any interval. Per-part, neither part
    # overlaps this mid-genome junction, so the gene is correctly dropped.
    parts = [(4992021, 4992088), (0, 1667)]
    assert merged_relative_coords(parts, 2782338, 2941580, 4992088) is None


def test_origin_wrapping_gene_inside_wrapping_interval_is_contiguous():
    # Origin-spanning gene inside an origin-spanning interval [9000,10000)+[0,1000):
    # both parts land adjacent across the wrap point (rel pos 1000), so the merged
    # span is contiguous and fully contained (not partial).
    parts = [(9500, 10000), (0, 500)]
    assert merged_relative_coords(parts, 9000, 1000, 10000) == (500, 1500, False)
