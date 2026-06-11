import pandas as pd
import pytest

from defensefinder_df_preformat import system_beg_end

# system_beg_end derives each system's [beg, end] span in bp from the coordinates
# of its boundary genes (sys_beg / sys_end, which are gene hit_ids). It consumes:
#   sdf -- systems table indexed by sys_id, columns sys_beg / sys_end (hit_ids)
#   gdf -- genes table indexed by (hit_id, sys_id), columns gene_beg / gene_end


def build_frames(systems, genes):
    """systems: [(sys_id, sys_beg_hit, sys_end_hit)]; genes: [(hit_id, sys_id, beg, end)]."""
    sdf = pd.DataFrame(
        [{"sys_beg": b, "sys_end": e} for _, b, e in systems],
        index=[s for s, _, _ in systems],
    )
    sdf.index.name = "sys_id"
    gdf = pd.DataFrame(
        [{"gene_beg": gb, "gene_end": ge} for _, _, gb, ge in genes],
        index=pd.MultiIndex.from_tuples(
            [(h, s) for h, s, _, _ in genes], names=["hit_id", "sys_id"]
        ),
    )
    return sdf, gdf


def test_single_gene_system():
    # sys_beg == sys_end: the span is just that gene's coordinates
    sdf, gdf = build_frames([("S1", "g1", "g1")], [("g1", "S1", 100, 200)])
    beg, end = system_beg_end(sdf, gdf)
    assert beg["S1"] == 100
    assert end["S1"] == 200


def test_normal_gene_order():
    # end gene starts after begin gene ends (es > be): span = begin.beg .. end.end
    sdf, gdf = build_frames(
        [("S2", "g1", "g2")], [("g1", "S2", 100, 200), ("g2", "S2", 300, 400)]
    )
    beg, end = system_beg_end(sdf, gdf)
    assert beg["S2"] == 100
    assert end["S2"] == 400


def test_inverse_gene_order():
    # end gene ends before begin gene ends (ee < be): genes are reversed, so the
    # span runs end.beg .. begin.end
    sdf, gdf = build_frames(
        [("S3", "g1", "g2")], [("g1", "S3", 300, 400), ("g2", "S3", 100, 200)]
    )
    beg, end = system_beg_end(sdf, gdf)
    assert beg["S3"] == 100
    assert end["S3"] == 400


def test_overlapping_genes():
    # genes overlap (bs < es and be < ee, but es <= be): span = begin.beg .. end.end
    sdf, gdf = build_frames(
        [("S4", "g1", "g2")], [("g1", "S4", 100, 300), ("g2", "S4", 250, 400)]
    )
    beg, end = system_beg_end(sdf, gdf)
    assert beg["S4"] == 100
    assert end["S4"] == 400


def test_multiple_systems():
    # each system is resolved independently and keyed by sys_id
    sdf, gdf = build_frames(
        [("A", "a1", "a1"), ("B", "b1", "b2")],
        [("a1", "A", 10, 20), ("b1", "B", 100, 200), ("b2", "B", 300, 400)],
    )
    beg, end = system_beg_end(sdf, gdf)
    assert (beg["A"], end["A"]) == (10, 20)
    assert (beg["B"], end["B"]) == (100, 400)


def test_unassigned_layout_raises():
    # end gene fully contains begin gene -> no branch matches, S/E stay -1
    sdf, gdf = build_frames(
        [("S5", "g1", "g2")], [("g1", "S5", 200, 300), ("g2", "S5", 100, 400)]
    )
    with pytest.raises(AssertionError):
        system_beg_end(sdf, gdf)


def test_begin_gene_must_be_forward():
    # begin gene with beg >= end (e.g. origin-wrapping) trips the bs < be guard
    sdf, gdf = build_frames(
        [("S6", "g1", "g2")], [("g1", "S6", 300, 100), ("g2", "S6", 400, 500)]
    )
    with pytest.raises(AssertionError):
        system_beg_end(sdf, gdf)


def test_end_gene_must_be_forward():
    # end gene with beg >= end trips the es < ee guard
    sdf, gdf = build_frames(
        [("S7", "g1", "g2")], [("g1", "S7", 100, 200), ("g2", "S7", 500, 400)]
    )
    with pytest.raises(AssertionError):
        system_beg_end(sdf, gdf)


def test_system_too_long_raises():
    # span >= 1e5 bp signals a likely wrap-around and is rejected
    sdf, gdf = build_frames(
        [("S8", "g1", "g2")], [("g1", "S8", 0, 10), ("g2", "S8", 200000, 200010)]
    )
    with pytest.raises(AssertionError):
        system_beg_end(sdf, gdf)
