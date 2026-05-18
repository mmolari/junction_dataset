import subprocess
import sys
from pathlib import Path

import numpy as np
from recombination_filter import (
    apply_window_mask,
    consensus_per_column,
    mask_to_intervals,
    triggering_positions,
)


# --- consensus_per_column -------------------------------------------------


def test_consensus_majority_per_column():
    # Column 0: AAA -> A; col 1: CCA -> C; col 2: GCG -> G; col 3: TTA -> T.
    A = np.array(
        [
            ["A", "C", "G", "T"],
            ["A", "C", "C", "T"],
            ["A", "A", "G", "A"],
        ]
    )
    assert "".join(consensus_per_column(A)) == "ACGT"


def test_consensus_tie_breaks_lexicographically():
    # Two A's and two T's: np.unique returns ['A', 'T'] sorted, so
    # argmax of [2, 2] picks index 0 -> 'A'.
    A = np.array([["A"], ["T"], ["A"], ["T"]])
    assert consensus_per_column(A).tolist() == ["A"]


def test_consensus_single_row_is_identity():
    A = np.array([["A", "C", "G", "T"]])
    assert "".join(consensus_per_column(A)) == "ACGT"


# --- triggering_positions -------------------------------------------------


def test_triggering_linear_cluster():
    # SNPs at 2,3,5 fall within a radius-3 window of each other on a
    # linear length-20 alignment; SNP at 10 is isolated.
    pos = np.array([2, 3, 5, 10])
    out = triggering_positions(pos, window=3, max_snps=2, length=20)
    assert out.tolist() == [2, 3, 5]


def test_triggering_circular_wrap_around_origin():
    # Positions 1, 18, 19 are all circular neighbours on L=20 within
    # radius 3 (e.g. 1 and 19 are 2 bp apart going through 0).
    pos = np.array([1, 18, 19])
    out = triggering_positions(pos, window=3, max_snps=1, length=20)
    assert out.tolist() == [1, 18, 19]


def test_triggering_empty_input():
    out = triggering_positions(np.array([], dtype=int), window=3, max_snps=2, length=20)
    assert out.size == 0


def test_triggering_sparse_no_trigger():
    # Two isolated SNPs far apart -> each window contains only itself.
    pos = np.array([5, 50])
    out = triggering_positions(pos, window=3, max_snps=1, length=100)
    assert out.size == 0


def test_triggering_boundary_inclusive():
    # Distance exactly == window must count as "within" the window.
    pos = np.array([0, 3])
    out = triggering_positions(pos, window=3, max_snps=1, length=100)
    assert out.tolist() == [0, 3]


# --- apply_window_mask ----------------------------------------------------


def test_apply_window_mask_no_wrap():
    keep = np.ones(20, dtype=bool)
    apply_window_mask(keep, np.array([10]), window=2, length=20)
    # focal=10, window=2 -> mask [10-2, 10+2+1) = [8, 13).
    assert keep[8:13].sum() == 0
    assert keep[:8].all() and keep[13:].all()


def test_apply_window_mask_wraps_origin():
    keep = np.ones(20, dtype=bool)
    apply_window_mask(keep, np.array([19]), window=3, length=20)
    # [(19-3) % 20, (19+3+1) % 20) = [16, 3) -> mask 16..19 and 0..2.
    assert not keep[16:20].any()
    assert not keep[0:3].any()
    assert keep[3:16].all()


def test_apply_window_mask_overlapping_focals_union():
    keep = np.ones(20, dtype=bool)
    apply_window_mask(keep, np.array([5, 7]), window=2, length=20)
    # [3, 8) U [5, 10) = [3, 10).
    assert not keep[3:10].any()
    assert keep[:3].all() and keep[10:].all()


# --- mask_to_intervals ----------------------------------------------------


def test_mask_to_intervals_empty_when_nothing_masked():
    keep = np.ones(10, dtype=bool)
    assert mask_to_intervals(keep) == []


def test_mask_to_intervals_single_middle_run():
    keep = np.ones(10, dtype=bool)
    keep[3:7] = False
    assert mask_to_intervals(keep) == [(3, 7)]


def test_mask_to_intervals_wrap_around_appears_as_two_intervals():
    # A region that wraps the origin in circular space shows up as a
    # head interval starting at 0 and a tail interval ending at L.
    keep = np.ones(10, dtype=bool)
    keep[0:2] = False
    keep[8:10] = False
    assert mask_to_intervals(keep) == [(0, 2), (8, 10)]


def test_mask_to_intervals_all_masked():
    keep = np.zeros(10, dtype=bool)
    assert mask_to_intervals(keep) == [(0, 10)]


# --- end-to-end smoke test ------------------------------------------------


def _parse_fasta(path):
    strains, seqs, current = [], [], ""
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if current:
                seqs.append(current)
            strains.append(line[1:].strip())
            current = ""
        else:
            current += line.strip()
    if current:
        seqs.append(current)
    return strains, seqs


def test_smoke_end_to_end_small_alignment(tmp_path):
    # 10 genomes x 100 sites, all 'A' except:
    #   s0: 4 clustered 'T' SNPs at positions 20..23 (recombination signal)
    #   s1: 1 isolated 'T' SNP at position 70 (real point mutation)
    # With window=10, max_snps=2 each cluster focal sees 4 neighbours in
    # its window so all four trigger; merged mask is [10, 34).
    n_strains, length = 10, 100
    matrix = np.full((n_strains, length), "A", dtype="U1")
    matrix[0, 20:24] = "T"
    matrix[1, 70] = "T"

    fa_in = tmp_path / "in.fa"
    with fa_in.open("w") as fh:
        for i, row in enumerate(matrix):
            fh.write(f">s{i}\n{''.join(row)}\n")

    fa_out = tmp_path / "out.fa"
    csv_out = tmp_path / "intervals.csv"
    script = Path(__file__).resolve().parent.parent / "scripts" / "recombination_filter.py"

    subprocess.run(
        [
            sys.executable, str(script),
            "--alignment", str(fa_in),
            "--window", "10",
            "--max_snps", "2",
            "--out_alignment", str(fa_out),
            "--out_intervals", str(csv_out),
        ],
        check=True, capture_output=True, text=True,
    )

    # Single merged interval covering the cluster window.
    assert csv_out.read_text().splitlines() == ["start,end,length", "10,34,24"]

    out_strains, out_seqs = _parse_fasta(fa_out)
    assert out_strains == [f"s{i}" for i in range(n_strains)]
    assert all(len(s) == length - 24 for s in out_seqs)
    # s0's cluster fell entirely inside the masked region -> no SNPs survive.
    assert out_seqs[0] == "A" * 76
    # s1's isolated SNP at original position 70 -> position 70-24=46 after filtering.
    assert out_seqs[1][46] == "T"
    assert out_seqs[1][:46] + out_seqs[1][47:] == "A" * 75
    # Untouched strains stay all-A.
    for seq in out_seqs[2:]:
        assert seq == "A" * 76
