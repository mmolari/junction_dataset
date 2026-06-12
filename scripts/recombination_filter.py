import argparse
import numpy as np
import pandas as pd
from Bio import AlignIO, Seq, SeqIO, SeqRecord


def parse_args():
    parser = argparse.ArgumentParser(
        description="""
        Mask recombined regions in an ungapped core-genome alignment.
        For every non-consensus site on every genome, count how many
        non-consensus sites lie within a circular window of radius
        `--window` (modulo alignment length) on the same genome,
        including the focal site itself. If that count exceeds
        `--max_snps`, the full [focal-window, focal+window] interval
        (2*window+1 bp) is masked from the alignment. Overlapping
        masked intervals are merged. Emits the column-filtered
        alignment as FASTA and the merged masked intervals as CSV.
        Wrap-around regions are reported as two linear intervals (one
        starting at 0 and one ending at L).
        """
    )
    parser.add_argument("--alignment", type=str, required=True)
    parser.add_argument(
        "--window",
        type=int,
        required=True,
        help="half-window radius in bp; masked region per focal SNP is 2*window+1 bp",
    )
    parser.add_argument(
        "--max_snps",
        type=int,
        required=True,
        help="trigger masking when SNP count (incl. focal) in the window is strictly greater than this",
    )
    parser.add_argument("--out_alignment", type=str, required=True)
    parser.add_argument("--out_intervals", type=str, required=True)
    return parser.parse_args()


def load_alignment(path):
    aln = AlignIO.read(path, "fasta")
    strains = [rec.id for rec in aln]
    matrix = np.array([list(str(rec.seq)) for rec in aln], dtype="U1")
    return strains, matrix


def consensus_per_column(matrix):
    # Most frequent character per column (ties broken by lexicographic
    # order via np.unique). Mirrors consensus() in the reference script
    # at scripts/pangraph/corealn_without_recombination.py.
    def site_consensus(col):
        letters, counts = np.unique(col, return_counts=True)
        return letters[np.argmax(counts)]

    return np.apply_along_axis(site_consensus, 0, matrix)


def triggering_positions(snp_pos, window, max_snps, length):
    # snp_pos: sorted positions of non-consensus sites on one genome.
    # Return the subset of snp_pos whose circular window of radius
    # `window` contains strictly more than `max_snps` SNPs (focal
    # included). Vectorised via np.searchsorted on a tripled array so
    # wrap-around is handled naturally.
    if snp_pos.size == 0:
        return snp_pos
    tripled = np.concatenate([snp_pos - length, snp_pos, snp_pos + length])
    lo = np.searchsorted(tripled, snp_pos - window, side="left")
    hi = np.searchsorted(tripled, snp_pos + window, side="right")
    counts = hi - lo
    return snp_pos[counts > max_snps]


def apply_window_mask(keep, focals, window, length):
    # In-place: clear `keep[(p-window) mod L : (p+window+1) mod L)` for
    # every focal p, splitting at the origin if the interval wraps.
    for p in focals:
        s = (p - window) % length
        e = (p + window + 1) % length
        if s < e:
            keep[s:e] = False
        else:
            keep[s:] = False
            keep[:e] = False


def mask_to_intervals(keep):
    # Linear runs of False as (start, end) half-open. Wrap-around
    # regions show up as a (0, e) and (s, L) pair.
    masked = ~keep
    if not masked.any():
        return []
    diff = np.diff(masked.astype(np.int8))
    starts = np.flatnonzero(diff == 1) + 1
    ends = np.flatnonzero(diff == -1) + 1
    if masked[0]:
        starts = np.concatenate([[0], starts])
    if masked[-1]:
        ends = np.concatenate([ends, [len(masked)]])
    return list(zip(starts.tolist(), ends.tolist()))


def write_alignment(strains, matrix, path):
    records = [
        SeqRecord.SeqRecord(Seq.Seq("".join(row)), id=strain, description="")
        for strain, row in zip(strains, matrix)
    ]
    SeqIO.write(records, path, "fasta")


def write_intervals(intervals, path):
    df = pd.DataFrame(intervals, columns=["start", "end"])
    df["length"] = df["end"] - df["start"]
    df.to_csv(path, index=False)


if __name__ == "__main__":
    args = parse_args()
    strains, matrix = load_alignment(args.alignment)
    n_strains, length = matrix.shape
    consensus = consensus_per_column(matrix)
    is_snp = matrix != consensus

    keep = np.ones(length, dtype=bool)
    for row in range(n_strains):
        snp_pos = np.flatnonzero(is_snp[row])
        focals = triggering_positions(snp_pos, args.window, args.max_snps, length)
        apply_window_mask(keep, focals, args.window, length)

    filtered = matrix[:, keep]
    intervals = mask_to_intervals(keep)
    write_alignment(strains, filtered, args.out_alignment)
    write_intervals(intervals, args.out_intervals)

    print(
        f"alignment {length} -> {filtered.shape[1]} bp "
        f"(masked {length - filtered.shape[1]} bp across "
        f"{len(intervals)} intervals)"
    )
