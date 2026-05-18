import argparse
import pathlib
import numpy as np
import pandas as pd
import pypangraph as pp
from Bio import AlignIO, Seq, SeqIO, SeqRecord


def parse_args():
    parser = argparse.ArgumentParser(
        description="""
        Build a concatenated core-genome alignment by walking a guide
        strain's path through the pangenome graph, selecting core blocks
        above a minimum-length threshold, and gluing their per-block
        MAFFT alignments together. Blocks the guide traverses in reverse
        are reverse-complemented so the concatenation stays in guide
        orientation. Emits the joined FASTA plus a per-block coordinate
        table.
        """
    )
    parser.add_argument("--pangraph", type=str, required=True)
    parser.add_argument(
        "--aln_folder",
        type=str,
        required=True,
        help="folder containing per-block MAFFT alignments named block_<id>.fa",
    )
    parser.add_argument("--guide_strain", type=str, required=True)
    parser.add_argument(
        "--min_length",
        type=int,
        required=True,
        help="minimum block consensus length to include in the alignment",
    )
    parser.add_argument("--out_alignment", type=str, required=True)
    parser.add_argument("--out_coordinates", type=str, required=True)
    parser.add_argument(
        "--ungapped",
        action="store_true",
        help="drop any column containing a gap from all strains",
    )
    return parser.parse_args()


def anchor_block_order(pan, guide_strain, min_length):
    # Walk the guide path, then keep only core blocks longer than the minimum length.
    # Order = guide traversal order, so the resulting alignment reads in
    # the guide's chromosomal direction.
    if guide_strain not in pan.paths:
        raise ValueError(f"guide strain {guide_strain!r} not in pangraph paths")
    bdf = pan.to_blockstats_df()
    rows = []
    for node_id in pan.paths[guide_strain].nodes:
        node = pan.nodes[node_id]
        rows.append({"block_id": node.block_id, "strand": node.strand})
    order = pd.DataFrame(rows).merge(
        bdf[["core", "len"]],
        left_on="block_id",
        right_index=True,
        how="left",
    )
    mask = order["core"] & (order["len"] >= min_length)
    return order[mask].drop(columns="core").reset_index(drop=True)


def block_alignment_matrix(pan, block_id, aln_path, strains, flip):
    # The per-block FASTA is keyed by pangraph *node ids*, not strain
    # names. Build the node->strain map for this block and re-key the
    # alignment by strain. If the guide traverses the block in reverse,
    # flip every sequence so the concatenation stays in guide
    # orientation. Returns a (n_strains, aln_length) char matrix with
    # rows ordered as `strains`.
    aln = AlignIO.read(aln_path, "fasta")
    path_id_to_name = {path.id: path.name for path in pan.paths}
    block = pan.blocks[block_id]
    node_to_strain = {
        nid: path_id_to_name[pan.nodes[nid].path_id]
        for nid in block.alignment.node_ids()
    }
    seqs = {node_to_strain[rec.id]: str(rec.seq) for rec in aln}
    if flip:
        seqs = {s: str(Seq.Seq(v).reverse_complement()) for s, v in seqs.items()}
    return np.array([list(seqs[s]) for s in strains], dtype="U1")


def concatenate_alignments(pan, anchor_order, aln_folder, strains, ungapped):
    # Materialise each anchor block as a char matrix, optionally drop
    # any column with a gap, and concat horizontally. Coordinates track
    # the post-strip width when `ungapped` is set.
    mats = []
    coords = []
    cursor = 0
    for _, row in anchor_order.iterrows():
        block_id = row["block_id"]
        aln_path = aln_folder / f"block_{block_id}.fa"
        mat = block_alignment_matrix(
            pan, block_id, aln_path, strains, flip=not row["strand"]
        )
        if ungapped:
            keep = ~np.any(mat == "-", axis=0)
            mat = mat[:, keep]
        aln_length = mat.shape[1]
        mats.append(mat)
        coords.append(
            {
                "block_id": block_id,
                "strand": row["strand"],
                "len": int(row["len"]),
                "start": cursor,
                "end": cursor + aln_length,
                "aln_length": aln_length,
            }
        )
        cursor += aln_length
    full = np.concatenate(mats, axis=1)
    alignment = {strain: "".join(full[i]) for i, strain in enumerate(strains)}
    return alignment, pd.DataFrame(coords)


def write_fasta(alignment, out_path):
    records = [
        SeqRecord.SeqRecord(Seq.Seq(seq), id=strain, description="")
        for strain, seq in alignment.items()
    ]
    SeqIO.write(records, out_path, "fasta")


if __name__ == "__main__":
    args = parse_args()
    pan = pp.Pangraph.from_json(args.pangraph)
    strains = sorted(pan.strains())
    anchor_order = anchor_block_order(pan, args.guide_strain, args.min_length)
    aln_folder = pathlib.Path(args.aln_folder)
    alignment, coords = concatenate_alignments(
        pan, anchor_order, aln_folder, strains, args.ungapped
    )
    write_fasta(alignment, args.out_alignment)
    coords.to_csv(args.out_coordinates, index=False)
