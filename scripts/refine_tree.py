import argparse

import numpy as np
import treetime
from Bio import AlignIO, Phylo


def parse_args():
    parser = argparse.ArgumentParser(
        description="Refine branch lengths of a core-genome tree using TreeTime."
    )
    parser.add_argument("--tree_in", type=str, required=True, help="input newick tree")
    parser.add_argument(
        "--aln",
        type=str,
        required=True,
        help="alignment FASTA: full filtered alignment",
    )
    parser.add_argument(
        "--tree_out", type=str, required=True, help="output newick tree"
    )
    return parser.parse_args()


def aln_length(aln_file):
    """Number of columns in a gap-free alignment."""
    aln = AlignIO.read(aln_file, format="fasta")
    A = np.array(aln)
    assert np.all(A != "-"), f"{aln_file} contains gaps"
    return A.shape[1]


if __name__ == "__main__":
    args = parse_args()

    # Midpoint root + ladderize before optimization so TreeTime sees a
    # stable topology and the optimization is reproducible.
    tree = Phylo.read(args.tree_in, format="newick")
    tree.root_at_midpoint()
    tree.ladderize()

    myTree = treetime.TreeAnc(
        gtr="Jukes-Cantor",
        tree=tree,
        aln=args.aln,
        verbose=0,
    )
    myTree.tree.root.branch_length = 0.0

    myTree.optimize_tree(prune_short=True)

    myTree.tree.root_at_midpoint()
    myTree.tree.ladderize()
    Phylo.write(
        myTree.tree,
        args.tree_out,
        format="newick",
        format_branch_length="%.5e",
    )
