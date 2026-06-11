import argparse

import pandas as pd
from Bio import SeqIO


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="""
        Attach genomic coordinates to a DefenseFinder genes table by parsing the
        per-gene begin/end/strand from the prodigal protein FASTA (.prt) headers.
        """,
    )
    parser.add_argument(
        "--input_gene_df",
        required=True,
        help="DefenseFinder genes TSV (must contain a hit_id column).",
    )
    parser.add_argument(
        "--proteins",
        required=True,
        help="Prodigal protein FASTA (.prt) whose ids match the gene hit_id values.",
    )
    parser.add_argument(
        "--output_gene_df",
        required=True,
        help="Destination TSV with added gene_beg/gene_end/gene_strand columns.",
    )
    return parser.parse_args()


def protein_coordinates(prot_file: str) -> dict[str, tuple[int, int, int]]:
    """Map each protein id to (begin, end, strand) from its prodigal FASTA header."""
    coords = {}
    for record in SeqIO.parse(prot_file, "fasta"):
        # prodigal/pyrodigal headers: ">id # begin # end # strand # ...extra"
        _, beg, end, strand, *_ = record.description.split("#")
        coords[record.id] = (int(beg), int(end), int(strand))
    return coords


def add_gene_locations(
    genes: pd.DataFrame, coords: dict[str, tuple[int, int, int]]
) -> pd.DataFrame:
    """Join begin/end/strand onto each gene row via its protein hit_id."""
    # lookup raises KeyError on any hit_id absent from the protein FASTA
    locs = pd.DataFrame(
        [coords[hit_id] for hit_id in genes["hit_id"]],
        columns=["gene_beg", "gene_end", "gene_strand"],
        index=genes.index,
    )
    return genes.join(locs)


if __name__ == "__main__":
    args = parse_args()

    genes = pd.read_csv(args.input_gene_df, sep="\t")
    coords = protein_coordinates(args.proteins)
    genes = add_gene_locations(genes, coords)

    genes.to_csv(args.output_gene_df, sep="\t", index=False)
