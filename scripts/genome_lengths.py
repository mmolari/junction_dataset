import argparse
import csv
from pathlib import Path
import pandas as pd

from Bio import SeqIO


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute genome lengths from GenBank files and write a CSV summary.",
    )
    parser.add_argument(
        "--gbk_files",
        nargs="+",
        help="Input GenBank files to process.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Destination CSV filename.",
    )
    return parser.parse_args()


def genome_length(path: Path) -> int:
    """Return the total sequence length recorded in a GenBank file."""
    record = SeqIO.read(path, "genbank")
    return len(record.seq)


if __name__ == "__main__":
    args = parse_args()
    rows = []

    for gbk in args.gbk_files:
        gbk_path = Path(gbk)
        if not gbk_path.exists():
            raise FileNotFoundError(f"GenBank file not found: {gbk}")
        acc = gbk_path.stem
        length = genome_length(gbk_path)
        rows.append((acc, length))

    df = pd.DataFrame(rows, columns=["accession", "length"])
    df.to_csv(args.output, index=False)
