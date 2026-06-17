import argparse

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="""
        Concatenate per-isolate abricate result tables and reshape them to the
        unified MGE schema (id, iso, beg, end, type).
        """,
    )
    parser.add_argument(
        "--input_tsvs",
        required=True,
        nargs="+",
        help="Input abricate result TSV files to concatenate.",
    )
    parser.add_argument(
        "--output_df",
        required=True,
        help="Destination CSV in the unified MGE schema, indexed by id.",
    )
    return parser.parse_args()


def preformat(dfs: list[pd.DataFrame]) -> pd.DataFrame:
    """Reshape abricate rows to the unified MGE schema.

    abricate reports one row per resistance-gene hit with SEQUENCE (the isolate),
    1-based START/END (START <= END), and GENE. Hits are labelled with the constant
    type "AMR"; the gene identity is preserved in the id, composed from iso, gene and
    the global row index to stay unique across the concatenated tables.
    """
    df = pd.concat(dfs, ignore_index=True)
    df["iso"] = df["SEQUENCE"]
    df["beg"] = df["START"]
    df["end"] = df["END"]
    df["id"] = df["iso"] + "|" + df["GENE"] + "|" + df.index.astype(str)
    df["type"] = "AMR"
    df = df[["id", "iso", "beg", "end", "type"]]
    return df.set_index("id", verify_integrity=True)


if __name__ == "__main__":
    args = parse_args()
    dfs = [pd.read_csv(tsv, sep="\t") for tsv in args.input_tsvs]
    df = preformat(dfs)
    df.to_csv(args.output_df)
