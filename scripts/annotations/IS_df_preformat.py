import argparse

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="""
        Concatenate per-isolate ISEScan summary tables and reshape them to the
        unified MGE schema (id, iso, beg, end, type).
        """,
    )
    parser.add_argument(
        "--input_tsvs",
        required=True,
        nargs="+",
        help="Input ISEScan summary TSV files to concatenate.",
    )
    parser.add_argument(
        "--output_df",
        required=True,
        help="Destination CSV in the unified MGE schema, indexed by id.",
    )
    return parser.parse_args()


def preformat(dfs: list[pd.DataFrame]) -> pd.DataFrame:
    """Reshape ISEScan rows to the unified MGE schema.

    ISEScan reports one row per IS copy with seqID (the isolate), family, and
    1-based isBegin/isEnd. The id is composed from iso, family and the global row
    index to stay unique across the concatenated tables.
    """
    df = pd.concat(dfs, ignore_index=True)
    df["iso"] = df["seqID"]
    df["beg"] = df["isBegin"]
    df["end"] = df["isEnd"]
    df["id"] = df["iso"] + "|" + df["family"] + "|" + df.index.astype(str)
    df["type"] = "IS"
    df = df[["id", "iso", "beg", "end", "type"]]
    return df.set_index("id", verify_integrity=True)


if __name__ == "__main__":
    args = parse_args()
    dfs = [pd.read_csv(tsv, sep="\t") for tsv in args.input_tsvs]
    df = preformat(dfs)
    df.to_csv(args.output_df)
