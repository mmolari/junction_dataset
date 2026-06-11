import argparse

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="""
        Concatenate per-isolate geNomad virus_summary tables and reshape them to
        the unified MGE schema (id, iso, beg, end, type).
        """,
    )
    parser.add_argument(
        "--input_tsvs",
        required=True,
        nargs="+",
        help="Input geNomad virus_summary TSV files to concatenate.",
    )
    parser.add_argument(
        "--output_df",
        required=True,
        help="Destination CSV in the unified MGE schema, indexed by id.",
    )
    return parser.parse_args()


def preformat(dfs: list[pd.DataFrame]) -> pd.DataFrame:
    """Reshape geNomad virus_summary rows to the unified MGE schema.

    Proviruses are named "{iso}|provirus_{beg}_{end}" with a "{beg}-{end}"
    coordinates field; iso is the seq_name prefix and beg/end come from
    coordinates. The id index must be unique across the concatenated tables.
    """
    df = pd.concat(dfs, ignore_index=True)
    df["iso"] = df["seq_name"].str.split("|").str[0]
    df["beg"] = df["coordinates"].str.split("-").str[0].astype(int)
    df["end"] = df["coordinates"].str.split("-").str[1].astype(int)
    df["type"] = "prophage"
    df = df.rename(columns={"seq_name": "id"})
    df = df[["id", "iso", "beg", "end", "type"]]
    return df.set_index("id", verify_integrity=True)


if __name__ == "__main__":
    args = parse_args()
    dfs = [pd.read_csv(tsv, sep="\t") for tsv in args.input_tsvs]
    df = preformat(dfs)
    df.to_csv(args.output_df)
