import argparse

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="""
        Concatenate per-isolate Integron_Finder .integrons tables and reshape them
        to the unified MGE schema (id, iso, beg, end, type).
        """,
    )
    parser.add_argument(
        "--input_tsvs",
        required=True,
        nargs="+",
        help="Input Integron_Finder .integrons files to concatenate.",
    )
    parser.add_argument(
        "--output_df",
        required=True,
        help="Destination CSV in the unified MGE schema, indexed by id.",
    )
    return parser.parse_args()


def read_integrons(path: str) -> pd.DataFrame:
    """Read one .integrons file, skipping the leading '#' comment line(s).

    Integron_Finder prefixes the table with a version comment ('# Integron Finder
    version ...'); the header row ('ID_integron\\tID_replicon\\t...') follows."""
    with open(path) as fh:
        n_comment = 0
        for line in fh:
            if line.startswith("#"):
                n_comment += 1
            else:
                break
    return pd.read_csv(path, sep="\t", skiprows=n_comment)


def preformat(dfs: list[pd.DataFrame]) -> pd.DataFrame:
    """Reshape Integron_Finder element rows to the unified MGE schema.

    Each .integrons row is one element (attC site, integrase or protein) belonging
    to an integron; the integron's span is min(pos_beg)..max(pos_end) over its
    elements (1-based). The per-integron `type` (complete/CALIN/In0) is constant
    within a group. The id is "{iso}|{n}|{type}", with n parsed from ID_integron
    ("integron_01" -> 1), to stay unique across the concatenated tables.
    """
    cols = ["id", "iso", "beg", "end", "type"]
    df = pd.concat(dfs, ignore_index=True)
    # rows without an integron id are 'no integron found' placeholders
    df = df.dropna(subset=["ID_integron"])
    if df.empty:
        return pd.DataFrame(columns=cols).set_index("id")

    rows = []
    for (iso, id_integron), g in df.groupby(["ID_replicon", "ID_integron"]):
        types = g["type"].unique()
        assert len(types) == 1, f"integron {iso}|{id_integron} has mixed types {types}"
        itype = types[0]
        n = int(id_integron.split("_")[-1])
        rows.append(
            {
                "id": f"{iso}|{n}|{itype}",
                "iso": iso,
                "beg": int(g["pos_beg"].min()),
                "end": int(g["pos_end"].max()),
                "type": itype,
            }
        )
    return pd.DataFrame(rows, columns=cols).set_index("id", verify_integrity=True)


if __name__ == "__main__":
    args = parse_args()
    dfs = [read_integrons(f) for f in args.input_tsvs]
    df = preformat(dfs)
    df.to_csv(args.output_df)
