import argparse

import pypangraph as pp


def parse_args():
    parser = argparse.ArgumentParser(
        description="""
        Detect backbone junctions in the all-isolates pangraph and write their
        flanking-core-block genomic positions to a CSV. Columns:
        edge, iso, left_start, left_end, right_start, right_end, strand.
        """
    )
    parser.add_argument(
        "--pangraph", type=str, required=True, help="Input pangraph JSON file."
    )
    parser.add_argument(
        "--min_length",
        type=int,
        default=500,
        help="Minimum block length to be considered backbone (L_thr).",
    )
    parser.add_argument(
        "--out_csv", type=str, required=True, help="Output CSV file."
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    pan = pp.Pangraph.from_json(args.pangraph)
    bj = pp.junctions.BackboneJunctions(pan, L_thr=args.min_length)
    # MultiIndex (edge, iso) -> flatten to columns for a readable CSV
    df = bj.positions().reset_index()
    df.to_csv(args.out_csv, index=False)
    print(f"Wrote {len(df)} junction positions to {args.out_csv}")
