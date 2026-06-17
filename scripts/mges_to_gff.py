import pandas as pd
import argparse
import utils as ut


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build a per-junction GFF3 of MGE annotations with coordinates "
        "relative to the extracted junction region."
    )
    parser.add_argument(
        "--junc_id", required=True, help="Junction edge id to extract."
    )
    parser.add_argument(
        "--mges_to_junctions",
        nargs="+",
        required=True,
        help="Per-tool junction-assignment CSVs (results/mges_to_junctions/*.csv).",
    )
    parser.add_argument(
        "--mges",
        nargs="+",
        required=True,
        help="Per-tool raw element CSVs (results/mges/*.csv), used to recover 'type'.",
    )
    parser.add_argument(
        "--iso_len", required=True, help="CSV with isolate genome lengths."
    )
    parser.add_argument("--out_gff", required=True, help="Output GFF3 file.")
    return parser.parse_args()


def load_type_map(element_csvs):
    """Map element id -> type from the raw per-tool element tables.

    ids are unique across tools (genomad ids carry '|', defensefinder/ISEScan use
    distinct suffixes), so a single flat dict is unambiguous."""
    type_map = {}
    for f in element_csvs:
        df = pd.read_csv(f)
        type_map.update(dict(zip(df["id"], df["type"])))
    return type_map


def load_assignments(assigned_csvs, junc_id):
    """Concatenate the per-tool assignment tables and keep one junction's rows.

    A tool with zero assignments writes an empty file -- skip it so concat doesn't
    raise."""
    frames = []
    for f in assigned_csvs:
        try:
            frames.append(pd.read_csv(f))
        except pd.errors.EmptyDataError:
            continue
    if not frames:
        return pd.DataFrame()
    df = pd.concat(frames, ignore_index=True)
    return df[df["junction"] == junc_id]


def junction_relative_coords(ib, ie, jcb, jce, strand, L):
    """Map an element interval [ib, ie] (0-based) to coordinates relative to the
    junction region [jcb, jce] on a length-L circular genome, co-oriented with the
    canonical edge direction. The junction wraps the origin when jcb >= jce.
    Returns (start, end, interval_len, is_partial)."""
    interval_len = jce + L - jcb if jcb >= jce else jce - jcb
    start, end, is_partial = ut.transform_coordinates(ib, ie, jcb, jce, L)
    # reverse-orient when the edge is stored in the non-canonical direction
    if not strand:
        start, end = interval_len - end, interval_len - start
    return start, end, interval_len, is_partial


def to_gff_entries(annotations):
    gff_entries = []
    for ann in annotations:
        attr = {
            "ID": ann["ann_id"],
            "is_partial": ann["is_partial"],
            "tool": ann["type"],
        }
        gff_entries.append(
            {
                "seqid": ann["iso"],
                "source": ann["type"],
                "type": ann["type"],
                # GFF3 is 1-based, end-inclusive: a 0-based half-open [s, e) -> [s+1, e]
                "start": ann["start"] + 1,
                "end": ann["end"],
                "score": ".",
                "strand": ".",
                "phase": ".",
                "attributes": attr,
            }
        )
    return gff_entries


if __name__ == "__main__":
    args = parse_args()

    type_map = load_type_map(args.mges)
    iso_L = pd.read_csv(args.iso_len, index_col=0)["length"].to_dict()
    assigned = load_assignments(args.mges_to_junctions, args.junc_id)

    annotations = []
    region_lengths = {}
    for row in assigned.itertuples(index=False):
        iso = row.iso
        # element coords (already 0-based) transformed relative to the junction region
        start, end, interval_len, is_partial = junction_relative_coords(
            row.ib, row.ie, row.jcb, row.jce, bool(row.j_strand), iso_L[iso]
        )
        region_lengths[iso] = interval_len

        annotations.append(
            {
                "iso": iso,
                "ann_id": row.id,
                "start": start,
                "end": end,
                "strand": bool(row.j_strand),
                "is_partial": is_partial,
                "type": type_map[row.id],
            }
        )

    gff_entries = to_gff_entries(annotations)
    ut.write_gff3(gff_entries, args.out_gff, region_lengths)
    print(f"Wrote {len(gff_entries)} MGE annotations to {args.out_gff}")
