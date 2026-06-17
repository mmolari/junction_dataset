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
                "start": ann["start"],
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
        L = iso_L[iso]
        # full junction region [left_start, right_end]; origin-wrapping when start >= end
        start, end = row.jcb, row.jce
        strand = bool(row.j_strand)
        interval_len = end + L - start if start >= end else end - start
        region_lengths[iso] = interval_len

        # element coords (already 0-based) transformed relative to the junction interval
        new_start, new_end, is_partial = ut.transform_coordinates(
            row.ib, row.ie, start, end, L
        )

        # co-orient with the canonical (reverse) edge direction
        if not strand:
            new_start, new_end = interval_len - new_end, interval_len - new_start

        annotations.append(
            {
                "iso": iso,
                "ann_id": row.id,
                "start": new_start,
                "end": new_end,
                "strand": strand,
                "is_partial": is_partial,
                "type": type_map[row.id],
            }
        )

    gff_entries = to_gff_entries(annotations)
    ut.write_gff3(gff_entries, args.out_gff, region_lengths)
    print(f"Wrote {len(gff_entries)} MGE annotations to {args.out_gff}")
