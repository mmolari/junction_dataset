from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import json
import pathlib
import utils as ut
import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract junction sequences and annotations from GenBank files."
    )
    parser.add_argument(
        "--gbk-fld", type=str, required=True, help="Folder containing GenBank files."
    )
    parser.add_argument(
        "--junc-id", type=str, required=True, help="Junction ID to extract."
    )
    parser.add_argument(
        "--junc-pos-file",
        type=str,
        required=True,
        help="JSON file with junction positions.",
    )
    parser.add_argument(
        "--out-fa",
        type=str,
        required=True,
        help="Output FASTA file for extracted sequences.",
    )
    parser.add_argument(
        "--out-ann",
        type=str,
        required=True,
        help="Output GFF3 file for extracted annotations.",
    )
    return parser.parse_args()


def load_junction_positions(json_file):
    with open(json_file, "r") as f:
        data = json.load(f)
    return data


if __name__ == "__main__":
    args = parse_args()

    j_pos = load_junction_positions(args.junc_pos_file)
    gbk_fld = pathlib.Path(args.gbk_fld)
    junct_id = args.junc_id

    records = []
    annotations = []
    region_lengths = {}
    for iso, positions in j_pos[junct_id].items():
        gbk_file = gbk_fld / f"{iso}.gbk"
        print(f"Processing {gbk_file}...")
        start, _, _, end, strand = positions
        seq = ut.extract_genbank_region(
            gbk_file, start, end, forward_orientation=strand
        )
        region_lengths[iso] = len(seq)
        ann = ut.extract_genbank_annotations(
            gbk_file, start, end, forward_orientation=strand
        )
        records.append(
            SeqIO.SeqRecord(
                Seq(seq),
                id=iso,
                name=junct_id,
                description=f"{start}-{end} ({'+' if strand else '-'})",
            )
        )
        annotations.extend(ann)
        print(f"  Extracted sequence length: {len(seq)}")
        print(f"  Found {len(ann)} overlapping annotations.")

    SeqIO.write(records, args.out_fa, "fasta")
    ut.write_gff3(annotations, args.out_ann, sequence_region_lengths=region_lengths)
    print(f"Wrote {len(records)} sequences to {args.out_fa}")
    print(f"Wrote {len(annotations)} annotations to {args.out_ann}")
