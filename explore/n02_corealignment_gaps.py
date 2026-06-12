# %%
import pathlib

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import pypangraph as pp
from Bio import AlignIO

fig_fld = pathlib.Path("figs/n02_corealignment_gaps")
fig_fld.mkdir(exist_ok=True, parents=True)

threshold_length = 500
aln_folder = pathlib.Path("../results/pangraph_core_block_alignments")

# %%
pan = pp.Pangraph.from_json("../results/pangraph.json")
bdf = pan.to_blockstats_df()

# anchor blocks = same set the core-genome alignment is built from
mask = bdf["core"] & (bdf["len"] >= threshold_length)
anchor_block_ids = bdf.index[mask].tolist()
print(f"Selected {len(anchor_block_ids)} core blocks (>= {threshold_length} bp)")


# %%
def trim_flanking_gaps(arr):
    # cumprod stays 1 only inside a row's leading-gap run, so sum(cumprod)
    # per row = leading-gap count; max across rows = how many columns must
    # be dropped to consume the longest such run. Same trick on the
    # reversed array for trailing gaps.
    is_gap = arr == "-"
    n_cols = arr.shape[1]
    left = int(np.cumprod(is_gap, axis=1).sum(axis=1).max())
    right = int(np.cumprod(is_gap[:, ::-1], axis=1).sum(axis=1).max())
    if left + right >= n_cols:
        return arr[:, 0:0], left, right
    return arr[:, left : n_cols - right], left, right


def interior_gap_frequencies(arr):
    # For each column with at least one gap, return the fraction of rows
    # (isolates) carrying that gap. Gap-free columns are dropped.
    is_gap = arr == "-"
    col_counts = is_gap.sum(axis=0)
    gapped = col_counts > 0
    return col_counts[gapped] / arr.shape[0], int(gapped.sum())


# %%
freqs_all = []
per_block = []
for block_id in anchor_block_ids:
    aln = AlignIO.read(aln_folder / f"block_{block_id}.fa", "fasta")
    arr = np.array(aln)
    trimmed, left, right = trim_flanking_gaps(arr)
    if trimmed.shape[1] == 0:
        continue
    freqs, n_gapped = interior_gap_frequencies(trimmed)
    freqs_all.append(freqs)
    per_block.append(
        {
            "block_id": block_id,
            "n_iso": arr.shape[0],
            "n_cols_orig": arr.shape[1],
            "n_cols_trimmed": trimmed.shape[1],
            "left_trim": left,
            "right_trim": right,
            "n_gapped_cols": n_gapped,
        }
    )

freqs_all = np.concatenate(freqs_all) if freqs_all else np.array([])
per_block = pd.DataFrame(per_block)
total_trimmed = int(per_block["n_cols_trimmed"].sum())
total_gapped = int(per_block["n_gapped_cols"].sum())
print(f"Trimmed core-alignment total length: {total_trimmed:,} bp")
print(
    f"Columns with interior gaps: {total_gapped:,} "
    f"({total_gapped / total_trimmed * 100:.2f}% of trimmed)"
)

# %%
# Cumulative distribution of interior gap frequencies across all core blocks.
# x = fraction of isolates with a gap at a given column;
# y = fraction of gapped columns at or below that frequency.
fig, ax = plt.subplots(figsize=(6, 4))
sns.ecdfplot(freqs_all, ax=ax)
ax.set_xlabel("gap frequency (fraction of isolates with a gap)")
ax.set_ylabel("cumulative fraction of gapped columns")
ax.set_title(
    f"{total_gapped:,} gapped cols / {total_trimmed / 1e6:.2f} Mbp "
    f"({total_gapped / total_trimmed * 100:.2f}%)"
)
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(fig_fld / "interior_gap_frequency_ecdf.png", dpi=300)
plt.show()

# %%
