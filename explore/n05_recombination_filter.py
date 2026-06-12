# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pathlib
from matplotlib.colors import LogNorm
from Bio import AlignIO, Phylo


fig_fld = pathlib.Path("figs/n05_recombination_filter")
fig_fld.mkdir(exist_ok=True, parents=True)

aln_path = "../results/core_genome_alignment/ungapped.fa"
tree_path = "../config/polished_tree.nwk"
intervals_path = "../results/recombination_filter/masked_intervals.csv"
bin_size = 10_000
win_size = 1_000

# %%
# Load the ungapped core-genome alignment as a (N_iso, L) uint8 array.
# This is the same input fed to scripts/recombination_filter.py, so
# masked-interval positions index directly into arr.
aln = AlignIO.read(aln_path, "fasta")
strain_order = [rec.id for rec in aln]
n_iso, aln_length = len(aln), aln.get_alignment_length()
arr = np.empty((n_iso, aln_length), dtype=np.uint8)
for i, rec in enumerate(aln):
    arr[i] = np.frombuffer(str(rec.seq).encode("ascii"), dtype=np.uint8)
print(f"Alignment: {arr.shape[0]} isolates x {arr.shape[1]:,} columns")

# %%
# Per-column consensus + (N, L) non-consensus boolean matrix. Anything
# outside ACGT is flagged as non-consensus (matches n01).
bases_u8 = np.array([ord(b) for b in "acgt"], dtype=np.uint8)
counts = np.stack([(arr == b).sum(axis=0) for b in bases_u8])  # (4, L)
consensus = bases_u8[counts.argmax(axis=0)]  # (L,)
non_consensus = arr != consensus[None, :]  # (N, L) bool
polymorphic_col = non_consensus.any(axis=0)  # (L,)
poly_positions = np.where(polymorphic_col)[0]
print(
    f"Polymorphic columns: {polymorphic_col.sum():,} "
    f"({polymorphic_col.mean() * 100:.3f}% of total)"
)

# %%
# Load masked intervals and derive a (L,) "kept" bool array. Intervals
# are half-open [start, end) in gap-free coords (see write_intervals in
# scripts/recombination_filter.py).
intervals = pd.read_csv(intervals_path)
kept = np.ones(arr.shape[1], dtype=bool)
for s, e in zip(intervals["start"], intervals["end"]):
    kept[s:e] = False
poly_kept = poly_positions[kept[poly_positions]]
kept_length = int(kept.sum())
n_poly_kept = int(poly_kept.size)
n_poly_total = int(poly_positions.size)
print(
    f"masked {len(intervals)} intervals, "
    f"alignment {arr.shape[1]:,} -> {kept_length:,} bp, "
    f"polymorphic {n_poly_total:,} -> {n_poly_kept:,}"
)

# %%
# Figure 1: 10 kb SNP-density histogram, all-polymorphic in light gray
# with surviving-polymorphic overlaid in blue. y-log so rare hot bins
# stay visible. Title quotes kept counts.
bins = np.arange(0, arr.shape[1] + bin_size, bin_size)
fig, ax = plt.subplots(figsize=(10, 3))
ax.hist(poly_positions, bins=bins, color="lightgray", label="all polymorphic")
ax.hist(poly_kept, bins=bins, color="C0", label="kept polymorphic")
ax.set_yscale("log")
ax.set_xlabel("ungapped alignment position (bp)")
ax.set_ylabel("polymorphic columns per 10 kb (log)")
ax.set_title(
    f"polymorphic cols / 10 kb — "
    f"kept {n_poly_kept:}/{n_poly_total:} polymorphic, "
    f"{kept_length / 1_000_000:.2f}/{arr.shape[1] / 1_000_000:.2f} Mbp"
)
ax.legend()
plt.tight_layout()
plt.savefig(fig_fld / "snp_density_10kb_filtered.png", dpi=300)
plt.show()

# %%
# Reorder rows by polished_tree leaf order (same idiom as n01).
tree = Phylo.read(tree_path, "newick")
leaf_order = [leaf.name for leaf in tree.get_terminals()]
strain_to_idx = {s: i for i, s in enumerate(strain_order)}
missing = [s for s in leaf_order if s not in strain_to_idx]
if missing:
    print(f"WARNING: {len(missing)} tree leaves not in alignment: {missing[:5]}...")
row_order = [strain_to_idx[s] for s in leaf_order if s in strain_to_idx]

# %%
# Figure 2: per-isolate non-consensus density in 1 kb windows,
# tree-ordered, with masked intervals overlaid in transparent red.
L = arr.shape[1]
L_trim = (L // win_size) * win_size
nc_windows = (
    non_consensus[:, :L_trim]
    .reshape(arr.shape[0], L_trim // win_size, win_size)
    .sum(axis=2)
)
nc_windows = nc_windows[row_order]
print(
    f"Heatmap: {nc_windows.shape[0]} isolates x "
    f"{nc_windows.shape[1]} windows of {win_size} bp"
)

fig, (tree_ax, ax) = plt.subplots(
    1,
    2,
    figsize=(16, 8),
    gridspec_kw={"width_ratios": [1, 10], "wspace": 0.05},
)

# Left panel: phylogeny. Suppress leaf labels (222 names is unreadable)
# and pin y-limits so the tree spans the same N slots as the heatmap.
Phylo.draw(
    tree,
    axes=tree_ax,
    do_show=False,
    show_confidence=False,
    label_func=lambda c: "",
)
n_leaves = len(leaf_order)
tree_ax.set_ylim(n_leaves + 0.5, 0.5)
tree_ax.set_xlabel("branch length")
tree_ax.set_ylabel("")
tree_ax.set_title("")
tree_ax.set_yticks([])
for spine in ("top", "right", "left"):
    tree_ax.spines[spine].set_visible(False)

# Right panel: heatmap + masked-region red overlay. axvspan accepts
# floats so we convert interval bp to window-index units; clipping to
# axis limits handles intervals overlapping the trimmed last window.
vmax = max(2, int(nc_windows.max()))
im = ax.imshow(
    nc_windows,
    aspect="auto",
    interpolation="nearest",
    cmap="Greys",
    norm=LogNorm(vmin=0.4, vmax=vmax),
)
for s, e in zip(intervals["start"], intervals["end"]):
    ax.axvspan(s / win_size, e / win_size, color="red", alpha=0.25, zorder=3)
ax.set_xlabel(f"ungapped position ({win_size // 1000} kb windows)")
ax.set_yticks([])
ax.set_title(
    f"non-consensus sites per {win_size // 1000} kb window (masked regions in red)"
)
fig.colorbar(im, ax=ax, label=f"non-consensus sites / {win_size // 1000} kb")
plt.savefig(
    fig_fld / "non_consensus_heatmap_1kb_filtered.png",
    dpi=300,
    bbox_inches="tight",
    pad_inches=0.1,
)
plt.show()

# %%
