# %%
import pypangraph as pp
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pathlib
import numpy as np
from Bio import AlignIO, Phylo


fig_fld = pathlib.Path("figs/n01_core_alignment")
fig_fld.mkdir(exist_ok=True, parents=True)

# %%
pan = pp.Pangraph.from_json("../results/pangraph.json")
bdf = pan.to_blockstats_df()

# %%
# total core block length
core_length = bdf[bdf["core"]]["len"].sum()
print(f"Total core block length: {core_length / 1e6:.2f} Mb")

# distributon of core-block lengths
sns.histplot(
    data=bdf,
    x="len",
    weights="len",
    bins=500,
    hue="core",
    element="step",
    fill=False,
    cumulative=True,
)
plt.axhline(
    core_length,
    color="k",
    linestyle="--",
)
# annotate total core-genome size above the horizontal line
plt.text(
    bdf["len"].max(),
    core_length * 1.02,
    f"core: {core_length / 1e6:.2f} Mb",
    ha="right",
    va="bottom",
    fontsize=8,
)
plt.xscale("log")
plt.ylabel("cumulative block length (bp)")
plt.xlabel("block length (log scale)")
plt.xlim(bdf["len"].min(), bdf["len"].max() * 1.01)
plt.tight_layout()
plt.savefig(fig_fld / "block_length_distribution.png", dpi=300)
plt.show()


# %%


aln_path = "../results/core_genome_alignment.fa"
tree_path = "../config/polished_tree.nwk"

# %%
# Load the concatenated core-genome alignment as a (N_iso, L) uint8
# array. uint8 instead of <U1 (UCS-4) means 4x less memory and much
# faster vectorised comparisons downstream.
aln = AlignIO.read(aln_path, "fasta")
strain_order = [rec.id for rec in aln]
n_iso, aln_length = len(aln), aln.get_alignment_length()
arr = np.empty((n_iso, aln_length), dtype=np.uint8)
for i, rec in enumerate(aln):
    arr[i] = np.frombuffer(str(rec.seq).encode("ascii"), dtype=np.uint8)
print(f"Alignment: {arr.shape[0]} isolates x {arr.shape[1]:,} columns")

# %%
# Drop every column that contains at least one gap. Coordinates after
# this step are "gap-free coordinates" — used for both panels below.
GAP = ord("-")  # 45
gap_free_mask = ~(arr == GAP).any(axis=0)
arr_gf = arr[:, gap_free_mask]
print(
    f"Gap-free alignment: {arr_gf.shape[1]:,} columns "
    f"({arr_gf.shape[1] / arr.shape[1] * 100:.2f}% of original)"
)

# %%
# Per-column consensus (majority base) and a (N, L) boolean matrix of
# non-consensus cells. Anything outside ACGT (e.g. N) gets flagged as
# non-consensus, which is what we want here.
bases_u8 = np.array([ord(b) for b in "acgt"], dtype=np.uint8)
counts = np.stack([(arr_gf == b).sum(axis=0) for b in bases_u8])  # (4, L)
consensus = bases_u8[counts.argmax(axis=0)]  # (L,)
non_consensus = arr_gf != consensus[None, :]  # (N, L) bool
polymorphic_col = non_consensus.any(axis=0)  # (L,)
print(
    f"Polymorphic columns: {polymorphic_col.sum():,} "
    f"({polymorphic_col.mean() * 100:.3f}% of gap-free)"
)

# %%
# SNP density along the gap-free genome: histogram of polymorphic-column
# positions in 10 kb bins. y-log so the rare hot bins are visible.
poly_positions = np.where(polymorphic_col)[0]
bin_size = 10_000
bins = np.arange(0, arr_gf.shape[1] + bin_size, bin_size)
fig, ax = plt.subplots(figsize=(10, 3))
ax.hist(poly_positions, bins=bins)
ax.set_yscale("log")
ax.set_xlabel("gap-free alignment position (bp)")
ax.set_ylabel("polymorphic columns per 10 kb (log)")
ax.set_title(
    f"{polymorphic_col.sum():,} polymorphic / {arr_gf.shape[1]:,} gap-free cols"
)
plt.tight_layout()
plt.savefig(fig_fld / "snp_density_10kb.png", dpi=300)
plt.show()

# %%
# Reorder rows by polished_tree leaf order.
tree = Phylo.read(tree_path, "newick")
leaf_order = [leaf.name for leaf in tree.get_terminals()]
strain_to_idx = {s: i for i, s in enumerate(strain_order)}
missing = [s for s in leaf_order if s not in strain_to_idx]
if missing:
    print(f"WARNING: {len(missing)} tree leaves not in alignment: {missing[:5]}...")
row_order = [strain_to_idx[s] for s in leaf_order if s in strain_to_idx]

# %%


# Per-isolate non-consensus density in 1 kb windows, tree-ordered.
# Drop trailing cols so length is divisible by win_size, then reshape
# and sum within each window. Result is (N_iso, N_windows) small-int
# counts — cheap to imshow with a sequential colormap.
win_size = 1_000
L = arr_gf.shape[1]
L_trim = (L // win_size) * win_size
nc_windows = (
    non_consensus[:, :L_trim]
    .reshape(arr_gf.shape[0], L_trim // win_size, win_size)
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

# Left panel: phylogeny. label_func returns "" to suppress per-leaf
# text (222 names would be unreadable). y-limits below force the tree
# to span exactly the same N "slots" as the heatmap so rows align.
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

# Right panel: heatmap. Rows are tree-ordered via row_order above, so
# the top row of the heatmap corresponds to the top leaf of the tree.
vmax = max(2, int(nc_windows.max()))
im = ax.imshow(
    nc_windows,
    aspect="auto",
    interpolation="nearest",
    cmap="Greys",
    norm=LogNorm(vmin=0.4, vmax=vmax),
)
ax.set_xlabel(f"gap-free position ({win_size // 1000} kb windows)")
ax.set_yticks([])
ax.set_title(f"non-consensus sites per {win_size // 1000} kb window")
fig.colorbar(im, ax=ax, label=f"non-consensus sites / {win_size // 1000} kb")
plt.savefig(
    fig_fld / "non_consensus_heatmap_1kb.png",
    dpi=300,
    bbox_inches="tight",
    pad_inches=0.1,
)
plt.show()

# %%
