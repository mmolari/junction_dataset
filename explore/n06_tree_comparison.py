# %%
import itertools
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Phylo
from matplotlib.ticker import MaxNLocator
from scipy.stats import pearsonr


fig_fld = pathlib.Path("figs/n06_tree_comparison")
fig_fld.mkdir(exist_ok=True, parents=True)

tree_paths = {
    "original": "../config/polished_tree.nwk",
    "filter": "../results/core_trees/polished_from_filtering.nwk",
    "gubbins": "../results/gubbins/gubbins.final_tree.tre",
}

# %%
# Load the three trees and intersect their tip sets. Trees from different
# pipelines occasionally diverge by a leaf or two; we restrict every
# downstream comparison (tip distances + RF) to the common set.
trees = {name: Phylo.read(p, "newick") for name, p in tree_paths.items()}
tip_sets = {
    name: {leaf.name for leaf in t.get_terminals()} for name, t in trees.items()
}
for name, s in tip_sets.items():
    print(f"{name}: {len(s)} tips")
common_tips = sorted(set.intersection(*tip_sets.values()))
print(f"common: {len(common_tips)} tips")
for name, s in tip_sets.items():
    diff = s - set(common_tips)
    if diff:
        print(f"  {name} drops {len(diff)} tips: {sorted(diff)[:5]}...")


# %%
# Patristic distance matrix per tree, restricted to common_tips and saved
# square-form so we can re-load without re-walking the tree.
def pairwise_distances(tree, ordered_tips):
    tip_by_name = {leaf.name: leaf for leaf in tree.get_terminals()}
    n = len(ordered_tips)
    D = np.zeros((n, n), dtype=float)
    for i, j in itertools.combinations(range(n), 2):
        d = tree.distance(tip_by_name[ordered_tips[i]], tip_by_name[ordered_tips[j]])
        D[i, j] = d
        D[j, i] = d
    return D


D = {}
for name, tree in trees.items():
    print(f"computing pairwise distances: {name}")
    D[name] = pairwise_distances(tree, common_tips)
    # pd.DataFrame(D[name], index=common_tips, columns=common_tips).to_csv(
    #     fig_fld / f"tip_distances_{name}.csv"
    # )

# %%
# Scatter pairwise distances tree-A vs tree-B for the three pairs.
# Branch-length units differ across trees (substitutions/site vs SNPs),
# so a y=x reference is not meaningful — annotate Pearson + Spearman r
# in each title instead.
pairs = [("original", "filter"), ("original", "gubbins"), ("filter", "gubbins")]
iu = np.triu_indices(len(common_tips), k=1)

fig, axes = plt.subplots(1, 3, figsize=(12, 4))
for ax, (a, b) in zip(axes, pairs):
    xa = D[a][iu]
    xb = D[b][iu]
    r_p, _ = pearsonr(xa, xb)
    # least-squares fit y = m*x (through origin): m = sum(xy) / sum(x^2)
    slope = float(np.dot(xa, xb) / np.dot(xa, xa))
    ax.scatter(xa, xb, s=4, alpha=0.2, color="C0", rasterized=True)
    x_line = np.array([0.0, xa.max()])
    ax.plot(x_line, slope * x_line, color="C3", lw=1, label=f"y = {slope:.3g} x")
    ax.legend(loc="upper left")
    ax.set_xlabel(f"{a} tip distance")
    ax.set_ylabel(f"{b} tip distance")
    ax.set_title(f"{a} vs {b}\nPearson r={r_p:.3f}")
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.yaxis.set_major_locator(MaxNLocator(4))
    ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(fig_fld / "tip_distance_scatter.png", dpi=300)
plt.show()


# %%
# Robinson–Foulds via bipartition sets.
# Each non-terminal clade defines a bipartition (the leaf set under it, vs
# the rest of the common tips). We canonicalize by the smaller side and
# drop trivial bipartitions (empty / full / single leaf).
def bipartitions(tree, all_tips):
    all_tips = set(all_tips)
    bps = set()
    for clade in tree.get_nonterminals():
        side = frozenset(leaf.name for leaf in clade.get_terminals()) & all_tips
        if len(side) < 2 or len(side) >= len(all_tips) - 1:
            continue
        other = frozenset(all_tips - side)
        bps.add(side if len(side) < len(other) else other)
    return bps


bps = {name: bipartitions(t, common_tips) for name, t in trees.items()}
for name, s in bps.items():
    print(f"{name}: {len(s)} non-trivial bipartitions")

rows = []
for a, b in pairs:
    rf = len(bps[a].symmetric_difference(bps[b]))
    rows.append(
        {
            "pair": f"{a}_vs_{b}",
            "rf": rf,
            f"n_bipartitions_{a}": len(bps[a]),
            f"n_bipartitions_{b}": len(bps[b]),
        }
    )
rf_df = pd.DataFrame(rows).set_index("pair")
print(rf_df)
rf_df.to_csv(fig_fld / "rf_distances.csv")

# %%
# RF matrix as a heatmap for a quick visual.
names = list(trees.keys())
M = np.zeros((len(names), len(names)), dtype=int)
for i, a in enumerate(names):
    for j, b in enumerate(names):
        M[i, j] = len(bps[a].symmetric_difference(bps[b]))

fig, ax = plt.subplots(figsize=(4.5, 4))
im = ax.imshow(M, cmap="viridis")
ax.set_xticks(range(len(names)), names)
ax.set_yticks(range(len(names)), names)
for i in range(len(names)):
    for j in range(len(names)):
        ax.text(
            j,
            i,
            str(M[i, j]),
            ha="center",
            va="center",
            color="white" if M[i, j] < M.max() / 2 else "black",
        )
ax.set_title("Robinson–Foulds distance")
fig.colorbar(im, ax=ax, label="RF")
plt.tight_layout()
plt.savefig(fig_fld / "rf_matrix.png", dpi=300)
plt.show()

# %%
