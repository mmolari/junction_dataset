# %%
import pathlib

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


fig_fld = pathlib.Path("figs/n07_junction_scatter")
fig_fld.mkdir(exist_ok=True, parents=True)

# uniform +/- jitter on the (integer) category count, applied before the
# log transform so points with the same n_categories don't fully overlap.
JITTER = 0.25
np.random.seed(0)  # reproducible jitter

# %%
df = pd.read_csv("../results/junction_stats.csv")

# log x-axis can't show tot_acc_len == 0 (transitive / core-only junctions)
n_zero = int((df["tot_acc_len"] == 0).sum())
df = df[df["tot_acc_len"] > 0].copy()
print(f"Dropped {n_zero} junctions with zero accessory length; {len(df)} remain")

# %%
y = df["n_categories"] + np.random.uniform(-JITTER, JITTER, size=len(df))

fig, ax = plt.subplots(figsize=(7, 5))
sc = ax.scatter(
    df["tot_acc_len"],
    y,
    c=df["nonempty_freq"],
    cmap="coolwarm",
    vmin=0,
    vmax=1,
    s=15,
    alpha=0.8,
    edgecolor="none",
)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("total accessory length (bp)")
ax.set_ylabel("n. path categories")
fig.colorbar(sc, ax=ax, label="fraction of isolates with accessory content")
plt.tight_layout()
plt.savefig(fig_fld / "acc_len_vs_n_categories.png", dpi=300)
plt.show()

# %%
