# %%
import pathlib
import subprocess

import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO

GUIDE = "NZ_CP096110.1"
ALN_PATH = pathlib.Path("../results/core_genome_alignment.fa")
GENOME_PATH = pathlib.Path(f"../data/fasta/{GUIDE}.fa")

fig_fld = pathlib.Path("figs/n03_guide_remap")
fig_fld.mkdir(exist_ok=True, parents=True)

# %%
# Stream-pull the guide row out of the 771 MB core-genome alignment so
# we don't load all 222 rows. Strip alignment gaps to recover what
# should be a verbatim, ordered concatenation of guide-genome substrings.
guide_seq = None
for rec in SeqIO.parse(ALN_PATH, "fasta"):
    if rec.id == GUIDE:
        guide_seq = str(rec.seq).replace("-", "").upper()
        break
assert guide_seq is not None, f"{GUIDE} not found in {ALN_PATH}"
print(f"Gap-stripped guide row: {len(guide_seq):,} bp")

# %%
# Map the gap-stripped row back onto the full guide genome with
# minimap2. asm5 = lowest-divergence asm preset (query and target are
# the same DNA up to the chop-and-glue), --secondary=no keeps the PAF
# to one row per region.
query_fa = fig_fld / "guide_coreseq.fa"
paf_path = fig_fld / "guide_coreseq.paf"
with open(query_fa, "w") as fh:
    fh.write(f">{GUIDE}_coreseq\n{guide_seq}\n")

subprocess.run(
    [
        "conda",
        "run",
        "-n",
        "minimap2",
        "minimap2",
        # "-x",
        # "asm5",
        "-r",
        "50",
        "-g",
        "50",
        "-c",
        "--secondary=no",
        str(GENOME_PATH),
        str(query_fa),
        "-o",
        str(paf_path),
    ],
    check=True,
)
print(f"PAF written: {paf_path}")

# %%
PAF_COLS = [
    "q_name",
    "q_len",
    "q_start",
    "q_end",
    "strand",
    "t_name",
    "t_len",
    "t_start",
    "t_end",
    "n_match",
    "aln_len",
    "mapq",
]
paf = pd.read_csv(paf_path, sep="\t", header=None, usecols=range(12), names=PAF_COLS)
paf = paf.sort_values("q_start").reset_index(drop=True)
print(paf.head())
print(
    f"\n{len(paf)} hits, total query covered: "
    f"{(paf['q_end'] - paf['q_start']).sum():,} / {len(guide_seq):,}"
)

# %%
# Verification. Minimap2 is a heuristic mapper, so even on identical
# DNA short or repetitive stretches near block boundaries can be left
# unanchored — we report query coverage % separately rather than
# demanding q_gap == 0. The structural checks are:
#   - every hit on + strand,
#   - query hits non-overlapping (q_gap >= 0),
#   - target moves forward, with one allowed wrap for the (at most one)
#     core block that straddles the chromosome origin: prev hit ends at
#     t_len, next hit starts at 0.
all_forward = (paf["strand"] == "+").all()
q_gaps = paf["q_start"].values[1:] - paf["q_end"].values[:-1]
t_gaps = paf["t_start"].values[1:] - paf["t_end"].values[:-1]
t_len = int(paf["t_len"].iloc[0])
q_covered = int((paf["q_end"] - paf["q_start"]).sum())
non_overlapping = bool((q_gaps >= 0).all())
q_uncovered = int(q_gaps[q_gaps > 0].sum())

neg_t = (t_gaps < 0).nonzero()[0]
if len(neg_t) == 0:
    target_ok = bool((t_gaps >= 0).all())
    wrap_msg = "no origin wrap"
elif len(neg_t) == 1:
    i = int(neg_t[0])
    wrap_ok = paf["t_end"].iloc[i] == t_len and paf["t_start"].iloc[i + 1] == 0
    target_ok = wrap_ok
    wrap_msg = (
        f"origin wrap between hits {i} and {i + 1}: "
        f"t_end={paf['t_end'].iloc[i]} (t_len={t_len}), "
        f"t_start={paf['t_start'].iloc[i + 1]}"
    )
else:
    target_ok = False
    wrap_msg = f"{len(neg_t)} negative t_gaps — too many for a single origin wrap"

print(f"all forward strand?               {all_forward}")
print(
    f"query non-overlapping (q_gap>=0)? {non_overlapping}  "
    f"(min/max q_gap = {q_gaps.min()}/{q_gaps.max()})"
)
print(f"target forward (origin wrap OK)?  {target_ok}  ({wrap_msg})")
print(
    f"query coverage:                   {q_covered:,} / {len(guide_seq):,} "
    f"({q_covered / len(guide_seq) * 100:.3f}%, "
    f"{q_uncovered:,} bp in inter-hit gaps)"
)

# %%
# Dot plot: each hit as a segment (t_start, q_start) -> (t_end, q_end).
# On a same-strand consecutive remap this is a strict staircase along
# the diagonal; reverse-strand hits would have negative slope and are
# drawn in red so they pop visually.
fig, ax = plt.subplots(figsize=(7, 7))
for _, row in paf.iterrows():
    if row["strand"] == "+":
        ax.plot(
            [row["t_start"], row["t_end"]],
            [row["q_start"], row["q_end"]],
            color="C0",
            lw=1.5,
        )
    else:
        ax.plot(
            [row["t_start"], row["t_end"]],
            [row["q_end"], row["q_start"]],
            color="C3",
            lw=1.5,
        )
ax.set_xlabel(f"target position on {GUIDE} (bp)")
ax.set_ylabel("query position on gap-stripped core row (bp)")
ax.set_title(f"{len(paf)} minimap2 hits  |  asm5  |  guide core row vs full genome")
ax.set_xlim(0, paf["t_len"].iloc[0])
ax.set_ylim(0, len(guide_seq))
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(fig_fld / "dotplot.png", dpi=300)
plt.show()

# %%
