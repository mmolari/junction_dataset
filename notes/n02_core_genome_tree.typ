#import "@preview/cheq:0.3.1": checklist
#show: checklist

#show link: underline
#show link: set text(fill: blue)
#set table(stroke: none, align: left)
#show figure.caption: set par(justify: true)
#show figure.caption: it => align(center, box(width: 80%, align(left, it)))

= Core-genome trees

I build three core-genome trees for the 222 ST131 _E. coli_ isolates and plan to compare them. All trees use the same isolate set as the alignment described in note `n01`.

== Tree from the published paper

The reference tree from @molari2025quantifying is included in the repository at `config/polished_tree.nwk` and is used throughout the pipeline as the canonical phylogeny (e.g. as leaf order in the non-consensus heatmaps). It was inferred on a recombination-filtered core-genome alignment of the same 222 isolates, with branch lengths polished by `treetime`. I treat this tree as the reference when comparing the trees below.

== Tree from a custom recombination-filtering procedure

I implemented a lightweight SNP-density recombination filter in `scripts/recombination_filter.py` (rule `filter_core_alignment`). The intent is to reproduce — in a simpler form — the masking step used in @molari2025quantifying.

*Procedure.* Starting from the ungapped core-genome alignment (3.56 Mbp, see note `n01`), I compute the column-wise consensus and classify every allele that differs from the consensus as a "non-consensus" site. For each isolate independently, I slide a _circular_ window of half-width $w = 1000$ bp along the genome: any focal non-consensus site whose $2 w + 1$ bp window contains more than `max_snps = 3` non-consensus sites (focal included) triggers masking of its full $[p - w, p + w]$ interval (2001 bp). The filter is *per-isolate* in detection but produces a *single global* mask: any column flagged for any isolate is removed from the alignment for every isolate.

*Effect of the filter* is visible in @snp_density_10kb_filtered: the per-10 kbp histogram of polymorphic sites collapses cleanly onto a background rate, with the hot regions removed. @non_consensus_heatmap_1kb_filtered overlays the masked intervals (red) on a per-isolate non-consensus heatmap, showing that the filter targets the recombination hotspots identified visually in note `n01`.

#figure(
  image("assets/n02/snp_density_10kb_filtered.png", width: 100%),
  caption: [Polymorphic columns per 10 kbp window (log scale). Grey: all polymorphic columns in the ungapped alignment. Blue: only those surviving the filter. The surviving signal reflects the background mutational density.],
) <snp_density_10kb_filtered>

#figure(
  image("assets/n02/non_consensus_heatmap_1kb_filtered.png", width: 100%),
  caption: [Per-isolate density of non-consensus sites in 1 kbp windows, rows ordered by `polished_tree.nwk`. Masked intervals are overlaid in red.],
) <non_consensus_heatmap_1kb_filtered>

*Tree inference.* From the filtered alignment I build a quick ML tree with FastTree (`-gtr -nt` options) and then refine its branch lengths with `treetime` against the same alignment.

== Tree from Gubbins

As an independent recombination-aware baseline, I also run Gubbins @croucher2015rapid on the *unfiltered* ungapped core alignment (`run_gubbins`). Gubbins iterates between ML tree inference (RAxML on the polymorphic-site sub-alignment) and per-branch detection of recombination tracts, masking the latter at each iteration until convergence. The final tree (`results/gubbins/gubbins.final_tree.tre`) has branch lengths in substitutions, since RAxML sees only the `filtered_polymorphic_sites.fasta` alignment. This makes it also hard to rescale branches by alignment length, since different branches have different unfiltered support length in the alignment.

== Tree comparison

After this procedure I have these three trees:
- The _"original"_ one from @molari2025quantifying.
- The _"filtered"_ from the adapted filtering procedure.
- the _"gubbins"_ one.

I compare them in two ways: pairwise tip-to-tip distances (a quantitative, branch-length view) and tanglegrams (a topological view).

*Tip-to-tip distances.* For every pair of trees I compute tip-to-tip distances between all unique pairs of leaves, and scatter them one tree against the other (see @tip_distance_scatter). The two recombination-filtering procedures yield essentially the same metric tree. The Gubbins branch lengths are in substitution (SNP) counts on the polymorphic-site sub-alignment rather than substitutions per site, so direct comparisons are harder, but the fact that the cloud is not well-captured by a single slope indicates that gubbin's masking is not uniform across all branches.

*Topology.* @tanglegrams shows TreeKnit tanglegrams, with colors indicating MCCs (Maximally Compatible Clades).

#figure(
  image("assets/n02/tip_distance_scatter.png", width: 100%),
  caption: [Pairwise tip-to-tip distances for the three tree pairs.],
) <tip_distance_scatter>

#figure(
  grid(
    columns: 2,
    gutter: 3pt,
    [#align(center)[*filtering vs original*] #v(-3mm) #image("assets/n02/tanglegram_original_vs_filter.png")],
    [#align(center)[*Gubbins vs original*] #v(-3mm) #image("assets/n02/tanglegram_original_vs_gubbins.png")],

    [#align(center)[*Gubbins vs filtered*] #v(-3mm) #image("assets/n02/tanglegram_filter_vs_gubbins.png")],
  ),
  caption: [TreeKnit tanglegrams for the three tree pairs. Colors indicate MCCs (Maximally Compatible Clades).],
) <tanglegrams>

#bibliography("bibliography.bib")
