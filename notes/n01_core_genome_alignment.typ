#show link: underline
#show link: set text(fill: blue)
#set table(stroke: none, align: left)
#show figure.caption: set par(justify: true)
#show figure.caption: it => align(center, box(width: 80%, align(left, it)))

= Core genome alignment construction and analysis

== Pangenome graph construction

I reconstructed the pangenome graph of #link("https://github.com/mmolari/ecoliST131-structural-evo/blob/main/config/datasets/ST131_ABC/chromosomes.txt", [the pre-selected 222 _E.coli_ ST131 strains]) using pangraph, with the command and options:
```sh
pangraph --circular -k minimap2 -s 20 -a 100 -b 5 -l 100 genomes.fa -o graph.json
```

== Core genome alignment construction

In the original paper @molari2025quantifying, the core genome alignment had length of 3.59 bp. We used isolate `NZ_CP096110.1` as a guide strain to order and orient core-blocks.

With the new pangraph version, the total length of core blocks is the same: 3.59 Mbp.

For simplicity, I decided to *exclude core blocks shorter than 500bp* (too short to be considered valid flanking regions for junctions). This only slightly reduces the total alignments length to 3.57 Mbp.

#figure(image("assets/n01/block_length_distribution.png", width: 70%), caption: [
  Cumulative sum of block length distribution, stratified by core (orange) and accessory (blue).
])

For all the remaining _anchor blocks_ (i.e. core-blocks with length $>= 500$ bp), I export the full sequences using the `pangraph export block-sequences --unaligned` command, and then re-align each block using `mafft --auto`. This is done to ensure that the alignments are as accurate as possible.

I then concatenate all the core-block alignments in the order (and strandedness) of the guide strain (isolate `NZ_CP096110.1`, the same used in @molari2025quantifying), keeping track of the block boundaries. Note that this alignment also includes all gaps.

== Diversity along the alignment

As done for @molari2025quantifying, I compute the density of polymorphic sites along the core genome alignment using 10kb windows, see @snp_density_10kb. This confirms the background density of 10 polymorphic sites per 10kb window.

I also visualize the density of non-consensus alleles across the core genome alignment, see @non_consensus_heatmap_1kb. This highlights putative recombination events.

#figure(image("assets/n01/snp_density_10kb.png", width: 100%), caption: [
  density of polymorphic sites per 10kb window.
]) <snp_density_10kb>

#figure(image("assets/n01/non_consensus_heatmap_1kb.png", width: 100%), caption: [
  Non-consensus allele frequency across the core genome alignment, computed in 1kb windows.
]) <non_consensus_heatmap_1kb>

== Short deletions are more common than short insertions

I also tried to quantify the relative frequency of short insertions and deletions in the core block alignments.

To do this, I take all core-block alignments and first trim flanking gaps. This is done by selecting the longest starting/ending gap and removing all columns up to that gap. Notice that this procedure is not perfect, as it may fail to trim some flanking gaps if by any reason there are non-alignable regions at the start/end of the block. But visual inspection of the full core alignment suggests that these should be rare.

I then take all columns that contain a gap in the alignment, and compute the gap frequency, i.e. the fraction of isolates that have a gap in that column. This is reported in @gap_freq. The rationale is that "recent" short insertion should be detected as high-frequency gaps, while "recent" short deletions should be detected as low-frequency gaps. We observe an excess of low-frequency gaps (roughly 60% vs 40%) suggesting an excess of recent "small" deletions, indicating a possible bias towards *sequence erosion*.

#figure(image("assets/n01/interior_gap_frequency_ecdf.png", width: 70%), caption: [
  Cumulative distribution of gap frequencies across all columns of core-block alignments (flanking gap regions have been excluded). We observe an excess of columns with low-gap frequencies (60% vs 40%) suggesting an excess of recent "small" deletions.
]) <gap_freq>

#bibliography("bibliography.bib")
