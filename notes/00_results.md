# `results/` contents

One-line description of each top-level file and folder produced by the pipeline.
(`results/` is gitignored — these are pipeline outputs, not source.)

## Files

| Item | Description |
|------|-------------|
| `pangraph.json` | All-isolates PanGraph built over every isolate FASTA (`--circular`); the master graph all downstream core/junction analyses derive from. |
| `pangraph_block_stats.csv` | Per-block stats of the all-isolates pangraph (`id, count, n_strains, duplicated, core, len`). |
| `junction_positions.csv` | Junction coordinates per isolate from the `junction_positions` checkpoint (`edge, iso, left_start, left_end, right_start, right_end, strand`). |
| `junction_stats.csv` | Per-junction summary statistics (path categories, entropy, block/length metrics, core flank lengths, accessory content). |
| `genome_lengths.csv` | Total sequence length per isolate accession (`accession, length`); used for MGE position assignment. |

## Folders

| Item | Description |
|------|-------------|
| `pangraph_block_sequences/` | Per-block unaligned FASTAs exported from the all-isolates pangraph (~6.3k `block_<id>.fa` files). |
| `pangraph_core_block_alignments/` | Per-core-block MAFFT (`--auto`) alignments (~1.3k `block_<id>.fa` files). |
| `core_genome_alignment/` | Concatenated core-genome alignment in guide-strain order — `gapped.fa`/`ungapped.fa` plus matching per-block coordinate CSVs. |
| `recombination_filter/` | SNP-density recombination-masked alignment (`alignment.fa`), the `masked_intervals.csv` it dropped, and a quick FastTree (`fasttree.nwk`). |
| `gubbins/` | Full Gubbins recombination-aware ML run (final tree, recombination predictions, filtered polymorphic sites, per-branch stats, logs). |
| `core_trees/` | TreeTime branch-length-refined core-genome tree(s) (`polished_from_filtering.nwk`). |
| `junction_sequences/` | Per-junction extracted region FASTA, one file per junction (~549 `<junc_id>.fa`). |
| `junction_annotations/` | Per-junction feature annotations as GFF, one file per junction (~549 `<junc_id>.gff`). |
| `junction_pangraphs/` | Per-junction PanGraph JSON built on the extracted region, one file per junction (~549 `<junc_id>.json`). |
| `mges/` | Collapsed per-isolate MGE calls into single CSVs — `genomad.csv` and `defensefinder.csv` (`id, iso, beg, end, type`). |
| `mges_to_junctions/` | MGE hits mapped onto junctions (`genomad.csv`, `defensefinder.csv`), with overlap classification and junction coordinates. |
| `plasmids/` | Per-isolate plasmid GenBank copies, one subfolder per host accession (~212 dirs of `<plasmid_acc>.gbk`). |
