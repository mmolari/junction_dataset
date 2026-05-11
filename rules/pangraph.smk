# This snakemake file contains rules to build the pangenome graph for all the input genomes,
# extract a core-genome alignment and build a phylogenetic tree


rule build_pangraph:
    input:
        expand(rules.gbk_to_fa.output, acc=acc_nums),
    output:
        "results/pangraph.json.gz",
    log:
        "logs/build_pangraph.log",
    threads: 8
    params:
        opt="--circular -k minimap2 -s 20 -a 100 -b 5 -l 100",
    shell:
        """
        pangraph build {input} {params.opt} -j {threads} -o {output} -v &> {log}
        """


checkpoint block_stats:
    input:
        graph="results/pangraph.json.gz",
    output:
        stats="results/pangraph_block_stats.csv",
    log:
        "logs/block_stats.log",
    conda:
        "../config/conda_envs/bioinfo.yaml"
    shell:
        """
        python -c "import pypangraph as pp; pan = pp.Pangraph.from_json('{input.graph}'); bdf = pan.to_blockstats_df(); bdf.to_csv('{output.stats}', index_label='id')" &> {log}
        """


checkpoint export_block_sequences:
    input:
        "results/pangraph.json.gz",
    output:
        directory("results/pangraph_block_sequences"),
    log:
        "logs/export_block_sequences.log",
    shell:
        """
        pangraph export block-sequences {input} -o {output} --unaligned -v &> {log}
        """


rule align_core_block:
    input:
        "results/pangraph_block_sequences/block_{block_id}.fa",
    output:
        "results/pangraph_core_block_alignments/block_{block_id}.fa",
    log:
        "logs/align_core_block/{block_id}.log",
    conda:
        "../config/conda_envs/mafft.yaml"
    shell:
        """
        mafft --auto {input} > {output} 2> {log}
        """


def core_block_alignments(wildcards):
    checkpoints.export_block_sequences.get()
    stats_file = checkpoints.block_stats.get().output.stats
    import pandas as pd

    stats = pd.read_csv(stats_file)
    min_len = config["core_block_min_length"]
    mask = stats["core"] & (stats["len"] >= min_len)
    block_ids = stats.loc[mask, "id"].astype(str).tolist()
    return expand(
        "results/pangraph_core_block_alignments/block_{block_id}.fa",
        block_id=block_ids,
    )


rule align_all_core_blocks:
    input:
        core_block_alignments,


rule pangraph_all:
    input:
        rules.export_block_sequences.output,
        rules.block_stats.output,
        rules.align_all_core_blocks.input,
