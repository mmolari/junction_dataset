# This snakemake file contains rules to build the pangenome graph for all the input genomes,
# extract a core-genome alignment and build a phylogenetic tree


rule build_pangraph:
    input:
        expand(rules.gbk_to_fa.output, acc=acc_nums),
    output:
        "results/pangraph.json",
    log:
        "logs/build_pangraph.log",
    threads: 8
    params:
        opt="--circular -k minimap2 -s 20 -a 100 -b 5 -l 100",
    shell:
        """
        pangraph build {input} {params.opt} -j {threads} -o {output} -v &> {log}
        """


rule pangraph_all:
    input:
        rules.build_pangraph.output,
