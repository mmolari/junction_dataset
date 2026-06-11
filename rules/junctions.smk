# This snakemake file contains rules to extract variable junctions from the full pangenome graph,
# and build secondary junction graphs


checkpoint junction_positions:
    input:
        graph=rules.build_pangraph.output,
    output:
        csv="results/junction_positions.csv",
    log:
        "logs/junction_positions.log",
    conda:
        "../config/conda_envs/bioinfo.yaml"
    params:
        min_length=config["core_block_min_length"],
    shell:
        """
        python scripts/junction_positions.py \
            --pangraph {input.graph} \
            --min_length {params.min_length} \
            --out_csv {output.csv} \
            &>{log}
        """


def junction_ids(wildcards):
    # junction IDs are derived dynamically from the pangraph; keep only junctions
    # present in more than one isolate (single-isolate junctions are degenerate).
    csv = checkpoints.junction_positions.get().output.csv
    df = pd.read_csv(csv)
    counts = df.groupby("edge")["iso"].nunique()
    return counts.index[counts > 1].tolist()


def junction_pangraphs(wildcards):
    return expand(rules.build_junction_pangraph.output, junc=junction_ids(wildcards))


rule extract_junction_sequences:
    input:
        gbk=expand(rules.download_gbk.output, acc=acc_nums),
        j_pos=rules.junction_positions.output.csv,
    output:
        fa="results/junction_sequences/{junc}.fa",
        gff="results/junction_annotations/{junc}.gff",
    log:
        "logs/extract_junctions/{junc}.log",
    conda:
        "../config/conda_envs/bioinfo.yaml"
    shell:
        """
        python scripts/extract_junctions.py \
            --gbk-fld "data/gbk" \
            --junc-id {wildcards.junc} \
            --junc-pos-file {input.j_pos} \
            --out-fa {output.fa} \
            --out-ann {output.gff} \
            &>{log}
        """


rule build_junction_pangraph:
    input:
        fa=rules.extract_junction_sequences.output.fa,
    output:
        "results/junction_pangraphs/{junc}.json",
    params:
        opt="-s 20 -a 100 -b 5 -l 100",
    shell:
        """
        pangraph build {input.fa} {params.opt} -o {output}
        """


rule junction_stats:
    input:
        pangraph=junction_pangraphs,
    output:
        "results/junction_stats.csv",
    conda:
        "../config/conda_envs/bioinfo.yaml"
    shell:
        """
        python scripts/junction_stats.py \
            --junct_pangraphs {input.pangraph} \
            --df_csv {output}
        """
