# This snakemake file contains rules to extract variable junctions from the full pangenome graph,
# and build secondary junction graphs

# load the table of junction positions
with open(config["junction_positions_file"]) as f:
    junc_pos = json.load(f)
# exclude single-isolate junctions
junc_pos = {j: p for j, p in junc_pos.items() if len(p) > 1}
# list of junction IDs
junc_ids = list(junc_pos.keys())
print(f"N. junctions: {len(junc_ids)}")


rule extract_junction_sequences:
    input:
        gbk=expand(rules.download_gbk.output, acc=acc_nums),
        j_pos=config["junction_positions_file"],
    output:
        fa="results/junction_sequences/{junc}.fa",
        gff="results/junction_annotations/{junc}.gff",
    log:
        "logs/extract_junctions/{junc}.log",
    conda:
        "config/conda_envs/bioinfo.yaml"
    shell:
        """
        python scripts/extract_junctions.py \
            --gbk-fld "data/gbk" \
            --junc-id {wildcards.junc} \
            --junc-pos-file {input.j_pos} \
            --out-fa {output.fa} \
            --out-ann {output.gff} \
            &> {log}
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
        pangraph=expand(rules.build_junction_pangraph.output, junc=junc_ids),
    output:
        "results/junction_stats.csv",
    conda:
        "config/conda_envs/bioinfo.yaml"
    shell:
        """
        python scripts/junction_stats.py \
            --junct_pangraphs {input.pangraph} \
            --df_csv {output}
        """
