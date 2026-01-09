rule genomad_download_db:
    output:
        db=directory("data/genomad_db"),
    conda:
        "../config/conda_envs/genomad.yaml"
    shell:
        """
        genomad download-database data/
        """


rule genomad_run:
    input:
        db=rules.genomad_download_db.output.db,
        fa=rules.gbk_to_fa.output.fa,
    output:
        d=directory("data/genomad/{acc}"),
        s="data/genomad/{acc}/{acc}_summary/{acc}_virus_summary.tsv",
    log:
        "logs/genomad/{acc}.log",
    conda:
        "../config/conda_envs/genomad.yaml"
    shell:
        """
        genomad end-to-end {input.fa} {output.d} {input.db} \
            --cleanup \
            --threads 4 \
            &> {log}
        """


rule genomad_preformat:
    input:
        lambda w: expand(rules.genomad_run.output.s, acc=acc_nums),
    output:
        "results/mges/genomad.csv",
    conda:
        "../config/conda_envs/bioinfo.yaml"
    shell:
        """
        python3 scripts/genomad_df_preformat.py \
            --input_tsvs {input} \
            --output_df {output}
        """


rule mge_assign_positions:
    input:
        el="results/mges/{tool}.csv",
        j_pos=config["junction_positions_file"],
        iso_len=rules.genome_lengths.output,
    output:
        "results/mges_to_junctions/{tool}.csv",
    conda:
        "../config/conda_envs/bioinfo.yaml"
    params:
        zero_based="--zero_based",
        random="",
    shell:
        """
        python3 scripts/assign_junctions.py \
            --iso_len {input.iso_len} \
            --junction_pos_json {input.j_pos} \
            --element_pos_df {input.el} \
            --output_pos {output} \
            {params.zero_based} \
            {params.random}
        """


rule mge_all:
    input:
        expand(rules.mge_assign_positions.output, tool=["genomad"]),
