rule GM_download_db:
    output:
        db=directory("data/genomad_db"),
    conda:
        "../config/conda_envs/genomad.yaml"
    shell:
        """
        genomad download-database data/
        """


rule GM_run:
    input:
        db=rules.GM_download_db.output.db,
        fa=rules.gbk_to_fa.output.fa,
    output:
        d=directory("data/genomad/{acc}"),
        s="data/genomad/{acc}/{acc}_summary/{acc}_virus_summary.tsv",
    conda:
        "../config/conda_envs/genomad.yaml"
    shell:
        """
        genomad end-to-end {input.fa} {output.d} {input.db} \
            --cleanup \
            --threads 4
        """


rule GM_preformat:
    input:
        lambda w: expand(rules.GM_run.output.s, acc=acc_nums),
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
