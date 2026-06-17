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
            &>{log}
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
        python3 scripts/annotations/genomad_df_preformat.py \
            --input_tsvs {input} \
            --output_df {output}
        """


rule defensefinder_models_download:
    output:
        directory("data/defensefinder_models"),
    conda:
        "../config/conda_envs/defensefinder.yaml"
    shell:
        """
        TMPDIR=$(mktemp -d -t defensefinder_model_download_XXXXXXXXX)
        echo "created temporary directory $TMPDIR"
        defense-finder update --models-dir {output}
        rm -r $TMPDIR
        """


rule defensefinder_find:
    input:
        fa=rules.gbk_to_fa.output.fa,
        mod=rules.defensefinder_models_download.output,
    output:
        a=directory("data/defense_finder/{acc}"),
        g="data/defense_finder/{acc}/{acc}_defense_finder_genes.tsv",
        s="data/defense_finder/{acc}/{acc}_defense_finder_systems.tsv",
        p="data/defense_finder/{acc}/{acc}.prt",
    conda:
        "../config/conda_envs/defensefinder.yaml"
    shell:
        """
        defense-finder run \
            -o {output.a} \
            --models-dir {input.mod} \
            --skip-model-version-check \
            {input.fa}
        """


rule defensefinder_gene_location:
    input:
        g=rules.defensefinder_find.output.g,
        p=rules.defensefinder_find.output.p,
    output:
        temp("data/defense_finder/{acc}/{acc}_genes_loc.tsv"),
    conda:
        "../config/conda_envs/bioinfo.yaml"
    shell:
        """
        python3 scripts/annotations/defensefinder_gene_location.py \
            --input_gene_df {input.g} \
            --proteins {input.p} \
            --output_gene_df {output}
        """


rule defensefinder_preformat:
    input:
        s=expand(rules.defensefinder_find.output.s, acc=acc_nums),
        g=expand(rules.defensefinder_gene_location.output, acc=acc_nums),
    output:
        "results/mges/defensefinder.csv",
    conda:
        "../config/conda_envs/bioinfo.yaml"
    shell:
        """
        python3 scripts/annotations/defensefinder_df_preformat.py \
            --input_genes {input.g} \
            --input_systems {input.s} \
            --output_df {output}
        """


rule antidefense_find:
    input:
        fa=rules.gbk_to_fa.output.fa,
        mod=rules.defensefinder_models_download.output,
    output:
        a=directory("data/antidefense_finder/{acc}"),
        g="data/antidefense_finder/{acc}/{acc}_defense_finder_genes.tsv",
        s="data/antidefense_finder/{acc}/{acc}_defense_finder_systems.tsv",
        p="data/antidefense_finder/{acc}/{acc}.prt",
    conda:
        "../config/conda_envs/defensefinder.yaml"
    shell:
        """
        defense-finder run \
            -o {output.a} \
            --models-dir {input.mod} \
            --skip-model-version-check \
            --antidefensefinder-only \
            {input.fa}
        """


rule antidefense_gene_location:
    input:
        g=rules.antidefense_find.output.g,
        p=rules.antidefense_find.output.p,
    output:
        temp("data/antidefense_finder/{acc}/{acc}_genes_loc.tsv"),
    conda:
        "../config/conda_envs/bioinfo.yaml"
    shell:
        """
        python3 scripts/annotations/defensefinder_gene_location.py \
            --input_gene_df {input.g} \
            --proteins {input.p} \
            --output_gene_df {output}
        """


rule antidefense_preformat:
    input:
        s=expand(rules.antidefense_find.output.s, acc=acc_nums),
        g=expand(rules.antidefense_gene_location.output, acc=acc_nums),
    output:
        "results/mges/antidefense.csv",
    conda:
        "../config/conda_envs/bioinfo.yaml"
    shell:
        """
        python3 scripts/annotations/defensefinder_df_preformat.py \
            --input_genes {input.g} \
            --input_systems {input.s} \
            --type antidefense_system \
            --output_df {output}
        """


rule ISEScan_run:
    input:
        fa=rules.gbk_to_fa.output.fa,
    output:
        d=directory("data/ISEScan/{acc}"),
        s="data/ISEScan/{acc}/fasta/{acc}.fa.tsv",
    log:
        "logs/ISEScan/{acc}.log",
    conda:
        "../config/conda_envs/isescan.yaml"
    threads: 4
    shell:
        """
        isescan.py --seqfile {input.fa} --output {output.d} --nthread {threads} &>{log}
        """


rule ISEScan_preformat:
    input:
        lambda w: expand(rules.ISEScan_run.output.s, acc=acc_nums),
    output:
        "results/mges/ISEScan.csv",
    conda:
        "../config/conda_envs/bioinfo.yaml"
    shell:
        """
        python3 scripts/annotations/IS_df_preformat.py \
            --input_tsvs {input} \
            --output_df {output}
        """


rule abricate_run:
    input:
        fa=rules.gbk_to_fa.output.fa,
    output:
        "data/abricate/{acc}/{acc}.tsv",
    log:
        "logs/abricate/{acc}.log",
    conda:
        "../config/conda_envs/abricate.yaml"
    threads: 4
    shell:
        """
        abricate --db card --threads {threads} {input.fa} >{output} 2>{log}
        """


rule abricate_preformat:
    input:
        lambda w: expand(rules.abricate_run.output, acc=acc_nums),
    output:
        "results/mges/abricate.csv",
    conda:
        "../config/conda_envs/bioinfo.yaml"
    shell:
        """
        python3 scripts/annotations/abricate_df_preformat.py \
            --input_tsvs {input} \
            --output_df {output}
        """


rule integron_finder_run:
    input:
        fa=rules.gbk_to_fa.output.fa,
    output:
        i="data/integron_finder/{acc}/{acc}.integrons",
    log:
        "logs/integron_finder/{acc}.log",
    conda:
        "../config/conda_envs/integron_finder.yaml"
    threads: 4
    params:
        outdir="data/integron_finder/{acc}",
    shell:
        """
        integron_finder --local-max --func-annot --circ --cpu {threads} \
            --outdir {params.outdir} {input.fa} &>{log}
        # IF nests output under Results_Integron_Finder_<name>/; flatten the .integrons
        cp {params.outdir}/Results_Integron_Finder_*/*.integrons {output.i}
        """


rule integron_finder_preformat:
    input:
        lambda w: expand(rules.integron_finder_run.output.i, acc=acc_nums),
    output:
        "results/mges/integron_finder.csv",
    conda:
        "../config/conda_envs/bioinfo.yaml"
    shell:
        """
        python3 scripts/annotations/integron_finder_df_preformat.py \
            --input_tsvs {input} \
            --output_df {output}
        """


# whether each tool reports element coordinates 0-based;
# ISEScan/abricate/integron_finder are 1-based; antidefense mirrors defensefinder
ZERO_BASED = {
    "genomad": True,
    "defensefinder": True,
    "antidefense": True,
    "ISEScan": False,
    "abricate": False,
    "integron_finder": False,
}

MGE_TOOLS = list(ZERO_BASED)


rule mge_assign_positions:
    input:
        el="results/mges/{tool}.csv",
        j_pos=rules.junction_positions.output.csv,
        iso_len=rules.genome_lengths.output,
    output:
        "results/mges_to_junctions/{tool}.csv",
    conda:
        "../config/conda_envs/bioinfo.yaml"
    params:
        zero_based=lambda w: "--zero_based" if ZERO_BASED[w.tool] else "",
    shell:
        """
        python3 scripts/assign_junctions.py \
            --iso_len {input.iso_len} \
            --junction_pos_csv {input.j_pos} \
            --element_pos_df {input.el} \
            --output_pos {output} \
            {params.zero_based}
        """


rule mge_to_junction_gff:
    input:
        assigned=expand(rules.mge_assign_positions.output, tool=MGE_TOOLS),
        elements=expand("results/mges/{tool}.csv", tool=MGE_TOOLS),
        iso_len=rules.genome_lengths.output,
    output:
        "results/junction_mges/{junc}.gff",
    log:
        "logs/mge_to_junction_gff/{junc}.log",
    conda:
        "../config/conda_envs/bioinfo.yaml"
    shell:
        """
        python3 scripts/mges_to_gff.py \
            --junc_id {wildcards.junc} \
            --mges_to_junctions {input.assigned} \
            --mges {input.elements} \
            --iso_len {input.iso_len} \
            --out_gff {output} \
            &>{log}
        """


def junction_mge_gffs(wildcards):
    return expand(rules.mge_to_junction_gff.output, junc=junction_ids(wildcards))


rule mges_all:
    input:
        expand(rules.mge_assign_positions.output, tool=MGE_TOOLS),
        junction_mge_gffs,
