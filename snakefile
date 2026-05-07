import pandas as pd
import json
from pathlib import Path


configfile: "config/config.yaml"


# load the list of accession numbers
acc_nums = pd.read_csv(config["acc_nums_file"], header=None)[0].tolist()
print(f"N. isolates: {len(acc_nums)}")
# list of plasmids
with open(config["plasmids_file"]) as f:
    plasmids = json.load(f)

# optional NCBI API key (gitignored). Raises efetch rate limit from 3 to 10 req/s.
ncbi_api_key_file = Path("config/ncbi_api_key.txt")
NCBI_API_KEY = (
    ncbi_api_key_file.read_text().strip() if ncbi_api_key_file.exists() else ""
)
# seconds to sleep after each NCBI download to stay under rate limits


wildcard_constraints:
    acc="[^/]+",
    junc="[^/]+",
    acc_iso="[^/]+",


rule download_gbk:
    output:
        "data/gbk/{acc}.gbk",
    conda:
        "config/conda_envs/ncbi_acc_download.yaml"
    params:
        api_key_arg=(f"--api-key {NCBI_API_KEY}" if NCBI_API_KEY else ""),
    shell:
        """
        ncbi-acc-download --out {output} {params.api_key_arg} {wildcards.acc}
        """


rule gbk_to_fa:
    input:
        gbk="data/gbk/{acc}.gbk",
    output:
        fa="data/fasta/{acc}.fa",
    conda:
        "config/conda_envs/bioinfo.yaml"
    shell:
        """
        python3 scripts/gbk_to_fa.py --gbk {input.gbk} --fa {output.fa}
        """


rule genome_lengths:
    input:
        expand(rules.download_gbk.output, acc=acc_nums),
    output:
        "results/genome_lengths.csv",
    conda:
        "config/conda_envs/bioinfo.yaml"
    shell:
        """
        python scripts/genome_lengths.py \
            --output {output} \
            --gbk_files {input}
        """


rule plasmids:
    input:
        rules.download_gbk.output,
    output:
        "results/plasmids/{acc_iso}/{acc}.gbk",
    shell:
        """
        cp {input} {output}
        """


def all_plasmid_outputs(wildcards):
    outs = []
    for iso, plasmid_list in plasmids.items():
        outs.extend(expand(rules.plasmids.output, acc_iso=iso, acc=plasmid_list))
    return outs


include: "rules/pangraph.smk"
include: "rules/junctions.smk"
include: "rules/mges.smk"


rule all:
    input:
        expand(rules.build_junction_pangraph.output, junc=junc_ids),
        rules.genome_lengths.output,
        all_plasmid_outputs,
        rules.junction_stats.output,


localrules:
    download_gbk,
    plasmids,
    genomad_download_db,
    defensefinder_models_download,
