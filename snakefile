import pandas as pd


configfile: "config/config.yaml"


# load the list of accession numbers
acc_nums = pd.read_csv(config["acc_nums_file"], header=None)[0].tolist()
print(f"N. isolates: {len(acc_nums)}")
# load the table of junction positions


rule download_gbk:
    output:
        "data/gbk/{acc}.gbk",
    conda:
        "entrez-direct"
    shell:
        """
        esearch -db nucleotide -query {wildcards.acc}[Accession] | efetch -format gbwithparts > {output}
        """


rule all:
    input:
        expand(rules.download_gbk.output, acc=acc_nums),
