# Junction dataset

This repository contains a simple pipeline to prepare and format the dataset of genomic junctions from the _E. coli_ ST131 collection described in [our paper](https://academic.oup.com/mbe/article/42/1/msae272/7942412).

## Setup

To run the pipeline you will need to have [Snakemake](https://snakemake.readthedocs.io/en/stable/) and [Conda](https://docs.conda.io/en/latest/) installed.

## Usage

Run the pipeline simply with:

```sh
snakemake --use-conda --cores <num_cores> all
```