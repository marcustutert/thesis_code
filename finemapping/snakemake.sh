#!/bin/bash

#$ -N mtest

#$ -q short.qc

#$ -P mcvean.prjc

#$ -cwd -V

#$ -pe shmem 2
source ~/.bashrc
conda activate snakemake
snakemake --snakefile pop_split_imputation_pipeline.smk -j 100000000 --max-status-checks-per-second 0.01 --profile /well/mcvean/mtutert/snakemake/profile -f --rerun-incomplete --use-conda 
