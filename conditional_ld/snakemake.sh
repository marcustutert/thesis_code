#!/bin/bash

#$ -N mtest

#$ -q short.qc

#$ -P mcvean.prjc

#$ -cwd -V

#$ -pe shmem 2
source ~/.bashrc
#export PATH="/well/mcvean/mtutert/anaconda3/bin/:$PATH"
#conda init bash
conda activate snakemake
snakemake --snakefile pop_split_pipeline.smk -j 100000000 --max-status-checks-per-second 0.01 --profile /well/mcvean/mtutert/snakemake/profile -f --rerun-incomplete --use-conda 
