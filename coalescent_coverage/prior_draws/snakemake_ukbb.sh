

#!/bin/bash

#$ -N mtest

#$ -q long.qc

#$ -P mcvean.prjc

#$ -cwd -V

#$ -pe shmem 2

#Sumbit OOA pipeline script
snakemake --snakefile UKBB_pipeline.smk -j 5000 --max-status-checks-per-second 0.01 --profile /well/mcvean/mtutert/snakemake/profile -f --rerun-incomplete
