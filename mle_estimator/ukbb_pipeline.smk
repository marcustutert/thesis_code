"""
#########README##########
Running the inference of UKBB coalescent simulations on the cluster
Assume we have generated the coalescent simulations locally with files located in: /Users/marcustutert/Desktop/Oxford_Dphil/InferLD_Validations/coalescent_simulations/prior_draws
We run this file on the directory:  /well/mcvean/mtutert/thesis/coalescent_coverage/prior_draws
With the command:
snakemake --snakefile PATH/TO/SNAKEMAKE-j 5000 --max-status-checks-per-second 0.01 --profile /well/mcvean/mtutert/snakemake/profile -f --rerun-incomplete -n
Note that we may need to add specific things to the profile as we go on the fly--but unlikely, besides changing around default
"""

import itertools as itertools
regions = list(range(1,38)) #Regions in UKBB GWAS data
#Reference panel used will be all 1000G super populations
populations = ["EUR","AFR","AMR","SAS","EAS"]

rule all:
    input:
        #expand("pop_split_msprime/target_1_replicate_{replicates}_dt_{divergence_time}.csv", replicates = replicates, divergence_time = divergence_time)
        expand("ukbb_mle_results/noise_{populations}_panel_replicate_{regions}", populations = populations, regions = regions)

#Generate the msprime simulations

#Do the MLE estimates
rule mle_estimate:
    output:
        "ukbb_mle_results/noise_{populations}_panel_replicate_{regions}"
    conda:
        "/well/mcvean/mtutert/snakemake/envs/r.yaml"
    params:
        populations              = "{populations}",
        regions                  = "{regions}"
    script:
        "ukbb_mle_estimate.R"

