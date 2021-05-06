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
populations   = ["EUR","AFR","EUR_AFR"] #Populations used in the reference panel
regions       = list(range(1,221))                     #This will be hard-coded in for now, and correspond to the result of the genome_segmentations
fst           = [0.1,0.01,0.001]
rule all:
    input:
        #expand("pop_split_msprime/target_1_replicate_{replicates}_dt_{divergence_time}.csv", replicates = replicates, divergence_time = divergence_time)
        expand("results/inference_{populations}_panel_region_{regions}_fst_{fst}", populations = populations, regions = regions, fst = fst)


#This only needs to be done for the different populations (or combinations thereof) (so no need to do over params)
#Will also do the genome_segmentation that is required as well, note that this will hardcode in the regions parameter above
#Can parallelis later since Im being lazy now, takes about 15min

######UNCOMMENT IF CHANGING AROUND THE WINDOWING
# rule prep_data:
#     conda:
#         "/well/mcvean/mtutert/snakemake/envs/r.yaml"
#     script:
#         "prep_ld_data.R"


#Run the Gibbs Sampling
rule run_gibbs:
    output:
        "results/inference_{populations}_panel_region_{regions}_fst_{fst}"
    conda:
        "/well/mcvean/mtutert/snakemake/envs/r.yaml"
    params:
        populations                  = "{populations}",
        regions                      = "{regions}",
        fst                          = "{fst}"
    script:
        "ld_inference_across_populations_parameters.R"
