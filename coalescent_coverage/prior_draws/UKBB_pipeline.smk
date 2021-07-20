"""
#########README##########
Running the inference of UKBB coalescent simulations on the cluster
Assume we have generated the coalescent simulations locally with files located in: /Users/marcustutert/Desktop/Oxford_Dphil/InferLD_Validations/coalescent_simulations/prior_draws
We run this file on the directory:  /well/mcvean/mtutert/thesis/coalescent_coverage/prior_draws
With the command:
snakemake --snakefile PATH/TO/SNAKEMAKE-j 5000 --max-status-checks-per-second 0.01 --profile /well/mcvean/mtutert/snakemake/profile -f --rerun-incomplete -n
Note that we may need to add specific things to the profile as we go on the fly--but unlikely, besides changing around default
"""
###########################
#Define snakemakeparameters
###########################

import itertools as itertools
regions = list(range(6,10)) #Regions in Chr21 UKBB (Change size of regions w in Rscript)
regions = [str(i) for i in regions]
#Number of chunks, tells us how many total nSamples to perform (Change nSamples w in Rscript)
chunks     = list(range(1,2))
chunks     = [str(i) for i in chunks]

#Reference panel used will be all 1000G super populations while using UKBB WBA populations
populations = ["EUR","AFR","AMR","SAS","EAS"]
fst         = ["0.0001","0.05","0.1","0.25"]           #Manually set this!!!!!
#Again we want to form a grid of 1000G reference panel, UKBB population as the GWAS panel, and different Fst values

#Calculate all pairs of these values (note that the diagonal is what we want in the end)
paired_values    = list(itertools.product(populations,fst))
#Convert this to a string that looks appropriate
def join_tuple_string(strings_tuple) -> str:
   return '_'.join(strings_tuple)

# joining all the tuples
result = map(join_tuple_string, paired_values)

# converting and printing the result
paired_values = (list(result))
print(paired_values) #This has a format of: EAS_0.08 (ie: ref_pop_fst)

####Prepopulate the UKBB data and reference panels!
##############################################
#Rules

rule all:
    input:
        #expand("results/UKBB/{paired_values}_AF_Region_{regions}_Chunk_{chunks}.RData", paired_values = paired_values, regions = regions, chunks = chunks),
        #expand("results/UKBB/{paired_values}_AF_Region_{regions}_Chunk_{chunks}.RData", paired_values = paired_values, regions = regions, chunks = chunks)
        #expand("prior_draw_analysis_panels/regions/matched_panels/1000G_{populations}_region_{regions}_matched_ukbb_ref_panel.haps", populations = populations, regions = regions),
        #expand("prior_draw_analysis_panels/regions/matched_panels/ukbb_region_{regions}_matched_1000G_{populations}.haps", populations = populations, regions = regions)
        expand("results/UKBB/{paired_values}_AF_pooled_quantile_counts.RData", paired_values = paired_values),
        expand("results/UKBB/{paired_values}_LD_pooled_quantile_counts.RData", paired_values = paired_values)

#rule match_panels:
#Generate the UKBB and 1000G Panels (On the cluster, this is a new functionality)
    #output:
    #The output will be the matched panels in the matched_panels dir
        #"prior_draw_analysis_panels/regions/matched_panels/1000G_{populations}_region_{regions}_matched_ukbb_ref_panel.haps",
        #"prior_draw_analysis_panels/regions/matched_panels/ukbb_region_{regions}_matched_1000G_{populations}.haps"
    #params:
        #populations = "{populations}",
        #regions     = "{regions}"
    #conda:
        #"/well/mcvean/mtutert/snakemake/envs/r.yaml"
    #script:
        #"match_ukbb.R"

rule draw_weights:
    input:
        expand("prior_draw_analysis_panels/regions/matched_panels/1000G_{populations}_region_{regions}_matched_ukbb_ref_panel.haps", regions = regions, populations = populations),
        expand("prior_draw_analysis_panels/regions/matched_panels/ukbb_region_{regions}_matched_1000G_{populations}.haps", regions = regions, populations = populations)
    output:
        "results/UKBB/{paired_values}_AF_Region_{regions}_Chunk_{chunks}.RData",
        "results/UKBB/{paired_values}_LD_Region_{regions}_Chunk_{chunks}.RData"
    params:
        paired_values    = "{paired_values}",
        regions          = "{regions}",
        chunks           = "{chunks}"
    conda:
        "/well/mcvean/mtutert/snakemake/envs/r.yaml"
    script:
        "UKBB_draw_weights.R"


rule calculate_quantiles:
    input:
        expand("results/UKBB/{paired_values}_AF_Region_{regions}_Chunk_{chunks}.RData", regions = regions, paired_values = paired_values, chunks = chunks),
        expand("results/UKBB/{paired_values}_LD_Region_{regions}_Chunk_{chunks}.RData", regions = regions, paired_values = paired_values, chunks = chunks)
    output:
        "results/UKBB/{paired_values}_LD_quantile_counts_region_{regions}.RData",
        "results/UKBB/{paired_values}_AF_quantile_counts_region_{regions}.RData"
    params:
        regions           = "{regions}",
        paired_values     = "{paired_values}"
    conda:
        "/well/mcvean/mtutert/snakemake/envs/r.yaml"
    script:
        "UKBB_calculate_quantiles.R"
#
rule pool_quantiles:
    input:
        expand("results/UKBB/{paired_values}_LD_quantile_counts_region_{regions}.RData", regions = regions, paired_values = paired_values),
        expand("results/UKBB/{paired_values}_AF_quantile_counts_region_{regions}.RData", regions = regions, paired_values = paired_values)
    output:
        "results/UKBB/{paired_values}_AF_pooled_quantile_counts.RData",
        "results/UKBB/{paired_values}_LD_pooled_quantile_counts.RData"
    params:
        paired_values  = "{paired_values}"
    conda:
        "/well/mcvean/mtutert/snakemake/envs/r.yaml"
    script:
        "UKBB_pool_quantiles.R"


