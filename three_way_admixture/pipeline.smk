"""
#########README##########
Running the simulations of the finemapping across different patterns of LD
Assume we have generated the coalescent simulations locally with files located in: /Users/marcustutert/Desktop/Oxford_Dphil/InferLD_Validations/Chapter_2_Simulations/finemapping
We run this file on the directory:  /well/mcvean/mtutert/thesis/finemapping
With the command:
snakemake --snakefile PATH/TO/SNAKEMAKE-j 100 --max-status-checks-per-second 0.01 --profile /well/mcvean/mtutert/snakemake/profile -f --rerun-incomplete -n
Note that we may need to add specific things to the profile as we go on the fly--but unlikely, besides changing around default
"""

import itertools as itertools
print("Executing snakefile")
#Run the simulation across different number of reference panel regions
regions             = list(range(1,185)) #
rule all:
    input:
        #expand("hapgen2/pop_ADMIX_region_{regions}.map", regions = regions)
        #expand("hapgen2/gwas_region_{regions}.summary", regions = regions),
        #expand("snptest/sumstats_gwas_region_{regions}", regions = regions),
        #expand("finemap/ref_admixed_ld_true_sumstats_region_{regions}.config", regions = regions)
        #expand("finemap/ref_admixed_ld_true_sumstats_region_{regions}_cred_set_fdr.RData", regions = regions)
        "results/pooled_ref_admixed_ld_true_sumstats_region_indicator.RData"

# #Generate msprime simulations (GWAS panel and the associated references)

# rule msprime_sims:
#     output:
#         "pop_split/msprime/gwas_99_dt_{divergence_time}.csv" #...ref_MAX... in replicates (HARD CODED)
#     conda:
#         "/well/mcvean/mtutert/snakemake/envs/python_3.7.yaml"
#     params:
#         divergence_time = "{divergence_time}",
#     script:
#         "pop_split/msprime_pop_split_sim.py"

#Prepare files for HAPGEN2
rule hapgen2_prep:
    input:
        "msprime_data/ref_pop_EAS_region_186.csv"
    output:
        "hapgen2/pop_ADMIX_region_{regions}.map"
    params:
        regions = "{regions}",
    conda:
        "/well/mcvean/mtutert/snakemake/envs/r.yaml"
    script:
        "HAPGEN2_Prep.R"

#Run Hapgen2 to generate GWAS case and control data
rule runhapgen2:
    input:
        "hapgen2/pop_ADMIX_region_{regions}.map"
    output:
        "hapgen2/gwas_region_{regions}.summary"
    params:
        regions = "{regions}"
    conda:
        "/well/mcvean/mtutert/snakemake/envs/r.yaml"
    script:
        "run_HAPGEN2.R"

#Run SNPTEST to generate the GWAS summary statisticsrule runsnptest:
rule runsnptest:
    input:
        "hapgen2/gwas_region_{regions}.summary"
    output:
        "snptest/sumstats_gwas_region_{regions}"
    params:
        regions = "{regions}"
    conda:
        "/well/mcvean/mtutert/snakemake/envs/r.yaml"
    script:
        "run_SNPTEST2.R"

rule run_weight_inference:
    input:
        "snptest/sumstats_gwas_region_{regions}"
    params:
        regions = "{regions}"
    conda:
        "/well/mcvean/mtutert/snakemake/envs/r.yaml"
    output:
        "weight_inference_results/weights_region_{regions}"
    script:
        "weight_inference.R"

rule generate_ld_distributions:
    input:
        "weight_inference_results/weights_region_{regions}"
    params:
        regions = "{regions}"
    conda:
        "/well/mcvean/mtutert/snakemake/envs/r.yaml"
    output:
        "distribution_ld_results/inferred_LD_sample_100_region_{regions}"
    script:
        "generate_ld_distributions.R"

rule finemap_distributions_ld:
    input:
        "distribution_ld_results/inferred_LD_sample_100_region_{regions}"
    params:
        regions = "{regions}"
    conda:
        "/well/mcvean/mtutert/snakemake/envs/r.yaml"
    output:
        "finemap_distributions/inferred_ld_true_sumstats_region_{regions}_sample_100.snp"
    script:
        "finemap_distributions_ld.R"

#Prepare FINEMAP script across all the reference panels and then run it
rule prep_and_run_FINEMAP:
    input:
        #Use SNPTEST output to run this
        #"snptest/sumstats_gwas_region_{regions}"
        #"distribution_ld_results/inferred_LD_sample_100_region_{regions}"
        #"weight_inference_results/weights_region_{regions}"
        "finemap_distributions/inferred_ld_true_sumstats_region_{regions}_sample_100.snp"
        #expand("pop_split/msprime/ref_{replicates}_dt_{divergence_time}.csv", replicates = replicates, divergence_time = divergence_time)
    output:
        #FINEMAP log files
        "finemap/ref_admixed_ld_true_sumstats_region_{regions}.config"
    conda:
        "/well/mcvean/mtutert/snakemake/envs/r.yaml"
    params:
        regions = "{regions}"
    script:
        "prep_run_FINEMAP.R"

#Check the results of the credible set across FINEMAP runs
rule calculate_credible_set:
    input:
        "finemap/ref_admixed_ld_true_sumstats_region_{regions}.config"
    output:
        "finemap/ref_admixed_ld_true_sumstats_region_{regions}_cred_set_fdr.RData"
    conda:
        "/well/mcvean/mtutert/snakemake/envs/r.yaml"
    params:
        regions = "{regions}"
    script:
        "calculate_credible_set.R"

#Check the results of the credible set across FINEMAP runs
rule pool_credible_sets:
    input:
        expand("finemap/ref_admixed_ld_true_sumstats_region_{regions}_cred_set_fdr.RData", regions = regions)
    output:
        "results/pooled_ref_admixed_ld_true_sumstats_region_indicator.RData"
    conda:
        "/well/mcvean/mtutert/snakemake/envs/r.yaml"
    script:
        "pool_credible_sets.R"
