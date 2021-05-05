"""
#########README##########
Running the simulations of the finemapping across different patterns of LD (with summary statistics imputation!)
We run this file on the directory:  /well/mcvean/mtutert/thesis/finemapping
snakemake --snakefile PATH/TO/SNAKEMAKE-j 10000 --max-status-checks-per-second 0.01 --profile /well/mcvean/mtutert/snakemake/profile -f --rerun-incomplete -n

####README#####
If signal is FALSE: only run rules up to pool pool_sumstat_imputation
If signal is TRUE run all rules

"""

#Note we uncomment out second last line if we are just looking at sumstat imputation


import itertools as itertools
print("Executing snakefile")

#Run the simulation across different number of reference panel regions
replicates             = list(range(1,1000))         #Ask Rao how to not have this be hard-coded in
divergence_time        = ["0","10","50","125"]
model                  = "pop_split"              #Change this between OOA model and pop_split
signal                 = "FALSE"                   #True means causal variants, False is under the null model

rule all:
    input:
        #expand("pop_split/msprime/gwas_4_dt_{divergence_time}.csv", divergence_time = divergence_time), #...ref_MAX... in replicates (HARD CODED)
        #expand("pop_split/hapgen2/gwas_{replicates}_dt_{divergence_time}.map", divergence_time = divergence_time, replicates = replicates, signal = signal),
        #expand("pop_split/hapgen2/gwas_{signal}_{replicates}_dt_{divergence_time}.summary", divergence_time = divergence_time, replicates = replicates, signal = signal),
        #expand("pop_split/snptest/sumstats_gwas_{signal}_{replicates}_dt_{divergence_time}", replicates = replicates, divergence_time = divergence_time, signal = signal),
        #expand("pop_split/snptest/sumstats_imputed_results_gwas_LD_gwas_{signal}_{replicates}_dt_{divergence_time}", signal = signal,replicates = replicates, divergence_time = divergence_time),
        #expand("pop_split/sumstat_imputation/results/pooled_rare_variants_{signal}_gwas_LD_dt_{divergence_time}", signal = signal, divergence_time = divergence_time)
        #expand("pop_split/finemap/ref_ld_imputed_sumstats_by_gwas_ld_region_{replicates}_dt_{divergence_time}_{signal}.log_sss", replicates = replicates, divergence_time = divergence_time, signal = signal),
        #expand("pop_split/finemap/gwas_ld_imputed_sumstats_by_gwas_ld_region_{replicates}_dt_{divergence_time}_cred_set_{signal}_indicator.RData", replicates = replicates, divergence_time = divergence_time, signal = signal),
        expand("pop_split/sumstat_imputation/results/pooled_rare_variants_{signal}_gwas_LD_dt_{divergence_time}", signal = signal, divergence_time = divergence_time)
        #expand("pop_split/results/pooled_gwas_ld_imputed_sumstats_by_ref_ld_region_fdr_dt_{divergence_time}.RData", divergence_time = divergence_time)
#Generate msprime simulations (GWAS panel and the associated references)
rule msprime_sims:
    output:
        "pop_split/msprime/gwas_999_dt_{divergence_time}.csv"     ######HARD CODED IN
    conda:
        "/well/mcvean/mtutert/snakemake/envs/python_3.7.yaml"
    params:
        divergence_time = "{divergence_time}"
    script:
        #Change number of replicates within the msprime_pop_split_sim.py script itself!
        "pop_split/msprime_pop_split_sim.py"

#Prepare files for HAPGEN2
rule hapgen2_prep:
    input:
        "pop_split/msprime/gwas_999_dt_{divergence_time}.csv"      ######HARD CODED IN
    output:
        "pop_split/hapgen2/gwas_{replicates}_dt_{divergence_time}.map"
    params:
        divergence_time = "{divergence_time}",
        replicates      = "{replicates}",
    conda:
        "/well/mcvean/mtutert/snakemake/envs/r.yaml"
    script:
        "pop_split/hapgen2_prep.R"

#Run Hapgen2 to generate GWAS case and control data
rule runhapgen2:
    input:
        "pop_split/hapgen2/gwas_{replicates}_dt_{divergence_time}.map"
    output:
        "pop_split/hapgen2/gwas_{signal}_{replicates}_dt_{divergence_time}.summary"
    params:
        divergence_time = "{divergence_time}",
        replicates      = "{replicates}",
        signal          = "{signal}"
    conda:
        "/well/mcvean/mtutert/snakemake/envs/r.yaml"
    script:
        "pop_split/runhapgen2.R"

#Run SNPTEST to generate the GWAS summary statisticsrule runsnptest:
rule runsnptest:
    input:
        "pop_split/hapgen2/gwas_{signal}_{replicates}_dt_{divergence_time}.summary"
    output:
        "pop_split/snptest/sumstats_gwas_{signal}_{replicates}_dt_{divergence_time}"
    params:
        divergence_time = "{divergence_time}",
        replicates      = "{replicates}",
        signal          = "{signal}"
    conda:
        "/well/mcvean/mtutert/snakemake/envs/r.yaml"
    script:
        "pop_split/runsnptest.R"

#Perform the sumstat imputation, also get the r2's for rare and common variants across each region
rule sumstat_imputation:
    input:
        #Use SNPTEST output to run this
        "pop_split/snptest/sumstats_gwas_{signal}_{replicates}_dt_{divergence_time}"
    output:
        #Files with same format as SNPTEST output (for purposes of using FINEMAP) but with summary stat imputation
        "pop_split/snptest/sumstats_imputed_results_gwas_LD_gwas_{signal}_{replicates}_dt_{divergence_time}"
    conda:
        "/well/mcvean/mtutert/snakemake/envs/r.yaml"
    params:
        divergence_time = "{divergence_time}",
        replicates      = "{replicates}",
        signal          = "{signal}"
    script:
        "pop_split/sumstat_imputation.R"

#Perform the sumstat imputation, also get the r2's for rare and common variants across each region
rule pool_sumstat_imputation:
    input:
        #Use SNPTEST output to run this
        expand("pop_split/snptest/sumstats_imputed_results_gwas_LD_gwas_{signal}_{replicates}_dt_{divergence_time}", signal = signal, replicates = replicates, divergence_time = divergence_time)
    output:
        #Files with same format as SNPTEST output (for purposes of using FINEMAP) but with summary stat imputation
        "pop_split/sumstat_imputation/results/pooled_rare_variants_{signal}_gwas_LD_dt_{divergence_time}"
    conda:
        "/well/mcvean/mtutert/snakemake/envs/r.yaml"
    params:
        divergence_time = "{divergence_time}",
        signal          = "{signal}"
    script:
        "pop_split/pool_sumstat_imputation.R"

#Prepare FINEMAP script across all the reference panels and then run it
rule prep_and_run_FINEMAP:
    input:
        #Use SNPTEST output to run this
        "pop_split/snptest/sumstats_imputed_results_gwas_LD_gwas_{signal}_{replicates}_dt_{divergence_time}"
        #expand("pop_split/msprime/ref_{replicates}_dt_{divergence_time}.csv", replicates = replicates, divergence_time = divergence_time)
    output:
        #FINEMAP log files
        "pop_split/finemap/ref_ld_imputed_sumstats_by_gwas_ld_region_{replicates}_dt_{divergence_time}_{signal}.log_sss"

    conda:
        "/well/mcvean/mtutert/snakemake/envs/r.yaml"
    params:
        divergence_time = "{divergence_time}",
        replicates      = "{replicates}",
        signal          = "{signal}"
    script:
        "pop_split/prep_and_run_FINEMAP.R"

#Check the results of the credible set across FINEMAP runs
rule calculate_credible_set:
    input:
        "pop_split/finemap/ref_ld_imputed_sumstats_by_gwas_ld_region_{replicates}_dt_{divergence_time}_{signal}.log_sss"
    output:
        "pop_split/finemap/gwas_ld_imputed_sumstats_by_gwas_ld_region_{replicates}_dt_{divergence_time}_cred_set_{signal}_indicator.RData"
    conda:
        "/well/mcvean/mtutert/snakemake/envs/r.yaml"
    params:
        replicates      = "{replicates}",
        divergence_time = "{divergence_time}",
        signal          = "{signal}"
    script:
        "pop_split/calculate_credible_set.R"

#Pool the credible set metrics across sumstats
rule pool_credible_sets:
    input:
        expand("pop_split/finemap/gwas_ld_imputed_sumstats_by_gwas_ld_region_{replicates}_dt_{divergence_time}_cred_set_{signal}_indicator.RData", replicates = replicates, divergence_time = divergence_time, signal = signal)
    output:
        "pop_split/results/pooled_gwas_ld_imputed_sumstats_by_ref_ld_region_fdr_dt_{divergence_time}.RData"
    conda:
        "/well/mcvean/mtutert/snakemake/envs/r.yaml"
    params:
        divergence_time = "{divergence_time}"
    script:
        "pop_split/pool_credible_sets.R"

