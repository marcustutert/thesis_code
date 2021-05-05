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
replicates             = list(range(1,100)) #
divergence_time        = ["0","10","50","125"]
rule all:
    input:
        #expand("pop_split/msprime/gwas_4_dt_{divergence_time}.csv", replicates = replicates,divergence_time = divergence_time) #...ref_MAX... in replicates (HARD CODED)
        #expand("pop_split/hapgen2/gwas_{replicates}_dt_{divergence_time}.map", divergence_time = divergence_time, replicates = replicates)
        #expand("pop_split/hapgen2/gwas__{replicates}_dt_{divergence_time}.summary", divergence_time = divergence_time, replicates = replicates)
        #expand("pop_split/snptest/sumstats_gwas_{replicates}_dt_{divergence_time}", replicates = replicates, divergence_time = divergence_time)
        #expand("pop_split/finemap/gwas_{replicates}_dt_{divergence_time}.config", replicates = replicates, divergence_time = divergence_time)
        #expand("pop_split/results/gwas_{replicates}_dt_{divergence_time}_cred_set_fdr.RData", replicates = replicates, divergence_time = divergence_time)
        expand("pop_split/results/pooled_fdr_dt_{divergence_time}.RData", divergence_time = divergence_time)
#Generate msprime simulations (GWAS panel and the associated references)
rule msprime_sims:
    output:
        "pop_split/msprime/gwas_99_dt_{divergence_time}.csv" #...ref_MAX... in replicates (HARD CODED)
    conda:
        "/well/mcvean/mtutert/snakemake/envs/python_3.7.yaml"
    params:
        divergence_time = "{divergence_time}",
    script:
        "pop_split/msprime_pop_split_sim.py"

#Prepare files for HAPGEN2
rule hapgen2_prep:
    input:
        "pop_split/msprime/gwas_99_dt_{divergence_time}.csv"
    output:
        "pop_split/hapgen2/gwas_{replicates}_dt_{divergence_time}.map"
    params:
        divergence_time = "{divergence_time}",
        replicates      = "{replicates}"
    conda:
        "/well/mcvean/mtutert/snakemake/envs/r.yaml"
    script:
        "pop_split/hapgen2_prep.R"

#Run Hapgen2 to generate GWAS case and control data
rule runhapgen2:
    input:
        "pop_split/hapgen2/gwas_{replicates}_dt_{divergence_time}.map"
    output:
        "pop_split/hapgen2/gwas_{replicates}_dt_{divergence_time}.summary"
    params:
        divergence_time = "{divergence_time}",
        replicates      = "{replicates}"
    conda:
        "/well/mcvean/mtutert/snakemake/envs/r.yaml"
    script:
        "pop_split/runhapgen2.R"

#Run SNPTEST to generate the GWAS summary statisticsrule runsnptest:
rule runsnptest:
    input:
        "pop_split/hapgen2/gwas_{replicates}_dt_{divergence_time}.summary"
    output:
        "pop_split/snptest/sumstats_gwas_{replicates}_dt_{divergence_time}"
    params:
        divergence_time = "{divergence_time}",
        replicates      = "{replicates}"
    conda:
        "/well/mcvean/mtutert/snakemake/envs/r.yaml"
    script:
        "pop_split/runsnptest.R"

#Prepare FINEMAP script across all the reference panels and then run it
rule prep_and_run_FINEMAP:
    input:
        #Use SNPTEST output to run this
        "pop_split/snptest/sumstats_gwas_{replicates}_dt_{divergence_time}"
        #expand("pop_split/msprime/ref_{replicates}_dt_{divergence_time}.csv", replicates = replicates, divergence_time = divergence_time)
    output:
        #FINEMAP log files
        "pop_split/finemap/gwas_{replicates}_dt_{divergence_time}.config"
    conda:
        "/well/mcvean/mtutert/snakemake/envs/r.yaml"
    params:
        divergence_time = "{divergence_time}",
        replicates      = "{replicates}"
    script:
        "pop_split/prep_and_run_FINEMAP.R"

#Check the results of the credible set across FINEMAP runs
rule calculate_credible_set:
    input:
        "pop_split/finemap/gwas_{replicates}_dt_{divergence_time}.config"
    output:
        "pop_split/results/gwas_{replicates}_dt_{divergence_time}_cred_set_fdr.RData"
    conda:
        "/well/mcvean/mtutert/snakemake/envs/r.yaml"
    params:
        replicates      = "{replicates}",
        divergence_time = "{divergence_time}"
    script:
        "pop_split/calculate_credible_set.R"

#Check the results of the credible set across FINEMAP runs
rule pool_credible_sets:
    input:
        expand("pop_split/results/gwas_{replicates}_dt_{divergence_time}_cred_set_fdr.RData", replicates = replicates, divergence_time = divergence_time)
    output:
        "pop_split/results/pooled_fdr_dt_{divergence_time}.RData"
    conda:
        "/well/mcvean/mtutert/snakemake/envs/r.yaml"
    params:
        divergence_time = "{divergence_time}"
    script:
        "pop_split/pool_credible_sets.R"

