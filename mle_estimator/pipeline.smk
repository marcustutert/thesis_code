"""
#########README##########
Running the simulations of the finemapping across different patterns of LD (with summary statistics imputation!)
We run this file on the directory:  /well/mcvean/mtutert/thesis/finemapping
snakemake --snakefile PATH/TO/SNAKEMAKE-j 10000 --max-status-checks-per-second 0.01 --profile /well/mcvean/mtutert/snakemake/profile -f --rerun-incomplete -n

####README#####
If signal is FALSE: only run rules up to pool pool_sumstat_imputation
If signal is TRUE run all rules

"""
import itertools as itertools
print("Executing snakefile")

#Run the simulation across different number of replicates
divergence_time             = ["0","10","50","100"]
noise                       = ["0.001","0.005","0.01"]
replicates                  = list(range(1,1000))

rule all:
    input:
        #expand("pop_split_msprime/target_1_replicate_{replicates}_dt_{divergence_time}.csv", replicates = replicates, divergence_time = divergence_time)
        expand("mle_results/glm_se/fst_replicate_{replicates}_dt_{divergence_time}_noise_{noise}", replicates = replicates, divergence_time = divergence_time, noise = noise)

#Generate msprime simulations (GWAS panel and the associated references)
rule msprime_sims:
    output:
        "msprime_data/gwas_999_dt_{divergence_time}.csv"
    conda:
        "/well/mcvean/mtutert/snakemake/envs/python_3.7.yaml"
    params:
        divergence_time         = "{divergence_time}",
    script:
        "msprime_simulate.py"

#Do the MLE estimates
rule mle_estimate:
    input:
        expand("msprime_data/gwas_999_dt_{divergence_time}.csv", divergence_time = divergence_time)
    output:
        "mle_results/glm_se/fst_replicate_{replicates}_dt_{divergence_time}_noise_{noise}"
    conda:
        "/well/mcvean/mtutert/snakemake/envs/r.yaml"
    params:
        divergence_time         = "{divergence_time}",
        replicates              = "{replicates}",
        noise                   = "{noise}"
    script:
        "mle_estimate.R"

