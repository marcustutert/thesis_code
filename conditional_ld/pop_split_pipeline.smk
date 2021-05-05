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
divergence_time             = ["10","20","50"]         #Ask Rao how to not have this be hard-coded in
recombination_rate          = ["1e-3","1e-4","1e-5"]

rule all:
    input:
        #expand("pop_split_msprime/target_1_replicate_{replicates}_dt_{divergence_time}.csv", replicates = replicates, divergence_time = divergence_time)
        expand("pop_split_results/results_replicate_999_dt_{divergence_time}_recomb_{recombination_rate}", divergence_time = divergence_time, recombination_rate = recombination_rate)
#Generate msprime simulations (GWAS panel and the associated references)
rule msprime_sims:
    output:
        "pop_split_msprime/target_1_replicate_999_dt_{divergence_time}_recomb_{recombination_rate}.csv"
    conda:
        "/well/mcvean/mtutert/snakemake/envs/python_3.7.yaml"
    params:
        divergence_time         = "{divergence_time}",
        recombination_rate      = "{recombination_rate}"
    script:
        "msprime_pop_split_generation.py"

rule conditional_calculation:
    input:
        expand("pop_split_msprime/target_1_replicate_999_dt_{divergence_time}_recomb_{recombination_rate}.csv", recombination_rate = recombination_rate, divergence_time = divergence_time)
    output:
        "pop_split_results/results_replicate_999_dt_{divergence_time}_recomb_{recombination_rate}"
    conda:
        "/well/mcvean/mtutert/snakemake/envs/r.yaml"
    params:
        divergence_time         = "{divergence_time}",
        recombination_rate      = "{recombination_rate}"
    script:
        "pop_split_conditional_calculation.R"

