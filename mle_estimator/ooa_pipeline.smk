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
replicates = list(range(1,1000))


#Need to write out all possible combo's for the reference and gwas panel
#Three populations are YRI, CEU & CHB
populations = ["CEU","CHB","YRI"]

#Get all permutations of the GWAS and reference population_pairs1
#Will sort out the matching of them later on
pairs = itertools.combinations(populations, 2)
population_pairs = []
for x in pairs:
    population_pairs.append('%s_GWAS_%s_Ref' %(x[0], x[1]))

rule all:
    input:
        #expand("pop_split_msprime/target_1_replicate_{replicates}_dt_{divergence_time}.csv", replicates = replicates, divergence_time = divergence_time)
        expand("ooa_mle_results/noise_{population_pairs}_replicate_{replicates}", population_pairs = population_pairs, replicates = replicates)

#Generate the msprime simulations

#Do the MLE estimates
rule mle_estimate:
    output:
        "ooa_mle_results/noise_{population_pairs}_replicate_{replicates}"
    conda:
        "/well/mcvean/mtutert/snakemake/envs/r.yaml"
    params:
        population_pairs         = "{population_pairs}",
        replicates               = "{replicates}"
    script:
        "ooa_mle_estimate.R"

