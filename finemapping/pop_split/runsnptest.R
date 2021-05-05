#Code to run SNPTEST on the output of HAPGEN2
#No need (usually) to prep data since they work off each other
divergence_time = snakemake@params$divergence_time
replicates      = snakemake@params$replicates
signal          = snakemake@params$signal

system(sprintf("/apps/well/snptest/2.5/snptest -data pop_split/hapgen2/gwas_%s_%s_dt_%s.controls.gen pop_split/hapgen2/gwas_%s_%s_dt_%s.controls.sample pop_split/hapgen2/gwas_%s_%s_dt_%s.cases.gen pop_split/hapgen2/gwas_%s_%s_dt_%s.cases.sample -o pop_split/snptest/sumstats_gwas_%s_%s_dt_%s -frequentist 1 -method score -pheno pheno",
               signal,replicates,divergence_time,signal, replicates, divergence_time,signal,replicates,divergence_time,signal,replicates,divergence_time,signal,replicates,divergence_time))

