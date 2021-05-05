#Code to run SNPTEST on the output of HAPGEN2
#No need (usually) to prep data since they work off each other
regions = snakemake@params$regions


system(sprintf("/apps/well/snptest/2.5/snptest -data hapgen2/gwas_region_%s.controls.gen hapgen2/gwas_region_%s.controls.sample hapgen2/gwas_region_%s.cases.gen hapgen2/gwas_region_%s.cases.sample -o snptest/sumstats_gwas_region_%s -frequentist 1 -method score -pheno pheno",
               regions,regions,regions,regions,regions))
