#####File to run PER-SAMPLE LD calculation
sample_number = snakemake@params[[1]][1]
replicate     = snakemake@params[[2]][1]
#Read in the reference and GWAS panel:
library(data.table)
source("/well/mcvean/mtutert/thesis/coalescent_coverage/helper_functions.R")
gwas = as.matrix(fread(sprintf("./msprime_data/population_split/GWAS_panel_replicate_%s_split_0.csv",replicate), header = T))#Perform LD inference and store as RData object
ref  = as.matrix(fread(sprintf("./msprime_data/population_split/GWAS_panel_replicate_%s_split_0.csv",replicate), header = T))

#Generate the matched ref/gwas/sumstats data so that everything is the same size
res           = msprime_gwas_sumstats(gwas_haplotypes = gwas, reference_haplotypes = ref)
se            = res[[1]]
gwas          = res[[2]]
ref           = res[[3]]
#Read in the RData
inference_results = readRDS(sprintf("./inference_results/inference_results_rep_%s", replicate))
sample_weights    =inference_results$Gibbs_Array[,,as.numeric(sample_number)]
cor.wt<-cov.wt(ref, sample_weights[,1], method="ML", cor=TRUE)
saveRDS(c(cor.wt$cor),sprintf("./inference_results/ld_results/ld_rep_%s_sample_%s",replicate,sample_number))

