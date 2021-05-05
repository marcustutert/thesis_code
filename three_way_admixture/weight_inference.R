#Code to take in the SNPTEST data and infer the weights
#Do this with a reference panel that has 1000G proportions (but change this around if need be!)
regions      = snakemake@params$regions
library(data.table)
#Load in the InferLD functions I need (double check its up to date!)
source("/well/mcvean/mtutert/thesis_code/thesis_code/ld_inference/InferLD/R/HMM_Functions.R")
source("/well/mcvean/mtutert/thesis_code/thesis_code/ld_inference/InferLD/R/validation_pipeline_functions.R")

#Read in reference panels that we will be using, for now it will be all three populations (EAS/EUR/AFR)
afr_reference_panel    = as.matrix(fread(sprintf("msprime_data/joint_filtered/ref_pop_AFR_region_%s", regions), header = T))
eur_reference_panel    = as.matrix(fread(sprintf("msprime_data/joint_filtered/ref_pop_EUR_region_%s", regions), header = T))
eas_reference_panel    = as.matrix(fread(sprintf("msprime_data/joint_filtered/ref_pop_EAS_region_%s", regions), header = T))

#Get 1000G proportions: 1332 AFRs, 1006 EUR and 1004 EAS
KG_reference_panel = rbind(afr_reference_panel[1:1332,],eur_reference_panel[1:1006,])

#Read in the "SNPTEST" summary statistics (note that we NO LONGER have to change scales---I dont think?)
summary_statistics = fread(sprintf("snptest/sumstats_gwas_region_%s", regions), skip = 10)
nsnps = ncol(KG_reference_panel)

#Run inference
mbar = 10000 #Change around if necessary! 
results                             = LD_from_GSHMM(ref_allele_matrix     = KG_reference_panel,
                                                    fst                   = 0.1,
                                                    alpha                 = 1e4,
                                                    nSamples              = 5,
                                                    weights_resolution    = 1000,
                                                    likelihood_toggle     = TRUE,
                                                    gwas_variance         = summary_statistics$frequentist_add_se_1^2,
                                                    LD_Infer              = FALSE,
                                                    recombination         = FALSE,
                                                    case_control_constant = mbar/2,
                                                    BurnIn                = 0.9,
                                                    debug                 = TRUE)

#Get metrics that we probably care about?
#Convert AF to MINOR ALLELE FREQUENCY
pi = c()
for (i in 1:nsnps) {
  if (results$inferred_af_given_weights[i,ncol(results$inferred_af_given_weights)] > 0.5) {
    pi[i] = 1-results$inferred_af_given_weights[i,ncol(results$inferred_af_given_weights)]
  }
  else{
    pi[i] = results$inferred_af_given_weights[i,ncol(results$inferred_af_given_weights)]
  }
}

#Get the posterior weights mean
wts          = rowMeans(results$Gibbs_Array[,floor(.9*ncol(results$Gibbs_Array)):ncol(results$Gibbs_Array)])
inferred_LD  = c(cov.wt(KG_reference_panel, wt = wts, method="ML", cor=TRUE)$cor)

#Write out the LD to use for finemapping
saveRDS(inferred_LD,sprintf("weight_inference_results/LD_inferred_region_%s", regions))
saveRDS(results,sprintf("weight_inference_results/weights_region_%s", regions))


