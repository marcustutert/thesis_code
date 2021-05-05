regions = snakemake@params$regions
library(data.table)
afr_reference_panel    = as.matrix(fread(sprintf("msprime_data/joint_filtered/ref_pop_AFR_region_%s", regions), header = T))
#Code to do the FINEMAPPING across the distributions of LD
#Read in reference panel:
afr_reference_panel    = as.matrix(fread(sprintf("msprime_data/joint_filtered/ref_pop_AFR_region_%s", regions), header = T))
eur_reference_panel    = as.matrix(fread(sprintf("msprime_data/joint_filtered/ref_pop_EUR_region_%s", regions), header = T))
eas_reference_panel    = as.matrix(fread(sprintf("msprime_data/joint_filtered/ref_pop_EAS_region_%s", regions), header = T))
#Get 1000G proportions: 1332 AFRs, 1006 EUR and 1004 EAS
KG_reference_panel = rbind(afr_reference_panel[1:1332,],eur_reference_panel[1:1006,])
nSamples           = 100
inferred_wts       = readRDS(sprintf("weight_inference_results/weights_region_%s",regions))
nsnps              = nrow(inferred_wts$inferred_af_given_weights)
burnin_region      = floor(0.5*length(inferred_wts$log_likelihood))
indexes_sampling   = sample(floor(burnin_region):length(inferred_wts$log_likelihood), size = nSamples)

#Store the array
ld_array = array(data = NA, dim = c(nsnps,nsnps,nSamples))
for (i in 1:length(indexes_sampling)) {
  #Calculate the LD on 100 samples in the post-burnin period, store as an array?
  #Get the LD across the samples
  wts                = inferred_wts$Gibbs_Array[,indexes_sampling[i]]
  inferred_LD        = cov.wt(KG_reference_panel, wt = wts, method="ML", cor=TRUE)$cor
  saveRDS(inferred_LD,sprintf("distribution_ld_results/inferred_LD_sample_%s_region_%s",i,regions))
  ld_array[,,i]           = inferred_LD
}

#Run the finemapping for each 