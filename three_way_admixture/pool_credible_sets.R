#Pool credible sets together

results_fdr_files = c("in_sample_ld_true_sumstats_region",
                      "ref_admixed_ld_true_sumstats_region",
                      "EAS_ld_true_sumstats_region",
                      "EUR_ld_true_sumstats_region",
                      "AFR_ld_true_sumstats_region",
                      "1000G_proportions_ld_true_sumstats_region",
                      "1000G_weighted_ld_true_sumstats_region")

results_indicator_files = c("in_sample_ld_true_sumstats_region",
                            "ref_admixed_ld_true_sumstats_region",
                            "EAS_ld_true_sumstats_region",
                            "EUR_ld_true_sumstats_region",
                            "AFR_ld_true_sumstats_region",
                            "1000G_proportions_ld_true_sumstats_region",
                            "1000G_weighted_ld_true_sumstats_region")

for (i in 1:length(results_indicator_files)) {
  files_1           = list.files("finemap", pattern = sprintf("%s_\\w+_cred_set_fdr.RData",results_fdr_files[i]),full.names = T)
  files_2           = list.files("finemap", pattern = sprintf("%s_\\w+_cred_set_indicator.RData",results_indicator_files[i]),full.names = T)
  results_1         = lapply(files_1,readRDS,.GlobalEnv)
  results_2         = lapply(files_2,readRDS,.GlobalEnv)
  
  #Write out the pooled quantiles
  saveRDS(results_2,sprintf("results/pooled_%s_indicator.RData",results_indicator_files[i]), version = 2)
  saveRDS(results_1,sprintf("results/pooled_%s_fdr.RData", results_fdr_files[i]), version = 2)
}

