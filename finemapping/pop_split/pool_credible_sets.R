#Pool credible sets together
divergence_time = snakemake@params$divergence_time
signal          = "TRUE"

results_fdr_files = c("gwas_ld_true_sumstats_region",
                            "gwas_ld_imputed_sumstats_by_gwas_ld_region",
                            "gwas_ld_imputed_sumstats_by_ref_ld_region",
                            "ref_ld_true_sumstats_region",
                            "ref_ld_imputed_sumstats_by_gwas_ld_region",
                            "ref_ld_imputed_sumstats_by_ref_ld_region")

results_indicator_files = c("gwas_ld_true_sumstats_region",
                      "gwas_ld_imputed_sumstats_by_gwas_ld_region",
                      "gwas_ld_imputed_sumstats_by_ref_ld_region",
                      "ref_ld_true_sumstats_region",
                      "ref_ld_imputed_sumstats_by_gwas_ld_region",
                      "ref_ld_imputed_sumstats_by_ref_ld_region")

for (i in 1:length(results_indicator_files)) {
  print(sprintf("%s_\\w+_dt_%s_cred_set_%s_fdr.RData",results_fdr_files[i],divergence_time,signal))
  files_1           = list.files("pop_split/finemap", pattern = sprintf("%s_\\w+_dt_%s_cred_set_%s_fdr.RData",results_fdr_files[i],divergence_time,signal),full.names = T)
  files_2           = list.files("pop_split/finemap", pattern = sprintf("%s_\\w+_dt_%s_cred_set_%s_indicator.RData",results_indicator_files[i],divergence_time,signal),full.names = T)
  results_1         = lapply(files_1,readRDS,.GlobalEnv)
  results_2         = lapply(files_2,readRDS,.GlobalEnv)

  #Write out the pooled quantiles
  saveRDS(results_2,sprintf("pop_split/results/pooled_%s_indicator_dt_%s.RData",results_indicator_files[i], divergence_time), version = 2)
  saveRDS(results_1,sprintf("pop_split/results/pooled_%s_fdr_dt_%s.RData", results_fdr_files[i], divergence_time), version = 2)
}


