#Calculate size of credible set given the results from FINEMAP
#Do this across all 6 FINEMAP configurations outlined in the previous rule
#True causal SNPs (position)

#Read in snakemake params
replicates      = snakemake@params$replicates
divergence_time = snakemake@params$divergence_time
signal          = snakemake@params$signal

#Load in table that contains the positions of the causal variants from HAPGEN2
causal_snps            = read.table(sprintf("pop_split/hapgen2/snps_gwas_%s_causal_dt_%s",replicates,divergence_time), header = F)
snp_1_causal_position  = causal_snps[1]
snp_2_causal_position  = causal_snps[2]
snp_3_causal_position  = causal_snps[3]

#Loop through all possible configurations for the GWAS FINEMAPPING (See previous rule for details)
finemapped_file_configs = c("pop_split/finemap/gwas_ld_true_sumstats_region_%s_dt_%s",
                            "pop_split/finemap/gwas_ld_imputed_sumstats_by_gwas_ld_region_%s_dt_%s",
                            "pop_split/finemap/gwas_ld_imputed_sumstats_by_ref_ld_region_%s_dt_%s",
                            "pop_split/finemap/ref_ld_true_sumstats_region_%s_dt_%s",
                            "pop_split/finemap/ref_ld_imputed_sumstats_by_gwas_ld_region_%s_dt_%s",
                            "pop_split/finemap/ref_ld_imputed_sumstats_by_ref_ld_region_%s_dt_%s")

for (i in 1:length(finemapped_file_configs)) {
  #Read in files
  finemap_snp_results    = read.table(sprintf(paste(finemapped_file_configs[i],".snp",sep =""), replicates,divergence_time), header = T, stringsAsFactors = F)
  finemap_config_results = read.table(sprintf(paste(finemapped_file_configs[i],".config",sep =""), replicates,divergence_time), header = T, stringsAsFactors = F)
  #Go through the config results until we have reached the 95% credible set
  culmative_prob = cumsum(finemap_config_results$prob)
  index_cred_set = min(which(culmative_prob > 0.95))

  #Find RSID of the causal SNP index
  causal_snp_rsids = finemap_snp_results$rsid[which(finemap_snp_results$position %in% c(snp_1_causal_position,snp_2_causal_position,snp_3_causal_position))]

  #Loop through the rows in the credible set configs
  pp_weighted_tp = c()
  cred_set_snps  = c()
  for (j in 1:(index_cred_set)) {
    #Check how many of the causal variants exist
    tp_in_set         = length(which(unlist(strsplit(finemap_config_results$config[j],",")) %in% as.character(causal_snp_rsids)))
    frac_tp           = tp_in_set/(length(unlist(strsplit(finemap_config_results$config[j],","))))
    pp_weighted_tp[j] = frac_tp * finemap_config_results$prob[j]
    cred_set_snps     = append(cred_set_snps,unlist(strsplit(finemap_config_results$config[j],",")))
  }
  sum_tp_pp = sum(pp_weighted_tp)
  indicator = length(which(unique(cred_set_snps) %in% causal_snp_rsids))

  #Save the results
  saveRDS(sum_tp_pp,sprintf(paste(finemapped_file_configs[i],"_cred_set_%s_fdr.RData",sep=""),replicates,divergence_time,signal), version = 2)
  saveRDS(indicator,sprintf(paste(finemapped_file_configs[i],"_cred_set_%s_indicator.RData",sep=""),replicates,divergence_time,signal), version = 2)

}

