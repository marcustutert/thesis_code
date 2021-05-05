#pool_summary_stat_imputation.R

#Pool summary stat results together
divergence_time = snakemake@params$divergence_time
signal          = snakemake@params$signal

#We need to do this across all common variants, rare variants and within gwas LD/ref LD (4 ways)

ref_ld_common         = list.files("pop_split/sumstat_imputation/results", pattern = sprintf("common_variants_ref_LD_gwas_%s_\\d+\\_dt_%s.RData",signal,divergence_time), full.names = T)
results = lapply(ref_ld_common,readRDS,.GlobalEnv)
#Write out the pooled quantiles
saveRDS(results,sprintf("pop_split/sumstat_imputation/results/pooled_common_variants_%s_ref_LD_dt_%s",signal,divergence_time), version = 2)

ref_ld_rare           = list.files("pop_split/sumstat_imputation/results", pattern = sprintf("rare_variants_ref_LD_gwas_%s_\\d+\\_dt_%s.RData",signal,divergence_time), full.names = T)
results = lapply(ref_ld_rare,readRDS,.GlobalEnv)
#Write out the pooled quantiles
saveRDS(results,sprintf("pop_split/sumstat_imputation/results/pooled_rare_variants_%s_ref_LD_dt_%s",signal,divergence_time), version = 2)

gwas_ld_common        = list.files("pop_split/sumstat_imputation/results", pattern = sprintf("common_variants_gwas_LD_gwas_%s_\\d+\\_dt_%s.RData",signal,divergence_time), full.names = T)
results = lapply(gwas_ld_common,readRDS,.GlobalEnv)
#Write out the pooled quantiles
saveRDS(results,sprintf("pop_split/sumstat_imputation/results/pooled_common_variants_%s_gwas_LD_dt_%s",signal,divergence_time), version = 2)

gwas_ld_rare          = list.files("pop_split/sumstat_imputation/results", pattern = sprintf("rare_variants_gwas_LD_gwas_%s_\\d+\\_dt_%s.RData",signal,divergence_time), full.names = T)
results = lapply(gwas_ld_rare,readRDS,.GlobalEnv)
#Write out the pooled quantiles
saveRDS(results,sprintf("pop_split/sumstat_imputation/results/pooled_rare_variants_%s_gwas_LD_dt_%s",signal,divergence_time), version = 2)


