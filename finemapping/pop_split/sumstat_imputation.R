#This script will perform sumstat imputation on the SNPTEST files we will basically be just replacing the betas and standard errors
#Source the sumstat imputation function and LD matrix function
source("pop_split/helper_functions.R")
replicates       = snakemake@params$replicates
divergence_time  = snakemake@params$divergence_time
signal           = snakemake@params$signal
#We need to read in the sumstats file from SNPTEST
sumstats         = read.table(sprintf("pop_split/snptest/sumstats_gwas_%s_%s_dt_%s",signal,replicates, divergence_time), header = T, skip = 10, stringsAsFactors = F)
#Read in the reference panel
ref_panel        = read.table(sprintf("pop_split/msprime/filtered_panels/ref_%s_dt_%s", replicates,divergence_time), header = T)
#Read in the GWAS panel from the HAPGEN2 haplotypes (case & control) note we will need to transpose
HAPGEN2_case_haplotypes    = t(read.table(sprintf("pop_split/hapgen2/gwas_%s_%s_dt_%s.cases.haps", signal, replicates,divergence_time), header = F))
HAPGEN2_control_haplotypes = t(read.table(sprintf("pop_split/hapgen2/gwas_%s_%s_dt_%s.controls.haps", signal,replicates,divergence_time), header = F))
gwas_panel       = rbind(HAPGEN2_case_haplotypes,HAPGEN2_control_haplotypes)

#Create list of imputed SNPs to use, 50% of the data
nsnps               = nrow(sumstats)
imputed_snps_index  = sort(sample(1:nsnps, 0.5*nsnps, replace = F)) #Note that this refers to the INDEX and NOT the RSID
typed_snps_index    = setdiff(1:nsnps,imputed_snps_index)

#####Bit hacky here fix this later #####
#Check if any of our imputed SNPs are ones which have NA (means allele frequency will be zero)
if (any(is.na(sumstats[typed_snps_index,]$frequentist_add_beta_1))) {
  #Remove these SNPs from the typed list
  typed_snps_index   = typed_snps_index[-which(is.na(sumstats[typed_snps_index,]$frequentist_add_beta_1))]
  imputed_snps_index = setdiff(1:nsnps,typed_snps_index)
}

#Extract the LD from reference panel
ref_LD                    = LD_Matrix(ref_panel)
gwas_LD                   = LD_Matrix(gwas_panel)

#Perform the sumstat imputation
imputed_zs_ref_LD            = sumstat_impute(typed_snps         = typed_snps_index,
                                              untyped_snps_index = imputed_snps_index,
                                              genotyped_sumstats = sumstats[-imputed_snps_index,],
                                              imputed_sumstats   = sumstats[imputed_snps_index,],
                                               LD                 = ref_LD)

imputed_zs_gwas_LD            = sumstat_impute(typed_snps         = typed_snps_index,
                                              untyped_snps_index = imputed_snps_index,
                                              genotyped_sumstats = sumstats[-imputed_snps_index,],
                                              imputed_sumstats   = sumstats[imputed_snps_index,],
                                              LD                 = gwas_LD)

#Write out the results to /results directory as an Rdata object (Just for book marking purposes!)
saveRDS(sumstats[imputed_snps_index,2],sprintf("pop_split/sumstat_imputation/sumstat_imputed_rsids_gwas_%s_%s_dt_%s",signal,replicates,divergence_time), version =  2)
saveRDS(imputed_zs_ref_LD,sprintf("pop_split/sumstat_imputation/sumstat_imputed_results_ref_ld_%s_%s_dt_%s",signal,replicates, divergence_time), version =  2)
saveRDS(imputed_zs_gwas_LD,sprintf("pop_split/sumstat_imputation/sumstat_imputed_results_gwas_ld_%s_%s_dt_%s",signal,replicates, divergence_time), version =  2)

############################### Calculate the accuracy of these results for downstream plotting ###########################
#We want to get the r2 compared to the truth across each of the divergence times and replicates, then plot it LOCALLY

#First, we want to compare data to the truth data
true_sumstats           = read.table(sprintf("pop_split/snptest/sumstats_gwas_%s_%s_dt_%s",signal, replicates,divergence_time), header = T,skip = 1)
#Find the z-scores
true_zscores            = true_sumstats$frequentist_add_beta_1/true_sumstats$frequentist_add_se_1
#Find the low frequency variants
low_freq_variants       = which(colMeans(gwas_panel[,imputed_snps_index]) < 0.05 | colMeans(gwas_panel[,imputed_snps_index]) > 0.95)
#Now filter to the SNPs that were actually imputed
true_zscores            = true_zscores[imputed_snps_index]
#Compare this with the truth data across low freq and common variants

#Do this for the reference LD
true_imputed_zscores_r2_rare_variants_ref_LD = summary(lm(true_zscores[low_freq_variants]~imputed_zs_ref_LD[low_freq_variants], na.action = na.exclude))$r.squared
true_imputed_zscores_r2_common_variants_ref_LD = summary(lm(true_zscores[-low_freq_variants]~imputed_zs_ref_LD[-low_freq_variants], na.action = na.exclude))$r.squared
#Save these results to pool locally later
print(true_imputed_zscores_r2_rare_variants_ref_LD)
saveRDS(true_imputed_zscores_r2_rare_variants_ref_LD,  sprintf("pop_split/sumstat_imputation/results/true_imputed_zscores_r2_rare_variants_ref_LD_gwas_%s_%s_dt_%s.RData", signal, replicates, divergence_time), version = 2)
saveRDS(true_imputed_zscores_r2_common_variants_ref_LD,  sprintf("pop_split/sumstat_imputation/results/true_imputed_zscores_r2_common_variants_ref_LD_gwas_%s_%s_dt_%s.RData", signal, replicates, divergence_time), version = 2)

#Do this for the GWAS LD
true_imputed_zscores_r2_rare_variants_gwas_LD = summary(lm(true_zscores[low_freq_variants]~imputed_zs_gwas_LD[low_freq_variants], na.action = na.exclude))$r.squared
true_imputed_zscores_r2_common_variants_gwas_LD = summary(lm(true_zscores[-low_freq_variants]~imputed_zs_gwas_LD[-low_freq_variants], na.action = na.exclude))$r.squared
#Save these results to pool locally later
saveRDS(true_imputed_zscores_r2_rare_variants_gwas_LD,  sprintf("pop_split/sumstat_imputation/results/true_imputed_zscores_r2_rare_variants_gwas_LD_gwas_%s_%s_dt_%s.RData", signal, replicates, divergence_time), version = 2)
saveRDS(true_imputed_zscores_r2_common_variants_gwas_LD,  sprintf("pop_split/sumstat_imputation/results/true_imputed_zscores_r2_common_variants_gwas_LD_gwas_%s_%s_dt_%s.RData", signal, replicates, divergence_time), version = 2)


###########################################################################################################################


#Generate a new SNPTEST results file in snptest/ which contains the IMPUTED Zs in the Beta column and 1's in the SE column
#Do this across the GWAS results & the reference LD results
sumstats_ref_LD = sumstats
sumstats_gwas_LD = sumstats

sumstats_ref_LD$frequentist_add_beta_1[imputed_snps_index] = imputed_zs_ref_LD
sumstats_ref_LD$frequentist_add_se_1[imputed_snps_index]   = 1
#Write out this file
write.table(sumstats_ref_LD,sprintf("pop_split/snptest/sumstats_imputed_results_ref_LD_gwas_%s_%s_dt_%s",signal, replicates, divergence_time), quote = F, row.names = F, col.names = T)

sumstats_gwas_LD$frequentist_add_beta_1[imputed_snps_index] = imputed_zs_gwas_LD
sumstats_gwas_LD$frequentist_add_se_1[imputed_snps_index]   = 1
#Write out this file
write.table(sumstats_gwas_LD,sprintf("pop_split/snptest/sumstats_imputed_results_gwas_LD_gwas_%s_%s_dt_%s",signal, replicates, divergence_time), quote = F, row.names = F, col.names = T)


