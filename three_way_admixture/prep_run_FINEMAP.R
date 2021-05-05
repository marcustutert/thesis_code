#Prep FINEMAP files and run FINEMAP
#We will do this (eventually) with summary statistic imputation added in, but for now just using the true sumstats data
#We need to calculate the FINEMAP statistics across the true GWAS haplotypes, the different reference panels, including those combined- and my weighted reference panel;

#Read in snakemake params
regions = snakemake@params$regions
source("/well/mcvean/mtutert/thesis_code/thesis_code/finemapping/helper_functions.R") #Get LD function to calculate from haplotypes including from the weights!

#Read in GWAS panel
gwas_case_haplotypes     = t(as.matrix(read.table(sprintf("hapgen2/gwas_region_%s.cases.haps",regions, header = T)))) #Read in the GWAS haplotypes (cases and control)
gwas_control_haplotypes  = t(as.matrix(read.table(sprintf("hapgen2/gwas_region_%s.controls.haps",regions, header = T)))) #Read in the GWAS haplotypes (cases and control)
#Merge together the haplotypes
gwas_haplotypes = rbind(gwas_case_haplotypes,gwas_control_haplotypes)

#Read in Ref panel
ref_pop_AFR_haplotypes    = as.matrix(read.table(sprintf("msprime_data/joint_filtered/ref_pop_AFR_region_%s", regions), header = T))
ref_pop_EUR_haplotypes    = as.matrix(read.table(sprintf("msprime_data/joint_filtered/ref_pop_EUR_region_%s", regions), header = T))
ref_pop_EAS_haplotypes    = as.matrix(read.table(sprintf("msprime_data/joint_filtered/ref_pop_EAS_region_%s", regions), header = T))


#Read in the TRUE sumstats from SNPTEST
true_sumstats              = read.table(sprintf("snptest/sumstats_gwas_region_%s",regions), header = T, skip = 10, stringsAsFactors = F)

#### Create master file ####
master_colnames = paste("z","ld","snp","config","cred","log","n_samples",sep = ";")

#Go across all possible configurations for the reference panels: 1) in-sample LD, 2) 1000G proportions, 3)Each reference population (x3), 5)Reference LD admixed proportions 6)1000G WEIGHTED version (TBD)

true_sumstats_in_sample_LD    = paste(sprintf("finemap/true_sumstats_region_%s.z", regions),
                                      sprintf("finemap/in_sample_ld_region_%s.ld", regions),
                                      sprintf("finemap/in_sample_ld_true_sumstats_region_%s.snp", regions),
                                      sprintf("finemap/in_sample_ld_true_sumstats_region_%s.config", regions),
                                      sprintf("finemap/in_sample_ld_true_sumstats_region_%s.cred",regions),
                                      sprintf("finemap/in_sample_true_sumstats_region_%s.log",regions),
                                      10000,
                                      sep = ";")


true_sumstats_1000G_proportions_LD    = paste(sprintf("finemap/true_sumstats_region_%s.z", regions),
                                      sprintf("finemap/1000G_proportions_ld_region_%s.ld", regions),
                                      sprintf("finemap/1000G_proportions_ld_true_sumstats_region_%s.snp", regions),
                                      sprintf("finemap/1000G_proportions_ld_true_sumstats_region_%s.config", regions),
                                      sprintf("finemap/1000G_proportions_ld_true_sumstats_region_%s.cred",regions),
                                      sprintf("finemap/1000G_proportions_ld_true_sumstats_region_%s.log",regions),
                                      10000,
                                      sep = ";")

true_sumstats_1000G_EUR_LD    = paste(sprintf("finemap/true_sumstats_region_%s.z", regions),
                                              sprintf("finemap/EUR_ld_region_%s.ld", regions),
                                              sprintf("finemap/EUR_ld_true_sumstats_region_%s.snp", regions),
                                              sprintf("finemap/EUR_ld_true_sumstats_region_%s.config", regions),
                                              sprintf("finemap/EUR_ld_true_sumstats_region_%s.cred",regions),
                                              sprintf("finemap/EUR_ld_true_sumstats_region_%s.log",regions),
                                              10000,
                                              sep = ";")

true_sumstats_1000G_AFR_LD    = paste(sprintf("finemap/true_sumstats_region_%s.z", regions),
                                      sprintf("finemap/AFR_ld_region_%s.ld", regions),
                                      sprintf("finemap/AFR_ld_true_sumstats_region_%s.snp", regions),
                                      sprintf("finemap/AFR_ld_true_sumstats_region_%s.config", regions),
                                      sprintf("finemap/AFR_ld_true_sumstats_region_%s.cred",regions),
                                      sprintf("finemap/AFR_ld_true_sumstats_region_%s.log",regions),
                                      10000,
                                      sep = ";")

true_sumstats_1000G_EAS_LD    = paste(sprintf("finemap/true_sumstats_region_%s.z", regions),
                                      sprintf("finemap/EAS_ld_region_%s.ld", regions),
                                      sprintf("finemap/EAS_ld_true_sumstats_region_%s.snp", regions),
                                      sprintf("finemap/EAS_ld_true_sumstats_region_%s.config", regions),
                                      sprintf("finemap/EAS_ld_true_sumstats_region_%s.cred",regions),
                                      sprintf("finemap/EAS_ld_true_sumstats_region_%s.log",regions),
                                      10000,
                                      sep = ";")

true_sumstats_1000G_Ref_ADMIXED_LD    = paste(sprintf("finemap/true_sumstats_region_%s.z", regions),
                                      sprintf("finemap/ref_admixed_ld_region_%s.ld", regions),
                                      sprintf("finemap/ref_admixed_ld_true_sumstats_region_%s.snp", regions),
                                      sprintf("finemap/ref_admixed_ld_true_sumstats_region_%s.config", regions),
                                      sprintf("finemap/ref_admixed_ld_true_sumstats_region_%s.cred",regions),
                                      sprintf("finemap/ref_admixed_ld_true_sumstats_region_%s.log",regions),
                                      10000,
                                      sep = ";")


true_sumstats_1000G_weighted_LD    = paste(sprintf("finemap/true_sumstats_region_%s.z", regions),
                                              sprintf("finemap/1000G_weighted_ld_region_%s.ld", regions),
                                              sprintf("finemap/1000G_weighted_ld_true_sumstats_region_%s.snp", regions),
                                              sprintf("finemap/1000G_weighted_ld_true_sumstats_region_%s.config", regions),
                                              sprintf("finemap/1000G_weighted_ld_true_sumstats_region_%s.cred",regions),
                                              sprintf("finemap/1000G_weighted_ld_true_sumstats_region_%s.log",regions),
                                              10000,
                                              sep = ";")

#Rbind all the data together
write.table(rbind(master_colnames,
                  true_sumstats_in_sample_LD,
                  true_sumstats_1000G_proportions_LD,
                  true_sumstats_1000G_EUR_LD,
                  true_sumstats_1000G_AFR_LD,
                  true_sumstats_1000G_EAS_LD,
                  true_sumstats_1000G_Ref_ADMIXED_LD,
                  true_sumstats_1000G_weighted_LD),
            sprintf("finemap/gwas_region_%s_master",regions),col.names = F, row.names = F, quote = F)

#### Create Z file ####
#To do this each of the 3 elements in these vectors will correspond to:
rsid       = replicate(1,true_sumstats$rsid)
chromosome = replicate(1,true_sumstats$chromosome)
position   = replicate(1,true_sumstats$position)
allele1    = replicate(1,true_sumstats$alleleA)
allele2    = replicate(1,true_sumstats$alleleB)
maf        = replicate(1,rep(0.420,length(true_sumstats$rsid))) #Nice
beta       = cbind(true_sumstats$frequentist_add_beta_1)
se         = cbind(true_sumstats$frequentist_add_se_1)

#Bind this data all together into a matrix for all three cases
true_sumstats_z = cbind(rsid[,1],chromosome[,1],position[,1],allele1[,1],allele2[,1],maf[,1],beta[,1],se[,1])
colnames(true_sumstats_z) = c("rsid","chromosome","position","allele1","allele2","maf","beta","se")

#Note that we will need to restrict any variants that are NA in the sumstats (because HAPGEN create far too low freq variants to pick up)
#Find a way to fix this!
low_freq_snps = which(is.na(beta[,1]))
#Remove this SNPs from all 3 z files
if (length(low_freq_snps)>0) {
  true_sumstats_z               = true_sumstats_z[-low_freq_snps,]
}

#Extract each row for each of the three cases
write.table(true_sumstats_z, file = sprintf("finemap/true_sumstats_region_%s.z",regions),quote = F, row.names = F, col.names = T)

#### Create LD #######
gwas_LD        = LD_Matrix(gwas_haplotypes)
write.table(gwas_LD, file = sprintf("finemap/in_sample_ld_region_%s.ld", regions), quote = F, col.names = F, row.names = F)
KG_Proportions_LD = LD_Matrix(rbind(ref_pop_AFR_haplotypes[1:1332,],ref_pop_EUR_haplotypes[1:1006,]))
write.table(KG_Proportions_LD, file = sprintf("finemap/1000G_proportions_ld_region_%s.ld", regions), quote = F, col.names = F, row.names = F)
EUR_LD         = LD_Matrix(ref_pop_EUR_haplotypes)
write.table(EUR_LD, file = sprintf("finemap/EUR_ld_region_%s.ld", regions), quote = F, col.names = F, row.names = F)
AFR_LD         = LD_Matrix(ref_pop_AFR_haplotypes)
write.table(AFR_LD, file = sprintf("finemap/AFR_ld_region_%s.ld", regions), quote = F, col.names = F, row.names = F)
EAS_LD         = LD_Matrix(ref_pop_EAS_haplotypes)
write.table(EAS_LD, file = sprintf("finemap/EAS_ld_region_%s.ld", regions), quote = F, col.names = F, row.names = F)
KG_Admixed_LD  = LD_Matrix(rbind(ref_pop_AFR_haplotypes[1:583,],ref_pop_EUR_haplotypes[1:1750,]))
write.table(KG_Admixed_LD, file = sprintf("finemap/ref_admixed_ld_region_%s.ld", regions), quote = F, col.names = F, row.names = F)

#Read in the inferred LD panel:
inferred_LD = as.matrix(readRDS(sprintf("distribution_ld_results/inferred_LD_sample_50_region_%s",regions)))
inferred_LD = matrix(data = inferred_LD, nrow = length(inferred_LD)^0.5, ncol = length(inferred_LD)^0.5)
write.table(inferred_LD, file = sprintf("finemap/1000G_weighted_ld_region_%s.ld", regions), quote = F, col.names = F, row.names = F)

# #Remove the variants that are low freq
# if (length(low_freq_snps)>0) {
#   gwas_LD = gwas_LD[,-low_freq_snps]
#   gwas_LD = gwas_LD[-low_freq_snps,]
#   ref_LD  = ref_LD[,-low_freq_snps]
#   ref_LD  = ref_LD[-low_freq_snps,]
# }
# 
# #Fix the LD bounds if necessary
# gwas_LD[gwas_LD < -1] = -1
# gwas_LD[gwas_LD > 1]  =  1
# ref_LD[ref_LD < -1]   = -1
# ref_LD[ref_LD > 1]    =  1

#####RUN FINEMAP######
#Do this across different configurations of reference panel LD and my inferred LD
system(sprintf("/well/mcvean/mtutert/software/FINEMAP/finemap_v1.4_x86_64 --sss --in-files finemap/gwas_region_%s_master --n-causal-snps 3 --log",regions))


