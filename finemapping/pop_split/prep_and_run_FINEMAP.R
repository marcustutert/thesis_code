#Prep FINEMAP files and run FINEMAP
#We want to run this across both the summary stat imputed data, as well as with the GWAS LD and the Reference LD
#We want to do this across all reference replicates and across all divergence times for a given TRUE signal
#For a given replicate (and divergence time) there are 6 different fine-mappings that will need to be done (3 choose 2):

# 1) GWAS LD & True sumstats
# 2) GWAS LD & Imputed (GWAS LD) sumstats
# 3) GWAS LD & Imputed (Ref LD) sumstats
# 4) Ref LD & Imputed (Ref LD) sumstats
# 5) Ref LD & True sumstats
# 6) Ref LD & Imputed (GWAS LD) sumstats

#Read in snakemake params
divergence_time = snakemake@params$divergence_time
replicates      = snakemake@params$replicates
signal          = snakemake@params$signal
print("here")
source("helper_functions.R") #Get LD function

#Read in the reference panel and the gwas panel (necessary for the LD input into FINEMAP)
gwas_panel = as.matrix(read.table(sprintf("pop_split/msprime/filtered_panels/gwas_%s_dt_%s",replicates,divergence_time), header = T)) #Read in the GWAS panel
ref_panel  = as.matrix(read.table(sprintf("pop_split/msprime/filtered_panels/ref_%s_dt_%s",replicates,divergence_time), header = T))

#Read in the TRUE sumstats from SNPTEST
true_sumstats              = read.table(sprintf("pop_split/snptest/sumstats_gwas_%s_%s_dt_%s",signal,replicates,divergence_time), header = T, skip = 10, stringsAsFactors = F)
#Read in the IMPUTED sumstats calculated from the previous step (w both the GWAS LD and the Reference LD)
imputed_sumstats_gwas_ld   = read.table(sprintf("pop_split/snptest/sumstats_imputed_results_gwas_LD_gwas_%s_%s_dt_%s",signal,replicates,divergence_time), header = T,stringsAsFactors = F)
imputed_sumstats_ref_ld    = read.table(sprintf("pop_split/snptest/sumstats_imputed_results_ref_LD_gwas_%s_%s_dt_%s",signal,replicates,divergence_time), header = T,stringsAsFactors = F)
print("here2")

#### Create master file ####
master_colnames = paste("z","ld","snp","config","cred","log","n_samples",sep = ";")

# 1) GWAS LD & True sumstats
gwas_ld_true_sumstats     = paste(sprintf("pop_split/finemap/true_sumstats_region_%s_dt_%s.z", replicates, divergence_time),
                                  sprintf("pop_split/finemap/gwas_%s_dt_%s.ld", replicates,divergence_time),
                                  sprintf("pop_split/finemap/gwas_ld_true_sumstats_region_%s_dt_%s.snp", replicates,divergence_time),
                                  sprintf("pop_split/finemap/gwas_ld_true_sumstats_region_%s_dt_%s.config", replicates,divergence_time),
                                  sprintf("pop_split/finemap/gwas_ld_true_sumstats_region_%s_dt_%s.cred",replicates,divergence_time),
                                  sprintf("pop_split/finemap/gwas_ld_true_sumstats_region_%s_dt_%s_%s.log",replicates,divergence_time,signal),
                                  10000,
                                  sep = ";")

# 2) GWAS LD & Imputed (GWAS LD) sumstats
gwas_ld_imputed_sumstats_by_gwas_ld     = paste(sprintf("pop_split/finemap/imputed_sumstats_gwas_ld_region_%s_dt_%s.z", replicates, divergence_time),
                                                sprintf("pop_split/finemap/gwas_%s_dt_%s.ld", replicates,divergence_time),
                                                sprintf("pop_split/finemap/gwas_ld_imputed_sumstats_by_gwas_ld_region_%s_dt_%s.snp", replicates,divergence_time),
                                                sprintf("pop_split/finemap/gwas_ld_imputed_sumstats_by_gwas_ld_region_%s_dt_%s.config", replicates,divergence_time),
                                                sprintf("pop_split/finemap/gwas_ld_imputed_sumstats_by_gwas_ld_region_%s_dt_%s.cred",replicates,divergence_time),
                                                sprintf("pop_split/finemap/gwas_ld_imputed_sumstats_by_gwas_ld_region_%s_dt_%s_%s.log",replicates,divergence_time,signal),
                                                10000,
                                                sep = ";")

# 3) GWAS LD & Imputed (Ref LD) sumstats
gwas_ld_imputed_sumstats_by_ref_ld  = paste(sprintf("pop_split/finemap/imputed_sumstats_ref_ld_region_%s_dt_%s.z", replicates, divergence_time),
                                            sprintf("pop_split/finemap/gwas_%s_dt_%s.ld", replicates,divergence_time),
                                            sprintf("pop_split/finemap/gwas_ld_imputed_sumstats_by_ref_ld_region_%s_dt_%s.snp", replicates,divergence_time),
                                            sprintf("pop_split/finemap/gwas_ld_imputed_sumstats_by_ref_ld_region_%s_dt_%s.config", replicates,divergence_time),
                                            sprintf("pop_split/finemap/gwas_ld_imputed_sumstats_by_ref_ld_region_%s_dt_%s.cred",replicates,divergence_time),
                                            sprintf("pop_split/finemap/gwas_ld_imputed_sumstats_by_ref_ld_region_%s_dt_%s_%s.log",replicates,divergence_time,signal),
                                            10000,
                                            sep = ";")
# 4) Ref LD & Imputed (Ref LD) sumstats
ref_ld_imputed_sumstats_by_ref_ld    = paste(sprintf("pop_split/finemap/imputed_sumstats_ref_ld_region_%s_dt_%s.z", replicates, divergence_time),
                                             sprintf("pop_split/finemap/ref_%s_dt_%s.ld", replicates,divergence_time),
                                             sprintf("pop_split/finemap/ref_ld_imputed_sumstats_by_ref_ld_region_%s_dt_%s.snp", replicates,divergence_time),
                                             sprintf("pop_split/finemap/ref_ld_imputed_sumstats_by_ref_ld_region_%s_dt_%s.config", replicates,divergence_time),
                                             sprintf("pop_split/finemap/ref_ld_imputed_sumstats_by_ref_ld_region_%s_dt_%s.cred",replicates,divergence_time),
                                             sprintf("pop_split/finemap/ref_ld_imputed_sumstats_by_ref_ld_region_%s_dt_%s_%s.log",replicates,divergence_time,signal),
                                             10000,
                                             sep = ";")

# 5) Ref LD & True sumstats
ref_ld_true_sumstats     =  paste(sprintf("pop_split/finemap/true_sumstats_region_%s_dt_%s.z", replicates, divergence_time),
                                  sprintf("pop_split/finemap/ref_%s_dt_%s.ld", replicates,divergence_time),
                                  sprintf("pop_split/finemap/ref_ld_true_sumstats_region_%s_dt_%s.snp", replicates,divergence_time),
                                  sprintf("pop_split/finemap/ref_ld_true_sumstats_region_%s_dt_%s.config", replicates,divergence_time),
                                  sprintf("pop_split/finemap/ref_ld_true_sumstats_region_%s_dt_%s.cred",replicates,divergence_time),
                                  sprintf("pop_split/finemap/ref_ld_true_sumstats_region_%s_dt_%s_%s.log",replicates,divergence_time,signal),
                                  10000,
                                  sep = ";")

# 6) Ref LD & Imputed (GWAS LD) sumstats
ref_ld_imputed_sumstats_by_gwas_ld     = paste(sprintf("pop_split/finemap/imputed_sumstats_gwas_ld_region_%s_dt_%s.z", replicates, divergence_time),
                                               sprintf("pop_split/finemap/ref_%s_dt_%s.ld", replicates,divergence_time),
                                               sprintf("pop_split/finemap/ref_ld_imputed_sumstats_by_gwas_ld_region_%s_dt_%s.snp", replicates,divergence_time),
                                               sprintf("pop_split/finemap/ref_ld_imputed_sumstats_by_gwas_ld_region_%s_dt_%s.config", replicates,divergence_time),
                                               sprintf("pop_split/finemap/ref_ld_imputed_sumstats_by_gwas_ld_region_%s_dt_%s.cred",replicates,divergence_time),
                                               sprintf("pop_split/finemap/ref_ld_imputed_sumstats_by_gwas_ld_region_%s_dt_%s_%s.log",replicates,divergence_time,signal),
                                               10000,
                                               sep = ";")

#Rbind all the data together
write.table(rbind(master_colnames,
                  gwas_ld_true_sumstats,
                  gwas_ld_imputed_sumstats_by_gwas_ld,
                  gwas_ld_imputed_sumstats_by_ref_ld,
                  ref_ld_imputed_sumstats_by_ref_ld,
                  ref_ld_true_sumstats,
                  ref_ld_imputed_sumstats_by_gwas_ld
),sprintf("pop_split/finemap/region_%s_dt_%s_master",replicates, divergence_time),col.names = F, row.names = F, quote = F)

#### Create Z file ####
#To do this each of the 3 elements in these vectors will correspond to:
# 1) TRUE SUMSTATS 2) IMPUTED SUMSTATS (GWAS LD) 3) IMPUTED SUMSTATS (REF LD)
#EVERYTHING WILL STAY THE SAME EXCEPT THE FOR THE BETA AND SEs
rsid       = replicate(3,true_sumstats$rsid)
chromosome = replicate(3,true_sumstats$chromosome)
position   = replicate(3,true_sumstats$position)
allele1    = replicate(3,true_sumstats$alleleA)
allele2    = replicate(3,true_sumstats$alleleB)
maf        = replicate(3,rep(0.420,length(true_sumstats$rsid))) #Nice
beta       = cbind(true_sumstats$frequentist_add_beta_1,imputed_sumstats_gwas_ld$frequentist_add_beta_1,imputed_sumstats_ref_ld$frequentist_add_beta_1)
se         = cbind(true_sumstats$frequentist_add_se_1,imputed_sumstats_gwas_ld$frequentist_add_se_1,imputed_sumstats_ref_ld$frequentist_add_se_1)

#Bind this data all together into a matrix for all three cases
true_sumstats_z = cbind(rsid[,1],chromosome[,1],position[,1],allele1[,1],allele2[,1],maf[,1],beta[,1],se[,1])
colnames(true_sumstats_z) = c("rsid","chromosome","position","allele1","allele2","maf","beta","se")
imputed_sumstats_by_gwas_ld_z = cbind(rsid[,2],chromosome[,2],position[,2],allele1[,2],allele2[,2],maf[,2],beta[,2],se[,2])
colnames(imputed_sumstats_by_gwas_ld_z) = c("rsid","chromosome","position","allele1","allele2","maf","beta","se")
imputed_sumstats_by_ref_ld_z = cbind(rsid[,3],chromosome[,3],position[,3],allele1[,3],allele2[,3],maf[,3],beta[,3],se[,3])
colnames(imputed_sumstats_by_ref_ld_z) = c("rsid","chromosome","position","allele1","allele2","maf","beta","se")

#Note that we will need to restrict any variants that are NA in the sumstats (because HAPGEN create far too low freq variants to pick up)
#Find a way to fix this!
low_freq_snps = which(is.na(beta[,1]))
#Remove this SNPs from all 3 z files
if (length(low_freq_snps)>0) {
  true_sumstats_z               = true_sumstats_z[-low_freq_snps,]
  imputed_sumstats_by_gwas_ld_z = imputed_sumstats_by_gwas_ld_z[-low_freq_snps,]
  imputed_sumstats_by_ref_ld_z  = imputed_sumstats_by_ref_ld_z[-low_freq_snps,]
}
#No need to extract signal (need this to be TRUE)
#Extract each row for each of the three cases
write.table(true_sumstats_z, file = sprintf("pop_split/finemap/true_sumstats_region_%s_dt_%s.z",replicates, divergence_time),quote = F, row.names = F, col.names = T)
write.table(imputed_sumstats_by_gwas_ld_z, file = sprintf("pop_split/finemap/imputed_sumstats_gwas_ld_region_%s_dt_%s.z",replicates, divergence_time),quote = F, row.names = F, col.names = T)
write.table(imputed_sumstats_by_ref_ld_z, file = sprintf("pop_split/finemap/imputed_sumstats_ref_ld_region_%s_dt_%s.z",replicates, divergence_time),quote = F, row.names = F, col.names = T)

#### Create LD #######
#This is done across the reference and GWAS panel
gwas_LD = LD_Matrix(gwas_panel)
ref_LD  = LD_Matrix(ref_panel)

#Remove the variants that are low freq
if (length(low_freq_snps)>0) {
  gwas_LD = gwas_LD[,-low_freq_snps]
  gwas_LD = gwas_LD[-low_freq_snps,]
  ref_LD  = ref_LD[,-low_freq_snps]
  ref_LD  = ref_LD[-low_freq_snps,]
}

gwas_LD[gwas_LD < -1] = -1
gwas_LD[gwas_LD > 1]  =  1
ref_LD[ref_LD < -1]   = -1
ref_LD[ref_LD > 1]    =  1

write.table(gwas_LD, file = sprintf("pop_split/finemap/gwas_%s_dt_%s.ld", replicates, divergence_time), quote = F, col.names = F, row.names = F)
write.table(ref_LD, file = sprintf("pop_split/finemap/ref_%s_dt_%s.ld", replicates, divergence_time), quote = F, col.names = F, row.names = F)

#####RUN FINEMAP######
system(sprintf("/well/mcvean/mtutert/software/FINEMAP/finemap_v1.4_x86_64 --sss --in-files pop_split/finemap/region_%s_dt_%s_master --n-causal-snps 3 --log",replicates, divergence_time))

