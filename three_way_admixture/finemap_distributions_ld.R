#Run FINEMAP across the distributions of LD
#We will do this (eventually) with summary statistic imputation added in, but for now just using the true sumstats data
#We need to calculate the FINEMAP statistics across the true GWAS haplotypes, the different reference panels, including those combined- and my weighted reference panel;

#Read in snakemake params
regions = snakemake@params$regions

#Read in GWAS panel (this stays constant through the samples)
gwas_case_haplotypes     = t(as.matrix(read.table(sprintf("hapgen2/gwas_region_%s.cases.haps",regions, header = T)))) #Read in the GWAS haplotypes (cases and control)
gwas_control_haplotypes  = t(as.matrix(read.table(sprintf("hapgen2/gwas_region_%s.controls.haps",regions, header = T)))) #Read in the GWAS haplotypes (cases and control)
#Merge together the haplotypes
gwas_haplotypes = rbind(gwas_case_haplotypes,gwas_control_haplotypes)

#Read in the TRUE sumstats from SNPTEST
true_sumstats              = read.table(sprintf("snptest/sumstats_gwas_region_%s",regions), header = T, skip = 10, stringsAsFactors = F)


#Loop across the LD distribution samples
nSamples = 100
for (i in 1:nSamples) {
  #### Create master file ####
  master_colnames               = paste("z","ld","snp","config","cred","log","n_samples",sep = ";") #Get the colnames
  true_sumstats_inferred_ld     = paste(sprintf("finemap_distributions/inferred_ld_true_sumstats_region_%s_sample_%s.z", regions,i),
                                        sprintf("finemap_distributions/inferred_ld_true_sumstats_region_%s_sample_%s.ld", regions,i),
                                        sprintf("finemap_distributions/inferred_ld_true_sumstats_region_%s_sample_%s.snp", regions,i),
                                        sprintf("finemap_distributions/inferred_ld_true_sumstats_region_%s_sample_%s.config", regions,i),
                                        sprintf("finemap_distributions/inferred_ld_true_sumstats_region_%s_sample_%s.cred",regions,i),
                                        sprintf("finemap_distributions/inferred_ld_true_sumstats_region_%s_sample_%s.log",regions,i),
                                        10000,
                                        sep = ";")
  #Rbind all the data together
  write.table(rbind(master_colnames,
                    true_sumstats_inferred_ld),
              sprintf("finemap_distributions/gwas_region_%s_sample_%s_master",regions,i),col.names = F, row.names = F, quote = F)
  
  
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
  write.table(true_sumstats_z, file = sprintf("finemap_distributions/inferred_ld_true_sumstats_region_%s_sample_%s.z",regions,i),quote = F, row.names = F, col.names = T)
  
  #Read in the inferred LD panel:
  inferred_LD = readRDS(sprintf("distribution_ld_results/inferred_LD_sample_%s_region_%s",i,regions))
  write.table(inferred_LD, file = sprintf("finemap_distributions/inferred_ld_true_sumstats_region_%s_sample_%s.ld", regions,i), quote = F, col.names = F, row.names = F)

  #####RUN FINEMAP######
  #Do this across different configurations of reference panel LD and my inferred LD
  system(sprintf("/well/mcvean/mtutert/software/FINEMAP/finemap_v1.4_x86_64 --sss --in-files finemap_distributions/gwas_region_%s_sample_%s_master --n-causal-snps 3 --log",regions,i))
}

