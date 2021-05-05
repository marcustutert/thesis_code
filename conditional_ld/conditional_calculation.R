library(data.table)

#Do this across replicates and pool the data later (will be quicker this way)
#SNAKEMAKE PARAMS
recombination_rate = snakemake@params$recombination_rate
divergence_time = 0

for (i in 1:1000) {
  print(i)
  #Read in target pop file
  target_pop = as.matrix(fread(sprintf("msprime/target_1_replicate_%s_dt_%s_recomb_%s.csv",i,divergence_time,recombination_rate),header = T))
  #Get constants
  nhaps_target   = nrow(target_pop) #Typically 1000
  nsnps_target   = ncol(target_pop) #This should only be two SNPs
  
  #Find the LD between the two SNPs in the target population
  target_pop_ld = cov.wt(target_pop,method = "ML", cor = TRUE)$cor[1,2]
  #Get list of the peripheral populations
  peripheral_list_names  = list.files(path = "msprime", pattern = sprintf("peripheral_\\d+_replicate_%s_dt_%s_recomb_%s.csv",i,divergence_time, recombination_rate), full.names = T)
  #Read out ONLY the SNPs that are in the target haplotype configuration
  peripheral_results     = lapply(peripheral_list_names, fread, header = T)
  peripheral_pop         = length(peripheral_results) 
  #Transform this into an array
  peripheral_array       = array(unlist(peripheral_results),dim = c(nhaps_target,2,peripheral_pop)) 

  #Remove times when the SNP is fixed (need to do this across both SNPs)
  pops_removed_snp_1 = which(snps_af_peripheral[1,] == 0)
  pops_removed_snp_2 = which(snps_af_peripheral[2,] == 0)
  
  #Unionize the indexes here
  pops_removed = c(pops_removed_snp_1,pops_removed_snp_2)
  
  #Deal with edge case where nothing matches
  if (length(pops_removed) > (peripheral_pop-1)){
    results = list(-4,-4,-4,-4)
    saveRDS(results,sprintf("results/results_replicate_%s_dt_%s_recomb_%s",i,divergence_time, recombination_rate))
  }else{
    if (length(pops_removed) == 0) {
      #Calculate the geometric mean of the SNPs in target pop
      gm_snps            = exp(mean(log(colMeans(target_pop))))
      #Calculate the LD between the SNPs
      ld_peripheral             = apply(peripheral_array, c(3), cov.wt,method = "ML", cor = TRUE)
      ld_cor_peripheral_matrix  = lapply(ld_peripheral, `[[`, 4)
      ld_peripheral_values      = c()
      #Loop through list and extract the values of ld we need
      for (j in 1:length(ld_cor_peripheral_matrix)) {
        ld_peripheral_values[j] = ld_cor_peripheral_matrix[[j]][1,2]
      }
      #Get the variance and mean of these values
      var_ld_peripheral  = var(ld_peripheral_values)
      mean_ld_peripheral = mean(ld_peripheral_values)
    }
    
    if (length(pops_removed) > 0) {
      #Calculate the geometric mean of the SNPs
      gm_snps            = exp(mean(log(colMeans(target_pop))))
      
      #Calculate the LD between the SNPs
      ld_peripheral             = apply(peripheral_array[,,-pops_removed], c(3), cov.wt,method = "ML", cor = TRUE)
      ld_cor_peripheral_matrix  = lapply(ld_peripheral, `[[`, 4)
      ld_peripheral_values      = c()
      #Loop through list and extract the values of ld we need
      for (j in 1:length(ld_cor_peripheral_matrix)) {
        ld_peripheral_values[j] = ld_cor_peripheral_matrix[[j]][1,2]
      }
      
      #Get the variance and mean of these values
      var_ld_peripheral  = var(ld_peripheral_values)
      mean_ld_peripheral = mean(ld_peripheral_values)
    }
    
    #Save results
    results = list(gm_snps,target_pop_ld,mean_ld_peripheral,var_ld_peripheral)
    saveRDS(results,sprintf("results/results_replicate_%s_dt_%s_recomb_%s",i,divergence_time,recombination_rate))
  }
}
