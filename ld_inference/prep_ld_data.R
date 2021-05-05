#Prep LD inference data
#Source helper script
source("helper_functions.R")
source("genome_segmentation.R")
populations = list("EUR","AFR",c("EUR","AFR")) #List of 1000G populations used as a reference panel (this we will need ot change)

#Read in summary statistics from BOLT
sumstats = fread("/well/mcvean/ukbb12788/mtutert/test_genotyped_chr1", header = T)
#Loop through the different populations (or pairs thereof) in the list
for (i in 1:length(populations)) {
  #Run the LD prep function
  print(populations[[i]])
  matched_data = prep_ld_inference_sumstats(populations = populations[[i]],
                                            sumstats    = sumstats)
  #Write out results to data folder
  
  #Deal with the multiple panels in a smart way
  if (length(populations[[i]] > 1)){
    populations[[i]] = paste(populations[[i]], collapse = '_')
  }
  saveRDS(matched_data[[1]],sprintf("data/reference_panel_%s", populations[[i]]))
  saveRDS(matched_data[[2]],sprintf("data/summary_stats_%s", populations[[i]]))
  }

#Match SNPs between the different populations to get unique set between all populations
#Read in all the reference panels to find the SNPs to remove
reference_panels_file_names = list.files(pattern = "reference_panel", path = "data",all.files = T,full.names = T)
snp_index_to_remove = c()
for (i in 1:length(reference_panels_file_names)) {
  haps = readRDS(file = reference_panels_file_names[i])
  freq = colMeans(haps)
  snp_index_to_remove = c(snp_index_to_remove,which(freq < 0.01 | freq > 0.99))
}

snp_index_to_remove = unique(snp_index_to_remove)
#Go through each of the reference panels and summary statistics and remove these SNPs then re-save the data
for (i in 1:length(reference_panels_file_names)) {
  haps           = readRDS(file = reference_panels_file_names[i])
  haps           = haps[,-c(snp_index_to_remove)]
  #Add the column names as the rsids
  colnames(haps) = matched_data[[2]]$SNP[-c(snp_index_to_remove)]
  saveRDS(haps,reference_panels_file_names[i])
}

summary_statistic_file_names = list.files(pattern = "summary_stats", path = "data",all.files = T,full.names = T)
for (i in 1:length(summary_statistic_file_names)) {
  sumstats = readRDS(file = summary_statistic_file_names[i])
  sumstats = sumstats[-c(snp_index_to_remove),]
  saveRDS(sumstats,summary_statistic_file_names[i])
}

#All reference panels and summary statistics should now be filtered and contain the same SNPs (across panels)

#Now perform the genome segmentation as well across all the populations
for (i in 1:length(populations)) {
  genome_segmentation(pop_panel = populations[[i]])
}


