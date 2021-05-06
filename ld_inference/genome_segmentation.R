#Code to do the genome chunking (including overlapping and merging of the SNPs)
#Outputs the regions and associated reference panel for each population and each window
genome_segmentation = function(window_size  = 2e6,   #How big is the actual window (Default 2Mb chunks)
                               overlap_size = 1.5e6, #Overapping regions (Default 0.5Mb chunks) can be 0!
                               min_SNPs     = 50,    #Minimum number of SNPs in a region pre-merging
                               pop_panel)            #This will tell let us know which pop_panel/matched sumstats are we loading in (BOLT format)
{
  
  #Load in the sumstats corresponding to the pop_panel of our choice
  sumstats = as.data.table(readRDS(sprintf("data/summary_stats_%s",pop_panel)))
  #Mb Chunks
  chunks                        = seq(1,sumstats$BP[nrow(sumstats)],window_size)
  offsetted_chunks              = seq(overlap_size,sumstats$BP[nrow(sumstats)],window_size) 
  
  #Find intervals closest to these positions in the summary stats vector
  start_positions    = findInterval(chunks,sumstats$BP)
  start_positions[1] = 1 #Replace first position
  end_positions      = findInterval(offsetted_chunks,sumstats$BP)
  #Add end row
  end_positions      = c(end_positions,nrow(sumstats))
  
  #Loop through the start and end positions vectors to generate list of regions
  regions_list = list()
  counter = 1
  for (i in 1:(length(start_positions)-1)) {
    regions_list[[counter]] = sumstats[start_positions[i]:start_positions[i+1],]
    counter = counter +1
    regions_list[[counter]] = sumstats[end_positions[i]:end_positions[i+1],]
    counter = counter +1
  }
  
  #Get SNPs per region
  nsnps_per_region   = unlist(lapply(regions_list,nrow))
  #Find problematic regions
  problematic_region = which(nsnps_per_region < min_SNPs)
  
  #Remove the problematic regions
  for (i in 1:length(problematic_region)) {
    regions_list[[problematic_region[i]]] <- NULL
  }
  
  #Write out these regions 
  for (i in 1:length(regions_list)) {
    #Write out the reference panel that corresponds with this as well
    #Loop through each of the regions to get the start and end position (note that the sumstats and reference panel will already be matched at this point in)
    #So hopefully no bugs will occur
    saveRDS(regions_list[[i]],sprintf("data/segemented_regions/%s_sumstats_region_%s", pop_panel, i))
    #Load in reference panel
    full_ref_panel    = readRDS(sprintf("data/reference_panel_%s",pop_panel))
    #Extract the SNP indexes, note that this means the rsids will be stored as the column names
    full_rsid         = colnames(full_ref_panel)
    #Find the rsids in each region
    indexes_kept      = which(full_rsid %in%  regions_list[[i]]$SNP)
    #Extract these SNPs from the reference panel
    subset_ref_panel  = full_ref_panel[,indexes_kept]
    #Save the panel
    saveRDS(subset_ref_panel,sprintf("data/segemented_regions/%s_ref_hap_panel_region_%s", pop_panel,i))
  }
}