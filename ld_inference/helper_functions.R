#Create a function which takes in possibly one or more multiple populations and performs the input filtering for our method
prep_ld_inference_sumstats = function(populations, #Vector (can be of length of more than one) -> Currently only works for the 1000G populations
                                      sumstats,    #Summary statistics used as input to the method--Assumes BOLT input (SNP == rsiD),
                                      chr){        #Only works for a single chr (currently just chr1)
  #Extract out the separate populations from the populations vector, we will wish to generate a table of the joint samples
  samples = c()
  for (i in 1:length(populations)) { #Loop through the 1000G populations
    #Read in the samples
    pop_samples = fread(sprintf("/well/mcvean/mtutert/thesis_code/thesis_code/ld_inference/sample_pops/%s_samples",populations[i]),header = F)
    samples = c(samples,pop_samples$V1)
  }
  #Write out samples file
  write.table(samples,"/well/mcvean/mtutert/thesis_code/thesis_code/ld_inference/temp/samples", quote = F,row.names = F, col.names = F)
  
  #Extract the SNPs from the sumstats file
  bolt_snps = sumstats$SNP
  #Write out SNPs
  write.table(bolt_snps,"/well/mcvean/mtutert/thesis_code/thesis_code/ld_inference/temp/sumstats_snps", quote = F,row.names = F, col.names = F)
  
  #Use PLINK2 to extract these samples and SNPs from the 1000G files
  system("/well/mcvean/mtutert/software/plink2 --pfile /well/mcvean/mtutert/1000_Genomes/plink/chr1_1000G --export haps --out /well/mcvean/mtutert/thesis_code/thesis_code/ld_inference/panels/population_panel --extract /well/mcvean/mtutert/thesis_code/thesis_code/ld_inference/temp/sumstats_snps --keep /well/mcvean/mtutert/thesis_code/thesis_code/ld_inference/temp/samples --max-alleles 2")
  
  #Read in haplotype file
  haps            = as.matrix(fread("/well/mcvean/mtutert/thesis_code/thesis_code/ld_inference/panels/population_panel.haps", header = T,stringsAsFactors = T))
  #Extract the rsids from the 1000G file
  ref_panel_rsids = haps[,2]
  #Cut away naming columns and transpose to get nhaps x nsnps
  haps            = t(haps[,6:ncol(haps)])
  haps            = `class<-`(haps, 'numeric')
  
  #Do the matching of SNPs w MAF filtering, note that these SNPs will have already been filtered
  #KG_filtered   = ref_panel_rsids[which(ref_panel_rsids %in% bolt_snps)]
  bolt_filtered  = sumstats[which(sumstats$SNP  %in% ref_panel_rsids)]
  
  #Convert the bolt SEs onto odds scale
  bolt_filtered$SE       = (bolt_filtered$SE/(0.2201582*(1-0.2201582)))^2
  
  #Ref allele freqs
  pi_ref = colMeans(haps)
  
  #Write out from this function the reference allele frequencies, the GWAS allele frequencies, the GWAS standard errors
  print(dim(haps))
  return(list(haps,
              bolt_filtered))
}