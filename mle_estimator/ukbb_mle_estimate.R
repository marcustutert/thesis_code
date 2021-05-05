#Code to take in UK-Biobank GWAS data and run the Fst estimate using the different 1000G populations
library(data.table)

population                     = snakemake@params$population
region                         = snakemake@params$region
source("/well/mcvean/mtutert/thesis_code/thesis_code/mle_estimator/likelihood.R")
##Debugging code
region     = as.numeric(region)

# #Debugging
# population = "EUR"
# region     = 5

#First, read in the UKBB BOLT GWAS data
bolt_gwas_data = fread("/well/mcvean/ukbb12788/mtutert/test_genotyped_chr1", header = T)

#Extract the SNPs from the specific region (38 regions)
range = seq(from = region*1000, to = (region+1)*1000, by = 1)

#Extract out the rsIDs in this region:
bolt_snp_rsid = bolt_gwas_data$SNP[range]
bolt_range    = bolt_gwas_data[range,]
#Read in KG freq file
KG_freq = fread(sprintf("sample_pops/%s_freq.afreq",population),header = T)

#Match the SNPs
KG_filtered   = KG_freq[which(KG_freq$ID %in% bolt_snp_rsid)]
bolt_filtered = bolt_range[which(bolt_range$SNP  %in% KG_filtered$ID)]

#Get the list of pi_ref freqs and SEs
pi_ref        = as.numeric(KG_filtered$ALT_FREQS)

bolt_se       = bolt_filtered$SE
#Remove rare variants in either population
pi_ref_low_freq = which(pi_ref > 0.99 | pi_ref < 0.01 | is.na(pi_ref))

if (length(pi_ref_low_freq) > 0 ) {
  pi_ref  = pi_ref[-pi_ref_low_freq]
  bolt_se = bolt_se[-pi_ref_low_freq] 
}
#Calculate mbar (do this manually for now, will change later)
mbar = (2*(361381*0.2201582*(1-0.2201582)))
#Convert SE^2 to quanitative scale (BOLT Quirk)
bolt_se       = (bolt_se/((0.2201582*(1-0.2201582))))^2

#Run the MLE
pars = optim(par = c( 0.2 ,0.001), fn = gamma_likelihood, se = bolt_se, pi_ref = pi_ref, mbar = mbar,lower = c( 0.00001,0.0001), upper = c( 0.5,0.2),method = "L-BFGS-B")$par

#Save Rdata
saveRDS(object = pars[1],file = sprintf("ukbb_mle_results/fst_%s_panel_replicate_%s",population,region))
saveRDS(object = pars[2],file = sprintf("ukbb_mle_results/noise_%s_panel_replicate_%s",population,region))





