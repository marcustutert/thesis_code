#Code to run the LD inference across different parameters and populations
#This will be done on the segmented genome! Which means that some SNPs will have two sets of weights (figure out how to deal with this later....)

### Snakemake params ####
populations    = snakemake@params$populations
chains         = snakemake@params$chains
regions        = snakemake@params$regions
###########################

#First we want to source the code within the InferLD package (Double check it's up to date!)
source("InferLD/R/HMM_Functions.R")
source("InferLD/R/validation_pipeline_functions.R")

#Read in reference panel
reference_panel    = readRDS(sprintf("data/segemented_regions/%s_ref_hap_panel_region_%s", populations, regions))
summary_statistics = readRDS(sprintf("data/segemented_regions/%s_sumstats_region_%s", populations, regions))

#Run inference
mbar = 124090 #Change around if necessary!
results                             = LD_from_GSHMM(ref_allele_matrix     = reference_panel,
                                                    fst                   = 0.01,
                                                    alpha                 = 1e5,
                                                    nSamples              = 10,
                                                    weights_resolution    = 2000,
                                                    likelihood_toggle     = TRUE,
                                                    gwas_variance         = summary_statistics$SE,
                                                    LD_Infer              = FALSE,
                                                    recombination         = FALSE,
                                                    case_control_constant = mbar/2,
                                                    BurnIn                = 0.9,
                                                    debug                 = TRUE)
#Save results across the regions
saveRDS(results,sprintf("results/inference_%s_panel_region_%s_chain_%s",populations,regions,chains))

