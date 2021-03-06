file = snakemake@params[[1]][1]
nsnps = 100
source("/well/mcvean/mtutert/thesis/coalescent_coverage/helper_functions.R")
library(data.table)
library(InferLD)
#Read in coalescent simulations
gwas = as.matrix(fread(sprintf("./msprime_data/population_split/GWAS_panel_replicate_%s_split_0.csv",file), header = T))
ref  = as.matrix(fread(sprintf("./msprime_data/population_split/GWAS_panel_replicate_%s_split_0.csv",file), header = T))

#Generate the sumstats with glm() in R and perform necessary filtering
res           = msprime_gwas_sumstats(gwas_haplotypes = gwas, reference_haplotypes = ref)
se            = res[[1]]
gwas          = res[[2]]
ref           = res[[3]]
case_control_constant = 2*(nrow(gwas)*0.5*(1-0.5))
#Run the inference
print("Starting Inference")
inference_results_1000G = LD_from_GSHMM(ref_panel_haplotypes   = ref[,1:nsnps],
                                        fst                    = 0.01,
                                        betas                  = FALSE,
                                        alpha                  = 100,
                                        nSamples               = 3,
                                        recomb_rate            = 1e-300,
                                        weights_resolution     = 5,
                                        likelihood_toggle      = FALSE,
                                        se_observed            = se[1:nsnps,3]^2,
                                        LD_Infer               = TRUE,
                                        genetic_map            = FALSE,
                                        chain_likelihood       = TRUE,
                                        nChains                = 1,
                                        recombination          = FALSE,
                                        case_control_constant  = case_control_constant,
                                        BurnIn                 = TRUE)
saveRDS(inference_results_1000G, sprintf("./inference_results/inference_results_rep_%s",file))
print("Inference Done")
#Find which quantile the true allele frequency is in
true_quantile_location = c()
nquantiles             = 100
nsnps                  = dim(inference_results_1000G$inferred_af_given_weights)[1]
for (i in 1:nsnps) {
  #Get the quantiled data
  data               = inference_results_1000G$inferred_af_given_weights[i,]
  quantiles = split(data, cut(data, quantile(data, prob = 0:nquantiles / nquantiles, names = FALSE), include = TRUE))
  intervals = names(quantiles)
  #Grep out the second element of each name to get the overlapping sets
  boundaries = c()
  for (j in 1:length(intervals)) {
    boundaries[j]= as.numeric(substr(sub(".*,", "", intervals[j]) ,1,nchar(sub(".*,", "", intervals[j]) )-1))
  }
  #Now look for where in the data does our true value lie
  pi_pop                 = (colMeans(gwas)/2)[i]
  true_quantile_location[i] = which.min(abs(boundaries - pi_pop))
  print(true_quantile_location[i])
}
saveRDS(true_quantile_location, sprintf("./inference_results/true_quantiles_rep_%s",file))


