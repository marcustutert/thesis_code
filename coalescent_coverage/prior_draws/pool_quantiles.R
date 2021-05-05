#Pool quantiles together and plot resulting histograms

#Input the snakemake parameters
divergence_time_fst_parameter      = snakemake@params$paired_values #Do this across all divergence/Fst values (grid of graphs in the end)

divergence_time = strsplit(divergence_time_fst_parameter, "_")[[1]][1]
fst_parameter   = strsplit(divergence_time_fst_parameter, "_")[[1]][2]

#Grab the results across the replicates
LD_Results_File      = list.files(path = "results/pop_split", pattern = sprintf("%s_split_LD_quantile", divergence_time_fst_parameter), full.names = T)
AF_Results_File      = list.files(path = "results/pop_split", pattern = sprintf("%s_split_AF_quantile", divergence_time_fst_parameter), full.names = T)

LD_Pooled = lapply(LD_Results_File,readRDS,.GlobalEnv)
AF_Pooled = lapply(AF_Results_File,readRDS,.GlobalEnv)

#Write out the pooled quantiles
saveRDS(LD_Pooled,sprintf("results/pop_split/%s_split_LD_pooled_quantile_counts.RData", divergence_time_fst_parameter), version = 2)
saveRDS(AF_Pooled,sprintf("results/pop_split/%s_split_AF_pooled_quantile_counts.RData", divergence_time_fst_parameter), version = 2)




