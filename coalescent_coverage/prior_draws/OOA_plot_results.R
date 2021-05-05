#Plot the results
#Plot results
library(plotly)
#Don't hard code this in, import it directly from snakemake
#divergence_times = c('0','10','50',"125")
population_pairs  = c("CEU_GWAS_CHB_Ref",
                     "CEU_GWAS_YRI_Ref",
                     "CHB_GWAS_CEU_Ref",
                     "CHB_GWAS_YRI_Ref",
                     "YRI_GWAS_CEU_Ref",
                     "YRI_GWAS_CHB_Ref")

#Read in results file (across all divergence times)
AF      = list()
AF_Plot = list()
LD      = list()

for (i in 1:length(population_pairs)) {
  AF[[i]] = readRDS(sprintf("results/OOA/%s_AF_pooled_quantile_counts.RData",population_pairs[i]))
  LD[[i]] = readRDS(sprintf("results/OOA/%s_LD_pooled_quantile_counts.RData",population_pairs[i]))
  #Plot results
  ld_plot  = unlist(LD[[i]])
  ld_plot  = ld_plot[!is.na(ld_plot)]
  LD[[i]] = plot_ly(x = ~ld_plot, type = "histogram")
  export(LD[[i]], file = sprintf("results/OOA/figures/%s_LD_result.jpeg", population_pairs[i]))
  af_plot  = unlist(AF[[i]])
  AF[[i]] = plot_ly(x = ~af_plot, type = "histogram")
  export(AF[[i]] , file = sprintf("results/OOA/figures/%s_AF_result.jpeg", population_pairs[i]))
}

