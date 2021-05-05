#Plot results
library(plotly)
#Don't hard code this in, import it directly from snakemake
divergence_times = c('0','10','50',"125")

#Read in results file (across all divergence times)
AF      = list()
AF_Plot = list()
LD      = list()

for (i in 1:4) {
  AF[[i]] = readRDS(sprintf("results/split_%s_AF_pooled_quantile_counts.RData",divergence_times[i]))
  LD[[i]] = readRDS(sprintf("results/split_%s_LD_pooled_quantile_counts.RData",divergence_times[i]))
  #Plot results
  ld_plot  = unlist(LD[[i]])
  ld_plot  = ld_plot[!is.na(ld_plot)]
  LD[[i]] = plot_ly(x = ~ld_plot, type = "histogram")
  export(LD[[i]], file = sprintf("results/%s_LD_result.jpeg", divergence_times[i]))
  export(unlist(AF[[i]]), file = sprintf("results/%s_AF_result.jpeg", divergence_times[i]))
  
}

