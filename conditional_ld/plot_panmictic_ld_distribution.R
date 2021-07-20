#Plot panmictic LD distribution across recombination rates
recombination_rate = c("1e-7","1e-8","1e-9")
results            = list()
var_ld             = list()
gm_snp             = list()
target_ld          = list()
mean_ld_peripheral = list()
peripheral_pop     = 1000
#Grab the results
for (i in 1:length(recombination_rate)) {
  results_file = list.files(path = "results", pattern = sprintf("_0_recomb_%s",recombination_rate[i]), full.names = T)
  results[[i]] = lapply(results_file,readRDS)
  
  #Extract the LD values
  var_ld[[i]]             = lapply(results[[i]], `[[`, 4)
  gm_snp[[i]]             = lapply(results[[i]], `[[`, 1)
  target_ld[[i]]          = lapply(results[[i]], `[[`, 2)
  mean_ld_peripheral[[i]] = lapply(results[[i]], `[[`, 3)
  
  #Remove peripheral populations which are non-segregating
  remove_peripheral_pop = which(gm_snp[[i]] == -4)

  if (length(remove_peripheral_pop) > 0 ) {
    var_ld[[i]]    = var_ld[[i]][-remove_peripheral_pop]
    gm_snp[[i]]    = gm_snp[[i]][-remove_peripheral_pop]
    target_ld[[i]] = target_ld[[i]][-remove_peripheral_pop]
  }
}

#Plot the results
fig_1 <- plot_ly(data = iris, x = ~target_ld[[1]], y = ~var_ld[[1]],color = ~unlist(gm_snp[[1]])) %>% layout(showlegend = FALSE) %>% hide_colorbar()
fig_2 <- plot_ly(data = iris, x = ~target_ld[[2]], y = ~var_ld[[2]],color = ~unlist(gm_snp[[2]])) %>% hide_colorbar()
fig_3 <- plot_ly(data = iris, x = ~target_ld[[3]], y = ~var_ld[[3]],color = ~unlist(gm_snp[[3]])) %>% layout(showlegend = FALSE)
fig = subplot(fig_1,fig_2,fig_3,shareY = T) %>% layout(title = "Distribution of Conditional Linkage Disequibrium Across Recombination Rates",
                                                       yaxis = list(title = "Variance LD Peripheral"))
fig %>% layout(annotations = list(
  list(x = 0.2 , y = 1.0, text = "rho = 1e-7", showarrow = F, xref='paper', yref='paper'),
  list(x = 0.6 , y = 1.0, text = "rho = 1e-8", showarrow = F, xref='paper', yref='paper'),
  list(x = 1.0 , y = 1.0, text = "rho = 1e-9", showarrow = F, xref='paper', yref='paper')))

#####Plot the conditional distributions (each point in the graph above)
values = c(25,20,50)
rep_1 = list.files(path = "msprime", pattern = sprintf("peripheral_\\d+_replicate_%s_dt_0_recomb_1e-8.csv",values[[1]]), full.names = T)
rep_2 = list.files(path = "msprime", pattern = sprintf("peripheral_\\d+_replicate_%s_dt_0_recomb_1e-8.csv",values[[2]]), full.names = T)
rep_3 = list.files(path = "msprime", pattern = sprintf("peripheral_\\d+_replicate_%s_dt_0_recomb_1e-8.csv",values[[3]]), full.names = T)

ld                       = lapply(rep_1,fread,header = T)
peripheral_array          = array(unlist(ld),dim = c(1000,2,peripheral_pop)) 
ld_peripheral             = apply(peripheral_array, c(3), cov.wt,method = "ML", cor = TRUE)
ld_cor_peripheral_matrix  = lapply(ld_peripheral, `[[`, 4)
ld_peripheral_values      = c()
#Loop through list and extract the values of ld we need
for (j in 1:length(ld_cor_peripheral_matrix)) {
  ld_peripheral_values[j] = ld_cor_peripheral_matrix[[j]][1,2]
}

density1 = density(ld_peripheral_values)

ld                        = lapply(rep_2,fread,header = T)
peripheral_array          = array(unlist(ld),dim = c(nhaps_target,2,peripheral_pop)) 
ld_peripheral             = apply(peripheral_array, c(3), cov.wt,method = "ML", cor = TRUE)
ld_cor_peripheral_matrix  = lapply(ld_peripheral, `[[`, 4)
ld_peripheral_values      = c()
#Loop through list and extract the values of ld we need
for (j in 1:length(ld_cor_peripheral_matrix)) {
  ld_peripheral_values[j] = ld_cor_peripheral_matrix[[j]][1,2]
}
density2 = density(ld_peripheral_values)


ld                        = lapply(rep_3,fread,header = T)
peripheral_array          = array(unlist(ld),dim = c(nhaps_target,2,peripheral_pop)) 
ld_peripheral             = apply(peripheral_array, c(3), cov.wt,method = "ML", cor = TRUE)
ld_cor_peripheral_matrix  = lapply(ld_peripheral, `[[`, 4)
ld_peripheral_values      = c()
#Loop through list and extract the values of ld we need
for (j in 1:length(ld_cor_peripheral_matrix)) {
  ld_peripheral_values[j] = ld_cor_peripheral_matrix[[j]][1,2]
}

density3 = density(ld_peripheral_values)

fig <- plot_ly(x = ~density1$x, y = ~density1$y, type = 'scatter', mode = 'lines', name = signif(unlist(gm_snp[[1]][values[1]]),3), fill = 'tozeroy')
fig <- fig %>% add_trace(x = ~density2$x, y = ~density2$y, name = signif(unlist(gm_snp[[2]][values[2]]),3), fill = 'tozeroy')
fig <- fig %>% add_trace(x = ~density3$x, y = ~density3$y, name = signif(unlist(gm_snp[[3]][values[1]]),3), fill = 'tozeroy')
fig  = fig %>% layout(title = "Distribution of LD Metrics in Peripheral Populations for Select Haplotype Configurations",
               xaxis = list(title = "LD r"), yaxis = list(title = "Density"))
fig
