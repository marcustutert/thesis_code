#Plot genetically drifted LD distribution across recombination rates and across genetic drifts
recombination_rate = c("1e-3","1e-4","1e-5")
divergence_time    = c("10","20","50")
#We will store our results as a 9 length list
#First three are all divergence time results for a single recombination rate
var_ld             = list()
gm_snp             = list()
target_ld          = list()
mean_ld_peripheral = list()
peripheral_pop     = 1000

#Start a counter
counter = 1
#Grab the results
results            = list()
for (i in 1:length(recombination_rate)) {
  for (j in 1:length(divergence_time)) {
    #Grab the results files
    results_file = list.files(path = "pop_split_results", pattern = sprintf("replicate_\\d+_dt_%s_recomb_%s",divergence_time[j],recombination_rate[i]), full.names = T)
    results[[counter]] = lapply(results_file,readRDS)
    #Extract the LD values
    var_ld[[counter]]             = lapply(results[[counter]], `[[`, 4)
    gm_snp[[counter]]             = lapply(results[[counter]], `[[`, 1)
    target_ld[[counter]]          = lapply(results[[counter]], `[[`, 2)
    mean_ld_peripheral[[counter]] = lapply(results[[counter]], `[[`, 3)
    
    #Remove peripheral populations which are non-segregating
    remove_peripheral_pop = which(gm_snp[[counter]] == -4)
    
    if (length(remove_peripheral_pop) > 0 ) {
      var_ld[[counter]]    = var_ld[[counter]][-remove_peripheral_pop]
      gm_snp[[counter]]    = gm_snp[[counter]][-remove_peripheral_pop]
      target_ld[[counter]] = target_ld[[counter]][-remove_peripheral_pop]
    }
    counter = counter + 1
  }
}

#######################################

#Plot the results (9 graphs in total, for a 3x3 grid)
fig_1 <- plot_ly(data = iris, x = ~target_ld[[1]], y = ~var_ld[[1]],color = ~unlist(gm_snp[[1]])) %>% layout(showlegend = FALSE)
fig_2 <- plot_ly(data = iris, x = ~target_ld[[2]], y = ~var_ld[[2]],color = ~unlist(gm_snp[[2]])) %>% hide_colorbar() %>% hide_colorbar()
fig_3 <- plot_ly(data = iris, x = ~target_ld[[3]], y = ~var_ld[[3]],color = ~unlist(gm_snp[[3]])) %>% layout(showlegend = FALSE) %>% hide_colorbar()
fig_4 <- plot_ly(data = iris, x = ~target_ld[[4]], y = ~var_ld[[4]],color = ~unlist(gm_snp[[4]])) %>% layout(showlegend = FALSE)  %>% hide_colorbar()
fig_5 <- plot_ly(data = iris, x = ~target_ld[[5]], y = ~var_ld[[5]],color = ~unlist(gm_snp[[5]])) %>% hide_colorbar() %>% hide_colorbar()
fig_6 <- plot_ly(data = iris, x = ~target_ld[[6]], y = ~var_ld[[6]],color = ~unlist(gm_snp[[6]])) %>% layout(showlegend = FALSE) %>% hide_colorbar()
fig_7 <- plot_ly(data = iris, x = ~target_ld[[7]], y = ~var_ld[[7]],color = ~unlist(gm_snp[[7]])) %>% layout(showlegend = FALSE) %>% hide_colorbar()
fig_8 <- plot_ly(data = iris, x = ~target_ld[[8]], y = ~var_ld[[8]],color = ~unlist(gm_snp[[8]])) %>% hide_colorbar() 
fig_9 <- plot_ly(data = iris, x = ~target_ld[[9]], y = ~var_ld[[9]],color = ~unlist(gm_snp[[9]])) %>% layout(showlegend = FALSE) %>% hide_colorbar()
fig = subplot(fig_1,fig_2,fig_3,fig_4,fig_5,fig_6,fig_7,fig_8,fig_9,nrows = 3, shareY = T, shareX = T,titleX = F,titleY = F) %>% 
      layout(xaxis = list(title = ""),yaxis = list(title = ""), title = "Distribution of Conditional Linkage Disequibrium Across Recombination Rates and Divergence Times")
#Columns are the divergence times and rows are the recombination rates
#Manually add the xaxis and yaxis annotations to this graphs
m <- list ( l = 50, r = 50, b = 100, t = 100, pad = 4 )
fig %>% layout(margin = m)
fig
#####Plot the conditional distributions (each point in the graph above)
#Do this across three divergence times (so again plotting 3 graphs with 3 traces per graph) and only one single recombination value (1e-8)


values = c(7,8,9)
#First Divergence time
rep_1_dt_10 = list.files(path = "pop_split_msprime", pattern = sprintf("peripheral_\\d+_replicate_%s_dt_10_recomb_1e-5.csv",values[[1]]), full.names = T)
rep_2_dt_10 = list.files(path = "pop_split_msprime", pattern = sprintf("peripheral_\\d+_replicate_%s_dt_10_recomb_1e-5.csv",values[[2]]), full.names = T)
rep_3_dt_10 = list.files(path = "pop_split_msprime", pattern = sprintf("peripheral_\\d+_replicate_%s_dt_10_recomb_1e-5.csv",values[[3]]), full.names = T)

ld                       = lapply(rep_1_dt_10,fread,header = T)
peripheral_array          = array(unlist(ld),dim = c(nhaps_target,2,peripheral_pop)) 
ld_peripheral             = apply(peripheral_array, c(3), cov.wt,method = "ML", cor = TRUE)
ld_cor_peripheral_matrix  = lapply(ld_peripheral, `[[`, 4)
ld_peripheral_values      = c()
#Loop through list and extract the values of ld we need
for (j in 1:length(ld_cor_peripheral_matrix)) {
  ld_peripheral_values[j] = ld_cor_peripheral_matrix[[j]][1,2]
}

density1 = density(ld_peripheral_values, na.rm = T)

ld                        = lapply(rep_2_dt_10,fread,header = T)
peripheral_array          = array(unlist(ld),dim = c(nhaps_target,2,peripheral_pop)) 
ld_peripheral             = apply(peripheral_array, c(3), cov.wt,method = "ML", cor = TRUE)
ld_cor_peripheral_matrix  = lapply(ld_peripheral, `[[`, 4)
ld_peripheral_values      = c()
#Loop through list and extract the values of ld we need
for (j in 1:length(ld_cor_peripheral_matrix)) {
  ld_peripheral_values[j] = ld_cor_peripheral_matrix[[j]][1,2]
}
density2 = density(ld_peripheral_values,na.rm = T)

ld                        = lapply(rep_3_dt_10,fread,header = T)
peripheral_array          = array(unlist(ld),dim = c(nhaps_target,2,peripheral_pop)) 
ld_peripheral             = apply(peripheral_array, c(3), cov.wt,method = "ML", cor = TRUE)
ld_cor_peripheral_matrix  = lapply(ld_peripheral, `[[`, 4)
ld_peripheral_values      = c()
#Loop through list and extract the values of ld we need
for (j in 1:length(ld_cor_peripheral_matrix)) {
  ld_peripheral_values[j] = ld_cor_peripheral_matrix[[j]][1,2]
}

density3 = density(ld_peripheral_values,na.rm = T)

fig_1 <- plot_ly(data = iris, x = ~density1$x, y = ~density1$y, type = 'scatter', mode = 'lines', name = signif(unlist(gm_snp[[7]][values[1]]),3), fill = 'tozeroy')
fig_1 <- fig_1 %>% add_trace(x = ~density2$x, y = ~density2$y, name = signif(unlist(gm_snp[[7]][values[2]]),3), fill = 'tozeroy')
fig_1 <- fig_1 %>% add_trace(x = ~density3$x, y = ~density3$y, name = signif(unlist(gm_snp[[7]][values[3]]),3), fill = 'tozeroy')

#Now repeat for divergence time 20
rep_1_dt_20 = list.files(path = "pop_split_msprime", pattern = sprintf("peripheral_\\d+_replicate_%s_dt_20_recomb_1e-5.csv",values[[1]]), full.names = T)
rep_2_dt_20 = list.files(path = "pop_split_msprime", pattern = sprintf("peripheral_\\d+_replicate_%s_dt_20_recomb_1e-5.csv",values[[2]]), full.names = T)
rep_3_dt_20 = list.files(path = "pop_split_msprime", pattern = sprintf("peripheral_\\d+_replicate_%s_dt_20_recomb_1e-5.csv",values[[3]]), full.names = T)

ld                       = lapply(rep_1_dt_20,fread,header = T)
peripheral_array          = array(unlist(ld),dim = c(nhaps_target,2,peripheral_pop)) 
ld_peripheral             = apply(peripheral_array, c(3), cov.wt,method = "ML", cor = TRUE)
ld_cor_peripheral_matrix  = lapply(ld_peripheral, `[[`, 4)
ld_peripheral_values      = c()
#Loop through list and extract the values of ld we need
for (j in 1:length(ld_cor_peripheral_matrix)) {
  ld_peripheral_values[j] = ld_cor_peripheral_matrix[[j]][1,2]
}

density1_20 = density(ld_peripheral_values,na.rm = T)

ld                        = lapply(rep_2_dt_20,fread,header = T)
peripheral_array          = array(unlist(ld),dim = c(nhaps_target,2,peripheral_pop)) 
ld_peripheral             = apply(peripheral_array, c(3), cov.wt,method = "ML", cor = TRUE)
ld_cor_peripheral_matrix  = lapply(ld_peripheral, `[[`, 4)
ld_peripheral_values      = c()
#Loop through list and extract the values of ld we need
for (j in 1:length(ld_cor_peripheral_matrix)) {
  ld_peripheral_values[j] = ld_cor_peripheral_matrix[[j]][1,2]
}
density2_20 = density(ld_peripheral_values,na.rm = T)


ld                        = lapply(rep_3_dt_20,fread,header = T)
peripheral_array          = array(unlist(ld),dim = c(nhaps_target,2,peripheral_pop)) 
ld_peripheral             = apply(peripheral_array, c(3), cov.wt,method = "ML", cor = TRUE)
ld_cor_peripheral_matrix  = lapply(ld_peripheral, `[[`, 4)
ld_peripheral_values      = c()
#Loop through list and extract the values of ld we need
for (j in 1:length(ld_cor_peripheral_matrix)) {
  ld_peripheral_values[j] = ld_cor_peripheral_matrix[[j]][1,2]
}

density3_20 = density(ld_peripheral_values,na.rm = T)

fig_2 <- plot_ly(data = iris,x = ~density1_20$x, y = ~density1_20$y, type = 'scatter', mode = 'lines', name = signif(unlist(gm_snp[[8]][values[1]]),3), fill = 'tozeroy')
fig_2 <- fig_2 %>% add_trace(x = ~density2_20$x, y = ~density2_20$y, name = signif(unlist(gm_snp[[8]][values[2]]),3), fill = 'tozeroy')
fig_2 <- fig_2 %>% add_trace(x = ~density3_20$x, y = ~density3_20$y, name = signif(unlist(gm_snp[[8]][values[3]]),3), fill = 'tozeroy')

#Third Divergence time
rep_1_dt_50 = list.files(path = "pop_split_msprime", pattern = sprintf("peripheral_\\d+_replicate_%s_dt_50_recomb_1e-5.csv",values[[1]]), full.names = T)
rep_2_dt_50 = list.files(path = "pop_split_msprime", pattern = sprintf("peripheral_\\d+_replicate_%s_dt_50_recomb_1e-5.csv",values[[2]]), full.names = T)
rep_3_dt_50 = list.files(path = "pop_split_msprime", pattern = sprintf("peripheral_\\d+_replicate_%s_dt_50_recomb_1e-5.csv",values[[3]]), full.names = T)

ld                       = lapply(rep_1_dt_50,fread,header = T)
peripheral_array          = array(unlist(ld),dim = c(nhaps_target,2,peripheral_pop)) 
ld_peripheral             = apply(peripheral_array, c(3), cov.wt,method = "ML", cor = TRUE)
ld_cor_peripheral_matrix  = lapply(ld_peripheral, `[[`, 4)
ld_peripheral_values      = c()
#Loop through list and extract the values of ld we need
for (j in 1:length(ld_cor_peripheral_matrix)) {
  ld_peripheral_values[j] = ld_cor_peripheral_matrix[[j]][1,2]
}

density1_50 = density(ld_peripheral_values,na.rm = T)

ld                        = lapply(rep_2_dt_50,fread,header = T)
peripheral_array          = array(unlist(ld),dim = c(nhaps_target,2,peripheral_pop)) 
ld_peripheral             = apply(peripheral_array, c(3), cov.wt,method = "ML", cor = TRUE)
ld_cor_peripheral_matrix  = lapply(ld_peripheral, `[[`, 4)
ld_peripheral_values      = c()
#Loop through list and extract the values of ld we need
for (j in 1:length(ld_cor_peripheral_matrix)) {
  ld_peripheral_values[j] = ld_cor_peripheral_matrix[[j]][1,2]
}
density2_50 = density(ld_peripheral_values,na.rm = T)


ld                        = lapply(rep_3_dt_50,fread,header = T)
peripheral_array          = array(unlist(ld),dim = c(nhaps_target,2,peripheral_pop)) 
ld_peripheral             = apply(peripheral_array, c(3), cov.wt,method = "ML", cor = TRUE)
ld_cor_peripheral_matrix  = lapply(ld_peripheral, `[[`, 4)
ld_peripheral_values      = c()
#Loop through list and extract the values of ld we need
for (j in 1:length(ld_cor_peripheral_matrix)) {
  ld_peripheral_values[j] = ld_cor_peripheral_matrix[[j]][1,2]
}

density3_50 = density(ld_peripheral_values,na.rm = T)
fig_3 <- plot_ly()
fig_3 <- plot_ly(x = ~density1_50$x, y = ~density1_50$y, type = 'scatter', mode = 'lines', name = signif(unlist(gm_snp[[9]][values[1]]),2), fill = 'tozeroy')
fig_3 <- fig_3 %>% add_trace(x = ~density2_50$x, y = ~density2_50$y, name = signif(unlist(gm_snp[[9]][values[2]]),2), fill = 'tozeroy')
fig_3 <- fig_3 %>% add_trace(x = ~density3_50$x, y = ~density3_50$y, name = signif(unlist(gm_snp[[9]][values[3]]),2), fill = 'tozeroy')


#Get the subplot together
fig = subplot(fig_1,fig_2,fig_3, shareY = T)
fig %>% layout(title = "Linkage Disequibrium Haplotype Configuration Distributions Across Divergence Times",
               yaxis = list(title = "Density"))
