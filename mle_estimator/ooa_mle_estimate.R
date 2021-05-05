#OOA mle estimate

#Snakemake variables, do this w OOA model pairs
replicates                     = snakemake@params$replicates
population_pairs               = snakemake@params$population_pairs

#Split out the string to pick the correct OOA populations as GWAS & Ref (YRI/CEU/CHB)
split_pops = strsplit(population_pairs, "_")
gwas_pop   = split_pops[[1]][1]
ref_pop    = split_pops[[1]][3]

library(data.table)
source("likelihood.R")
gwas = as.matrix(fread(sprintf("ooa_msprime_data/%s_replicate_%s.csv",gwas_pop, replicates), header = T))
ref  = as.matrix(fread(sprintf("ooa_msprime_data/%s_replicate_%s.csv",ref_pop,replicates), header = T))

#Filter out variants that are at less than 1% MAF in either population
#Add colnames to both panels
colnames(gwas) = paste("rsid", seq(1:ncol(gwas)), sep = "")
colnames(ref) = paste("rsid", seq(1:ncol(ref)), sep = "")

#Search across the gwas and reference panel for any non-segregating variants
non_seg_index_gwas = which(colMeans(gwas) < 0.01 | colMeans(gwas) > .99 )
non_seg_index_ref  = which(colMeans(ref) < 0.01 | colMeans(ref) > .99 )

#Get union of these sets
union_seg = c(non_seg_index_gwas,non_seg_index_ref)

#Remove these indexes from both panels
gwas = gwas[,-union_seg]
ref  = ref[,-union_seg]

#Calculate gwas and ref allele frequencies
pi_pop = colMeans(gwas)
pi_ref = colMeans(ref)

#Get parameters
nsnps = ncol(gwas)
nhaps = nrow(gwas)

#Get the standard errors from a GLM
op<-cbind(1:nsnps, rep(0,nsnps), rep(0,nsnps))
colnames(op)<-c("Pos", "Beta", "se.Beta")
#Randomly assign case/controls at 50% freq.
phi = 0.5
y<-rbinom(nhaps, 1, phi); 
for (i in 1:nsnps) {
  #Run linear model
  slm<-summary(glm(y~gwas[,i], family="binomial"))
  op[i,2:3]<-slm$coef[2,1:2]
}

#Extract standard errors from this
gwas_se_noisy = op[,3] #Don't forget to square!
#Calculate mbar
mbar = phi*(1-phi)*nhaps
#Run the optimization
pars = optim(par = c( 0.2 ,0.001), fn = gamma_likelihood, se = gwas_se_noisy^2, pi_ref = pi_ref, mbar = mbar,lower = c( 0.0001,0.0001), upper = c( 0.5,0.2),method = "L-BFGS-B")$par
# #Save Rdata
saveRDS(object = pars[1],file = sprintf("ooa_mle_results/fst_%s_GWAS_%s_Ref_replicate_%s",gwas_pop,ref_pop,replicates))
saveRDS(object = pars[2],file = sprintf("ooa_mle_results/noise_%s_GWAS_%s_Ref_replicate_%s",gwas_pop,ref_pop,replicates))


