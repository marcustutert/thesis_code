#Take as input the msprime data, match the SNPs within ref & gwas panels (based on some allele frequency threshold)
#Generate GWAS summary statistics either with 1)AF->SE transformation directly and add noise or two 2)with glm and some genotype missigness 4)AF-SE w missingness

#Snakemake variables
replicates                     = snakemake@params$replicates
divergence_time                = snakemake@params$divergence_time
noise                          = snakemake@params$noise

# # ###Debugging values###
# replicates      = "5"
# divergence_time = "10"
# noise           = "0.005"

library(data.table)
source("likelihood.R")
#Read in the msprime table
gwas = as.matrix(fread(sprintf("msprime_data/gwas_%s_dt_%s.csv",replicates, divergence_time), header = T))
ref  = as.matrix(fread(sprintf("msprime_data/ref_%s_dt_%s.csv",replicates,divergence_time), header = T))

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

#Run this two ways: once with generative gamma noise
#1) Gamma Noise (Generative Model), Direct SEs
gwas_se_directly = 1/(pi_pop*(1-pi_pop))
noise            = as.numeric(noise)
#Add in lots of gamma noise
noisy_se         = rgamma(length(pi_pop), shape = 1 / noise, scale = noise) * gwas_se_directly
#Run the optimization
pars = optim(par = c( 0.2 ,0.001), fn = gamma_likelihood, se = noisy_se, pi_ref = pi_ref, lower = c( 0.0001,0.0001), upper = c( 0.5,0.2),method = "L-BFGS-B")$par
# #Save Rdata
saveRDS(object = pars[1],file = sprintf("mle_results/gamma_noise/fst_replicate_%s_dt_%s_noise_%s",replicates,divergence_time,noise))
saveRDS(object = pars[2],file = sprintf("mle_results/gamma_noise/noise_replicate_%s_dt_%s_noise_%s",replicates,divergence_time,noise))

#2) Get the standard errors from a GLM
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
saveRDS(object = pars[1],file = sprintf("mle_results/glm_se/fst_replicate_%s_dt_%s_noise_%s",replicates,divergence_time,noise))
saveRDS(object = pars[2],file = sprintf("mle_results/glm_se/noise_replicate_%s_dt_%s_noise_%s",replicates,divergence_time,noise))


