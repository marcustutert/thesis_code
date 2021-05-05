#Generating weights from prior
#Split into jobs in order to improve speed of the inference
library(data.table)
library(abind)
source("/well/mcvean/mtutert/thesis/coalescent_coverage/helper_functions.R")
file                 = snakemake@params$replicates
chunk                = snakemake@params$chunk
population_pairs     = snakemake@params$paired_values
print(population_pairs)

#Get populations (GWAS & Reference)
split_pops = strsplit(population_pairs, "_")
gwas_pop   = split_pops[[1]][1]
ref_pop    = split_pops[[1]][3]
fst        = split_pops[[1]][5]

#Import the reference & gwas panels in according to the population_pairs string
ref      = as.matrix(fread(sprintf("msprime_data/OOA/%s_replicate_%s.csv", ref_pop, file), header = T))
gwas     = as.matrix(fread(sprintf("msprime_data/OOA/%s_replicate_%s.csv", gwas_pop, file), header = T))
print(dim(gwas))
#Perform filtering (removing non-segregating and low freq SNPs)
res      = msprime_gwas_sumstats(gwas_haplotypes = gwas, reference_haplotypes = ref)
gwas     = res[[1]]
ref      = res[[2]]

#Write out GWAS & Ref (matched) tables
#Note that this DOESN'T have to be done on a per Fst basis (since wont change data structure)
write.table(gwas,sprintf("msprime_data/OOA/matched_panels/%s_GWAS_matched_to_%s_replicate_%s.csv",gwas_pop,ref_pop,file),quote = F,col.names = T, row.names = F)
write.table(ref,sprintf("msprime_data/OOA/matched_panels/%s_Ref_matched_to_%s_replicate_%s.csv",ref_pop,gwas_pop,file),quote = F,col.names = T, row.names = F)
#Back out the correct Fst given which population we are looking at
nhaps_ref     = nrow(ref)
nhaps_gwas    = nrow(gwas)
nsnps         = ncol(ref)
effective_fst = as.numeric(fst)


#Loop across to get the draws I need
qc_matrix_count = 0
nSamples        = 500
tol             = 1
while(qc_matrix_count<nSamples) {
  print(qc_matrix_count)
  #Draw from gamma_quantiled_weights nhaps times
  gamma_draw               = rgamma(n = nhaps_ref, shape =  1/( nhaps_ref * ( effective_fst / (1-effective_fst))), scale = ( nhaps_ref * (effective_fst/(1-effective_fst))))
  #Extend into matrix
  weight_matrix            = matrix(rep(gamma_draw,nsnps), ncol = nsnps)
  #Normalize matrix
  norm_weight_matrix       = weight_matrix/colSums(weight_matrix)[col(weight_matrix)]

  ####### Remove Ascertainment Bias
  #Ask if the weight matrix we are generating will break our filtering step
  tol = 1e-100
  if (all(colSums(ref*norm_weight_matrix)) > tol & all(colSums(ref*norm_weight_matrix)) < 1-tol) {
    qc_matrix_count = qc_matrix_count + 1
    #cbind the matrix results for AF if we have not generated a good weight matrix yet
    if (qc_matrix_count == 1) {
      print(qc_matrix_count)
      AF_Inferred_Results = colSums(ref*norm_weight_matrix)
      cov                 = cov.wt(ref,norm_weight_matrix[,1],cor = TRUE, method = "ML")
      LD_Results          = as.array(cov$cor)
    }
    else{
    AF_Inferred_Results      = cbind(AF_Inferred_Results,colSums(ref*norm_weight_matrix))
    cov                      = cov.wt(ref,norm_weight_matrix[,1],cor = TRUE, method = "ML")
    LD_Results               = abind(LD_Results,cov$cor, along = 3)
    }
  }
}

#Save object in format Fst_#_Replicate_#_Chunk_#
saveRDS(object = LD_Results, file = sprintf("results/OOA/%s_GWAS_%s_Ref_%s_LD_Replicate_%s_Chunk_%s.RData",gwas_pop,ref_pop,fst,file,chunk), version = 2)
saveRDS(object = AF_Inferred_Results, file = sprintf("results/OOA/%s_GWAS_%s_Ref_%s_AF_Replicate_%s_Chunk_%s.RData",gwas_pop,ref_pop,fst,file,chunk), version = 2)
#



