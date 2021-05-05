#Function that will actually perform the sumstat imputation
sumstat_impute = function(typed_snps,
                          untyped_snps_index,
                          genotyped_sumstats,
                          imputed_sumstats,
                          LD)
{
  LD_typed_untyped   = LD[typed_snps,untyped_snps_index] #LD of typed - untyped pairs
  inv_LD_typed       = solve( LD[typed_snps,typed_snps] + 0.001*diag(length(typed_snps))) #inverse of LD of typed SNPs
  W                  = LD[untyped_snps_index, typed_snps] %*% inv_LD_typed #these are the weights that turn typed to imputed
  infos              = as.numeric(rowSums(W * LD[untyped_snps_index,typed_snps])) #info measures per each untyped
  z.imp              = (W %*% (genotyped_sumstats$frequentist_add_beta_1/genotyped_sumstats$frequentist_add_se_1^0.5))/sqrt(infos) #use scaling 1/sqrt(infos) to get var=1 for z-scores
  #true_z             = -imputed_sumstats[,2]/(imputed_sumstats[,3]) #Flip beta alleles
  return((z.imp))
}


#This function calculates LD between all SNP's
LD_Matrix = function(haplotypes){

  #Takes in haplotype matrix (dim of haps x snps) to calculate various LD metrics (r here)
  haplotypes = t(haplotypes)
  pAB = haplotypes %*% t( haplotypes ) / ncol( haplotypes)
  pA  = rowMeans( haplotypes )
  pB  = rowMeans( haplotypes )
  D  = pAB - pA %*% t( pB )
  denA = sqrt( 1 / pA / ( 1 - pA ) )
  denB = sqrt( 1 / pB / ( 1 - pB ) )
  r  = t( D * denA ) * denB
  return(r)
}

