
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




