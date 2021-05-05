#Funcion that reads in a *haplotype* panel (generated from msprime) and reads out GWAS sumstats (SEs) form a GLM
msprime_gwas_sumstats = function(gwas_haplotypes,reference_haplotypes){
  #Add rsids to make the filtering way easy
  rsids            = paste("rsid",seq(1:ncol(gwas_haplotypes)),sep = "")
  colnames(gwas_haplotypes)      = rsids
  colnames(reference_haplotypes) = rsids

  #Split GWAS population into cases & controls
  #First convert GWAS haplotypes -> genotypes
  #gwas_genotypes          = t(sapply(seq(1,nrow(gwas_haplotypes),by=2),function(i) colSums(gwas_haplotypes[i:(i+1),])))
  #Remove duplicate haplotypes!
  #gwas_genotypes          = gwas_genotypes[!duplicated(gwas_genotypes), ]

  #status_index            = rbinom(seq(1:nrow(gwas_genotypes)),1,0.5)

  #Remove everything in the GWAS panel that is outside of AF range
  tol = 1e-4
  qc_variants_ref_min              = which(colMeans(reference_haplotypes) == 0)
  qc_variants_ref_max              = which(colMeans(reference_haplotypes) == 1)
  qc_variants_gwas_max             = which(colMeans(gwas_haplotypes)      == 0)
  qc_variants_gwas_min             = which(colMeans(gwas_haplotypes)      == 1)


  filtered_snps            = unique(c(qc_variants_ref_min,
                                      qc_variants_ref_max,
                                      qc_variants_gwas_min,
                                      qc_variants_gwas_max))

  if (length(filtered_snps) > 0) {
    # #Remove these variants
    gwas_haplotypes           = gwas_haplotypes[,-filtered_snps]
    reference_haplotypes      = reference_haplotypes[,-filtered_snps]
  }
  # nsnps = ncol(gwas_genotypes)
  # op<-cbind(1:nsnps, rep(0, nsnps), rep(0, nsnps));
  # colnames(op)<-c("Pos", "Beta", "se.Beta");
  # for (i in 1:nsnps) {
  #   flush.console()
  #   slm<-summary(glm(status_index~gwas_genotypes[,i], family="binomial"));
  #   op[i,2:3]<-slm$coef[2,1:2];
  # }
  return(list(gwas_haplotypes,reference_haplotypes))
}



