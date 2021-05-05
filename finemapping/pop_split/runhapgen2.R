#Run hapgen2
#Example command:
#/apps/eb/skylake/software/HAPGEN2/2.2.0/hapgen2 -m test.map -l test.leg -h test.hap -o test -dl 3 1 1.5 2.25 11 0 2 4
#Pick two variants to be 'causal'
divergence_time = snakemake@params$divergence_time
replicates      = snakemake@params$replicates
signal          = snakemake@params$signal

#Read in GWAS table
positions       = as.vector(read.table(sprintf("pop_split/msprime/filtered_panels/positions_gwas_%s_dt_%s",replicates, divergence_time)))
ncausal_snps    = 3
snps_causal     = sample(positions,size = ncausal_snps)
ncases          = 5000
ncontrols       = 5000

if (signal == "TRUE") {
  #Write out the causal SNPs so we can check it with FINEMAP accuracy
  write.table(snps_causal,sprintf("pop_split/hapgen2/snps_gwas_%s_causal_dt_%s", replicates, divergence_time),quote = F, row.names = F, col.names = F)
  causal_snp_1    = snps_causal[1]
  causal_snp_2    = snps_causal[2]
  causal_snp_3    = snps_causal[3]
  hetero_risk_1   = 1.5 #Change this
  homo_risk_1     = 2 #Change this
  hetero_risk_2   = 1.5 #Change this
  homo_risk_2     = 2 #Change this
  hetero_risk_3   = 1.5 #Change this
  homo_risk_3     = 2 #Change this
  system(sprintf("/apps/eb/skylake/software/HAPGEN2/2.2.0/hapgen2 -m pop_split/hapgen2/gwas_%s_dt_%s.map -l pop_split/hapgen2/gwas_%s_dt_%s.leg -h pop_split/hapgen2/gwas_%s_dt_%s.hap -o pop_split/hapgen2/gwas_%s_%s_dt_%s -dl %s 1 %s %s %s 1 %s %s %s 1 %s %s -n %s %s",
                 replicates,divergence_time,replicates, divergence_time,replicates,divergence_time,signal,replicates,divergence_time,causal_snp_1,hetero_risk_1,homo_risk_1,causal_snp_2,hetero_risk_2,homo_risk_2, causal_snp_3, hetero_risk_3, homo_risk_3,ncases,ncontrols))
}

if (signal == "FALSE") {
  #Run HAPGEN2 under the null
  system(sprintf("/apps/eb/skylake/software/HAPGEN2/2.2.0/hapgen2 -m pop_split/hapgen2/gwas_%s_dt_%s.map -l pop_split/hapgen2/gwas_%s_dt_%s.leg -h pop_split/hapgen2/gwas_%s_dt_%s.hap -o pop_split/hapgen2/gwas_%s_%s_dt_%s -dl %s 1 1 1 -n %s %s",
                 replicates,divergence_time,replicates, divergence_time,replicates,divergence_time,signal,replicates,divergence_time,snps_causal[1],ncases,ncontrols))
}


