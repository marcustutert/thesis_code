#Run hapgen2
#Example command:
#/apps/eb/skylake/software/HAPGEN2/2.2.0/hapgen2 -m test.map -l test.leg -h test.hap -o test -dl 3 1 1.5 2.25 11 0 2 4
#Pick two variants to be 'causal'
regions      = snakemake@params$regions

#Read in GWAS table
positions       = as.vector(read.table(sprintf("msprime_data/joint_filtered/positions_region_1", regions)))
ncausal_snps    = 3
snps_causal     = sample(positions,size = ncausal_snps)
ncases          = 5000
ncontrols       = 5000

#Write out the causal SNPs so we can check it with FINEMAP accuracy
write.table(snps_causal,sprintf("hapgen2/snps_gwas_causal_region_%s", regions),quote = F, row.names = F, col.names = F)
causal_snp_1    = snps_causal[1]
causal_snp_2    = snps_causal[2]
causal_snp_3    = snps_causal[3]
hetero_risk_1   = 1.1 #Change this
homo_risk_1     = 1.4 #Change this
hetero_risk_2   = 1.1 #Change this
homo_risk_2     = 1.4 #Change this
hetero_risk_3   = 1.1 #Change this
homo_risk_3     = 1.4 #Change this
system(sprintf("/apps/eb/skylake/software/HAPGEN2/2.2.0/hapgen2 -m hapgen2/pop_ADMIX_region_%s.map -l hapgen2/pop_ADMIX_region_%s.leg -h hapgen2/pop_ADMIX_region_%s.hap -o hapgen2/gwas_region_%s -dl %s 1 %s %s %s 1 %s %s %s 1 %s %s -n %s %s",
                 regions,regions,regions,regions,causal_snp_1,hetero_risk_1,homo_risk_1,causal_snp_2,hetero_risk_2,homo_risk_2, causal_snp_3, hetero_risk_3, homo_risk_3,ncases,ncontrols))


