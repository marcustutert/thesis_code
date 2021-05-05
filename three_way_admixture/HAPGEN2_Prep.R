#Convert msprime data into hapgen2 format
#First thing we should do is add rsIDs to all the SNPs, this will just be from 1:ncol(haplotypes)
library(data.table)
#This will need to be run across all divergence times and replicates!
regions      = snakemake@params$regions

#Read in the GWAS population (ADMIX) and the reference populations (A & B)
pop_ADMIX = as.matrix(fread(sprintf("msprime_data/pop_ADMIXED_region_%s.csv",regions),header = T))
pop_AFR     = as.matrix(fread(sprintf("msprime_data/ref_pop_AFR_region_%s.csv",regions),header = T))
pop_EUR     = as.matrix(fread(sprintf("msprime_data/ref_pop_EUR_region_%s.csv",regions),header = T))
pop_EAS     = as.matrix(fread(sprintf("msprime_data/ref_pop_EAS_region_%s.csv",regions),header = T))

#Add colnames to both panels
colnames(pop_ADMIX)   = paste("rsid", seq(1:ncol(pop_ADMIX)), sep = "")
colnames(pop_AFR)     = paste("rsid", seq(1:ncol(pop_AFR)), sep = "")
colnames(pop_EUR)     = paste("rsid", seq(1:ncol(pop_EUR)), sep = "")
colnames(pop_EAS)     = paste("rsid", seq(1:ncol(pop_EAS)), sep = "")

#Jointly filter the data 
#Calculate the AFs across the three populations
pi_ref_ADMIX = colMeans(pop_ADMIX)
pi_ref_AFR   = colMeans(pop_AFR)
pi_ref_EUR   = colMeans(pop_EUR)
pi_ref_EAS   = colMeans(pop_EAS)

#Get the list of the non-segregating variants in each population
non_seg_index_ADMIX = which(pi_ref_ADMIX < 0.01 | pi_ref_ADMIX > 0.99 ) 
non_seg_index_AFR     = which(pi_ref_AFR < 0.01 | pi_ref_AFR > 0.99 ) 
non_seg_index_EUR     = which(pi_ref_EUR < 0.01 | pi_ref_EUR > 0.99 ) 
non_seg_index_EAS     = which(pi_ref_EAS < 0.01 | pi_ref_EAS > 0.99 ) 


#Get unique list of othe SNPs to remove
union_seg = unique(c(non_seg_index_AFR,non_seg_index_EUR,non_seg_index_EAS,non_seg_index_ADMIX))

#Check if we have zero position BP
positions = read.table(sprintf("msprime_data/positions_pop_ADMIXED_region_%s",regions))
if (positions[1] == 0) {
  #Add this to the removal of the allele freqs (first index)
  union_seg = c(1,union_seg)
}

#Remove these SNP indexes from all panels
pop_ADMIX   = pop_ADMIX[,-union_seg]
pop_AFR     = pop_AFR[,-union_seg]
pop_EUR     = pop_EUR[,-union_seg]
pop_EAS    = pop_EAS[,-union_seg]


#We also need to remove the non-segregating sites from the positions table
positions = positions[-union_seg]
#Write out the filtered position indexes
write.table(positions,sprintf("msprime_data/joint_filtered/positions_region_%s",regions), quote = F, row.names = F, col.names = F)

#Write out the filtered populations 
write.table(pop_ADMIX,sprintf("msprime_data/joint_filtered/pop_ADMIX_region_%s", regions), quote = F, row.names = F, col.names = T)
write.table(pop_AFR, sprintf("msprime_data/joint_filtered/ref_pop_AFR_region_%s", regions), quote = F, row.names = F, col.names = T)
write.table(pop_EUR, sprintf("msprime_data/joint_filtered/ref_pop_EUR_region_%s", regions), quote = F, row.names = F, col.names = T)
write.table(pop_EAS, sprintf("msprime_data/joint_filtered/ref_pop_EAS_region_%s", regions), quote = F, row.names = F, col.names = T)

#####Generate HAPGEN2 .hap file####
gwas = read.table(sprintf("msprime_data/joint_filtered/pop_ADMIX_region_%s",regions), header = T)
write.table(t(gwas),sprintf("hapgen2/pop_ADMIX_region_%s.hap",regions), quote = F, col.names = F, row.names = F)

#####Generate hapgen2 legend file####
#This will have columns: rs, position, X0, X1 (ref and alt base: A->T); the latter which is just mock data

#Read out the rsids from the column names in the GWAS panel
rsids = colnames(gwas)
#Get the positions from the filtered msprime positions table
snp_bp  = as.vector(read.table(sprintf("msprime_data/joint_filtered/positions_region_%s",regions), header = F))
#Generate the fake mutations
X0 = rep("A",length(snp_bp))
X1 = rep("G",length(snp_bp))

#Bind together to form .leg format
leg = t(rbind(rsids,snp_bp,X0,X1))
colnames(leg) = c("rs", "position","X0","X1")
#Write out to hapgen2 directory
write.table(leg,sprintf("hapgen2/pop_ADMIX_region_%s.leg",regions), quote = F, col.names = T, row.names = F)

####Generate hapgen2 recombination file####
#Using the recombination map provided by msprime
first_position        = snp_bp[1]
last_position         = snp_bp[length(snp_bp)]
#Read in recombination rate map:
genetic_map            = as.matrix(as.vector(fread(sprintf("msprime_data/recombination_map"), header = F)))[1,]*100 #Read in the culmative mass (genetic map)
#Filter to only the positions that have mutations on them
positions = unlist(snp_bp)
genetic_map_positions  = genetic_map[unlist(snp_bp)]
#Get the recombination rates
recomb_rates           = as.matrix(as.vector(fread(sprintf("msprime_data/recombination_rate"), header = F)))[1,] #Read in the culmative mass (genetic map)
#Convert these positions into recombination rates using the formula on RELATE page:
#Add in the 1st position
combined_rate_across_positions = c()
for (i in 1:(length(genetic_map_positions)-1)) {
  combined_rate_across_positions[i] = ((genetic_map_positions[i+1] - genetic_map_positions[i])/(positions[i+1] - positions[i])) * 1e6
}
#Add in the 1st position
combined_rate_across_positions = c(combined_rate_across_positions,recomb_rates[1]*1e6*100)
#Calculate the recombination rate 
recomb_file           = cbind(positions,combined_rate_across_positions,genetic_map_positions)
colnames(recomb_file) = c("position","COMBINED_rate(cM/Mb)","Genetic_Map(cM)")
write.table(recomb_file,sprintf("hapgen2/pop_ADMIX_region_%s.map",regions), quote = F, row.names = F, col.names = T)



