#Convert msprime data into hapgen2 format
#First thing we should do is add rsIDs to all the SNPs, this will just be from 1:ncol(haplotypes)
library(data.table)
#This will need to be run across all divergence times and replicates!
regions      = snakemake@params$regions

#Read in the GWAS population (ADMIX) and the reference populations (A & B)
pop_ADMIX = as.matrix(fread(sprintf("msprime_data/pop_ADMIXED_region_%s.csv",regions),header = T))
pop_A     = as.matrix(fread(sprintf("msprime_data/ref_pop_A_region_%s.csv",regions),header = T))
pop_B     = as.matrix(fread(sprintf("msprime_data/ref_pop_B_region_%s.csv",regions),header = T))

#Add colnames to both panels
colnames(pop_ADMIX) = paste("rsid", seq(1:ncol(pop_ADMIX)), sep = "")
colnames(pop_A)     = paste("rsid", seq(1:ncol(pop_A)), sep = "")
colnames(pop_B)     = paste("rsid", seq(1:ncol(pop_B)), sep = "")


#Jointly filter the data (A, B and AB)
#Calculate the AFs across the three populations
pi_ref_ADMIX = colMeans(pop_ADMIX)
pi_ref_A = colMeans(pop_A)
pi_ref_B = colMeans(pop_B)

#Get the list of the non-segregating variants in each population
non_seg_index_ADMIX = which(pi_ref_ADMIX < 0.01 | pi_ref_ADMIX > 0.99 ) 
non_seg_index_A     = which(pi_ref_A < 0.01 | pi_ref_A > 0.99 ) 
non_seg_index_B     = which(pi_ref_B < 0.01 | pi_ref_B > 0.99 ) 
  
#Get unique list of othe SNPs to remove
union_seg = unique(c(non_seg_index_A,non_seg_index_B,non_seg_index_ADMIX))

#Check if we have zero position BP
positions = read.table(sprintf("msprime_data/positions_pop_ADMIXED_region_%s",regions))
if (positions[1] == 0) {
  #Add this to the removal of the allele freqs (first index)
  union_seg = c(1,union_seg)
}

#Remove these SNP indexes from all panels
pop_ADMIX = pop_ADMIX[,-union_seg]
pop_A     = pop_A[,-union_seg]
pop_B     = pop_B[,-union_seg]

#We also need to remove the non-segregating sites from the positions table
positions = positions[-union_seg]
#Write out the filtered position indexes
write.table(positions,sprintf("msprime_data/joint_filtered/positions_region_%s",regions), quote = F, row.names = F, col.names = F)

#Write out the filtered populations 
write.table(pop_ADMIX,sprintf("msprime_data/joint_filtered/pop_ADMIX_region_%s", regions), quote = F, row.names = F, col.names = T)
write.table(pop_A, sprintf("msprime_data/joint_filtered/ref_pop_A_region_%s", regions), quote = F, row.names = F, col.names = T)
write.table(pop_B, sprintf("msprime_data/joint_filtered/ref_pop_B_region_%s", regions), quote = F, row.names = F, col.names = T)

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
#For now we will assume a flat recombination rate across the genome
first_position        = snp_bp[1]
last_position         = snp_bp[length(snp_bp)]
human_level_recomb    = 1 #This can be changed perhaps?
hapmap_recomb_row_one = c(first_position,1,0)
hapmap_recomb_row_two = c(last_position,0,last_position*1e-10)                   #Change this depending on the recombination rate!
recomb_file           = rbind(hapmap_recomb_row_one,hapmap_recomb_row_two)
colnames(recomb_file) = c("position","COMBINED_rate(cM/Mb)","Genetic_Map(cM)")
write.table(recomb_file,sprintf("hapgen2/pop_ADMIX_region_%s.map",regions), quote = F, row.names = F, col.names = T)



