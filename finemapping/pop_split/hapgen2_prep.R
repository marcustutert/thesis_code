#Convert msprime data into hapgen2 format
#First thing we should do is add rsIDs to all the SNPs, this will just be from 1:ncol(haplotypes)
library(data.table)
#This will need to be run across all divergence times and replicates!
divergence_time = snakemake@params$divergence_time
replicates      = snakemake@params$replicates

#Read in the GWAS data
gwas = as.matrix(fread(sprintf("pop_split/msprime/gwas_%s_dt_%s.csv",replicates,divergence_time),header = T))
ref  = as.matrix(fread(sprintf("pop_split/msprime/ref_%s_dt_%s.csv",replicates,divergence_time),header = T))

#Add colnames to both panels
colnames(gwas) = paste("rsid", seq(1:ncol(gwas)), sep = "")
colnames(ref) = paste("rsid", seq(1:ncol(ref)), sep = "")

#Search across the gwas and reference panel for any non-segregating variants
non_seg_index_gwas = which(colMeans(gwas) == 0 | colMeans(gwas) == 1 )
non_seg_index_ref  = which(colMeans(ref) == 0 | colMeans(ref) == 1 )

#Get union of these sets
union_seg = c(non_seg_index_gwas,non_seg_index_ref)

#Check if we have zero position BP
positions = read.table(sprintf("pop_split/msprime/positions_gwas_%s_dt_%s",replicates, divergence_time))
if (positions[1] == 0) {
  #Add this to the removal of the allele freqs (first index)
  union_seg = c(1,union_seg)
}

#Remove these indexes from both panels
gwas = gwas[,-union_seg]
ref  = ref[,-union_seg]

#We also need to remove the non-segregating sites from the positions table
positions = positions[-union_seg]
#Write out position indexes
write.table(positions,sprintf("pop_split/msprime/filtered_panels/positions_gwas_%s_dt_%s",replicates, divergence_time), quote = F, row.names = F, col.names = F)

#Write out tables
write.table(gwas,sprintf("pop_split/msprime/filtered_panels/gwas_%s_dt_%s", replicates, divergence_time), quote = F, row.names = F, col.names = T)
write.table(ref, sprintf("pop_split/msprime/filtered_panels/ref_%s_dt_%s", replicates, divergence_time), quote = F, row.names = F, col.names = T)

#####Generate hapgen2 hap file####
#One per divergence time
gwas = read.table(sprintf("pop_split/msprime/filtered_panels/gwas_%s_dt_%s",replicates,divergence_time), header = T)
write.table(t(gwas),sprintf("pop_split/hapgen2/gwas_%s_dt_%s.hap",replicates,divergence_time), quote = F, col.names = F, row.names = F)

#####Generate hapgen2 legend file####
#This will have columns: rs, position, X0, X1 (ref and alt base: A->T)

#Read out the rsids from the column names in the GWAS panel
rsids = colnames(gwas)
#Get the positions from the filtered msprime positions table
snp_bp  = as.vector(read.table(sprintf("pop_split/msprime/filtered_panels/positions_gwas_%s_dt_%s",replicates,divergence_time), header = F))
#Generate the fake mutations
X0 = rep("A",length(snp_bp))
X1 = rep("G",length(snp_bp))

#Bind together to form .leg format
leg = t(rbind(rsids,snp_bp,X0,X1))
colnames(leg) = c("rs", "position","X0","X1")
#Write out to hapgen2 directory
write.table(leg,sprintf("pop_split/hapgen2/gwas_%s_dt_%s.leg",replicates,divergence_time), quote = F, col.names = T, row.names = F)

####Generate hapgen2 recombination file####
#For now we will assume a flat recombination rate across the genome
first_position        = snp_bp[1]
last_position         = snp_bp[length(snp_bp)]
human_level_recomb    = 1 #This can be changed perhaps?
hapmap_recomb_row_one = c(first_position,1,0)
hapmap_recomb_row_two = c(last_position,0,last_position*1e-5)
recomb_file           = rbind(hapmap_recomb_row_one,hapmap_recomb_row_two)
colnames(recomb_file) = c("position","COMBINED_rate(cM/Mb)","Genetic_Map(cM)")
write.table(recomb_file,sprintf("pop_split/hapgen2/gwas_%s_dt_%s.map",replicates,divergence_time), quote = F, row.names = F, col.names = T)



