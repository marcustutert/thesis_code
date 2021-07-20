#Code to generate split populations (ancestrally) in msprime
#If running this on Rstudio on rescomp server
#library(reticulate)
#use_condaenv(condaenv = "python_3.7", conda = "/well/mcvean/mtutert/anaconda3/bin/conda")

#Import modules
import msprime #NOTE THIS IS THE DEVELOPPER VERSION (1.00.6A)
import allel; print('scikit-allel', allel.__version__)
import numpy as np
import pandas as pd
import math
from itertools import product
from random import sample 

#Snakemake params
divergence_time    = int(snakemake.params[0])
recombination_rate = snakemake.params[1]

####Variables#### 
sequence_length    = 1000000
Ne_A               = 10000
Ne_B               = Ne_A/2 
peripheral_pops    = 500
###################

for i in range(1,2):
  demography = msprime.Demography()
  descendant_pops = []
  demography.add_population(name="ancestral_population", initial_size=Ne_A)
  for j in range(0,peripheral_pops):
    name = "peripheral_%s" % j
    demography.add_population(name=name, initial_size=Ne_B)
    demography.add_mass_migration(time=divergence_time, source=name, dest="ancestral_population", proportion=1)
    descendant_pops.append("peripheral_%s" % j)
    
  #Create the dictionary of sample sizes python
  sample_sizes = { k : 500 for k in descendant_pops }
  
  #Do a recombination map with a single recombination hotspot in the middle
  #rate_map = msprime.RateMap.read_hapmap(position=[0, sequence_length / 2 - 1, sequence_length / 2, sequence_length],rate=[0,recombination_rate, 0])
  rate_map = msprime.RateMap.read_hapmap("recombination_map")
  
  tree       = msprime.sim_ancestry(recombination_rate = rate_map, 
                                    demography         = demography,
                                    samples            = sample_sizes,
                                    random_seed        = i)
                                    
  #Sprinkle on some mutations onto the tree (Salt Bae style)                       
  mutated_ts = msprime.sim_mutations(tree, rate=1e-9, model = "binary",random_seed = i)
  
  #Choose which SNP to select on either side of the recombination breakpoint
  #This breakpoint occurs at math.floor(sequence_length/2)
  #Get list of the SNP positions 
  tables                      = mutated_ts.tables
  breakpoint_1                = 500500
  breakpoint_2                = 499500
  positions                   = tables.sites.position
  
  #Remove any SNPs that are NOT segregating in the target population
  #Get list of all samples, filter to those in pop_1 (target)
  my_arr                                   = np.zeros(mutated_ts.num_nodes, dtype=bool)
  my_arr[mutated_ts.samples(population=1)] = True
  sample_in_pop1                           = my_arr[mutated_ts.samples()]
  
  #Loop through the variants to find all frequencies
  freq_results = []
  for v in mutated_ts.variants():
    genos = v.genotypes[sample_in_pop1]
    freq=np.sum(genos>0)/len(genos)
    freq_results.append(freq)
  
  #Find index of variants with non-segregating allele frequencies
  non_seg_variants   = np.where(np.asarray(freq_results) == 0)[0]
  #Remove this from position vector
  positions_filtered = np.delete(positions, non_seg_variants)
  
  before_breakpoint_filtered_positions = [x for x in list(positions_filtered) if x > breakpoint_1]
  after_breakpoint_filtered_positions  = [x for x in list(positions_filtered) if x < breakpoint_2]
  
  #Randomly choose a SNP in each region (before & after)
  snp_1_index = np.random.choice(before_breakpoint_filtered_positions,1)
  snp_2_index = np.random.choice(after_breakpoint_filtered_positions,1)
  
  #Look back in original positions array to see which this is
  snp_1_position = list(np.where(positions == snp_1_index)[0])[0]
  snp_2_position = list(np.where(positions == snp_2_index)[0])[0]
  
  #Now convert the mutated_ts to an haplotype array
  haps        = np.transpose(np.asarray(allel.HaplotypeArray(mutated_ts.genotype_matrix())))
  target      = pd.DataFrame(haps[0:1000,[snp_1_position,snp_2_position]])
  target.to_csv("msprime/target_1_replicate_%s_dt_%s_recomb_%s.csv" % (i,divergence_time,recombination_rate) , index = False)
  #Write out the files
  for j in range(1,peripheral_pops):
    #Output the GWAS population
    gwas = pd.DataFrame(haps[j*1000:(j+1)*1000,[snp_1_position,snp_2_position]])
    gwas.to_csv("msprime/peripheral_%s_replicate_%s_dt_%s_recomb_%s.csv" % (j,i,divergence_time,recombination_rate), index = False)
