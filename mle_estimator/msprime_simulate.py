#Uncomment if running locally
#library(reticulate)
#use_condaenv(condaenv = "python_3.7", conda = "/well/mcvean/mtutert/anaconda3/bin/conda")

import msprime
import allel
print('scikit-allel', allel.__version__)
import numpy as np
import pandas as pd

#This will be done under models of population splits
divergence_time   = snakemake.params[0]#This might be wrong...?
divergence_time   = int(divergence_time)
#divergence_time = 10

####Variables#### 
#These are the DIPLOID counts
ngwas_haps         = 1000
nref_haps          = 1000
ngwas_replicates   = 1000       #This needs to be the same number as whatever is in the top of the SNAKEMAKE file
###################

### Population Split ###
#Generate msprime populations, under model of an ancestral population split some time t ago in the past (generations)
#Use the 'legacy version of the demographic model'
#Now we want to define our population parameters
Ne_A = 10000
Ne_B = Ne_A/2 #Define the population split size (in terms of Ne's). This is the REFERENCE PANEL POP
Ne_C = Ne_A/2 #Define the population split size (in terms of Ne's). This is the GWAS PANEL POP.

#Define our population configurations in msprime
population_configurations = [
      msprime.PopulationConfiguration(
          #Source population (id = 0), goes extinct at time t
          initial_size = Ne_A),
      msprime.PopulationConfiguration(
          initial_size = Ne_B), #Reference panel size
      msprime.PopulationConfiguration(
          initial_size = Ne_C) #GWAS panel size
]
  
#Define the demographic events
demographic_events = [
  # Merging of pops B into A at some time t
  msprime.MassMigration(
      time         = divergence_time,
      source       = 2,
      destination  = 0,
      proportion   = 1.0),
  # Merging of pops C into A at some time t
  msprime.MassMigration(
      time         = divergence_time,
      source       = 1,
      destination  = 0,
      proportion   = 1.0),
  #No need to kill of Pop A, since we sample B & C; can keep it for future use
  # msprime.PopulationParametersChange(
  #     time=merge_time, initial_size=0, population_id=0)
]

#Run msprime simulations
demography = msprime.Demography.from_old_style(population_configurations=population_configurations,demographic_events=demographic_events)

sample_sizes = {
  #Note the conversion to diploid sample sizes
  "pop_1": ngwas_haps,
  "pop_2": nref_haps,
}


#### Loop Across the GWAS replicates (regions) ####

for i in range(1,ngwas_replicates):

  tree       = msprime.sim_ancestry(recombination_rate = 1e-5, 
                                    sequence_length    = 2000,
                                    demography         = demography,
                                    samples            = sample_sizes)
  #Rate here refers to MUTATION rate                                
  mutated_ts = msprime.sim_mutations(tree, rate=1e-5, model = "binary") 
  #Now convert the mutated_ts to an haplotype array
  haps        = np.transpose(np.asarray(allel.HaplotypeArray(mutated_ts.genotype_matrix())))
  print(np.shape(haps))
  positions   = mutated_ts.tables.sites.position
  #Write out the array of the sites to directory
  np.savetxt("msprime_data/positions_gwas_%s_dt_%s" % (i,divergence_time), positions, newline=" ")
  ref = pd.DataFrame(haps[ngwas_haps*2:,])
  print(np.shape(ref))
  #Write out the haplotypes file for both the GWAS population and for the Reference population
  ref.to_csv("msprime_data/ref_%s_dt_%s.csv" % (i,divergence_time), index = False)
  #Output the GWAS population
  gwas = pd.DataFrame(haps[:ngwas_haps*2,])
  print(np.shape(gwas))
  gwas.to_csv("msprime_data/gwas_%s_dt_%s.csv" % (i,divergence_time), index = False)
    

