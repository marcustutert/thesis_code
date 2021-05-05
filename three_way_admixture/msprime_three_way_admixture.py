#Running 3 way admixture (American Admixture; Browning 2018)
#If running this on Rstudio on rescomp server
library(reticulate)
use_condaenv(condaenv = "python_3.7", conda = "/well/mcvean/mtutert/anaconda3/bin/conda")

import msprime
import allel
print('scikit-allel', allel.__version__)
import numpy as np
import pandas as pd
import stdpopsim

#Build demographic history!
T_OOA = 920

demography = msprime.Demography()
demography.add_population(name="AFR", description="African", initial_size=14474)
demography.add_population(
    name="EUR",
    description="European",
    initial_size=34039,
    growth_rate=0.0038,
)
demography.add_population(
    name="EAS",
    description="East Asian",
    initial_size=45852,
    growth_rate=0.0048,
)
demography.add_population(
    name="ADMIX",
    description="Admixed America",
    initial_size=54664,
    growth_rate=0.05,
)
demography.add_population(
    name="OOA",
    description="Bottleneck out-of-Africa",
    initial_size=1861,
)
demography.add_population(
    name="AMH", description="Anatomically modern humans", initial_size=14474
)
demography.add_population(
    name="ANC",
    description="Ancestral equilibrium",
    initial_size=7310,
)
demography.set_symmetric_migration_rate(["AFR", "EUR"], 2.5e-5)
demography.set_symmetric_migration_rate(["AFR", "EAS"], 0.78e-5)
demography.set_symmetric_migration_rate(["EUR", "EAS"], 3.11e-5)

demography.add_admixture(
    12,
    derived="ADMIX",
    ancestral=["AFR", "EUR"],
    proportions=[0.25,0.75],
)

demography.add_population_split(T_OOA, derived=["EUR", "EAS"], ancestral="OOA")
demography.add_symmetric_migration_rate_change(
    time=T_OOA, populations=["AFR", "OOA"], rate=15e-5
)
demography.add_population_split(2040, derived=["OOA", "AFR"], ancestral="AMH")
demography.add_population_split(5920, derived=["AMH"], ancestral="ANC")


####Variables#### 
#These are the DIPLOID sample sizes
ngwas_haps         = 1000     #GWAS population here is the ADMIX population
nref_haps          = 1000     #Reference popualtion is either AFR, EUR, EAS
nreplicates        = 1000      #This needs to be the same number as whatever is in the top of the SNAKEMAKE file
##################

sample_sizes = {
  #Note the conversion to diploid sample sizes
  "ADMIX": ngwas_haps,
  "AFR":   nref_haps,
  "EUR":   nref_haps,
  "EAS":   nref_haps
}

sequence_length = 1000000
#Add in a recombination map with a hotspot!
rate_map = msprime.RateMap(position=[0, sequence_length / 2 - 1, sequence_length / 2, sequence_length],rate=[1e-8,1e-4,1e-8])
#Write out this recombination map across ALL POSITIONS NOT JUST THOSE WITH SNPS
recomb_map = rate_map.get_cumulative_mass(range(0,sequence_length + 1))
recomb_rate = rate_map.get_rate(range(0,sequence_length + 1))
#Write out this array 
np.savetxt("msprime_data/recombination_rate", recomb_rate, newline=" ")
np.savetxt("msprime_data/recombination_map", recomb_map, newline=" ")

#### Loop Across the replicates (regions) ####
for i in range(1,nreplicates):
  tree       = msprime.sim_ancestry(recombination_rate = rate_map,          
                                    sequence_length    = sequence_length,          
                                    demography         = demography,    #Use the demographic history from above
                                    samples            = sample_sizes)  #Use the sample sizes from the variables
                                    
  #Rate here refers to MUTATION rate                                
  mutated_ts = msprime.sim_mutations(tree, rate=1e-8, model = "binary") 
  #Now convert the mutated_ts to an haplotype array
  haps        = np.transpose(np.asarray(allel.HaplotypeArray(mutated_ts.genotype_matrix())))
  print(np.shape(haps))
  positions   = mutated_ts.tables.sites.position
  #Write out the array of the sites to directory
  np.savetxt("msprime_data/positions_pop_ADMIXED_region_%s" % (i), positions, newline=" ")
  #Extract the populations from the haplotype array (ADMIX,AFR,EUR,EAS)
  pop_ADMIXED = pd.DataFrame(haps[0:ngwas_haps*2,])
  pop_ADMIXED.to_csv("msprime_data/pop_ADMIXED_region_%s.csv" % (i), index = False)
  pop_AFR = pd.DataFrame(haps[ngwas_haps*2:(ngwas_haps*2 + nref_haps*2),])
  pop_AFR.to_csv("msprime_data/ref_pop_AFR_region_%s.csv" % (i), index = False)
  pop_EUR = pd.DataFrame(haps[(ngwas_haps*2 + nref_haps*2):(ngwas_haps*2 + nref_haps*4),])
  pop_EUR.to_csv("msprime_data/ref_pop_EUR_region_%s.csv" % (i), index = False)
  pop_EAS = pd.DataFrame(haps[(ngwas_haps*2 + nref_haps*4):,])
  pop_EAS.to_csv("msprime_data/ref_pop_EAS_region_%s.csv" % (i), index = False)


