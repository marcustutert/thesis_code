#Running 2 way admixture (AA-like)
#Code to generate split populations (ancestrally) in msprime
#If running this on Rstudio on rescomp server
use_condaenv(condaenv = "python_3.7", conda = "/well/mcvean/mtutert/anaconda3/bin/conda")
library(reticulate)
import msprime
import allel
print('scikit-allel', allel.__version__)
import numpy as np
import pandas as pd

#Run A-B admixture with msprime
#Create the appropriate demographic history w variable 'demography'
demography = msprime.Demography()
demography.add_population(name="A",      initial_size     = 14000)
demography.add_population(name="B",      initial_size     = 5000)
demography.add_population(name="ADMIX",  initial_size     = 30000)
demography.add_population(name="ANC",    initial_size     = 7000)
demography.add_admixture(time=12, derived="ADMIX", ancestral=["A", "B"], proportions=[0.2, 0.8])
demography.add_population_split(time=1000, derived=["A", "B"], ancestral="ANC")

####Variables#### 
#These are the DIPLOID sample sizes
ngwas_haps         = 1000     #GWAS population here is the ADMIXED population
nref_haps          = 1000     #Reference popualtion is either A,B or AB (OUT of sample)
nreplicates        = 1      #This needs to be the same number as whatever is in the top of the SNAKEMAKE file
###################

sample_sizes = {
  #Note the conversion to diploid sample sizes
  "ADMIX": ngwas_haps,
  "A":     nref_haps,
  "B":     nref_haps
}

#### Loop Across the replicates (regions) ####
for i in range(1,nreplicates):
  tree       = msprime.sim_ancestry(recombination_rate = 1e-10,          #Change this if appropriate!
                                    sequence_length    = 4000,          #Change this if appropriate!
                                    demography         = demography,    #Use the demographic history from above
                                    samples            = sample_sizes)  #Use the sample sizes from the variables
                                    
  #Rate here refers to MUTATION rate                                
  mutated_ts = msprime.sim_mutations(tree, rate=1e-2, model = "binary") 
  #Now convert the mutated_ts to an haplotype array
  haps        = np.transpose(np.asarray(allel.HaplotypeArray(mutated_ts.genotype_matrix())))
  positions   = mutated_ts.tables.sites.position
  #Write out the array of the sites to directory
  np.savetxt("msprime_data/positions_pop_ADMIXED_region_%s" % (i), positions, newline=" ")
  #Extract the populations from the haplotype array (Admixed, A, B)
  pop_ADMIXED = pd.DataFrame(haps[0:ngwas_haps*2,])
  pop_ADMIXED.to_csv("msprime_data/pop_ADMIXED_region_%s.csv" % (i), index = False)
  pop_A = pd.DataFrame(haps[ngwas_haps*2:(ngwas_haps*2 + nref_haps*2),])
  pop_A.to_csv("msprime_data/ref_pop_A_region_%s.csv" % (i), index = False)
  pop_B = pd.DataFrame(haps[(ngwas_haps*2 + nref_haps*2):,])
  pop_B.to_csv("msprime_data/ref_pop_B_region_%s.csv" % (i), index = False)


