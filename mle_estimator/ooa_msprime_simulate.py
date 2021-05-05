#Run OOA simulations with MLE estimator
#We want to use the existing code from the previous chapter here
#Uncomment if running locally
library(reticulate)
use_condaenv(condaenv = "python_3.7", conda = "/well/mcvean/mtutert/anaconda3/bin/conda")

import msprime
import tskit
import allel; print('scikit-allel', allel.__version__)
import numpy as np
import random
import numpy as np
import sys
import pandas as pd
import os
import matplotlib.pyplot as plt
import stdpopsim

#Set demographic simulation with stdpopsim
species = stdpopsim.get_species("HomSap")
model   = species.get_demographic_model("OutOfAfrica_3G09")
#Create the contig (ie: change the length of the simulation and the recomb rate params)
contig  = species.get_contig("chr22", length_multiplier=0.01)
#Choose sample sizes for each population (YRI/CEU/CHB)
YRI_sample_size = 1000
CEU_sample_size = 1000
CHB_sample_size = 1000
samples = model.get_samples(YRI_sample_size, CEU_sample_size, CHB_sample_size)
engine = stdpopsim.get_engine("msprime")

#Run simulations across replicates and save to cluster
for i in range(1,1000):
  ts = engine.simulate(model, contig, samples, seed = i )
  #Convert to haplotype array using scikit allele again noting population structure
  haplotype_array = np.asarray(allel.HaplotypeArray(ts.genotype_matrix()))
  #Split into the different panels
  haplotype_array = np.transpose(haplotype_array)
  #Read out as csv with pandas
  YRI  = pd.DataFrame((haplotype_array[0:YRI_sample_size:,]))
  CEU  = pd.DataFrame((haplotype_array[(YRI_sample_size):(YRI_sample_size+CEU_sample_size):,]))
  CHB  = pd.DataFrame((haplotype_array[(YRI_sample_size+CEU_sample_size):(YRI_sample_size+CEU_sample_size+CHB_sample_size):,]))
  
  YRI.to_csv("ooa_msprime_data/YRI_replicate_%d.csv" % (i), index=False)    
  CEU.to_csv("ooa_msprime_data/CEU_replicate_%d.csv" % (i), index=False)  
  CHB.to_csv("ooa_msprime_data/CHB_replicate_%d.csv" % (i), index=False)




