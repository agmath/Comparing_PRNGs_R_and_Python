import os
import numpy as np
import pandas as pd

#The file where the random seeds will be saved
KS_TEST_STATISTICS_NULL_DISTRIB_SEEDS_FILE = 'data/random_seeds_ks_statistic_null_distrib.csv'


def get_seeds_for_null_distributions(min_rand = 0, 
                                   max_rand = 2**31-1,  #R's biggest integer
                                   num_seeds = 10000, 
                                   overwrite_prior_seeds = False):
                                     
   
  #If the file is not there or the overwrite flag is set to true, 
  #regenerate the seeds
  if not(os.path.isfile(KS_TEST_STATISTICS_NULL_DISTRIB_SEEDS_FILE)) or (overwrite_prior_seeds):
    
    #reset random's seed to the default state
    np.random.seed()
    
    # Generate random seeds - high is exclusive.
    seeds = np.random.randint(low=min_rand, high=max_rand + 1, size=num_seeds)
     
    # Save to CSV
    pd.DataFrame({'seed': seeds}).to_csv(KS_TEST_STATISTICS_NULL_DISTRIB_SEEDS_FILE, index=False)
    
  else: #Just return the seeds
    
    seeds = pd.read_csv(KS_TEST_STATISTICS_NULL_DISTRIB_SEEDS_FILE)
    
  #Return the seeds too
  return seeds
