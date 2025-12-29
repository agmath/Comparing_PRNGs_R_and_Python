##################################
### Run KS-test Comparing      ###
### random sequences generated ###
### using Mersenne Twister in  ###
### R and Python               ###
##################################

import importlib.util
from scipy.stats import ks_2samp
import numpy as np
import pandas as pd
import generate_sequences as gen_seqs
import get_seeds_for_null_distributions as gen_seeds

def compare_dists_mt_mt_ks_test(seq_length = 2**14 - 1, 
                                min_rand = 0,
                                max_rand = 99999, 
                                seed = 1234):
  """
  This function generates a sequence in R and a sequence in Python using the 
  given parameters.  The function that performs a Kolmogorov–Smirnov Test and
  returns a list containing:
    
  * KS Test Statistic
  * P-Value
  * Whether or not the null Hypothesis was rejected
    True:   The null was rejected
    False:  The null was not rejected
  """                                  

  #Get the sequences
  sequence_r, sequence_py = gen_seqs.generate_sequences(seq_length, 
                                                        min_rand, 
                                                        max_rand, 
                                                        seed)
  #Run the KS Test
  ks_stat, p_value = ks_2samp(sequence_r, sequence_py)

  #Return the desired values
  return [ks_stat, p_value, p_value < 0.05]  #True indicates reject_H0 at a 0.05 conf level


  
def generate_mt_mt_ks_test_null_distribution(seq_length = 2**14 - 1, 
                                            min_rand = 0,
                                            max_rand = 99999):
  
  #By default, the seeds will be between 0-2**31 (R's biggest integer) and
  #10000 seeds will be generated.  If the seeds have been generated prior, 
  #they will not be overwritten but a file with the previously generated
  #seeds will be read in.
  seeds_df = gen_seeds.get_seeds_for_null_distributions()
                                        
  
  #use each seed to launch a compare_dists_mt_mt_ks_test                                      
  seeds_df['ks_test_results'] = seeds_df['seed'].progress_apply(lambda current_seed: 
                                        compare_dists_mt_mt_ks_test(seq_length, 
                                                                    min_rand,
                                                                    max_rand, 
                                                                    seed = current_seed))
                                                                    
  #Theks_test_results column has a list with: [ks test statistic, p-value, reject H0 or not].
  #Let's split the column containing a 3 element list into 3 distinct columns.
  seeds_df['ks_test_statistic'] = seeds_df['ks_test_results'].apply(lambda row: row[0])
  seeds_df['p_value'] = seeds_df['ks_test_results'].apply(lambda row: row[1])
  seeds_df['reject_H0'] = seeds_df['ks_test_results'].apply(lambda row: row[2])
  seeds_df.drop('ks_test_results', axis=1, inplace=True)
  
  return seeds_df

