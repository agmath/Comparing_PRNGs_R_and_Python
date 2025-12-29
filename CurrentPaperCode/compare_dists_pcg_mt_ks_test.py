##################################
### Run KS-test for comparing  ###
### random sequences generated ###
### using Mersenne Twister and ###
### Permuted Congruence        ### 
### Generator Algorithm        ###
##################################

import importlib.util
from scipy.stats import ks_2samp
from numpy.random import Generator, PCG64DXSM
import generate_sequences as gen_seqs_mt
import get_seeds_for_null_distributions as get_null_seeds

def compare_dists_pcg_mt_ks_test(seq_length = 2**14 - 1, 
                                 min_rand = 0,
                                 max_rand = 99999, 
                                 rand_seed = 1234):
                                   
  # Generate random sequence with Mersenne Twister
  # We will not use the sequence_r_mt here
  sequence_r_mt, sequence_py_mt = gen_seqs_mt.generate_sequences(seq_length, 
                                              min_rand, 
                                              max_rand, 
                                              rand_seed)

  ###############################################################################
  # Generate random sequence with PCG
  ###############################################################################
  # https://realpython.com/numpy-random-number-generator/ indicates that this PCG
  # will become the default version in future Python versions.
  #
  # Documentation:
  # https://numpy.org/doc/2.3/reference/random/index.html#module-numpy.random
  ###############################################################################
  rng = Generator(PCG64DXSM(rand_seed))
  sequence_py_pcg = rng.integers(low=min_rand, 
                                 high=max_rand+1, 
                                 size=seq_length)  #high is exclusive


  ks_stat, p_value = ks_2samp(sequence_py_pcg, sequence_py_mt)
  #Printing if desired
  #print(f"KS Statistic: {ks_stat:.4f}")
  #print(f"p-value: {p_value:.4f}")

  return [ks_stat, p_value, p_value < 0.05]  #True indicates reject_H0 at a 0.05 conf level



#Function to generate many pairs of PCG and MT sequences to compare.
#Run a ks test on each pair.
#Record the KS test statistic, p-value, and if the hypothess is rejected or not.
#     H0:  The Distributions of the PCG sequence and MT sequence are the same
#    HA:  The Distributions of the PCG sequence and MT sequence are not the same
def generate_pcg_mt_ks_test_null_distribution(seq_length = 2**14 - 1, 
                                              min_rand = 0,
                                              max_rand = 99999):
  
  #By default, the seeds will be between 0-2**31 (R's biggest integer) and
  #10000 seeds will be generated.  If the seeds have been generated prior, 
  #they will not be overwritten but a file with the previously generated
  #seeds will be read in.
  seeds_df = get_null_seeds.get_seeds_for_null_distributions()
                                        
  #use each seed to launch a compare_dists_mt_mt_ks_test                                      
  seeds_df['ks_test_results'] = seeds_df['seed'].progress_apply(
    lambda current_seed: compare_dists_pcg_mt_ks_test(
        seq_length, 
        min_rand,
        max_rand, 
        current_seed))
        
  #The ks_test_results column has a list with: [ks test statistic, p-value, reject H0 or not].
  #Let's split the column containing a 3 element list into 3 distinct columns.                              
  seeds_df['ks_test_statistic'] = seeds_df['ks_test_results'].apply(lambda row: row[0])
  seeds_df['p_value'] = seeds_df['ks_test_results'].apply(lambda row: row[1])
  seeds_df['reject_H0'] = seeds_df['ks_test_results'].apply(lambda row: row[2])
  seeds_df.drop('ks_test_results', axis=1, inplace  =True)
  
  return seeds_df



