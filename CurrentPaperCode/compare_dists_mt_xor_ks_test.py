##################################
### Run KS-test Comparing      ###
### random sequences generated ###
### using Mersenne Twister and ###
### the XORshift128+           ###
##################################

import importlib.util
from scipy.stats import ks_2samp
import numpy as np
import pandas as pd
import generate_sequences as gen_seqs_mt
import get_seeds_for_null_distributions as get_null_seeds
import generate_rands_xor_shift as gen_seqs_xor

def compare_mt_xor_ks_test(seq_length = 2**14 - 1, 
                        min_rand = 0,
                        max_rand = 99999, 
                        rand_seed_mt = 1234, 
                        rand_seed_xor1 = 1234, 
                        rand_seed_xor2 = 5678):
  
  # Generate random sequence with Mersenne Twister
  # We will not use the sequence_r_mt here
  sequence_r_mt, sequence_py_mt = gen_seqs_mt.generate_sequences(seq_length, 
                                              min_rand, 
                                              max_rand, 
                                              rand_seed_mt)
                                              
                                              
  # Scale the MT sequence                                            
  sequence_mt_scaled = sequence_py_mt/max_rand
  
  ###############################################################################
  # Generate random sequence with sequence_xor
  ###############################################################################
  sequence_xor = gen_seqs_xor.generate_xor_shift_128(seed1 = rand_seed_xor1, 
                                                      seed2 = rand_seed_xor2, 
                                                      n = seq_length)
  
  ks_stat, p_value = ks_2samp(sequence_mt_scaled, sequence_xor)
  #Printing if desired
  #print(f"KS Statistic: {ks_stat:.4f}")
  #print(f"p-value: {p_value:.4f}")
  
  return [ks_stat, p_value, p_value < 0.05]


#Helper function:  We use this to apply one ks test to a single row
#when generating many ks test results for a null distribution
def perform_single_ks_test_per_row(current_row, 
                                  seq_length, 
                                  min_rand,
                                  max_rand):
  
  #Without these next two lines, we would get an error that extends because
  #rpy2 doesn’t know how to convert Python numpy integers to/from R integers.
  seed1 = int(current_row['seed1'])
  seed2 = int(current_row['seed2'])
    
  #Perform a ks test between an MT sequence and an XOR sequence  
  return compare_mt_xor_ks_test(seq_length = seq_length, 
                                min_rand = min_rand,
                                max_rand = max_rand,
                                rand_seed_mt = seed1,
                                rand_seed_xor1 = seed1,
                                rand_seed_xor2 = seed2)
                            

#Function to generate many pairs of XOR and MT sequences to compare.
#Run a ks test on each pair.
#Record the KS test statistic, p-value, and if the hypothess is rejected or not.
#     H0:  The Distributions of the XOR sequence and MT sequence are the same
#    HA:  The Distributions of the XOR sequence and MT sequence are not the same
def generate_xor_mt_ks_test_null_distribution(seq_length = 2**14 - 1, 
                                              min_rand = 0,
                                              max_rand = 99999):
  
  #By default, the seeds will be between 0-2**31 (R's biggest integer) and
  #10000 seeds will be generated.  If the seeds have been generated prior, 
  #they will not be overwritten but a file with the previously generated
  #seeds will be read in.
  seeds_df = get_null_seeds.get_seeds_for_null_distributions()
  
  #For the xorshift algorithm, we need 2 random seeds.
  #So split the seeds into 2 columns. 
  #After doing so, the half_seeds_df will look like this
  #but be half as long as the seeds_df.
  ###############################
  #        seed1       seed2
  #0  1263166662    60250403
  #1   289859650   122664920
  #2  1883615152  1705704388
  #3  1283397658  1277261284
  #4  1405373200  1931352802
  #    etc
  ###############################
  num_seeds = len(seeds_df)//2
  
  #The only column in seeds_df is "seed" so this is the 0th column 
  seed1_col = seeds_df.iloc[:num_seeds, 0].reset_index(drop=True)
  seed2_col = seeds_df.iloc[num_seeds:, 0].reset_index(drop=True)
  
  half_seeds_df = pd.DataFrame({"seed1": seed1_col, 
                                "seed2": seed2_col})
                              
  #Run a ks test for each row
  half_seeds_df['ks_test_results'] = half_seeds_df.progress_apply(
    lambda row: perform_single_ks_test_per_row(row, seq_length, min_rand, max_rand), 
    axis=1)
  
  #The ks_test_results column has a list with: [ks test statistic, p-value, reject H0 or not].
  #Let's split the column containing a 3 element list into 3 distinct columns.                              
  half_seeds_df['ks_test_statistic'] = half_seeds_df['ks_test_results'].apply(lambda row: row[0])
  half_seeds_df['p_value'] = half_seeds_df['ks_test_results'].apply(lambda row: row[1])
  half_seeds_df['reject_H0'] = half_seeds_df['ks_test_results'].apply(lambda row: row[2])
  half_seeds_df.drop('ks_test_results', axis=1, inplace  =True)
  
  return half_seeds_df
