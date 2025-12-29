##################################
### Helper function for        ###
### generating random sequences###
### using Mersenne Twister --  ###
### calls generate_rands_*()   ###
##################################

import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import importlib.util
import numpy as np
import generate_rands_mersenne_twister as py_gen

def generate_sequences(seq_length = 2**14 - 1, 
                       min_rand = 0, 
                       max_rand = 99999, 
                       rand_seed = 1234):
  

  # Source the R script
  robjects.r['source']('generate_rands_mersenne_twister.R')
  
  # Generate the random numbers with R
  generate_rands_r = robjects.globalenv['generate_rands_mersenne_twister']
  my_rands_r = generate_rands_r(seq_length, min_rand, max_rand, rand_seed)


  # Generate the random numbers with Python
  my_rands_py = py_gen.generate_rands_mersenne_twister(seq_length, 
                                                        min_rand, 
                                                        max_rand, 
                                                        rand_seed)
  
  # Compare the generated sequences
  sequence_r = np.array(my_rands_r.rx2("sequence"))
  sequence_py = np.array(my_rands_py["sequence"])

  return sequence_r, sequence_py


