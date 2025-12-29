##################################
### Helper function for        ###
### generating random sequences###
### using Mersenne Twister in  ###
### Python                     ###
##################################

import numpy as np
import pandas as pd
  
# Python Code (Mersenne)
## generate random numbers in python
def generate_rands_mersenne_twister(seq_length = 2**14 - 1, 
                                    min_rand = 0, 
                                    max_rand = 99999, 
                                    rand_seed = 1234):
  
  np.random.seed(rand_seed)
  #randint: low is inclusive, high is exclusive
  my_rands = np.random.randint(low = min_rand, high = max_rand + 1, size = seq_length)
  df_rands = pd.DataFrame({"sequence" : my_rands})
  
  return df_rands
