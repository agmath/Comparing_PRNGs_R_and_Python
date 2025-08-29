##################################
### Helper function for        ###
### generating random sequences###
### using Mersenne Twister in  ###
### Python                     ###
##################################

# Python Code (Mersenne)
## generate random numbers in python
def generate_rands(seq_length = 2**14 - 1, max_rand = 99999, rand_seed = 1234):
  import numpy as np
  import pandas as pd
  
  np.random.seed(rand_seed)
  my_rands = np.random.randint(low = 0, high = max_rand + 1, size = seq_length)
  df_rands = pd.DataFrame({"sequence" : my_rands})
  
  return df_rands
