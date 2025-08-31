##################################
### Helper function for        ###
### generating random sequences###
### using XORshift128+         ###
### implemented in Python      ###
##################################

def generate_xor_shift_128(seed1, seed2, n):
  import numpy as np
  import pandas as pd

  s1 = seed1
  s2 = seed2
  results = []
  for _ in range(n):
    s1 ^= (s1 << 23) & 0xFFFFFFFFFFFFFFFF
    s1 ^= (s1 >> 17)
    s1 ^= (s1 << 26) & 0xFFFFFFFFFFFFFFFF
    t = s1
    s1 = s2
    s2 ^= (s2 >> 13)
    s2 ^= (s2 << 17) & 0xFFFFFFFFFFFFFFFF
    s2 ^= (t + s2) & 0xFFFFFFFFFFFFFFFF
    results.append(((s1 + s2) & 0xFFFFFFFFFFFFFFFF) / 2**64)
    
  return results
