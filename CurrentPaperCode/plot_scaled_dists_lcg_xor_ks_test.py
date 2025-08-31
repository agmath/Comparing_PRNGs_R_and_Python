##################################
### Run KS-test for comparing  ###
### random sequences generated ###
### using Mersenne Twister and ###
### Linear Congruence          ###
##################################

import importlib.util
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def plot_dists_xor_lcg(seq_length = 2**14 - 1, max_rand = 99999, rand_seed_lcg = 1234, a = 1664525, c = 1013904223, m = 2**32, rand_seed_xor1 = 1234, rand_seed_xor2 = 5678):
  
  spec_xor = importlib.util.spec_from_file_location("gen_seqs_xor", "CurrentPaperCode/generate_rands_xor_shift.py")
  gen_seqs_xor = importlib.util.module_from_spec(spec_xor)
  spec_xor.loader.exec_module(gen_seqs_xor)
  
  sequence_xor = gen_seqs_xor.generate_xor_shift_128(seed1 = rand_seed_xor1, seed2 = rand_seed_xor2, n = seq_length)
  
  # Generate random sequence with LCG
  spec_lcg = importlib.util.spec_from_file_location("gen_seqs", "CurrentPaperCode/generate_rands_lcg.py")
  gen_seqs_lcg = importlib.util.module_from_spec(spec_lcg)
  spec_lcg.loader.exec_module(gen_seqs_lcg)
  
  sequence_lcg = gen_seqs_lcg.generate_lcg(rand_seed_lcg, seq_length, a, c, m)
  sequence_lcg_scaled = np.array(sequence_lcg) / m
  
  plt.clf()
  sns.kdeplot(sequence_xor, label='Sequence XOR', color='blue')
  sns.kdeplot(sequence_lcg_scaled, label='Sequence LCG', color='red')
  plt.title('Density Plot Comparison')
  plt.legend()
  plt.show()

  return None


