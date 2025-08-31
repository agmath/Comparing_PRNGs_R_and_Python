##################################
### Run KS-test Comparing      ###
### random sequences generated ###
### using Mersenne Twister and ###
### the XORshift128+           ###
##################################

import importlib.util
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.stats as stats

def plot_qq_mt_xor(seq_length = 2**14 - 1, max_rand = 99999, rand_seed_mt = 1234, rand_seed_xor1 = 1234, rand_seed_xor2 = 5678):
  spec_mt = importlib.util.spec_from_file_location("gen_seqs_mt", "CurrentPaperCode/generate_sequences.py")
  gen_seqs_mt = importlib.util.module_from_spec(spec_mt)
  spec_mt.loader.exec_module(gen_seqs_mt)
  
  _, mt_sequence_py = gen_seqs_mt.generate_sequences(seq_length, max_rand, rand_seed_mt)
  sequence_mt_scaled = mt_sequence_py/max_rand
  
  spec_xor = importlib.util.spec_from_file_location("gen_seqs_xor", "CurrentPaperCode/generate_rands_xor_shift.py")
  gen_seqs_xor = importlib.util.module_from_spec(spec_xor)
  spec_xor.loader.exec_module(gen_seqs_xor)
  
  sequence_xor = gen_seqs_xor.generate_xor_shift_128(seed1 = rand_seed_xor1, seed2 = rand_seed_xor2, n = seq_length)
  
  plt.clf()
  plt.figure(figsize = (12, 3))
  plt.subplot(1, 3, 1)
  stats.probplot(sequence_mt_scaled, dist="uniform", plot=plt)
  plt.title("Uniform Q-Q Plot (MT)")
  
  plt.subplot(1, 3, 3)
  stats.probplot(sequence_xor, dist="uniform", plot=plt)
  plt.title("Uniform Q-Q Plot (XOR)")
  plt.show()

  return None


