##################################
### Produce QQ Plots for       ###
### random sequences generated ###
### using Mersenne Twister in  ###
### R and Python               ###
##################################

import importlib.util
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.stats as stats

def plot_qq(seq_length = 2**14 - 1, max_rand = 99999, rand_seed = 1234):
  spec = importlib.util.spec_from_file_location("gen_seqs", "CurrentPaperCode/generate_sequences.py")
  gen_seqs = importlib.util.module_from_spec(spec)
  spec.loader.exec_module(gen_seqs)
  
  sequence_r, sequence_py = gen_seqs.generate_sequences(seq_length, max_rand, rand_seed)
  
  plt.clf()
  plt.figure(figsize = (12, 3))
  plt.subplot(1, 3, 1)
  stats.probplot(sequence_r, dist="uniform", plot=plt)
  plt.title("Uniform Q-Q Plot (R)")
  
  plt.subplot(1, 3, 3)
  stats.probplot(sequence_py, dist="uniform", plot=plt)
  plt.title("Uniform Q-Q Plot (Python)")
  plt.show()

  return None


