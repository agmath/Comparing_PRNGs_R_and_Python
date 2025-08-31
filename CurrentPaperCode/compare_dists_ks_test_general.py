##################################
### Run KS-test Comparing      ###
### random sequences generated ###
### using Mersenne Twister in  ###
### R and Python               ###
##################################

import importlib.util
from scipy.stats import ks_2samp
import numpy as np

def compare_dists(seq_length = 2**14 - 1, max_rand = 99999, rand_seed = 1234, seq_1 = "mt_r", seq_2 = "mt_py"):
  spec = importlib.util.spec_from_file_location("gen_seqs", "CurrentPaperCode/generate_sequences.py")
  gen_seqs = importlib.util.module_from_spec(spec)
  spec.loader.exec_module(gen_seqs)
  
  sequence_r, sequence_py = gen_seqs.generate_sequences(seq_length, max_rand, rand_seed)
  
  ks_stat, p_value = ks_2samp(sequence_r, sequence_py)
  print(f"KS Statistic: {ks_stat:.4f}")
  print(f"p-value: {p_value:.4f}")

  return None


