##################################
### Run KS-test Comparing      ###
### random sequences generated ###
### using Mersenne Twister and ###
### the XORshift128+           ###
##################################

import importlib.util
from scipy.stats import ks_2samp
import numpy as np

def compare_mt_xor_dists(seq_length = 2**14 - 1, max_rand = 99999, rand_seed_mt = 1234, rand_seed_xor1 = 1234, rand_seed_xor2 = 5678):
  spec_mt = importlib.util.spec_from_file_location("gen_seqs_mt", "CurrentPaperCode/generate_sequences.py")
  gen_seqs_mt = importlib.util.module_from_spec(spec_mt)
  spec_mt.loader.exec_module(gen_seqs_mt)
  
  _, mt_sequence_py = gen_seqs_mt.generate_sequences(seq_length, max_rand, rand_seed_mt)
  sequence_mt_scaled = mt_sequence_py/max_rand
  
  spec_xor = importlib.util.spec_from_file_location("gen_seqs_xor", "CurrentPaperCode/generate_rands_xor_shift.py")
  gen_seqs_xor = importlib.util.module_from_spec(spec_xor)
  spec_xor.loader.exec_module(gen_seqs_xor)
  
  sequence_xor = gen_seqs_xor.generate_xor_shift_128(seed1 = rand_seed_xor1, seed2 = rand_seed_xor2, n = seq_length)
  
  ks_stat, p_value = ks_2samp(sequence_mt_scaled, sequence_xor)
  print(f"KS Statistic: {ks_stat:.4f}")
  print(f"p-value: {p_value:.4f}")

  return None


