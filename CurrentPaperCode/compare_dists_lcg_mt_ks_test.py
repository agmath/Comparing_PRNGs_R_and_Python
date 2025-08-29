import importlib.util
from scipy.stats import ks_2samp
import numpy as np

def compare_dists_mt_lcg(seq_length = 2**14 - 1, max_rand = 99999, rand_seed = 1234, a = 1664525, c = 1013904223, m = 2**32):
  # Generate random sequrnce with Mersenne Twister
  spec_mt = importlib.util.spec_from_file_location("gen_seqs", "CurrentPaperCode/generate_sequences.py")
  gen_seqs_mt = importlib.util.module_from_spec(spec_mt)
  spec_mt.loader.exec_module(gen_seqs_mt)
  
  _, sequence_mt = gen_seqs_mt.generate_sequences(seq_length, max_rand, rand_seed)
  sequence_mt_scaled = np.array(sequence_mt) / max_rand
  
  # Generate random sequence with LCG
  spec_lcg = importlib.util.spec_from_file_location("gen_seqs", "CurrentPaperCode/generate_rands_lcg.py")
  gen_seqs_lcg = importlib.util.module_from_spec(spec_lcg)
  spec_lcg.loader.exec_module(gen_seqs_lcg)
  
  sequence_lcg = gen_seqs_lcg.generate_lcg(rand_seed, seq_length, a, c, m)
  sequence_lcg_scaled = np.array(sequence_lcg) / m
  
  ks_stat, p_value = ks_2samp(sequence_mt_scaled, sequence_lcg_scaled)
  print(f"KS Statistic: {ks_stat:.4f}")
  print(f"p-value: {p_value:.4f}")

  return None


