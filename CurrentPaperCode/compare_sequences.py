import importlib.util
import numpy as np

def compare_seqs(seq_length = 2**14 - 1, max_rand = 99999, rand_seed = 1234):
  spec = importlib.util.spec_from_file_location("gen_seqs", "CurrentPaperCode/generate_sequences.py")
  gen_seqs = importlib.util.module_from_spec(spec)
  spec.loader.exec_module(gen_seqs)
  
  sequence_r, sequence_py = gen_seqs.generate_sequences(seq_length, max_rand, rand_seed)
  
  are_equal = np.array_equal(sequence_r, sequence_py)
  print(f"Sequences identical: {are_equal}")
  
  return None


