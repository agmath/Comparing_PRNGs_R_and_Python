import importlib.util
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def plot_dists(seq_length = 2**14 - 1, max_rand = 99999, rand_seed = 1234):
  spec = importlib.util.spec_from_file_location("gen_seqs", "CurrentPaperCode/generate_sequences.py")
  gen_seqs = importlib.util.module_from_spec(spec)
  spec.loader.exec_module(gen_seqs)
  
  sequence_r, sequence_py = gen_seqs.generate_sequences(seq_length, max_rand, rand_seed)
  
  plt.clf()
  sns.kdeplot(sequence_r, label='Sequence R', color='blue')
  sns.kdeplot(sequence_py, label='Sequence Python', color='red')
  plt.title('Density Plot Comparison')
  plt.legend()
  plt.show()


  return None


