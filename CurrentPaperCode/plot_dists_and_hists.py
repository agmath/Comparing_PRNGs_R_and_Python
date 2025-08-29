import importlib.util
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def plot_dists_and_hists(seq_length = 2**14 - 1, max_rand = 99999, rand_seed = 1234):
  spec = importlib.util.spec_from_file_location("gen_seqs", "CurrentPaperCode/generate_sequences.py")
  gen_seqs = importlib.util.module_from_spec(spec)
  spec.loader.exec_module(gen_seqs)
  
  sequence_r, sequence_py = gen_seqs.generate_sequences(seq_length, max_rand, rand_seed)
  
  plt.clf()
  plt.hist(sequence_r, bins=50, density=True, alpha=0.5, color='blue')
  plt.hist(sequence_py, bins=50, density=True, alpha=0.5, color='red')
  sns.kdeplot(sequence_r, label='Sequence R', color='blue')
  sns.kdeplot(sequence_py, label='Sequence Python', color='red')
  plt.title('Distribution Comparisons')
  plt.legend()
  plt.show()


  return None


