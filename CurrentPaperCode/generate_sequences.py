def generate_sequences(seq_length = 2**14 - 1, max_rand = 99999, rand_seed = 1234):
  import rpy2.robjects as robjects
  from rpy2.robjects.packages import importr
  import importlib.util
  import numpy as np

  # Source the R script
  robjects.r['source']('CurrentPaperCode/generate_rands_mersenne_twister.R')
  
  # Generate the random numbers with R
  generate_rands_r = robjects.globalenv['generate_rands']
  my_rands_r = generate_rands_r(seq_length, max_rand, rand_seed)

  # Source the Python Script
  spec = importlib.util.spec_from_file_location("py_gen", "CurrentPaperCode/generate_rands_mersenne_twister.py")
  py_gen = importlib.util.module_from_spec(spec)
  spec.loader.exec_module(py_gen)
  
  # Generate the random numbers with Python
  my_rands_py = py_gen.generate_rands(seq_length, max_rand, rand_seed)
  
  # Compare the generated sequences
  sequence_r = np.array(my_rands_r.rx2("sequence"))
  sequence_py = np.array(my_rands_py["sequence"])

  return sequence_r, sequence_py


