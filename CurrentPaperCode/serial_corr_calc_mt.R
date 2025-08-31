library(dplyr)

serial_correlation_calc_mt_seqs <- function(num_lag = 1, num_seqs = 10, seq_length = 2^14 - 1, max_rand = 1e5 - 1, seed = 1234){
  source("CurrentPaperCode/generate_multiple_rand_seqs_mt.R")
  
  #Generate sequences using Mersenne Twister
  pRNGseqs_mt <- generate_mt_rands(num_seqs, seq_length, max_rand, seed)
  
  results_df <- data.frame()
  
  for(i in 1:num_seqs){
    curr_seq <- pRNGseqs_mt[[i]]
    lagged <- lag(pRNGseqs_mt[[i]], n = num_lag)
    serial_cor <- cor(curr_seq, lagged, 
                      use = "pairwise.complete.obs")
    
    results_df <- results_df %>%
      bind_rows(
        tibble(
          sequence = names(pRNGseqs_mt)[i],
          serial_correlation = serial_cor
        )
      )
  }
  return(results_df)
}
