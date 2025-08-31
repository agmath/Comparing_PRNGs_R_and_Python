library(dplyr)

calc_descriptives_mt_seqs <- function(num_seqs = 10, seq_length = 2^14 - 1, max_rand = 1e5 - 1, seed = 1234){
  source("CurrentPaperCode/generate_multiple_rand_seqs_mt.R")
  
  #Generate sequences using Mersenne Twister
  pRNGseqs_mt <- generate_mt_rands(num_seqs, seq_length, max_rand, seed)
  
  results_df <- data.frame(
    seq_name = NA,
    mean = NA,
    sd = NA,
    cv_pct = NA,
    min = NA,
    max = NA
  )
  
  for(i in 1:num_seqs){
    seq_sum_df <- data.frame(
      seq_name = names(pRNGseqs_mt)[i],
      mean = mean(pRNGseqs_mt[[i]]),
      sd = sd(pRNGseqs_mt[[i]]),
      min = min(pRNGseqs_mt[[i]]),
      max = max(pRNGseqs_mt[[i]])
    ) %>%
      mutate(cv_pct = 100*(sd/mean))
    
    results_df <- results_df %>%
      bind_rows(seq_sum_df)
  }
  
  return(results_df)
}