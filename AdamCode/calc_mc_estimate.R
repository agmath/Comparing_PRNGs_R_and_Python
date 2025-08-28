library(dplyr)

source("AdamCode/generate_rands_mersenne_twister.R")

mc_sim_collector <- function(my_seqs){
  collection_df <- tibble()
  n_seqs <- my_seqs %>% ncol()
  
  for(i in 1:n_seqs){
    curr_seq <- my_seqs[[i]]
    scaled_seq <- curr_seq/max_rand
    seq_len <- length(scaled_seq)
    pi_est <- 4*(1/seq_len)*sum(sqrt(1 - scaled_seq^2))
    
    collection_df <- collection_df %>%
      bind_rows(
        tibble(
          sequence = names(my_seqs)[i],
          pi_estimate = pi_est
        )
      )
  }
  return(collection_df)
}
