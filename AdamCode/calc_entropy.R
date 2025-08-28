library(dplyr)

entropy_collector <- function(my_seqs){
  collector_df <- tibble()
  n_seqs <- my_seqs %>% ncol()
  
  for(i in 1:n_seqs){
    seq_length <- length(my_seqs[[i]])
    entropy <- 0
    curr_seq_counts <- table(my_seqs[[i]])
    
    curr_seq_probs <- curr_seq_counts/seq_length
    entropy <- entropy + sum(curr_seq_probs*log2(curr_seq_probs))
    
    collector_df <- collector_df %>%
      bind_rows(
        tibble(
          sequence = names(my_seqs)[i],
          entropy = -entropy
        )
      )
  }
  return(collector_df)
}