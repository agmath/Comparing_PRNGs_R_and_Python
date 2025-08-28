library(dplyr)

serial_corr_collector <- function(my_seqs){
  collector_df <- tibble()
  n_seqs <- my_seqs %>% ncol()
  
  for(i in 1:n_seqs){
    curr_seq <- my_seqs[[i]]
    lagged <- lag(my_seqs[[i]])
    serial_cor <- cor(curr_seq, lagged, 
                      use = "pairwise.complete.obs")
    
    collector_df <- collector_df %>%
      bind_rows(
        tibble(
          sequence = names(my_seqs)[i],
          serial_correlation = serial_cor
        )
      )
  }
  return(collector_df)
}

