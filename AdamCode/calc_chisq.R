library(dplyr)

chi_sq_collector <- function(my_seqs){
  collection_df <- tibble()
  n_seqs <- my_seqs %>% ncol()
  for(i in 1:n_seqs){
    curr_seq_counts <- table(my_seqs[[i]])
    chisq_results <- chisq.test(curr_seq_counts) %>%
      broom.mixed::tidy() %>%
      mutate(sequence = names(my_seqs[i]))
    collection_df <- collection_df %>%
      bind_rows(chisq_results)
  }
  collection_df <- collection_df %>%
    select(sequence, statistic, p.value)
  return(collection_df)
}
