library(dplyr)
library(tidyr)

calc_mean <- function(my_seqs){
  collector_df <- my_seqs %>%
  summarize_all(~ mean(.)) %>%
  pivot_longer(everything()) %>%
  rename(
    sequence = name,
    mean = value
  )
  
  return(collector_df)
}