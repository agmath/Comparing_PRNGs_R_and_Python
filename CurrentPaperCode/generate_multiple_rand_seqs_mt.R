##################################
### Helper function for        ###
### generating random sequences###
### using Mersenne Twister in  ###
### R                          ###
##################################

library(dplyr)


generate_mt_rands <- function(num_seqs = 10, 
                              seq_length = 2^14 - 1,
                              min_rand = 0,
                              max_rand = 1e5 - 1,
                              seed = 1234){
  set.seed(seed)
  
  #sample uses the Mersenne Twister by default
  pRNGseqs_mt <- replicate(
    num_seqs,
    sample(min_rand:max_rand, size = seq_length, replace = TRUE),
    simplify = FALSE
  ) %>%
    as_tibble(.name_repair = "minimal") %>%
    setNames(paste0("r_mt_", 1:10))

  set.seed(NULL)
  
  return(pRNGseqs_mt)
}
