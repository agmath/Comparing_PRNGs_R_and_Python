##################################
### Helper function for        ###
### generating random sequences###
### using Mersenne Twister in  ###
### R                          ###
##################################

# R code (Mersenne)
## generate random numbers in r

generate_rands <- function(seq_length = 2^14 - 1, max_rand = 1e5 - 1, rand_seed = 1234){
  set.seed(rand_seed)
  my_rands <- sample(0:max_rand, size = seq_length, replace = TRUE)
  
  df_rands <- data.frame(sequence = my_rands)
  
  return(df_rands)
}