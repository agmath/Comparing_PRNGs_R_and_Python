##################################
### Helper function for        ###
### generating random sequences###
### using Mersenne Twister in  ###
### R                          ###
##################################

# R code (Mersenne)
## generate a single sequence of random numbers in r
generate_rands_mersenne_twister <- function(seq_length = 2^14 - 1, 
                           min_rand = 0, 
                           max_rand = 1e5 - 1, 
                           rand_seed = 1234){
  set.seed(rand_seed)
  
  #Inclusive min/max
  my_rands <- sample(min_rand:max_rand, size = seq_length, replace = TRUE)
  
  df_rands <- data.frame(sequence = my_rands)
  
  return(df_rands)
}