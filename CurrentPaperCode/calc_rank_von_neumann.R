library(EnvStats)

calc_rank_von_neumann <- function(num_seqs = 10, 
                                seq_length = 2^14 - 1, 
                                min_rand = 0,
                                max_rand = 1e5 - 1, 
                                seed = 1234){
  
  source("generate_multiple_rand_seqs_mt.R")
  
  #Generate sequences using Mersenne Twister
  pRNGseqs_mt <- generate_mt_rands(num_seqs, seq_length, min_rand, max_rand, seed)
  
  results_df <- data.frame(test_statistic = c(), p_value = c(), reject_H0 = c())
  
  for(i in 1:num_seqs){
        
    rank_von_neumann_test = serialCorrelationTest(pRNGseqs_mt[[i]], 
                                 test = "rank.von.Neumann", 
                                 alternative = "two.sided", 
                                 conf.level = 0.95)
    
    
    ###Checking if we have 15% repeated values
    #sorted_vector_asc <- sort(pRNGseqs_mt[[i]])
    #repeated_vals <- sum(c(NA, sorted_vector_asc[-1] == sorted_vector_asc[-length(sorted_vector_asc)]), 
    #                     na.rm = TRUE)
    
    value_counts = table(pRNGseqs_mt[[i]])
    repeated_vals = sum(value_counts[value_counts > 1])

    print(glue("For sequence {i}, there are {repeated_vals} of {seq_length} ({repeated_vals*100/seq_length}%)."))
    
    
    stat = rank_von_neumann_test$statistic
        
    p = rank_von_neumann_test$p.value
    
    reject = (p < 0.05)
        
    results_df = rbind(results_df, data.frame(test_statistic = c(stat), 
                                              p_value = c(p), 
                                              reject_H0 = c(reject)))
  }
  return(results_df)
}