library(dplyr)


binned_chi_square_tests_mt <- function(num_seqs = 10, 
                                       seq_length = 2^14 - 1, 
                                       min_rand = 0,
                                       max_rand = 1e5 - 1, 
                                       seed = 1234, 
                                       num_bins = 15) {
  
  source("generate_multiple_rand_seqs_mt.R")
  
  
  #Generate sequences using Mersenne Twister
  pRNGseqs_mt <- generate_mt_rands(num_seqs, 
                                   seq_length, 
                                   min_rand, 
                                   max_rand, 
                                   seed)
  
  #Create bins
  bin_breaks <- seq(min_rand, max_rand, length.out = num_bins + 1)
  
  
  #Initialize vector to store results
  chi_sq_test_statistics <- c()
  chi_sq_p_values <- c()
  
  #Degrees of Freedom for the test
  df <- num_bins - 1
  
  #Critical Value at 95% significance
  critical_val <- qchisq(0.95, df)
  
  
  #Loop through each generated sequence
  for(i in 1:ncol(pRNGseqs_mt)){

    
    
    #Calculate observed counts
    obs_counts <- table(cut(pRNGseqs_mt[[i]], breaks = bin_breaks, include.lowest = TRUE))
    
    #Calculate expected counts, assuming uniformity
    exp_counts <- rep(sum(obs_counts)/num_bins, num_bins)
    
    #Compute Chi-Square statistic
    chi_sq_test_statistics[i] <- sum((obs_counts - exp_counts)^2/exp_counts)
    
    #pvalue for the results
    chi_sq_p_values[i] <- pchisq(chi_sq_test_statistics[i], df = df, lower.tail = FALSE)
    
  }
  

  

  #Combine Results
  results_df <- data.frame(
    Column = colnames(pRNGseqs_mt),
    Chi_Square_Statistic = chi_sq_test_statistics,
    Chi_Square_PValue = chi_sq_p_values,
    DF = df,
    Critical_Value_95 = critical_val
  ) %>%
    mutate(
      Reject_H0 = (Chi_Square_Statistic > Critical_Value_95)
    )

  return(results_df)
}
