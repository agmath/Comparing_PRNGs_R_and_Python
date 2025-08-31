library(dplyr)


binned_chisq_tests_mt_seqs <- function(num_bins = 15, num_seqs = 10, seq_length = 2^14 - 1, max_rand = 1e5 - 1, seed = 1234){
  source("CurrentPaperCode/generate_multiple_rand_seqs_mt.R")
  
  #Generate sequences using Mersenne Twister
  pRNGseqs_mt <- generate_mt_rands(num_seqs, seq_length, max_rand, seed)
  
  #Create bins
  bin_breaks <- seq(0, max_rand, length.out = num_bins + 1)
  
  #Initialize vector to store results
  chi_sq_results <- c()
  
  #Loop through each generated sequence
  for(i in 1:ncol(pRNGseqs_mt)){
    
    #Calculate observed counts
    obs_counts <- table(cut(pRNGseqs_mt[[i]], breaks = bin_breaks, include.lowest = TRUE))
    
    #Calculate expected counts, assuming uniformity
    exp_counts <- rep(sum(obs_counts)/num_bins, num_bins)
    
    #Compute Chi-Square statistic
    chi_sq_results[i] <- sum((obs_counts - exp_counts)^2/exp_counts)
  }
  
  #Degrees of Freedom
  df <- num_bins - 1
  
  #Critical Value at 95% significance
  critical_val <- qchisq(0.95, df)
  
  #Combine Results
  
  results_df <- data.frame(
    Column = colnames(pRNGseqs_mt),
    Chi_Square_Statistic = chi_sq_results,
    DF = df,
    Critical_Value_95 = critical_val
  ) %>%
    mutate(
      Reject_H0 = (Chi_Square_Statistic > Critical_Value_95)
    )
  
  return(results_df)
}
