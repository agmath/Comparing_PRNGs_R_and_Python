library(dgof)

##############################################################################
#  This function generates the desired number of sequences in R using the 
#  given parameters.  The function that performs a 1 sample Kolmogorov–Smirnov Test
#  to determine if each sequence follows a discrete uniform distribution on the
#  given range of values used to generate the sequences.  The function
#  returns a dataframe containing:
#  * KS Test Statistic
#  * P-Value
#  * Whether or not the null Hypothesis was rejected
#    True:   The null was rejected  (H0:  the sequence follows a uniform distribution)
#    False:  The null was not rejected
###############################################################################3   
compare_dists_mt_uniform_ks_test <- function(num_seqs = 10, 
                                       seq_length = 2^14 - 1, 
                                       min_rand = 0,
                                       max_rand = 1e5 - 1, 
                                       seed = 1234) {

    
  source("generate_multiple_rand_seqs_mt.R")
  
  #Generate sequences using Mersenne Twister
  pRNGseqs_mt <- generate_mt_rands(num_seqs, 
                                   seq_length, 
                                   min_rand, 
                                   max_rand, 
                                   seed)
  
  
  #Initialize vectors to store results
  ks_test_statistics <- c()
  ks_test_pvalues <- c()
  
  #Loop through each generated sequence
  for(i in 1:ncol(pRNGseqs_mt)){
  
    #Here we run a 1 sample ks test comparing our sequence to a uniform discrete 
    #distribution.  We followed the example code from Section 5 here.
    #https://journal.r-project.org/articles/RJ-2011-016/#citation
    
    #ecdf(min_rand:max_rand) is a step function that jumps at each integer,
    #increasing by 1/n each time, so it gives the empirical CDF of a discrete uniform 
    #distribution on the given range.
    test_results = ks.test(pRNGseqs_mt[[i]],  ecdf(min_rand:max_rand))
    
    #Get the D KS statistic and p-value
    D = test_results$statistic
    p_value = test_results$p.value
    
    #Add these values to our vectors
    ks_test_statistics[i] <- D 
    ks_test_pvalues[i] <- p_value
  
  }
  
  
  # Combine Results
  results_df <- data.frame(
    Column = colnames(pRNGseqs_mt),
    KS_Statistic = ks_test_statistics,
    KS_p_value = ks_test_pvalues
  ) %>%
    mutate(
      Decision = ifelse(KS_p_value < 0.05, "Reject H0", "Do Not Reject")
    )
  
  
  return(results_df)
}
