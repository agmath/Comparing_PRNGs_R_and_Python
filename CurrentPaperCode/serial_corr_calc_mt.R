##################################
### Helper function for        ###
### analyzing autocorrelation  ###
### with a given lag amongst   ###
### terms in a randomly        ###
### generated sequence.        ###
##################################

library(dplyr)

serial_corr_calc_mt <- function(num_lag = 1, 
                                     num_seqs = 10, 
                                     seq_length = 2^14 - 1, 
                                     min_rand = 0,
                                     max_rand = 1e5 - 1, 
                                     seed = 1234){
  
  source("generate_multiple_rand_seqs_mt.R")
  
  #Generate sequences using Mersenne Twister
  pRNGseqs_mt <- generate_mt_rands(num_seqs, 
                                   seq_length, 
                                   min_rand, 
                                   max_rand, 
                                   seed)
  
  results_df <- data.frame()
  
  for(i in 1:num_seqs){
    curr_seq <- pRNGseqs_mt[[i]]
    lagged <- lag(pRNGseqs_mt[[i]], n = num_lag)
    
    #Generate the correlation coefficient of the current and lagged sequence. 
    #Due to the shift, there will be NaNs in the lagged sequence. 
    #Ignore any rows with an NAN in the lagged sequence
    serial_cor_test <- cor.test(curr_seq, 
                                lagged, 
                                method = "pearson", 
                                na.action = na.omit)

    #Critical Value at 95% significance
    #Since it is a 2 sided test, we use p = .05/2
    t_2sided_critical_val_95 <- qt(p = .025, df = serial_cor_test$parameter)
    test_statistic = serial_cor_test$statistic
      
    #Create a df of the results
    results_df <- results_df %>%
      bind_rows(
        tibble(
          sequence = names(pRNGseqs_mt)[i],
          serial_correlation = serial_cor_test$estimate,
          df = serial_cor_test$parameter, 
          critical_value = t_2sided_critical_val_95,  #This cv will be negative
          test_statistic = test_statistic, 
          p_value = serial_cor_test$p.value,
          #reject_H0 = (test_statistic >t_2sided_critical_val_95)
          reject_H0 = (test_statistic < t_2sided_critical_val_95 | 
                         test_statistic > -1 * t_2sided_critical_val_95)
        )
      )
  }
  return(results_df)
}