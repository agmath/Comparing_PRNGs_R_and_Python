#Find and load functionality from individual scripts
file_names <- dir("AdamCode")
file_names <- setdiff(file_names, 
                      c("generate_test_and_summarize.R"))

for(name in file_names){
  source(paste0("AdamCode/", name))
}

#Generate random numbers
max_rand <- 1e5 - 1
seq_length <- 2^14 - 1
num_seqs <- 10

pRNGseqs <- generate_mt_rands(num_seqs, seq_length, max_rand)

#Run tests and build tables
boot_chisq_results <- pRNGseqs %>%
  bootstrap_chi_sq_collector()

chisq_results <- pRNGseqs %>%
  chi_sq_collector()

mean_results <- pRNGseqs %>%
  calc_mean()

entropy_results <- pRNGseqs %>%
  entropy_collector()

serial_corr_results <- pRNGseqs %>%
  serial_corr_collector()

mc_sim_results <- pRNGseqs %>%
  mc_sim_collector()

final_tbl <- entropy_results %>%
  left_join(mean_results, by = "sequence") %>%
  left_join(boot_chisq_results, by = "sequence") %>%
  left_join(chisq_results, by = "sequence") %>%
  left_join(mc_sim_results, by = "sequence") %>%
  left_join(serial_corr_results, by = "sequence") %>%
  rename(
    boot_chisq_stat = statistic.x,
    boot_chisq_pval = p.value.x,
    chisq_stat = statistic.y,
    chisq_pval = p.value.y
  )

final_tbl
