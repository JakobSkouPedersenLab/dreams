# Testing mutation calling on mutated data

source("~/GenomeDK/ctDNA/read_level_prediction/r_source_functions/mutation_call_functions.R")

# Parameters --------------------------------------------------------------

tf <- 0.02
n <- 10000

max_error <- 1e-2

# -------------------------------------------------------------------------

# Simulate error rate and read count for ref and mut
err_ref_to_mut <- runif(n, 0, max_error)
err_mut_to_ref <- runif(n, 0, max_error)

mut_reads <- rbinom(1, n, tf)
ref_reads <- n - mut_reads

# Simulate observation (introduce errors)
observed_ref_pos <- rbinom(ref_reads, 1, err_ref_to_mut[1:ref_reads])
observed_mut_pos <- 1 - rbinom(mut_reads, 1, err_mut_to_ref[-(1:ref_reads)])

sample_mutations <- c(observed_ref_pos, observed_mut_pos)

# Do simple test
EM_test(sample_mutations, err_ref_to_mut, err_mut_to_ref, use_warp_speed = TRUE)


# Simulating EM algorithm results -----------------------------------------

em_sim <- function(n, max_err, tf, use_warp_speed) {
  err_ref_to_mut <- runif(n, 0, max_err)
  err_mut_to_ref <- runif(n, 0, max_err)

  mut_reads <- rbinom(1, n, tf / 2)
  ref_reads <- n - mut_reads

  observed_ref_pos <- rbinom(ref_reads, 1, err_ref_to_mut[1:ref_reads])
  observed_mut_pos <- 1 - rbinom(mut_reads, 1, err_mut_to_ref[-(1:ref_reads)])

  sample_mutations <- c(observed_ref_pos, observed_mut_pos)

  em_res <- EM_test(
    sample_mutations = sample_mutations,
    err_ref_to_mut = err_ref_to_mut,
    err_mut_to_ref = err_mut_to_ref,
    use_warp_speed = use_warp_speed
  ) %>% data.frame()

  return(em_res)
}

# Simulate and estimate tf using fast algorithm
tf <- 0.002
max_err = 0.001
fast_em_sim_res <- bind_rows(replicate(500, em_sim(n = 1000, max_err = max_err, tf = tf, use_warp_speed = TRUE), simplify = FALSE))

hist(fast_em_sim_res$tumor_freq, breaks = 50, freq = FALSE)
lines(density(fast_em_sim_res$tumor_freq, from = 0), lwd = 2)
abline(v = tf, col = "red", lwd = 2)
abline(v = mean(fast_em_sim_res$tumor_freq), col = "green", lwd = 2, lty = 2)

# Simulate and estimate tf using slow algorithm
slow_em_sim_res <- bind_rows(replicate(500, em_sim(n = 1000, max_err = max_err, tf = tf, use_warp_speed = FALSE), simplify = FALSE))

hist(slow_em_sim_res$tumor_freq, breaks = 50, freq = FALSE)
lines(density(slow_em_sim_res$tumor_freq, from = 0), lwd = 2)
abline(v = tf, col = "red", lwd = 2)
abline(v = mean(slow_em_sim_res$tumor_freq), col = "green", lwd = 2, lty = 2)



# Compare fast and slow algorithm via simulation --------------------------

em_compare_sim <- function(n, max_err, tf = 0.001) {
  err_ref_to_mut <- runif(n, 0, 1e-2)
  err_mut_to_ref <- runif(n, 0, 1e-2)

  mut_reads <- rbinom(1, n, tf)
  ref_reads <- n - mut_reads

  observed_ref_pos <- rbinom(ref_reads, 1, err_ref_to_mut[1:ref_reads])
  observed_mut_pos <- 1 - rbinom(mut_reads, 1, err_mut_to_ref[-(1:ref_reads)])

  sample_mutations <- c(observed_ref_pos, observed_mut_pos)

  em_res_slow <- EM_test(sample_mutations, err_ref_to_mut, err_mut_to_ref, use_warp_speed = FALSE) %>% data.frame()
  em_res_fast <- EM_test(sample_mutations, err_ref_to_mut, err_mut_to_ref, use_warp_speed = TRUE) %>% data.frame()

  colnames(em_res_fast) <- paste0(colnames(em_res_fast), "_fast")
  colnames(em_res_slow) <- paste0(colnames(em_res_slow), "_slow")

  em_res <- bind_cols(em_res_fast, em_res_slow)

  return(em_res)
}


em_compare_sim_res <- bind_rows(replicate(50, em_compare_sim(n = 10000, max_err = 1e-2), simplify = FALSE))

par(mfrow = c(3, 2))

plot(em_compare_sim_res$ll_A_fast, em_compare_sim_res$ll_A_slow)
abline(a = 0, b = 1)

plot(em_compare_sim_res$tumor_freq_fast, em_compare_sim_res$tumor_freq_slow)
abline(a = 0, b = 1)

plot(em_compare_sim_res$p_val_fast, em_compare_sim_res$p_val_slow, log = "xy")
abline(a = 0, b = 1)

plot(em_compare_sim_res$Q_fast + 0.01, em_compare_sim_res$Q_slow + 0.01, log = "xy")
abline(a = 0, b = 1)

plot(em_compare_sim_res$fct_evals_fast, em_compare_sim_res$fct_evals_slow, log = "xy")
abline(a = 0, b = 1)

plot(em_compare_sim_res$EM_steps_fast, em_compare_sim_res$EM_steps_slow)
abline(a = 0, b = 1)

par(mfrow = c(1, 1))

