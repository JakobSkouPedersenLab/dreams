test_that("Small example", {
  three_mutations_df <-
    data.frame(
      CHROM = c("chr1", "chr2", "chr3"),
      POS = c(10, 13, 17),
      REF = "A",
      ALT = "T"
    )

  # Reads on Mut 1,2 - No reads on Mut 3 - And random read
  some_reads_df <-
    data.frame(
      chr =
        c(
          "chr1", "chr1",
          "chr2", "chr2", "chr2",
          "chr5"
        ),
      genomic_pos =
        c(
          10, 10,
          13, 13, 13,
          19
        ),
      ref = "A",
      obs =
        c(
          "A", "T",
          "A", "A", "A",
          "A"
        )
    )

  model_path <- system.file("extdata", "model_test.h5", package = "dreams")
  model <- keras::load_model_hdf5(model_path)

  res <- call_mutations(mutations_df = three_mutations_df, all_reads = some_reads_df, model = model, beta = 0.01, alpha = 0.05)

  # Cancer info
  expect_equal(res$cancer_info$mutations_tested, 3)
  expect_equal(res$cancer_info$total_coverage, 5)
  expect_true(res$cancer_info$EM_converged)

  # Mutations info
  expect_equal(nrow(res$mutation_info), 3)
  expect_equal(res$mutation_info$count, c(1, 0, 0))
  expect_equal(res$mutation_info$coverage, c(2, 3, 0))
})


test_that("Snapshot", {
  three_mutations_df <-
    data.frame(
      CHROM = c("chr1", "chr2", "chr3"),
      POS = c(10, 13, 17),
      REF = "A",
      ALT = "T"
    )

  # Reads on Mut 1,2 - No reads on Mut 3 - And random read
  cov_1 = 17
  count_1 = 1
  cov_2 = 43
  count_2 = 3

  some_reads_df <-
    data.frame(
      chr =
        c(
          rep("chr1", cov_1),
          rep("chr2", cov_2)
        ),
      genomic_pos =
        c(
          rep(10, cov_1),
          rep(13, cov_2)
        ),
      ref = "A",
      obs =
        c(
          rep("A", cov_1 - count_1), rep("T", count_1),
          rep("A", cov_2 - count_2), rep("T", count_2)
        )
    )

  model_path <- system.file("extdata", "model_test.h5", package = "dreams")
  model <- keras::load_model_hdf5(model_path)

  res <- call_mutations(mutations_df = three_mutations_df, all_reads = some_reads_df, model = model, beta = 0.01, alpha = 0.05)

  expect_snapshot(res)
})






# Old ---------------------------------------------------------------------



# Testing mutation calling on mutated data

# Parameters --------------------------------------------------------------

tf <- 0.01
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
em_test_mutation(sample_mutations, err_ref_to_mut, err_mut_to_ref, use_warp_speed = TRUE)


# Simulating EM algorithm results -----------------------------------------

em_sim <- function(n, max_err, tf, use_warp_speed) {
  err_ref_to_mut <- runif(n, 0, max_err)
  err_mut_to_ref <- runif(n, 0, max_err)

  mut_reads <- rbinom(1, n, tf / 2)
  ref_reads <- n - mut_reads

  observed_ref_pos <- rbinom(ref_reads, 1, err_ref_to_mut[1:ref_reads])
  observed_mut_pos <- 1 - rbinom(mut_reads, 1, err_mut_to_ref[-(1:ref_reads)])

  sample_mutations <- c(observed_ref_pos, observed_mut_pos)

  em_res <- em_test_mutation(
    sample_mutations = sample_mutations,
    err_ref_to_mut = err_ref_to_mut,
    err_mut_to_ref = err_mut_to_ref,
    use_warp_speed = use_warp_speed
  ) %>% data.frame()

  return(em_res)
}

# Simulate and estimate tf using fast algorithm
tf <- 0.002
max_err <- 0.001
fast_em_sim_res <- bind_rows(replicate(1, em_sim(n = 1000, max_err = max_err, tf = tf, use_warp_speed = TRUE), simplify = FALSE))

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

  em_res_slow <- em_test_mutation(sample_mutations, err_ref_to_mut, err_mut_to_ref, use_warp_speed = FALSE) %>% data.frame()
  em_res_fast <- em_test_mutation(sample_mutations, err_ref_to_mut, err_mut_to_ref, use_warp_speed = TRUE) %>% data.frame()

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
