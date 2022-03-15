if (file.exists("/faststorage/project/PolyA/BACKUP/ctDNA_analysis/read_level_prediction/r_source_functions/bam_handling_functions_wf0_4.R")) {
  suppressMessages(source("/faststorage/project/PolyA/BACKUP/ctDNA_analysis/read_level_prediction/r_source_functions/bam_handling_functions_wf0_4.R"))
  suppressMessages(source("/faststorage/project/PolyA/BACKUP/ctDNA_analysis/read_level_prediction/r_source_functions/data_preprocessing.R"))
} else {
  suppressMessages(source("~/GenomeDK/ctDNA/read_level_prediction/r_source_functions/bam_handling_functions_wf0_4.R"))
  suppressMessages(source("~/GenomeDK/ctDNA/read_level_prediction/r_source_functions/data_preprocessing.R"))
}
suppressMessages(library(bedr))
suppressMessages(library(keras))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(SQUAREM))

predict_error_rates <- function(reads, model, beta) {
  prediction <- model %>%
    predict(reads) %>%
    data.frame() %>%
    rename(
      A = X1,
      T = X2,
      C = X3,
      G = X4
    )

  prediction$ref <- reads$ref

  corrected_errors <- prediction %>%
    correct_errors_predictions(beta = beta) %>%
    select(-ref)

  predicted_errors <- bind_cols(reads, corrected_errors)

  return(predicted_errors)
}

log_lik <- function(tf, sample_mutations, err_ref_to_mut, err_mut_to_ref) {
  log_lik_mut <- log((1 - err_mut_to_ref) * tf / 2 + err_ref_to_mut * (1 - tf / 2))
  log_lik_ref <- log((1 - err_ref_to_mut) * (1 - tf / 2) + err_mut_to_ref * tf / 2)

  return(sum(log_lik_mut[sample_mutations == 1]) + sum(log_lik_ref[sample_mutations == 0]))
}

update_P_Y1 <- function(x, tf, err_ref_to_mut, err_mut_to_ref) {
  # Get scaled probabilities for fragment state(Y)
  Q_Y0 <- ifelse(x == 0,
                 (1 - tf / 2) * (1 - err_ref_to_mut),
                 (1 - tf / 2) * err_ref_to_mut
  )

  Q_Y1 <- ifelse(x == 0,
                 tf / 2 * err_mut_to_ref,
                 tf / 2 * (1 - err_mut_to_ref)
  )

  # Normalize probabilities
  P_Y1 <- Q_Y1 / (Q_Y0 + Q_Y1)

  return(P_Y1)
}

get_tumor_freq <- function(x, err_ref_to_mut, err_mut_to_ref, use_warp_speed, max_it = 5000) {

  # Get starting guess for tf
  expected_signal <- sum(err_ref_to_mut) / length(x)
  excess_signal <- mean(x) - expected_signal

  tf <- max(1e-8, excess_signal * 2)

  # Improve guess if value close to 0 is better
  ll_0 <- log_lik(0, x, err_ref_to_mut = err_ref_to_mut, err_mut_to_ref = err_mut_to_ref)
  ll_eps <- log_lik(1e-8, x, err_ref_to_mut = err_ref_to_mut, err_mut_to_ref = err_mut_to_ref)
  ll_start <- log_lik(tf, x, err_ref_to_mut = err_ref_to_mut, err_mut_to_ref = err_mut_to_ref)

  if (ll_eps >= ll_start) {
    tf <- 1e-8

    if (ll_0 >= ll_eps) {
      tf <- 0
    }
  }

  if (use_warp_speed) {
    em_update <- function(par, sample_mutations, err_ref_to_mut, err_mut_to_ref) {
      tf <- par[1]

      # Update P_Y1
      P_Y1 <- update_P_Y1(x = sample_mutations, tf = tf, err_ref_to_mut = err_ref_to_mut, err_mut_to_ref = err_mut_to_ref)

      # Update tf
      tf_new <- mean(P_Y1) * 2

      return(c(tf_new, P_Y1))
    }

    em_objective <- function(par, sample_mutations, err_ref_to_mut, err_mut_to_ref) {
      tf <- par[1]

      ll <- log_lik(tf, sample_mutations, err_ref_to_mut = err_ref_to_mut, err_mut_to_ref = err_mut_to_ref)

      return(ll)
    }

    # Start speed up version of EM algorithm

    squarem_res <- squarem(
      par = c(tf, rep(0, length(err_ref_to_mut))),
      fixptfn = em_update,
      objfn = em_objective,
      sample_mutations = x,
      err_ref_to_mut = err_ref_to_mut,
      err_mut_to_ref = err_mut_to_ref,
      control = list(
        minimize = FALSE
      )
    )

    EM_converged <- squarem_res$convergence
    EM_steps <- squarem_res$iter
    fct_evals <- squarem_res$fpevals + squarem_res$objfevals

    tf <- squarem_res$par[1]
  } else {
    # Start regular EM algorithm
    EM_converged <- FALSE
    EM_steps <- 0

    for (i in 1:max_it) {
      EM_steps <- EM_steps + 1
      tf_old <- tf

      # Update z
      z <- update_P_Y1(x = x, tf = tf_old, err_ref_to_mut = err_ref_to_mut, err_mut_to_ref = err_mut_to_ref)

      # Update tf
      tf <- mean(z) * 2

      if (abs(tf - tf_old) / (tf_old + 1e-8) < 1e-8 & EM_steps >= 5) {
        EM_converged <- TRUE
        break
      }
    }

    fct_evals <- EM_steps
  }


  # Final check of tf value
  ll_final <- log_lik(tf, x, err_ref_to_mut = err_ref_to_mut, err_mut_to_ref = err_mut_to_ref)
  if (ll_0 >= ll_final) {
    tf <- 0
  }

  res <- list(
    tf = tf,
    EM_converged = EM_converged,
    EM_steps = EM_steps,
    fct_evals = fct_evals
  )

  return(res)
}


EM_test <- function(sample_mutations, err_ref_to_mut, err_mut_to_ref, use_warp_speed = TRUE) {
  get_tumor_freq_res <-
    get_tumor_freq(
      sample_mutations,
      err_ref_to_mut = err_ref_to_mut,
      err_mut_to_ref = err_mut_to_ref,
      use_warp_speed = use_warp_speed
    )

  tumor_freq <- get_tumor_freq_res$tf
  EM_converged <- get_tumor_freq_res$EM_converged
  EM_steps <- get_tumor_freq_res$EM_steps
  fct_evals <- get_tumor_freq_res$fct_evals

  ll_0 <- log_lik(0, sample_mutations, err_ref_to_mut = err_ref_to_mut, err_mut_to_ref = err_mut_to_ref)
  ll_A <- log_lik(tumor_freq, sample_mutations, err_ref_to_mut = err_ref_to_mut, err_mut_to_ref = err_mut_to_ref)
  Q_val <- -2 * (ll_0 - ll_A)
  p_val <- pchisq(Q_val, 1, lower.tail = FALSE)

  return(
    list(
      tumor_freq = tumor_freq,
      Q = Q_val,
      p_val = p_val,
      ll_A = ll_A,
      ll_0 = ll_0,
      EM_converged = EM_converged,
      EM_steps = EM_steps,
      fct_evals = fct_evals
    )
  )
}



call_mutations <- function(mutations, all_reads, beta) {
  res_df <- NULL

  if (nrow(mutations) == 0) {
    empty_row <- data.frame(
      chr = NA, pos = NA, ref = NA, alt = NA,
      tumor_freq = 0, p_val = 1,
      EM_converged = FALSE, EM_steps = 0, fct_evals = 0,
      count = 0, expected_error_count = 0,
      coverage = 0, full_coverage = 0,
      Q = 0, ll_A = 0, ll_0 = 0
    )
    return(empty_row)
  } else {
    for (i in 1:nrow(mutations)) {
      cat("mutation", i, "\n")
      chr <- mutations$CHROM[i]
      genomic_pos <- as.numeric(as.character(mutations$POS[i]))
      ref <- mutations$REF[i]
      alt <- mutations$ALT[i]

      reads <- all_reads %>% filter(chr == !!chr, genomic_pos == !!genomic_pos)
      full_coverage <- nrow(reads)
      print ("COVERAGE")
      cat(chr, genomic_pos, "\n")

      print (full_coverage)


      # if no reads cover position, skip to next mutation
      if (nrow(reads) == 0) {
        new_row <- data.frame(
          chr = chr, pos = genomic_pos, ref = ref, alt = alt,
          tumor_freq = 0, p_val = 1,
          EM_converged = FALSE, EM_steps = 0, fct_evals = 0,
          count = 0, expected_error_count = 0,
          coverage = 0, full_coverage = full_coverage,
          Q = 0, ll_A = 0, ll_0 = 0
        )
        res_df <- rbind(res_df, new_row)
        next
      }



      ###### Filer reads to only have reads that obser
      reads_ref_alt <- reads %>% filter(obs == !!ref | obs == !!alt)
      n_ref <- sum(reads_ref_alt$obs == ref)
      n_mut <- sum(reads_ref_alt$obs == alt)


      # No reads reads supporting reference or mutation
      if (n_ref == 0 | n_mut == 0) {
        new_row <- data.frame(
          chr = chr, pos = genomic_pos, ref = ref, alt = alt,
          tumor_freq = 0, p_val = 1,
          EM_converged = FALSE, EM_steps = 0, fct_evals = 0,
          count = n_mut, expected_error_count = NA,
          coverage = nrow(reads_ref_alt), full_coverage = full_coverage,
          Q = 0, ll_A = 0, ll_0 = 0
        )
        res_df <- rbind(res_df, new_row)
        next
      }


      ##### PREDICT ERROR RATES

      errors_ref_df <- predict_error_rates(reads_ref_alt, model, beta)
      errors_mut_df <- predict_error_rates(reads_ref_alt %>% mutate(ref = !!alt), model, beta)

      ##### EXTRACT CORRECTED ERRORS

      error_ref_to_ref <- errors_ref_df[[paste0(ref, "_corrected")]]
      error_ref_to_mut <- errors_ref_df[[paste0(alt, "_corrected")]]

      error_ref_to_mut_norm <- error_ref_to_mut / (error_ref_to_ref + error_ref_to_mut)

      error_mut_to_ref <- errors_mut_df[[paste0(ref, "_corrected")]]
      error_mut_to_mut <- errors_mut_df[[paste0(alt, "_corrected")]]

      error_mut_to_ref_norm <- error_mut_to_ref / (error_mut_to_ref + error_mut_to_mut)


      ##### EM Log-lik
      sample_mutations <- reads_ref_alt$obs == alt

      EM_res <- EM_test(sample_mutations, error_ref_to_mut_norm, error_mut_to_ref_norm)

      p_val <- EM_res$p_val
      tumor_freq <- EM_res$tumor_freq
      Q <- EM_res$Q
      ll_A <- EM_res$ll_A
      ll_0 <- EM_res$ll_0
      EM_converged <- EM_res$EM_converged
      EM_steps <- EM_res$EM_steps
      fct_evals <- EM_res$fct_evals

      new_row <- data.frame(
        chr = chr, pos = genomic_pos, ref = ref, alt = alt,
        tumor_freq = tumor_freq, p_val = p_val,
        EM_converged = EM_converged, EM_steps = EM_steps, fct_evals = fct_evals,
        count = n_mut, expected_error_count = sum(error_ref_to_mut_norm),
        coverage = nrow(reads_ref_alt), full_coverage = full_coverage,
        Q = Q, ll_A = ll_A, ll_0 = ll_0
      )
      res_df <- rbind(res_df, new_row)
    }
    return(res_df)
  }
}
