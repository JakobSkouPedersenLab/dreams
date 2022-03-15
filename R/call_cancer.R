
# Template for empty output
empty_cancer_info <- function(mutations_df, X_list) {
  data.frame(
    tf_est = 0, tf_min = 0, tf_max = 1,
    r_est = 0, r_min = 0, r_max = 1,
    Q_val = 0,
    ll_A = NA,
    ll_0 = NA,
    mutations_tested = nrow(mutations_df),
    est_mutations_present = 0,
    total_coverage = sum(sapply(X_list, length)),
    total_count = 0,
    EM_converged = NA,
    EM_steps = NA,
    fpeval = NA,
    objfeval = NA,
    p_val = 1,
    cancer_detected = FALSE
  )
}

prepare_em_input <- function(mutations_df, reads_df, model, beta) {
  X_list <- list()
  error_ref_to_mut_list <- list()
  error_mut_to_ref_list <- list()

  if (nrow(mutations_df) >= 1) {
    for (i in 1:nrow(mutations_df)) {
      chr <- mutations_df[i, "chr"]
      genomic_pos <- mutations_df[i, "genomic_pos"]
      ref <- mutations_df[i, "ref"]
      alt <- mutations_df[i, "alt"]

      mut_reads_ref_alt <-
        reads_df %>%
        filter(
          # Filter reads positions to mutation
          .data$chr == !!chr, .data$genomic_pos == !!genomic_pos,
          # Filter mutation type (remove N and "other" than ref/alt alleles)
          .data$obs == !!ref | .data$obs == !!alt
        )

      # Add to list
      X <- mut_reads_ref_alt$obs == alt
      X_list <- append(X_list, list(X))

      # Error rates
      error_ref_df <- predict_error_rates(mut_reads_ref_alt, model, beta)
      error_mut_df <- predict_error_rates(mut_reads_ref_alt %>% mutate(ref = !!alt), model, beta)

      # Pick relevant error rates for ref and alt
      error_ref_to_ref_raw <- error_ref_df[[paste0(ref, "_corrected")]]
      error_ref_to_mut_raw <- error_ref_df[[paste0(alt, "_corrected")]]

      error_mut_to_ref_raw <- error_mut_df[[paste0(ref, "_corrected")]]
      error_mut_to_mut_raw <- error_mut_df[[paste0(alt, "_corrected")]]

      # Normalize
      error_ref_to_mut <- error_ref_to_mut_raw / (error_ref_to_ref_raw + error_ref_to_mut_raw)
      error_mut_to_ref <- error_mut_to_ref_raw / (error_mut_to_ref_raw + error_mut_to_mut_raw)

      # Add to list
      error_ref_to_mut_list <- append(error_ref_to_mut_list, list(error_ref_to_mut))
      error_mut_to_ref_list <- append(error_mut_to_ref_list, list(error_mut_to_ref))
    }
  }

  return(
    list(
      X_list = X_list,
      error_ref_to_mut_list = error_ref_to_mut_list,
      error_mut_to_ref_list = error_mut_to_ref_list
    )
  )
}

#' Title
#'
#' @param mutations_df
#' @param reads_df
#' @param beta
#' @param use_warp_speed
#' @param calculate_confidence_intervals
#' @param model
#'
#' @return
#' @export
call_cancer <- function(mutations_df, reads_df, model, beta, alpha = 0.05, use_warp_speed = TRUE, calculate_confidence_intervals = FALSE) {
  # Clean up mutations
  mutations_df <- mutations_df %>%
    select(
      "chr" = matches("chr|CHR|CHROM"),
      "genomic_pos" = matches("pos|POS"),
      "ref" = matches("ref|REF"),
      "alt" = matches("alt|ALT|obs|OBS")
    )

  # Stop if mutations do not have the expected columns
  mutations_expected_columns <- c("chr", "genomic_pos", "ref", "alt")
  if (!all(mutations_expected_columns %in% colnames(mutations_df))) {
    stop("mutations_df should have the columns ['chr', genomic_pos', 'ref, 'alt']")
  }

  # Stop if reads do not have the expected columns
  reads_expected_columns <- c("chr", "genomic_pos", "ref", "obs")
  if (!all(reads_expected_columns %in% colnames(reads_df))) {
    stop("reads_df should have the columns ['chr', genomic_pos', 'ref, 'obs']")
  }

  # Prepare inputs for algorithm
  em_input <- prepare_em_input(mutations_df = mutations_df, reads_df = reads_df, model = model, beta = beta)
  X_list <- em_input$X_list
  error_ref_to_mut_list <- em_input$error_ref_to_mut_list
  error_mut_to_ref_list <- em_input$error_mut_to_ref_list

  # If no mutations return empty result
  if (nrow(mutations_df) == 0) {
    return(
      list(
        cancer_info = empty_cancer_info(mutations_df, em_input),
        mutation_info = data.frame()
      )
    )
  }

  # If no positions have relevant coverage (ref or alt) or count of mu is 0 for every candidate
  no_mutations_has_relevant_count <- all(sapply(X_list, sum) == 0)
  if (no_mutations_has_relevant_count) {
    # Return output with no mutations
    no_count_mutation_info <- data.frame(
      mutations_df,
      P_mut_is_present = NA,
      obs_freq = 0,
      exp_count = sapply(error_ref_to_mut_list, sum),
      count = sapply(X_list, sum),
      coverage = sapply(X_list, length)
    )

    return(
      list(
        cancer_info = empty_cancer_info(mutations_df, X_list),
        mutation_info = no_count_mutation_info
      )
    )
  }

  # EM algorithm
  if (use_warp_speed) {
    EM_res <- run_EM(
      X_list = X_list,
      error_ref_to_mut_list = error_ref_to_mut_list,
      error_mut_to_ref_list = error_mut_to_ref_list
    )
  } else {
    EM_res <- run_EM_full(
      X_list = X_list,
      error_ref_to_mut_list = error_ref_to_mut_list,
      error_mut_to_ref_list = error_mut_to_ref_list
    )
  }


  # Confidence intervals
  if (calculate_confidence_intervals) {
    tf_CI <- get_tf_CI(
      X_list = X_list,
      error_mut_to_ref_list = error_mut_to_ref_list,
      error_ref_to_mut_list = error_ref_to_mut_list,
      r_est = EM_res$r, tf_est = EM_res$tf,
      alpha = alpha
    )
    r_CI <- get_r_CI(
      X_list = X_list,
      error_mut_to_ref_list = error_mut_to_ref_list,
      error_ref_to_mut_list = error_ref_to_mut_list,
      r_est = EM_res$r, tf_est = EM_res$tf,
      alpha = alpha
    )
  } else {
    tf_CI <- list(tf_min = NA, tf_max = NA)
    r_CI <- list(r_min = NA, r_max = NA)
  }

  # Test significance
  ll_0 <- log_likelihood(
    X_list = X_list,
    error_mut_to_ref_list = error_mut_to_ref_list,
    error_ref_to_mut_list = error_ref_to_mut_list,
    r = 0,
    tf = 0
  )

  ll_A <- log_likelihood(
    X_list = X_list,
    error_mut_to_ref_list = error_mut_to_ref_list,
    error_ref_to_mut_list = error_ref_to_mut_list,
    r = EM_res$r,
    tf = EM_res$tf
  )

  Q_val <- -2 * (ll_0 - ll_A)
  p_val <- stats::pchisq(Q_val, 2, lower.tail = FALSE)

  # Collect cancer information
  cancer_info <-
    data.frame(
      tf_est = EM_res$tf, tf_min = tf_CI$tf_min, tf_max = tf_CI$tf_max,
      r_est = EM_res$r, r_min = r_CI$r_min, r_max = r_CI$r_max,
      Q_val = Q_val,
      ll_A = ll_A,
      ll_0 = ll_0,
      mutations_tested = length(X_list),
      est_mutations_present = sum(EM_res$P_mut_is_present),
      total_coverage = sum(sapply(X_list, length)),
      total_count = sum(sapply(X_list, sum)),
      EM_converged = EM_res$EM_converged,
      EM_steps = EM_res$EM_steps,
      fpeval = EM_res$fpeval,
      objfeval = EM_res$objfeval,
      p_val = p_val,
      cancer_detected = p_val <= alpha
    )

  # Collect mutation information
  mutation_info <-
    dplyr::bind_cols(
      mutations_df,
      data.frame(
        P_mut_is_present = EM_res$P_mut_is_present,
        obs_freq = sapply(X_list, mean),
        exp_count = sapply(error_ref_to_mut_list, sum),
        count = sapply(X_list, sum),
        coverage = sapply(X_list, length)
      )
    )

  return(
    list(
      cancer_info = cancer_info,
      mutation_info = mutation_info
    )
  )

  res <- list(
    cancer_info = EM_res$cancer_info,
    mutation_info =
    )

  return(res)
}
