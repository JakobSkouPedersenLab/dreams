#' Call mutations from read positions
#'
#' @description This function evaluate the presence (calls) of individual mutations from a predefined list.
#' TODO: Add describeIn instead
#' @inheritParams call_cancer
#' @param calculate_confidence_interval Logical. Should confidence intervals be calculated? Default is FALSE.
#'
#' @return
#'
#'
#'
#'
#'
#' @seealso [call_cancer()], [train_dreams_model()]
#'
#' @export
call_mutations <- function(mutations_df, read_positions_df, model, beta, alpha = 0.05, use_turboem = TRUE, calculate_confidence_interval = FALSE) {
  # If no mutations return empty result
  if (nrow(mutations_df) == 0) {
    return(data.frame())
  }

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
  if (!all(reads_expected_columns %in% colnames(read_positions_df))) {
    stop("read_positions_df should have the columns ['chr', genomic_pos', 'ref, 'obs']")
  }

  # Prepare EM input
  em_input <- prepare_em_input(mutations_df = mutations_df, read_positions_df = read_positions_df, model = model, beta = beta)

  obs_is_mut_list <- em_input$obs_is_mut_list
  error_ref_to_mut_list <- em_input$error_ref_to_mut_list
  error_mut_to_ref_list <- em_input$error_mut_to_ref_list

  # Add full coverage to mutations
  full_coverage_df <- read_positions_df %>%
    count(chr, genomic_pos, name = "full_coverage")

  mutations_df <- mutations_df %>%
    left_join(full_coverage_df, by = c("chr", "genomic_pos"))

  # Call mutations from reads
  res_df <- NULL
  for (i in 1:nrow(mutations_df)) {
    # Mutation information
    mutation_line <- mutations_df[i, ]

    # Read information
    obs_is_mut <- obs_is_mut_list[[i]]
    error_ref_to_mut <- error_ref_to_mut_list[[i]]
    error_mut_to_ref <- error_mut_to_ref_list[[i]]

    # Run EM algorithm
    em_res <- get_tf_estimate_vc(
      obs_is_mut = obs_is_mut,
      error_ref_to_mut = error_ref_to_mut,
      error_mut_to_ref = error_mut_to_ref,
      use_turboem = use_turboem
    )

    # Confidence intervals
    if (calculate_confidence_interval) {
      tf_CI <- get_tf_CI(
        obs_is_mut_list = obs_is_mut_list,
        error_mut_to_ref_list = error_mut_to_ref_list,
        error_ref_to_mut_list = error_ref_to_mut_list,
        r_est = 1,
        tf_est = em_res$tf_est,
        alpha = alpha
      )
    } else {
      tf_CI <- list(tf_min = NA, tf_max = NA)
    }

    # Test mutation
    ll_0 <- log_likelihood(
      tf = 0,
      r = 1,
      obs_is_mut_list = obs_is_mut,
      error_ref_to_mut_list = error_ref_to_mut,
      error_mut_to_ref_list = error_mut_to_ref
    )
    ll_A <- log_likelihood(
      tf = em_res$tf_est,
      r = 1,
      obs_is_mut_list = obs_is_mut,
      error_ref_to_mut_list = error_ref_to_mut,
      error_mut_to_ref_list = error_mut_to_ref
    )
    Q_val <- -2 * (ll_0 - ll_A)
    df <- 1
    p_val <- pchisq(Q_val, df, lower.tail = FALSE)

    # Add result to output
    new_row <- data.frame(
      # Mutation annotations
      mutation_line,
      # EM results
      em_res,
      # Confidence interval
      tf_min = tf_CI$tf_min, tf_max = tf_CI$tf_max,
      # Misc.
      exp_count = sum(error_ref_to_mut),
      count = sum(obs_is_mut),
      coverage = length(obs_is_mut),
      obs_freq = mean(obs_is_mut),
      # Test results
      ll_0,
      ll_A,
      Q_val,
      df,
      p_val,
      mutation_detected = p_val <= alpha
    )

    res_df <- rbind(res_df, new_row)
  }
  return(res_df)
}
