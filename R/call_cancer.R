#' Call cancer from read positions
#'
#' @description This function evaluates the presence of cancer in a sample by combining the cancerous signal across a catalogue of candidate mutations.
#'
#' @param mutations_df A [data.frame()] with candidate mutations (SNVs) (chromosome, positions, reference and alternative)
#' @param reads_df A [data.frame()] with read-positions.
#' @param model A dreams model. See [train_model()].
#' @param beta Down sampling parameter from (TODO: Link)
#' @param alpha Alpha-level used for testing and confidence intervals. Default is 0.05.
#' @param use_turboem Logical. Should [turboEM::turboem()] be used for EM algorithm? Default is TRUE.
#' @param calculate_confidence_intervals Logical. Should confidence intervals be calculated? Default is FALSE.
#'
#' TODO: Make sub list of items in data.frames
#' @return
#' A \code{list()} with:
#'   \item{cancer_info}{A [data.frame()] with results for cancer calling across all mutations}
#'   \item{mutation_info}{A [data.frame()] with information about the individual mutations}
#'
#' @seealso [call_mutations()], [train_model()]
#'
#' @export
call_cancer <- function(mutations_df, reads_df, model, beta, alpha = 0.05, calculate_confidence_intervals = FALSE, use_turboem = TRUE) {
  # If no mutations return empty result
  if (nrow(mutations_df) == 0) {
    return(
      list(
        cancer_info = data.frame(),
        mutation_info = data.frame()
      )
    )
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
  if (!all(reads_expected_columns %in% colnames(reads_df))) {
    stop("reads_df should have the columns ['chr', genomic_pos', 'ref, 'obs']")
  }

  # Prepare inputs for algorithm
  em_input <- prepare_em_input(mutations_df = mutations_df, reads_df = reads_df, model = model, beta = beta)
  obs_is_mut_list <- em_input$obs_is_mut_list
  error_ref_to_mut_list <- em_input$error_ref_to_mut_list
  error_mut_to_ref_list <- em_input$error_mut_to_ref_list


  # EM algorithm
  em_res <- get_em_parameter_estimates(
    obs_is_mut_list = obs_is_mut_list,
    error_ref_to_mut_list = error_ref_to_mut_list,
    error_mut_to_ref_list = error_mut_to_ref_list,
    use_turboem = use_turboem
  )

  # Confidence intervals
  if (calculate_confidence_intervals) {
    tf_CI <- get_tf_CI(
      obs_is_mut_list = obs_is_mut_list,
      error_mut_to_ref_list = error_mut_to_ref_list,
      error_ref_to_mut_list = error_ref_to_mut_list,
      r_est = em_res$r_est,
      tf_est = em_res$tf_est,
      alpha = alpha
    )
    r_CI <- get_r_CI(
      obs_is_mut_list = obs_is_mut_list,
      error_mut_to_ref_list = error_mut_to_ref_list,
      error_ref_to_mut_list = error_ref_to_mut_list,
      r_est = em_res$r_est,
      tf_est = em_res$tf_est,
      alpha = alpha
    )
  } else {
    tf_CI <- list(tf_min = NA, tf_max = NA)
    r_CI <- list(r_min = NA, r_max = NA)
  }

  # Test significance
  ll_0 <- log_likelihood(
    obs_is_mut_list = obs_is_mut_list,
    error_mut_to_ref_list = error_mut_to_ref_list,
    error_ref_to_mut_list = error_ref_to_mut_list,
    r = 0,
    tf = 0
  )

  ll_A <- log_likelihood(
    obs_is_mut_list = obs_is_mut_list,
    error_mut_to_ref_list = error_mut_to_ref_list,
    error_ref_to_mut_list = error_ref_to_mut_list,
    r = em_res$r_est,
    tf = em_res$tf_est
  )

  Q_val <- -2 * (ll_0 - ll_A)
  df <- 2
  p_val <- stats::pchisq(Q_val, df, lower.tail = FALSE)

  # Collect cancer information
  cancer_info <-
    data.frame(
      tf_est = em_res$tf_est,
      tf_min = tf_CI$tf_min, tf_max = tf_CI$tf_max,
      r_est = em_res$r_est,
      r_min = r_CI$r_min, r_max = r_CI$r_max,
      Q_val = Q_val,
      ll_A = ll_A,
      ll_0 = ll_0,
      mutations_tested = length(obs_is_mut_list),
      est_mutations_present = sum(em_res$P_mut_is_present),
      total_coverage = sum(sapply(obs_is_mut_list, length)),
      total_count = sum(sapply(obs_is_mut_list, sum)),
      EM_converged = em_res$EM_converged,
      EM_steps = em_res$EM_steps,
      fpeval = em_res$fpeval,
      objfeval = em_res$objfeval,
      p_val = p_val,
      cancer_detected = p_val <= alpha
    )

  # Collect mutation information
  mutation_info <-
    dplyr::bind_cols(
      mutations_df,
      data.frame(
        P_mut_is_present = em_res$P_mut_is_present,
        exp_count = sapply(error_ref_to_mut_list, sum),
        count = sapply(obs_is_mut_list, sum),
        coverage = sapply(obs_is_mut_list, length),
        obs_freq = sapply(obs_is_mut_list, mean)
      )
    )

  return(
    list(
      cancer_info = cancer_info,
      mutation_info = mutation_info
    )
  )
}
