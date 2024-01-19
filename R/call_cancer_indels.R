#' Cancer Detection in Indels Using DREAMS
#'
#' This function performs cancer detection based on indel information from a BAM file using a predefined list of mutations.
#' It processes the mutations data and uses the provided model to make cancer calls.
#'
#' @param mutations_df A dataframe containing the list of mutations to be analyzed.
#' @param bam_file_path Path to the BAM file containing sequencing data.
#' @param reference_path Path to the reference genome file, typically in FASTA format.
#' @param model The model to be used for cancer detection.
#' @param model_indels The model to be used for calling mutations for indels.
#' @param alpha Significance level for statistical testing, default is 0.05.
#' @param calculate_confidence_intervals Logical flag indicating whether to calculate confidence intervals, default is FALSE.
#' @param use_turboem Logical flag indicating whether to use the turboEM algorithm, default is TRUE.
#'
#' @return
#' A \code{list()} with:
#' \itemize{
#'   \item \strong{cancer_info} A [data.frame()] with results for cancer calling across all mutations:
#'     \describe{
#'       \item{tf_est}{The estiamted tumor fraction (allele fraction).}
#'       \item{tf_min, tf_max}{The confidence interval of \code{tf_est}.}
#'       \item{r_est, est_mutations_present}{The estiamted fraction/number of candidate mutations present in the sample.}
#'       \item{r_min, r_max}{The confidence interval of \code{r_est}.}
#'       \item{mutations_tested}{Number of candidate mutations tested.}
#'       \item{total_coverage, total_count}{Total count and coverage across all mutations (only reference and alternative allele(s).}
#'       \item{mutations_tested}{Number of candidate mutations tested.}
#'       \item{EM_converged}{If the EM algorithm converged.}
#'       \item{EM_steps, fpeval, objfeval}{Number of steps and function evaluations by the EM algorithm.}
#'       \item{ll_0, ll_A}{The value of the log-likelihood function under the null (tf=0) and alternative (tf>0) hypothesis.}
#'       \item{Q_val, df, p_val}{The chisq test statistic, degrees of freedom and p-value of the statistical test.}
#'       \item{cancer_detected}{Whether cancer was detected at the supplied alpha level.}
#'     }
#'
#'   \item \strong{mutation_info} A [data.frame()] with information about the individual mutations:
#'       \describe{
#'         \item{chr, genomic_pos}{The genomic position of the mutation.}
#'         \item{ref, alt}{The reference and alternative allele.}
#'         \item{P_mut_is_present}{The estimated probability the mutation is present in the sample.}
#'         \item{exp_count}{The expected count of the alternative allele under the error (null) model.}
#'         \item{count}{The count of the alternative allele.}
#'         \item{coverage}{The coverage used by the model (only referenceredas with and alternative allele).}
#'         \item{obs_freq}{The observed frequency of the alternative allele.}
#'      }
#'
#' }
#' @seealso [call_mutations()], [train_dreams_model()]
#'
#' @export
dreams_cc_indels <- function(mutations_df, bam_file_path, reference_path, model, model_indels,
                      alpha = 0.05, calculate_confidence_intervals = FALSE, use_turboem = TRUE) {
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


  # Get read positions
  read_positions_df <- get_read_positions_from_BAM(
    bam_file_path = bam_file_path,
    chr = mutations_df$chr,
    genomic_pos = mutations_df$genomic_pos,
    reference_path
  )

  read_positions_df <- get_read_positions_from_BAM_indels(
    bam_file_path = bam_file_path,
    chr = mutations_df$chr,
    genomic_pos = mutations_df$genomic_pos,
    reference_path
  )


  beta <- get_training_data_from_bam_indel(bam_file_path, reference_path)$info$beta
  beta_indels <- get_training_data_from_bam_indel(bam_file_path, reference_path)$info$beta

  # Call cancer
  cancer_call <- call_cancer_indels(
      mutations_df = mutations_df,
      read_positions_df = read_positions_df,
      read_positions_df_indels = read_positions_df_indels,
      model = model,
      model_indels = model_indels,
      beta = beta,
      beta_indels = beta_indels,
      alpha = alpha,
      use_turboem = use_turboem,
      calculate_confidence_intervals = calculate_confidence_intervals
    )

  return(cancer_call)


}


#' Call cancer from read positions
#'
#' @description This function evaluates the presence of cancer in a sample by combining the cancerous signal across a catalogue of candidate mutations.
#'
#' @param mutations_df A [data.frame()] with candidate mutations (SNVs) (chromosome, positions, reference and alternative)
#' @param read_positions_df A [data.frame()] with read-positions. See [get_read_positions_from_BAM()]
#' @param model A dreams model. See [train_dreams_model()].
#' @param beta Down sampling parameter from [get_training_data()] for correcting the error-rates from the DREAMS model.
#' @param alpha Alpha-level used for testing and confidence intervals. Default is 0.05.
#' @param use_turboem Logical. Should [turboEM::turboem()] be used for EM algorithm? Default is TRUE.
#' @param calculate_confidence_intervals Logical. Should confidence intervals be calculated? Default is FALSE.
#'
#' @return
#' A \code{list()} with:
#' \itemize{
#'   \item \strong{cancer_info} A [data.frame()] with results for cancer calling across all mutations:
#'     \describe{
#'       \item{tf_est}{The estiamted tumor fraction (allele fraction).}
#'       \item{tf_min, tf_max}{The confidence interval of \code{tf_est}.}
#'       \item{r_est, est_mutations_present}{The estiamted fraction/number of candidate mutations present in the sample.}
#'       \item{r_min, r_max}{The confidence interval of \code{r_est}.}
#'       \item{mutations_tested}{Number of candidate mutations tested.}
#'       \item{total_coverage, total_count}{Total count and coverage across all mutations (only reference and alternative allele(s).}
#'       \item{mutations_tested}{Number of candidate mutations tested.}
#'       \item{EM_converged}{If the EM algorithm converged.}
#'       \item{EM_steps, fpeval, objfeval}{Number of steps and function evaluations by the EM algorithm.}
#'       \item{ll_0, ll_A}{The value of the log-likelihood function under the null (tf=0) and alternative (tf>0) hypothesis.}
#'       \item{Q_val, df, p_val}{The chisq test statistic, degrees of freedom and p-value of the statistical test.}
#'       \item{cancer_detected}{Whether cancer was detected at the supplied alpha level.}
#'     }
#'
#'   \item \strong{mutation_info} A [data.frame()] with information about the individual mutations:
#'       \describe{
#'         \item{chr, genomic_pos}{The genomic position of the mutation.}
#'         \item{ref, alt}{The reference and alternative allele.}
#'         \item{P_mut_is_present}{The estimated probability the mutation is present in the sample.}
#'         \item{exp_count}{The expected count of the alternative allele under the error (null) model.}
#'         \item{count}{The count of the alternative allele.}
#'         \item{coverage}{The coverage used by the model (only referenceredas with and alternative allele).}
#'         \item{obs_freq}{The observed frequency of the alternative allele.}
#'      }
#'
#' }
#' @seealso [call_mutations()], [train_dreams_model()]
#'
#' @export
call_cancer_indels <- function(mutations_df, read_positions_df, read_positions_df_indels, model, model_indels, beta, beta_indels, alpha = 0.05, calculate_confidence_intervals = FALSE, use_turboem = TRUE) {
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
  if (!all(reads_expected_columns %in% colnames(read_positions_df))) {
    stop("read_positions_df should have the columns ['chr', genomic_pos', 'ref, 'obs']")
  }

  # Prepare inputs for algorithm
  em_input <- prepare_em_input_indels(mutations_df = mutations_df, read_positions_df = read_positions_df,
                                      read_positions_df_indels = read_positions_df_indels, model = model,
                                      model_indels = model_indels, beta = beta, beta_indels = beta_indels)
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
