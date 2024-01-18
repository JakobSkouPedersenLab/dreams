

#' Variant Calling in Indels
#'
#' This function performs variant calling on indels from a given BAM file using a predefined list of mutations.
#' It processes the mutations data, batch processes the genomic positions, and calls mutations using a specific model.
#'
#' @param mutations_df A dataframe containing the list of mutations to be analyzed.
#' @param bam_file_path Path to the BAM file containing sequencing data.
#' @param reference_path Path to the reference genome file, typically in FASTA format.
#' @param model The model to be used for calling mutations.
#' @param alpha Significance level for statistical testing, default is 0.05.
#' @param use_turboem Logical flag indicating whether to use the turboEM algorithm, default is TRUE.
#' @param calculate_confidence_intervals Logical flag indicating whether to calculate confidence intervals, default is FALSE.
#' @param batch_size Number of positions to process in each batch; if NULL, it's determined based on the data.
#'
#' @return A [data.frame()] with information about the individual mutation calls, including:
#' \describe{
#'   \item{chr, genomic_pos}{The genomic position of the mutation.}
#'   \item{ref, alt}{The reference and alternative allele.}
#'   \item{EM_converged}{If the EM algorithm converged.}
#'   \item{EM_steps, fpeval, objfeval}{Number of steps and function evaluations by the EM algorithm.}
#'   \item{tf_est}{The estiamted tumor fraction (allele fraction).}
#'   \item{tf_min, tf_max}{The confidence interval of \code{tf_est}.}
#'   \item{exp_count}{The expected count of the alternative allele under the error (null) model.}
#'   \item{count}{The count of the alternative allele.}
#'   \item{coverage}{The coverage used by the model (only referenceredas with and alternative allele).}
#'   \item{full_coverage}{The total coverage of the position (for reference).}
#'   \item{obs_freq}{The observed frequency of the alternative allele.}
#'   \item{ll_0, ll_A}{The value of the log-likelihood function under the null (tf=0) and alternative (tf>0) hypothesis.}
#'   \item{Q_val, df, p_val}{The chisq test statistic, degrees of freedom and p-value of the statistical test.}
#'   \item{mutation_detected}{Whether the mutation was detected at the supplied alpha level.}
#' }
#'
#' @seealso [call_mutations()], [call_cancer()], [train_dreams_model()]
#'
#' @export

dreams_vc_indels <- function(mutations_df, bam_file_path, reference_path, model, model_indels, alpha = 0.05,
                             use_turboem = TRUE, calculate_confidence_intervals = FALSE,
                             batch_size = NULL) {


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

  positions <- mutations_df %>%
    select(.data$chr, .data$genomic_pos) %>%
    distinct()

  if (is.null(batch_size)) {
    batch_size <- nrow(positions) + 1
  }

  position_batches <- positions %>% mutate(batch_idx = (row_number() %/% batch_size))


  mutation_calls <- NULL

  n_batches <- length(unique(position_batches$batch_idx))

  print(paste0("Calling mutations in ", n_batches, " batches:"))

  count <- 1

  # Get beta
  beta <- get_training_data_from_bam(bam_file_path, reference_path)$info$beta
  beta_indels <- get_training_data_from_bam_indel(bam_file_path, reference_path)$info$beta

  for (batch in sort(unique(position_batches$batch_idx))) {
    print(paste0("Calling batch ", count, "/", n_batches))
    count <- count + 1

    q <- position_batches %>% filter(batch_idx == batch)

    # Get read positions
    read_positions_df <- get_read_positions_from_BAM(
      bam_file_path = bam_file_path,
      chr = q$chr,
      genomic_pos = q$genomic_pos,
      reference_path
    )


    read_positions_df_indels <- get_read_positions_from_BAM_indels(
      bam_file_path = bam_file_path,
      chr = q$chr,
      genomic_pos = q$genomic_pos,
      reference_path
    )

    current_mutations <- mutations_df %>% filter(
      .data$chr %in% q$chr,
      .data$genomic_pos %in% q$genomic_pos
    )


    # Call mutations
    calls <- call_mutations_indels(
      mutations_df = current_mutations,
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

    mutation_calls <- rbind(mutation_calls, calls)
  }

  return(mutation_calls)
}

#' Call mutations from read positions
#'
#' @description This function evaluate the presence (calls) of individual mutations from a predefined list.
#' @inheritParams call_cancer
#' @param batch_size Number of positions to process at a time

#'
#' @return A [data.frame()] with information about the individual mutation calls, including:
#' \describe{
#'   \item{chr, genomic_pos}{The genomic position of the mutation.}
#'   \item{ref, alt}{The reference and alternative allele.}
#'   \item{EM_converged}{If the EM algorithm converged.}
#'   \item{EM_steps, fpeval, objfeval}{Number of steps and function evaluations by the EM algorithm.}
#'   \item{tf_est}{The estiamted tumor fraction (allele fraction).}
#'   \item{tf_min, tf_max}{The confidence interval of \code{tf_est}.}
#'   \item{exp_count}{The expected count of the alternative allele under the error (null) model.}
#'   \item{count}{The count of the alternative allele.}
#'   \item{coverage}{The coverage used by the model (only referenceredas with and alternative allele).}
#'   \item{full_coverage}{The total coverage of the position (for reference).}
#'   \item{obs_freq}{The observed frequency of the alternative allele.}
#'   \item{ll_0, ll_A}{The value of the log-likelihood function under the null (tf=0) and alternative (tf>0) hypothesis.}
#'   \item{Q_val, df, p_val}{The chisq test statistic, degrees of freedom and p-value of the statistical test.}
#'   \item{mutation_detected}{Whether the mutation was detected at the supplied alpha level.}
#' }
#' @import foreach
#' @seealso [call_cancer()], [train_dreams_model()]
#'
#' @export
#'
call_mutations_indels <- function(mutations_df, read_positions_df, read_positions_df_indels, model, model_indels, beta,
                                  beta_indels, alpha = 0.05, use_turboem = TRUE, calculate_confidence_intervals = FALSE,
                                  batch_size = NULL) {
  # If no mutations return empty result
  if (nrow(mutations_df) == 0) {
    return(data.frame())
  }

  if (nrow(read_positions_df) == 0) {
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

  positions <- mutations_df %>%
    select(.data$chr, .data$genomic_pos) %>%
    distinct()

  if (is.null(batch_size)) {
    batch_size <- nrow(positions) + 1
  }

  position_batches <- positions %>% mutate(batch_idx = (row_number() %/% batch_size))

  n_batches <- length(unique(position_batches$batch_idx))

  print(paste0("Calling mutations in ", n_batches, " batches:"))

  mutation_calls <- NULL


  for (batch in sort(unique(position_batches$batch_idx))) {
    q <- position_batches %>% filter(batch_idx == batch)


    current_read_positions_df <- read_positions_df %>% filter(
      .data$chr %in% q$chr,
      .data$genomic_pos %in% q$genomic_pos
    )

    current_read_positions_df_indels <- read_positions_df_indels %>% filter(
      .data$chr %in% q$chr,
      .data$genomic_pos %in% q$genomic_pos
    )

    if (nrow(current_read_positions_df) == 0) {
      return(data.frame())
    }


    current_mutations <- mutations_df %>% filter(
      .data$chr %in% q$chr,
      .data$genomic_pos %in% q$genomic_pos
    )

    # Prepare EM input
    em_input <- prepare_em_input_indels(mutations_df = current_mutations, read_positions_df = current_read_positions_df,
                                        read_positions_df_indels = current_read_positions_df_indels, model = model,
                                        model_indels = model_indels, beta = beta, beta_indels = beta_indels)

    obs_is_mut_list <- em_input$obs_is_mut_list
    error_ref_to_mut_list <- em_input$error_ref_to_mut_list
    error_mut_to_ref_list <- em_input$error_mut_to_ref_list

    # Add full coverage to mutations
    full_coverage_df <- current_read_positions_df %>%
      count(.data$chr, .data$genomic_pos, name = "full_coverage")

    current_mutations_df <- current_mutations %>%
      left_join(full_coverage_df, by = c("chr", "genomic_pos"))

    # Call mutations from reads
    res_df <- NULL
    for (i in 1:nrow(current_mutations_df)) {
      # Mutation information
      mutation_line <- current_mutations_df[i, ]

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
      if (calculate_confidence_intervals) {
        tf_CI <- get_tf_CI(
          obs_is_mut_list = obs_is_mut,
          error_mut_to_ref_list = error_mut_to_ref,
          error_ref_to_mut_list = error_ref_to_mut,
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
      p_val <- stats::pchisq(Q_val, df, lower.tail = FALSE)

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

    mutation_calls <- rbind(mutation_calls, res_df)
  }
  return(mutation_calls)
}



