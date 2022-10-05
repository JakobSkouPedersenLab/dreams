

#' Call mutations in a bam-file
#'
#' @description This function evaluate the presence (calls) of individual mutations from a predefined list.
#' @inheritParams call_cancer
#' @param bam_file_path Path to .BAM-file
#' @param reference_path Path to reference genome e.g. FASTA-file.
#' @param ncores Number of processing cores
#' @param batch_size Number of positions to process at a time
#' @param log_file write log-file to this path



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
#' @import parallel
#' @import doParallel
#'
#' @export

dreams_vc_parallel <- function(mutations_df, bam_file_path, reference_path, model,
                               beta, alpha = 0.05, use_turboem = TRUE, calculate_confidence_intervals = FALSE,
                               batch_size = NULL, ncores = 1, log_file = NULL) {
  if (nrow(mutations_df) == 0) {
    return(data.frame())
  }

  if (nrow(mutations_df) < ncores) {
    ncores <- max(nrow(mutations_df), 1)
  }


  serial_model <- keras::serialize_model(model)

  mutations_df <- mutations_df %>%
    dplyr::select(
      "chr" = matches("chr|CHR|CHROM"),
      "genomic_pos" = matches("pos|POS"),
      "ref" = matches("ref|REF"),
      "alt" = matches("alt|ALT|obs|OBS")
    ) %>%
    mutate(idx = sort(row_number() %% ncores, decreasing = F))

  index_list <- unique(mutations_df$idx)


  cl <- makeCluster(ncores)
  doParallel::registerDoParallel(cl)

  mutation_calls <- foreach::foreach(
    i = index_list,
    .combine = rbind,
    .packages = c("keras", "tensorflow", "parallel", "doParallel"),
    .errorhandling = "pass",
    .export = "dreams_vc"
  ) %dopar% {
    unserial_model <- keras::unserialize_model(serial_model)

    if (!is.null(log_file)) {
      sink(paste0(log_file, "_", i))
    }


    mutations <- mutations_df %>%
      dplyr::filter(idx == i) %>%
      dplyr::select(-idx)

    print(mutations)


    current_calls <- dreams_vc(
      mutations_df = mutations,
      bam_file_path = bam_file_path,
      reference_path = reference_path,
      model = unserial_model,
      beta = beta,
      use_turboem = use_turboem,
      batch_size = batch_size,
      calculate_confidence_intervals = calculate_confidence_intervals,
      alpha = alpha
    )


    if (nrow(current_calls) > 0) {
      print(head(current_calls))
      return(current_calls)
    } else {
      print(current_calls)
      return(NULL)
    }
  }

  sink()
  stopCluster(cl)
  return(mutation_calls)
}

#' Call mutations in a bam-file
#'
#' @description This function evaluate the presence (calls) of individual mutations from a predefined list.
#' @inheritParams call_cancer
#' @param bam_file_path Path to .BAM-file
#' @param reference_path Path to reference genome e.g. FASTA-file.
#' @param batch_size Number of positions to process at a time


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

dreams_vc <- function(mutations_df, bam_file_path, reference_path, model,
                      beta, alpha = 0.05, use_turboem = TRUE, calculate_confidence_intervals = FALSE,
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

  for (batch in sort(unique(position_batches$batch_idx))) {
    print(paste0("Calling batch ", count, "/", n_batches))
    count <- count + 1

    q <- position_batches %>% filter(batch_idx == batch)

    print("DIM Q")
    print(dim(q))

    # Get read positions
    read_positions_df <- get_read_positions_from_BAM(
      bam_file_path = bam_file_path,
      chr = q$chr,
      genomic_pos = q$genomic_pos,
      reference_path
    )

    print("DIM RP_DF")
    print(dim(read_positions_df))
    print (unique(read_positions_df$genomic_pos))

    current_mutations <- mutations_df %>% filter(
      .data$chr %in% q$chr,
      .data$genomic_pos %in% q$genomic_pos
    )

    print("current_mutations")
    print(dim(current_mutations))
    print (current_mutations)


    # Call mutations
    calls <- call_mutations(
      mutations_df = current_mutations,
      read_positions_df = read_positions_df,
      model = model,
      beta = beta,
      alpha = alpha,
      use_turboem = use_turboem,
      calculate_confidence_intervals = calculate_confidence_intervals
    )

    print("calls")
    print(dim(calls))


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
#'
#' @seealso [call_cancer()], [train_dreams_model()]
#'
#' @export
call_mutations <- function(mutations_df, read_positions_df, model, beta,
                           alpha = 0.05, use_turboem = TRUE, calculate_confidence_intervals = FALSE, batch_size = NULL) {
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

    print ("Q")
    print (q)

    current_read_positions_df <- read_positions_df %>% filter(
      .data$chr %in% q$chr,
      .data$genomic_pos %in% q$genomic_pos
    )

    current_mutations <- mutations_df %>% filter(
      .data$chr %in% q$chr,
      .data$genomic_pos %in% q$genomic_pos
    )


    # Prepare EM input
    em_input <- prepare_em_input(mutations_df = current_mutations, read_positions_df = current_read_positions_df, model = model, beta = beta)

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

      print ("LINE")
      print (mutation_line)

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

      print ("EM")
      print (em_res)

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
