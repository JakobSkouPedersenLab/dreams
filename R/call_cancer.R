#' Title
#'
#' @param mutations
#' @param reads
#' @param beta
#' @param use_warp_speed
#' @param calc_ci
#'
#' @return
#' @export
call_cancer <- function(mutations, reads, beta, use_warp_speed = TRUE, calc_ci = FALSE) {
  # Clean up mutations
  mutations_df <- mutations %>%
    rename(
      chr = CHROM,
      genomic_pos = POS,
      ref = REF,
      alt = ALT
    )

  # Summaries
  total_coverage <- nrow(reads)
  mutations_tested <- nrow(mutations_df)

  # Template for empty output
  empty_cancer_info <-
    data.frame(
      tf_est = 0, tf_min = 0, tf_max = 1,
      r_est = 0, r_min = 0, r_max = 1,
      Q_val = 0,
      ll_A = NA,
      ll_0 = NA,
      mutations_tested = mutations_tested,
      est_mutations_present = 0,
      total_coverage = total_coverage,
      total_count = 0,
      p_val = 1,
      EM_converged = NA,
      EM_steps = NA,
      fpeval = NA,
      objfeval = NA
    )

  # If no mutations return empty result
  if (mutations_tested == 0 | total_coverage == 0) {
    return(
      list(
        cancer_info = empty_cancer_info,
        mutation_info = NULL
      )
    )
  }

  # Prepare inputs for algorithm
  X_list <- list()
  error_ref_to_mut_list <- list()
  error_mut_to_ref_list <- list()

  mutations_w_cov_df <- NULL
  for (i in 1:mutations_tested) {
    chr <- mutations_df[i, "chr"]
    genomic_pos <- mutations_df[i, "genomic_pos"]
    ref <- mutations_df[i, "ref"]
    alt <- mutations_df[i, "alt"]

    mut_reads_ref_alt <-
      reads %>%
      filter(
        # Filter position
        chr == !!chr, genomic_pos == !!genomic_pos,
        # Filter mutation type (remove N and "other" alleles)
        obs == !!ref | obs == !!alt
      )

    # If no relevant coverage - skip to next mutation
    if (nrow(mut_reads_ref_alt) == 0) {
      next
    }

    # Save mutation
    mutations_w_cov_df <- rbind(mutations_w_cov_df, mutations_df[i, ])

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

  # If no positions have relevant coverage (ref or alt) or count of mu is 0 for every candidate
  n_mutaitons_w_relevant_coverage <- length(X_list)
  mutations_counts <- sapply(X_list, sum)
  if (n_mutaitons_w_relevant_coverage == 0 | all(mutations_counts == 0)) {
    return(
      list(
        cancer_info = empty_cancer_info,
        mutation_info = NULL
      )
    )
  }

  # EM algorithm
  EM_res <- EM_test(
    X_list = X_list,
    error_ref_to_mut_list = error_ref_to_mut_list,
    error_mut_to_ref_list = error_mut_to_ref_list,
    use_warp_speed = use_warp_speed,
    calc_ci = calc_ci
  )

  res <- list(
    cancer_info = EM_res$cancer_info,
    mutation_info = bind_cols(mutations_w_cov_df, EM_res$mutation_info)
  )

  return(res)
}
