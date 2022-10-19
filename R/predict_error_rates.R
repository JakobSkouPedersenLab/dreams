limit <- function(x, min = 0, max = 1) {
  return(pmax(pmin(x, max), min))
}

#' @import dplyr
correct_errors_predictions <- function(error_df, beta) {
  error_df %>%
    mutate(
      error_prob_sample = case_when(
        .data$ref == "A" ~ 1 - A,
        .data$ref == "C" ~ 1 - C,
        .data$ref == "T" ~ 1 - T,
        .data$ref == "G" ~ 1 - G
      ),
      error_prob_corrected = beta * .data$error_prob_sample / (beta * .data$error_prob_sample - .data$error_prob_sample + 1),
      A_corrected = ifelse(.data$ref == "A",
        1 - .data$error_prob_corrected,
        .data$error_prob_corrected * .data$A / (.data$error_prob_sample)
      ),
      C_corrected = ifelse(.data$ref == "C",
        1 - .data$error_prob_corrected,
        .data$error_prob_corrected * .data$C / (.data$error_prob_sample)
      ),
      G_corrected = ifelse(.data$ref == "G",
        1 - .data$error_prob_corrected,
        .data$error_prob_corrected * .data$G / (.data$error_prob_sample)
      ),
      T_corrected = ifelse(.data$ref == "T",
        1 - .data$error_prob_corrected,
        .data$error_prob_corrected * .data$T / (.data$error_prob_sample)
      ),
      A_corrected = ifelse(is.nan(.data$A_corrected), 0, .data$A_corrected) %>% limit(),
      C_corrected = ifelse(is.nan(.data$C_corrected), 0, .data$C_corrected) %>% limit(),
      G_corrected = ifelse(is.nan(.data$G_corrected), 0, .data$G_corrected) %>% limit(),
      T_corrected = ifelse(is.nan(.data$T_corrected), 0, .data$T_corrected) %>% limit()
    )
}
#' Predict error rates for read positions
#'
#' @param read_positions_df Dataframe with read positions
#' @param model model file
#' @param beta beta value
#'
#' @importFrom stats predict
#' @export

predict_error_rates <- function(read_positions_df, model, beta) {

  # Predict error rates for read positions from trained DREAM model
  if (nrow(read_positions_df) == 0) {
    predictions_df <-
      data.frame(
        A = numeric(),
        T = numeric(),
        C = numeric(),
        G = numeric()
      )
  } else {
    predictions_df <- model %>%
      predict(read_positions_df) %>%
      data.frame() %>%
      rename(
        A = .data$X1,
        T = .data$X2,
        C = .data$X3,
        G = .data$X4
      )
  }

  # Add predictions to read positions and correct error rates with bea
  predicted_errors <-
    read_positions_df %>%
    bind_cols(predictions_df) %>%
    correct_errors_predictions(beta = beta)

  return(predicted_errors)
}





#' Parallel predict error rates
#'
#' @description This function evaluate the presence (calls) of individual mutations from a predefined list.
#' @inheritParams call_cancer
#' @param bam_file_path Path to .BAM-file
#'
#' @import parallel
#' @import doParallel
#' @import dplyr
#'
#' @export

predict_error_rates_parallel <- function(mutations_df, bam_file_path, reference_path, model,
                                         beta = NULL, factor = NULL, bed_file = NULL, ncores = 1, log_file = NULL, mm_rate_max = 0.05, batch_size = NULL) {

  # If no beta value

  if (is.null(beta) && is.null(factor)) {
    stop("Please provide beta factor of scaling factor")
  } else if (is.null(beta)) {
    beta <- calculate_beta_factor(
      bam_file_path = bam_file_path,
      factor = factor,
      mm_rate_max = mm_rate_max,
      bed_file = bed_file,
      reference_path = reference_path
    )
  }

  if (nrow(mutations_df) == 0) {
    return(data.frame())
  }

  if (nrow(mutations_df) < ncores) {
    ncores <- max(nrow(mutations_df), 1)
  }


  serial_model <- keras::serialize_model(model)

  mutations_df <- mutations_df %>%
    select(
      "chr" = matches("chr|CHR|CHROM"),
      "genomic_pos" = matches("pos|POS"),
      "ref" = matches("ref|REF"),
      "alt" = matches("alt|ALT|obs|OBS")
    ) %>%
    mutate(idx = sort(row_number() %% ncores, decreasing = F))

  index_list <- unique(mutations_df$idx)


  cl <- makeCluster(ncores)
  doParallel::registerDoParallel(cl)

  print(mutations_df)

  error_rates <- foreach::foreach(
    i = index_list,
    .combine = rbind,
    .packages = c("keras", "tensorflow", "parallel", "doParallel", "dplyr"),
    .errorhandling = "pass",
    .export = "dreams_vc"
  ) %dopar% {
    unserial_model <- keras::unserialize_model(serial_model)


    if (!is.null(log_file)) {
      sink(paste0(log_file, "_", i))
    }

    print("INSIDE")
    print("INSIDE2")

    mutations <- mutations_df %>%
      dplyr::filter(idx == i) %>%
      dplyr::select(-idx)


    print("MUTATIONS")
    print (head(mutations))


    current_error_rates <- predict_error_rates_batches(
      mutations_df = mutations,
      bam_file_path = bam_file_path,
      reference_path = reference_path,
      batch_size = batch_size,
      model = model,
      beta = beta
    ) %>% dplyr::mutate(
      beta = beta
    )

    print(head(current_error_rates))


    return(current_error_rates)
  }

  sink()
  stopCluster(cl)
  return(error_rates)
}



#' Parallel predict error rates
#'
#' @description This function evaluate the presence (calls) of individual mutations from a predefined list.
#' @inheritParams call_cancer
#' @param bam_file_path Path to .BAM-file
#'
#' @import parallel
#' @import doParallel
#' @import dplyr
#'
#' @export

predict_error_rates_batches <- function(mutations_df, bam_file_path, reference_path, model,
                                        beta = NULL, factor = NULL, bed_file = NULL,  mm_rate_max = 0.05, batch_size = NULL) {

  # If no beta value

  print ("INSIDE BATCHES")

  if (is.null(beta) && is.null(factor)) {
    stop("Please provide beta factor of scaling factor")
  } else if (is.null(beta)) {
    beta <- calculate_beta_factor(
      bam_file_path = bam_file_path,
      factor = factor,
      mm_rate_max = mm_rate_max,
      bed_file = bed_file,
      reference_path = reference_path
    )
  }

  print(beta)


  # Clean up mutations
  mutations_df <- mutations_df %>%
    select(
      "chr" = matches("chr|CHR|CHROM"),
      "genomic_pos" = matches("pos|POS"),
      "ref" = matches("ref|REF"),
      "alt" = matches("alt|ALT|obs|OBS")
    )

  print (mutations_df)

  positions <- mutations_df %>%
    select(chr, genomic_pos) %>%
    distinct()

  if (is.null(batch_size)) {
    batch_size <- nrow(positions) + 1
  }

  position_batches <- positions %>% dplyr::mutate(batch_idx = (row_number() %/% !!batch_size))


  error_rates <- NULL

  n_batches <- length(unique(position_batches$batch_idx))

  print(paste0("Calling mutations in ", n_batches, " batches:"))

  count <- 1

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

    current_error_rates <- predict_error_rates(
      read_positions_df = read_positions_df,
      model = model,
      beta = beta
    ) %>% dplyr::mutate(
      beta = beta
    )


    error_rates <- rbind(error_rates, current_error_rates)
  }

  return(error_rates)
}

