
#' Extract training data from BAM files (indel)
#'
#' @param bam_paths Vector of strings. Paths to \code{.bam} files to extract
#'   training data from.
#' @param reference_path String. Path to reference genome fasta file.
#' @param bed_include_path String. Path to bed-file with regions to include.
#'   Default is \code{NULL}.
#' @param factor Number between 0 and 1. Ratio between negative and positive
#'   data. Default is 1.
#' @param common_positions_to_exclude_paths Vector of strings. List of files
#'   with positions to exclude from all samples. Default is \code{NULL}.
#' @param positions_to_exclude_paths Vector of strings. List of files with
#'   positions to exclude from training with length equal to number of samples.
#'   Default is \code{NULL}.
#' @param mm_rate_max Number between 0 and 1. Maximum mismatch rate in position.
#'   Default is 1.
#' @param verbose TODO: Write this
#'
#' @return A list containing two elements:
#' \itemize{
#'   \item `data`: A `tbl_df`
#'   \item `info`: A `data.frame`
#' }
#' @export
#'
#' @seealso [train_dreams_model()] Function for training model.
get_training_data_indels <- function(bam_paths,
                              reference_path,
                              bed_include_path = NULL,
                              factor = 1,
                              common_positions_to_exclude_paths = NULL,
                              positions_to_exclude_paths = NULL,
                              mm_rate_max = 1,
                              verbose = F) {

  # Check if there is a position exclude path for each bam file
  if ((!is.null(positions_to_exclude_paths) &
    (length(bam_paths) != length(positions_to_exclude_paths)))) {
    stop("Wrong number of exclude paths")
  }

  training_data <- NULL
  info <- NULL

  n_bam_files = length(bam_paths)

  for (bam_idx in 1:n_bam_files) {
    bam_path <- bam_paths[[bam_idx]]

    if (verbose) {
      cat("file ", bam_idx, "/", n_bam_files, "\n")
    }

    # Combine sample specific position exclusion with common exclusion
    if (!is.null(positions_to_exclude_paths)) {
      current_positions_to_exclude_paths <- c(common_positions_to_exclude_paths, positions_to_exclude_paths[[bam_idx]])
    } else {
      current_positions_to_exclude_paths <- common_positions_to_exclude_paths
    }


    # Get training data for single bam file
    current_training_data <- get_training_data_from_bam_indel(
      bam_path = bam_path,
      reference_path = reference_path,
      bed_include_path = bed_include_path,
      positions_to_exclude_paths = current_positions_to_exclude_paths,
      factor = factor,
      mm_rate_max = mm_rate_max
    )

    training_data <- rbind(training_data, current_training_data$data)
    info <- rbind(info, current_training_data$info)
  }

  # Collect output info for beta calculation
  output_info <- data.frame(
    total_mismatches = sum(info$n_mismatches),
    total_matches = sum(info$n_matches),
    total_coverage = sum(info$total_coverage)
  ) %>% mutate(beta = .data$total_matches / (.data$total_coverage - .data$total_mismatches))

  return(list(
    data = training_data,
    info = output_info
  ))
}


#' Extract Training Data from BAM File (indels)
#'
#' Extracts training data from a BAM file by integrating information from reference and BED files.
#' It processes genomic positions of indels, extracts features, filters based on mismatch rates,
#' and combines positive and negative samples to form the training dataset.
#'
#' @param bam_path Path to the BAM file.
#' @param reference_path Path to the reference genome file.
#' @param bed_include_path Optional; BED file defining regions to include in the analysis.
#' @param factor The ratio of negative to positive data in the output.
#' @param positions_to_exclude_paths Optional; paths to files defining positions to exclude from training.
#' @param mm_rate_max Maximum mismatch rate allowed in a position.
#' @keywords internal
#'
#' @return A list with two elements: `data`, a `data.frame` containing the combined positive and negative training data,
#' and `info`, a `data.frame` containing metadata about the training set.
#'
get_training_data_from_bam_indel <- function(bam_path, reference_path, bed_include_path = NULL, factor = 1, positions_to_exclude_paths = NULL, mm_rate_max = 1) {
  bam_df <- load_BAM(bam_path)

  # Add genomic positions of indels
  indel_bam_df <- extract_indel_info(bam_df)

  # Add features
  indels_positions_df <-
    extract_features_from_bam_indels(
      bam_df = indel_bam_df,
      reference_path = reference_path
    )

  # Filter indels
  indels <-
    filter_mismatch_positions(
      read_positions = indels_positions_df,
      bam_file = bam_path,
      mm_rate_max = mm_rate_max,
      bed_include_path = bed_include_path,
      positions_to_exclude_paths = positions_to_exclude_paths
    )

  positive_samples <- indels$data
  info <- indels$info

  n_samples <- nrow(positive_samples) * factor


  # Generate negative samples
  negative_read_positions_df <-
    sample_negative_read_positions_indels(
      bam_df = bam_df,
      n_samples = n_samples
    )

  # Add features
  negative_samples <-
    extract_features_from_bam_indels_negatives(
      bam_df = negative_read_positions_df,
      reference_path = reference_path
    )

  info <- info %>% mutate(n_matches = nrow(negative_samples))
  info <- info %>% mutate(beta = nrow(negative_samples) / (info$total_coverage - nrow(positive_samples)))

  output_data <- rbind(positive_samples, negative_samples)

  output_list <- list(
    data = output_data,
    info = info
  )

  return(output_list)
}





