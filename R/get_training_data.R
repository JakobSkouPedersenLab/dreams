
#' Extract training data from BAM files
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
#'   \item `data`: A `tbl_df` with dimensions 2 x 22.
#'   \item `info`: A `data.frame` with dimensions 1 x 4.
#' }
#' @export
#'
#' @seealso [train_dreams_model()] Function for training model.
get_training_data <- function(bam_paths,
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
    current_training_data <- get_training_data_from_bam(
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


#' Title
#'
#' @param bam_path Path to BAM file
#' @param reference_path Path to reference file
#' @param bed_include_path BED regions to include
#' @param factor ratio between negative and positive data
#' @param positions_to_exclude_paths positions to exclude from training
#' @param mm_rate_max maximum mismatch rate in position
#' @keywords internal
#'
#' @return `data.frame` with training data for a \code{.bam} file

get_training_data_from_bam <- function(bam_path, reference_path, bed_include_path = NULL, factor = 1, positions_to_exclude_paths = NULL, mm_rate_max = 1) {
  bam_df <- load_BAM(bam_path)

  # Add genomic positions of mismatches
  mismatch_bam_df <- extract_mismatch_positions(bam_df)

  # Add features
  mismatch_positions_df <-
    extract_features_from_bam(
      bam_df = mismatch_bam_df,
      reference_path = reference_path
    )

  # Filter mismatches
  mismatches <-
    filter_mismatch_positions(
      read_positions = mismatch_positions_df,
      bam_file = bam_path,
      mm_rate_max = mm_rate_max,
      bed_include_path = bed_include_path,
      positions_to_exclude_paths = positions_to_exclude_paths
    )

  positive_samples <- mismatches$data
  info <- mismatches$info

  n_samples <- nrow(positive_samples) * factor


  # Generate negative samples
  negative_read_positions_df <-
    sample_negative_read_positions(
      bam_df = bam_df,
      n_samples = n_samples
    )

  # Add features
  negative_samples <-
    extract_features_from_bam(
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



#' Filter Mismatch Positions in Sequencing Data
#'
#' This function filters read positions from sequencing data based on mismatch rates and specific genomic regions. It is designed to process data from BAM files and can include or exclude specific regions and positions as needed.
#'
#' @param read_positions A dataframe containing positions of reads (fx df obtained from extract_features_from_bam()) . Each row should represent a unique read position with relevant information such as observed nucleotides.
#' @param bam_file A string specifying the path to a BAM file.
#' @param mm_rate_max A numeric value representing the maximum mismatch rate allowed for positions. Defaults to 1. Positions with a mismatch rate higher than this threshold will be excluded.
#' @param bed_include_path An optional string specifying the path to a BED file that contains genomic regions to include in the analysis. If `NULL`, all regions are included.
#' @param positions_to_exclude_paths An optional vector of strings specifying paths to files containing positions that should be excluded from the analysis. If `NULL`, no positions are excluded.
#' @keywords internal
#'
#' @return A list with two elements: `data` containing the filtered read positions dataframe, and `info` containing a dataframe with summary information such as the number of mismatches and total coverage.
#'
#' @importFrom Rsamtools pileup
#' @importFrom readr read_csv

filter_mismatch_positions <- function(read_positions, bam_file, mm_rate_max = 1, bed_include_path = NULL, positions_to_exclude_paths = NULL) {

  read_positions_filtered = read_positions %>% filter(obs != "N")

  # Load coverage data

  included_regions_granges <- bed_to_granges(bed_include_path)

  pp <- Rsamtools::PileupParam(
    max_depth = 250000000, min_base_quality = 13, min_mapq = 0,
    min_nucleotide_depth = 1, min_minor_allele_depth = 0,
    distinguish_strands = FALSE, distinguish_nucleotides = FALSE,
    ignore_query_Ns = TRUE, include_deletions = TRUE, include_insertions = FALSE,
    left_bins = NULL, query_bins = NULL, cycle_bins = NULL
  )

  coverage_data <- Rsamtools::pileup(bam_file, pileupParam = pp, scanBamParam = ScanBamParam(which = included_regions_granges)) %>%
    rename(chr = .data$seqnames, genomic_pos = .data$pos, coverage = .data$count)

  # Filter heterozygote positions

  read_positions_summarized <- read_positions %>%
    group_by(.data$chr, .data$genomic_pos) %>%
    summarize(n_mismatches = n()) %>%
    ungroup()

  # Join with coverage dataframe - all positions if included_regions is NULL

  # Remove positions with high mismatch rate in mismatch  and coverage data

  read_position_mm_rate <- read_positions_summarized %>%
    inner_join(coverage_data, by = c("chr", "genomic_pos")) %>%
    mutate(mm_rate = .data$n_mismatches / .data$coverage)

  read_position_filter <- read_position_mm_rate %>%
    filter(.data$mm_rate < mm_rate_max)

  read_positions_filtered <- read_positions_filtered %>%
    semi_join(read_position_filter, by = c("chr", "genomic_pos"))

  coverage_data_filtered <- coverage_data %>%
    anti_join(read_position_mm_rate %>% filter(.data$mm_rate > mm_rate_max), by = c("chr", "genomic_pos"))

  # Remove unwanted positions based on exclude files

  if (!is.null(positions_to_exclude_paths)) {
    for (p in positions_to_exclude_paths) {
      positions_to_exclude <- read_csv(p, show_col_types = FALSE)

      read_positions_filtered <- read_positions_filtered %>%
        anti_join(positions_to_exclude, by = c("chr", "genomic_pos"))

      coverage_data_filtered <- coverage_data_filtered %>%
        anti_join(positions_to_exclude, by = c("chr", "genomic_pos"))
    }
  }

  # Output beta info

  beta_info <- data.frame(
    n_mismatches = nrow(read_positions_filtered),
    total_coverage = sum(coverage_data_filtered$coverage)
  )


  return(list(
    data = read_positions_filtered,
    info = beta_info
  ))
}

#' Convert a BED File to a GRanges Object
#'
#' This function reads a BED file and converts it into a GRanges object. BED files typically contain genomic intervals
#' (like regions of the genome associated with certain features or annotations) and are formatted in a specific way.
#' The function adjusts for the zero-based indexing of BED format to the one-based indexing used in GRanges.
#'
#' @param bed_pathA A string specifying the path to the BED file to be converted.
#' @return A GRanges object representing the genomic intervals from the BED file.
#' @keywords internal
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom dplyr mutate
#'
bed_to_granges <- function(bed_path) {
  # Check for a null input and return an empty GRanges object if true
  if (is.null(bed_path)) {
    return(GRanges())
  }

  # Read the BED file with specified column names and delimiter
  # Automatically adjust for BED's zero-based start position

  df <- readr::read_delim(bed_path, delim = "\t", col_names = c("chrom", "start", "end"),
                          show_col_types = FALSE) %>%
    mutate(start = .data$start + 1) # Adjust start positions to 1-based indexing

  # Convert the dataframe to a GRanges object
  grange_from_bed <- makeGRangesFromDataFrame(df, start.field = "start",
                                              end.field = c("end", "stop"))
  # Return the constructed GRanges object
  return(grange_from_bed)
}
