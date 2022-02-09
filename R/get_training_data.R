#' Title
#'
#' @param bam_path Path to BAM file
#' @param reference_path Path to reference file
#' @param bed_include_path BED regions to include
#' @param positions_to_exclude_path positions to exclude from training
#' @param factor ratio between negative and position data
#' @param mm_rate_max maximum mismatch rate in position
#'
#' @return dataframe with training data for a bam file
#' @export
#'
#' @examples
generate_training_samples <- function(bam_path, reference_path, bed_include_path = NULL, factor = 1, positions_to_exclude_path = NULL, mm_rate_max = 1) {
  bam_df <- load_BAM(bam_path)

  # Add genomic positions of mismatches
  mismatch_bam_df <- extract_mismatch_positions(bam_df)

  # Add features
  mismatch_positions_df <-
    extract_features_from_bam(
      bam_df = mismatch_bam_df,
      reference_path = reference_path
    )

  positive_samples <-
    filter_mismatch_positions(
      read_positions = mismatch_positions_df,
      bam_file = bam_path,
      mm_rate_max = mm_rate_max,
      bed_include_path = bed_include_path,
      positions_to_exclude_path = positions_to_exclude_path
    )



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


  output_data <- rbind(positive_samples, negative_samples)

  return(output_data)
}






#' Title
#'
#' @param read_positions dataframe of read positions
#' @param bam_file bam file path
#' @param mm_rate_max maximum mm_rate for positions
#' @param bed_include_path bed regions to include in training data
#' @param positions_to_exclude_path positions to exclude from training
#'
#' @return filtered read position dataframe
#'
#' @importFrom readr read_csv

filter_mismatch_positions <- function(read_positions, bam_file, mm_rate_max = 1, bed_include_path = NULL, positions_to_exclude_path = NULL) {
  read_positions_filtered <-
    read_positions %>%
    filter(.data$obs != "N")

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

  read_position_filter <- read_positions_summarized %>%
    inner_join(coverage_data, by = c("chr", "genomic_pos")) %>%
    mutate(mm_rate = .data$n_mismatches / .data$coverage) %>%
    filter(.data$mm_rate < mm_rate_max)

  read_positions_filtered <- read_positions_filtered %>%
    semi_join(read_position_filter, by = c("chr", "genomic_pos"))

  # Remove unwanted positions

  if (!is.null(positions_to_exclude_path)) {
    positions_to_exclude <- read_csv(positions_to_exclude_path)

    read_positions_filtered %>%
      anti_join(positions_to_exclude, by = c("chr", "genomic_pos"))
  }

  return(read_positions_filtered)
}

#' Title
#'
#' @param bed_path
#'
#' @return
#' @importFrom GenomicRanges makeGRangesFromDataFrame

bed_to_granges <- function(bed_path) {
  if (is.null(bed_path)) {
    return(GRanges())
  }

  df <- readr::read_delim(bed_path, delim = "\t", col_names = c("chrom", "start", "end"), show_col_types = FALSE) %>% mutate(start = start + 1)

  grange_from_bed <- makeGRangesFromDataFrame(df, start.field = "start", end.field = c("end", "stop"))

  return(grange_from_bed)
}
