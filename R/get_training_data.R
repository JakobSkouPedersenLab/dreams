#' Extract training data from BAM files
#'
#' @param bam_paths Vector of strings. Paths to \code{.bam} files to extract training data from.
#' @param reference_path String. Path to reference genome fasta file.
#' @param bed_include_path String. Path to bed-file with regions to include. Default is \code{NULL}.
#' @param positions_to_exclude_paths Vector of strings. List of files with positions to exclude from training with length equal to number of samples. Default is \code{NULL}.
#' @param common_positions_to_exclude_paths Vector of strings. List of files with positions to exclude from all samples. Default is \code{NULL}.
#' @param factor Number between 0 and 1. Ratio between negative and positive data. Default is 1.
#' @param mm_rate_max Number between 0 and 1. Maximum mismatch rate in position. Default is 1.
#' @param verbose TODO: Write this
#'
#' @export
#' @return \code{data.frame} with training data for a bam file
#' @seealso [train_dreams_model()] Function for training model
get_training_data_chr_wise <- function(bam_paths,
                                       reference_path,
                                       bed_include_path = NULL,
                                       factor = 1,
                                       common_positions_to_exclude_paths = NULL,
                                       positions_to_exclude_paths = NULL,
                                       mm_rate_max = 1,
                                       verbose = F,
                                       chr = NULL,
                                       batch_size = 32000,
                                       n_reads = NULL) {
  # Check if there is a position exclude path for each bam file
  if ((!is.null(positions_to_exclude_paths) &
    (length(bam_paths) != length(positions_to_exclude_paths)))) {
    stop("Wrong number of exclude paths")
  }

  training_data <- NULL
  info <- NULL

  n_bam_files <- length(bam_paths)


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
      mm_rate_max = mm_rate_max,
      chr = chr,
      n_reads = n_reads
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


#' Extract training data from BAM files
#'
#' @param bam_paths Vector of strings. Paths to \code{.bam} files to extract training data from.
#' @param reference_path String. Path to reference genome fasta file.
#' @param bed_include_path String. Path to bed-file with regions to include. Default is \code{NULL}.
#' @param positions_to_exclude_paths Vector of strings. List of files with positions to exclude from training with length equal to number of samples. Default is \code{NULL}.
#' @param common_positions_to_exclude_paths Vector of strings. List of files with positions to exclude from all samples. Default is \code{NULL}.
#' @param factor Number between 0 and 1. Ratio between negative and positive data. Default is 1.
#' @param mm_rate_max Number between 0 and 1. Maximum mismatch rate in position. Default is 1.
#' @param verbose TODO: Write this
#'
#' @export
#' @return \code{data.frame} with training data for a bam file
#' @seealso [train_dreams_model()] Function for training model
get_training_data <- function(bam_paths,
                              reference_path,
                              bed_include_path = NULL,
                              factor = 1,
                              common_positions_to_exclude_paths = NULL,
                              positions_to_exclude_paths = NULL,
                              mm_rate_max = 1,
                              verbose = F,
                              n_reads = NULL) {
  # Check if there is a position exclude path for each bam file
  if ((!is.null(positions_to_exclude_paths) &
    (length(bam_paths) != length(positions_to_exclude_paths)))) {
    stop("Wrong number of exclude paths")
  }

  training_data <- NULL
  info <- NULL

  n_bam_files <- length(bam_paths)

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
      mm_rate_max = mm_rate_max,
      n_reads = n_reads
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
#' @param positions_to_exclude_paths positions to exclude from training
#' @param factor ratio between negative and positive data
#' @param mm_rate_max maximum mismatch rate in position
#' @keywords internal
#'
#' @return dataframe with training data for a bam file
get_training_data_from_bam <- function(bam_path, reference_path, bed_include_path = NULL,
                                       factor = 1, positions_to_exclude_paths = NULL, mm_rate_max = 1, chr = NULL, n_reads = NULL) {
  bam_df <- load_BAM(bam_path, chr = chr)

  if (!is.null(n_reads)) {
    bam_df <- bam_df %>% sample_n(n_reads)
  }
  print("BAM FILE LOADED")

  # Add genomic positions of mismatches
  mismatch_bam_df <- extract_mismatch_positions(bam_df)

  print("MISMATCHES extract")


  # Add features
  mismatch_positions_df <-
    extract_features_from_bam(
      bam_df = mismatch_bam_df,
      reference_path = reference_path
    )

  print("FEATURES EXTRACTED")
  print(nrow(mismatch_positions_df))

  # Filter mismatches
  mismatches <-
    filter_mismatch_positions(
      read_positions = mismatch_positions_df,
      bam_file = bam_path,
      mm_rate_max = mm_rate_max,
      bed_include_path = bed_include_path,
      positions_to_exclude_paths = positions_to_exclude_paths,
      chr = chr
    )

  print("mismatches")
  print(nrow(mismatches))

  positive_samples <- mismatches$data
  info <- mismatches$info

  n_samples <- nrow(positive_samples) * factor


  # Generate negative samples
  negative_read_positions_df <-
    sample_negative_read_positions(
      bam_df = bam_df,
      n_samples = n_samples
    )

  print("negative")
  print(nrow(negative_read_positions_df))

  # Add features
  negative_samples <-
    extract_features_from_bam(
      bam_df = negative_read_positions_df,
      reference_path = reference_path
    )

  print("negative samples")
  print(nrow(negative_samples))

  info <- info %>% mutate(n_matches = nrow(negative_samples))
  info <- info %>% mutate(beta = nrow(negative_samples) / (info$total_coverage - nrow(positive_samples)))

  output_data <- rbind(positive_samples, negative_samples)

  output_list <- list(
    data = output_data,
    info = info
  )

  return(output_list)
}






#' Title
#'
#' @param read_positions dataframe of read positions
#' @param bam_file bam file path
#' @param mm_rate_max maximum mm_rate for positions
#' @param bed_include_path bed regions to include in training data
#' @param positions_to_exclude_paths positions to exclude from training
#' @keywords internal
#'
#' @return filtered read position dataframe
#'
#' @importFrom readr read_csv

filter_mismatch_positions <- function(read_positions, bam_file, mm_rate_max = 1, bed_include_path = NULL, chr = NULL, positions_to_exclude_paths = NULL) {
  read_positions_filtered <- read_positions %>% filter(obs != "N")

  # Load coverage data

  if (!is.null(bed_include_path)) {
    included_regions_granges <- bed_to_granges(bed_include_path)
  } else if (!is.null(chr)) {
    included_regions_granges <- GRanges(chr, IRanges(1, 536870912))
  } else {
    included_regions_granges <- GRanges(NULL, NULL)
  }


  pp <- Rsamtools::PileupParam(
    max_depth = 250000000, min_base_quality = 13, min_mapq = 0,
    min_nucleotide_depth = 1, min_minor_allele_depth = 0,
    distinguish_strands = FALSE, distinguish_nucleotides = T,
    ignore_query_Ns = TRUE, include_deletions = TRUE, include_insertions = FALSE,
    left_bins = NULL, query_bins = NULL, cycle_bins = NULL
  )

  count_data <- Rsamtools::pileup(bam_file, pileupParam = pp, scanBamParam = ScanBamParam(which = included_regions_granges)) %>%
    dplyr::rename(chr = .data$seqnames, genomic_pos = .data$pos) %>%
    mutate(
      ref = get_reference_seq(
        chr = .data$chr,
        genomic_pos = .data$genomic_pos,
        buffer = 0,
        reference_path = reference_path
      ),
      is_mm = ref != nucleotide,
    ) %>%
    select(-which_label)



  mm_data <- count_data %>%
    group_by(chr, genomic_pos) %>%
    mutate(coverage = sum(count)) %>%
    filter(is_mm) %>%
    group_by(chr, genomic_pos, ref, obs, coverage) %>%
    summarize(n_mismatches = sum(count)) %>%
    mutate(mm_rate = mm_count / coverage)

  high_mismatch_positions <- mm_data %>% filter(mm_rate > mm_rate_max)

  read_positions_filtered <- read_positions_filtered %>% anti_join(high_mismatch_positions)

  coverage_data <- count_data %>%
    anti_join(high_mismatch_positions) %>%
    group_by(chr, genomic_pos) %>%
    summarize(total_coverage = sum(count))


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
    n_mismatches = mm_data$n_mismatches,
    total_coverage = coverage_data$total_coverage
  )


  return(list(
    data = read_positions_filtered,
    info = beta_info
  ))
}

#' Title
#' @param bed_path path to bed-file
#' @return granges object
#' @keywords internal
#' @importFrom GenomicRanges makeGRangesFromDataFrame

bed_to_granges <- function(bed_path) {
  if (is.null(bed_path)) {
    return(GRanges())
  }

  df <- readr::read_delim(bed_path, delim = "\t", col_names = c("chrom", "start", "end"), show_col_types = FALSE) %>% mutate(start = .data$start + 1)

  grange_from_bed <- makeGRangesFromDataFrame(df, start.field = "start", end.field = c("end", "stop"))

  return(grange_from_bed)
}
