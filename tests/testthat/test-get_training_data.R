# test_that("simple examples 2", {
#   # Example 1
#   read_example_bam_file <- system.file("extdata", "mini_example.bam", package = "dreams")
#   reference_path <- system.file("extdata", "ref.fasta", package = "dreams")
#
#   bam_df <- load_BAM(read_example_bam_file)
#
#   # Add genomic positions of mismatches
#   mismatch_bam_df <- extract_mismatch_positions(bam_df)
#
#   mismatch_positions_df <-
#     extract_features_from_bam(
#       bam_df = mismatch_bam_df,
#       reference_path = reference_path
#     )
#
#   # Filter data
#   filtered_mismatch_positions_df <-
#     filter_mismatch_positions(
#       read_positions = mismatch_positions_df,
#       bam_file = read_example_bam_file,
#       mm_rate_max = 0.49,
#       bed_include_path = NULL
#     )
# })
#
#
# test_that("Individual example", {
#   read_example_bam_file <- system.file("extdata", "mini_example.bam", package = "dreams")
#   reference_path <- system.file("extdata", "ref.fasta", package = "dreams")
#   bed_file_3 <- system.file("extdata", "bed_cov3.bed", package = "dreams")
#
#   bam_df <- load_BAM(read_example_bam_file)
#
#
#   # Add genomic positions of mismatches
#   mismatch_bam_df <- extract_mismatch_positions(bam_df)
#
#   # Add features
#   mismatch_positions_df <-
#     extract_features_from_bam(
#       bam_df = mismatch_bam_df,
#       reference_path = reference_path
#     )
#
#   filtered_mismatch_positions_df <-
#     filter_mismatch_positions(
#       read_positions = mismatch_positions_df,
#       bam_file = read_example_bam_file,
#       mm_rate_max = 1,
#       bed_include_path = bed_file_3
#     )
#
#   n_samples <- nrow(filtered_mismatch_positions_df$data) * 1
#
#   # Generate negative samples
#   negative_read_positions_df <-
#     sample_negative_read_positions(
#       bam_df = bam_df,
#       n_samples = n_samples
#     )
#
#   # Add features
#   negative_samples <-
#     extract_features_from_bam(
#       bam_df = negative_read_positions_df,
#       reference_path = reference_path
#     )
# })
#




test_that("full example - get_training_data_from_bam", {
  read_example_bam_file <- system.file("extdata", "mini_example.bam", package = "dreams")
  reference_path <- system.file("extdata", "ref.fasta", package = "dreams")

  n_errors <- 1
  factor <- 1

  samples_1 <- get_training_data_from_bam(
    bam_path = read_example_bam_file,
    reference_path = reference_path,
    bed_include_path = NULL,
    factor = factor,
    mm_rate_max = 0.51
  )

  expect_true(nrow(samples_1$data) == n_errors + n_errors * factor)

  n_errors <- 1
  factor <- 10

  samples_2 <- get_training_data_from_bam(
    bam_path = read_example_bam_file,
    reference_path = reference_path,
    bed_include_path = NULL,
    factor = factor,
    mm_rate_max = 0.51
  )
  expect_true(nrow(samples_2$data) == n_errors + n_errors * factor)


  n_errors <- 1
  factor <- 1
  positions_to_exclude_paths <- c(system.file("extdata", "positions_to_exclude_1.csv", package = "dreams"))


  samples_3 <- get_training_data_from_bam(
    bam_path = read_example_bam_file,
    reference_path = reference_path,
    bed_include_path = NULL,
    factor = factor,
    mm_rate_max = 0.51,
    positions_to_exclude_paths = positions_to_exclude_paths
  )

  n_errors <- 1
  factor <- 1
  positions_to_exclude_paths_4 <- c(system.file("extdata", "positions_to_exclude_2.csv", package = "dreams"))


  samples_4 <- get_training_data_from_bam(
    bam_path = read_example_bam_file,
    reference_path = reference_path,
    bed_include_path = NULL,
    factor = factor,
    mm_rate_max = 0.51,
    positions_to_exclude_paths = positions_to_exclude_paths_4
  )
})



test_that("pileup example", {
  # Example 1
  read_example_bam_file <- system.file("extdata", "mini_example.bam", package = "dreams")
  Rsamtools::pileup(read_example_bam_file)

  pp <- Rsamtools::PileupParam(
    max_depth = 250, min_base_quality = 13, min_mapq = 0,
    min_nucleotide_depth = 1, min_minor_allele_depth = 0,
    distinguish_strands = FALSE, distinguish_nucleotides = FALSE,
    ignore_query_Ns = TRUE, include_deletions = TRUE, include_insertions = FALSE,
    left_bins = NULL, query_bins = NULL, cycle_bins = NULL
  )

  coverage_all <- Rsamtools::pileup(read_example_bam_file, pileupParam = pp)
  expect_true(nrow(coverage_all) == 11)

  ranges_NULL <- bed_to_granges(NULL)
  coverage_NULL <- Rsamtools::pileup(read_example_bam_file, pileupParam = pp, scanBamParam = ScanBamParam(which = ranges_NULL))

  # FILTER WITH BED

  bed_file_0 <- system.file("extdata", "bed_cov0.bed", package = "dreams")
  bed_file_1 <- system.file("extdata", "bed_cov1.bed", package = "dreams")
  bed_file_2 <- system.file("extdata", "bed_cov2.bed", package = "dreams")
  bed_file_3 <- system.file("extdata", "bed_cov2.bed", package = "dreams")

  ranges_0 <- bed_to_granges(bed_file_0)
  coverage_0 <- Rsamtools::pileup(read_example_bam_file, pileupParam = pp, scanBamParam = ScanBamParam(which = ranges_0))

  ranges_1 <- bed_to_granges(bed_file_1)
  coverage_1 <- Rsamtools::pileup(read_example_bam_file, pileupParam = pp, scanBamParam = ScanBamParam(which = ranges_1))
  expect_true(coverage_1$count == 1)

  ranges_2 <- bed_to_granges(bed_file_2)
  coverage_2 <- Rsamtools::pileup(read_example_bam_file, pileupParam = pp, scanBamParam = ScanBamParam(which = ranges_2))
  expect_true(coverage_2$count == 2)

  ranges_3 <- bed_to_granges(bed_file_3)
  coverage_3 <- Rsamtools::pileup(read_example_bam_file, pileupParam = pp, scanBamParam = ScanBamParam(which = ranges_3))
  expect_true(coverage_3$count == 2)
})



test_that("full example - get_training_data_from_bam", {
  read_example_bam_file <- system.file("extdata", "mini_example.bam", package = "dreams")
  reference_path <- system.file("extdata", "ref.fasta", package = "dreams")

  factor <- 1

  # Basic

  samples_1 <- get_training_data(
    bam_paths = read_example_bam_file,
    reference_path = reference_path,
    bed_include_path = NULL,
    factor = factor,
    mm_rate_max = 0.51
  )

  samples_2 <- get_training_data(
    bam_paths = c(read_example_bam_file, read_example_bam_file),
    reference_path = reference_path,
    bed_include_path = NULL,
    factor = factor,
    mm_rate_max = 0.51
  )

  # Personal positions to exclude

  positions_to_exclude_paths <- c(
    system.file("extdata", "positions_to_exclude_1.csv", package = "dreams"),
    system.file("extdata", "positions_to_exclude_1.csv", package = "dreams"),
    system.file("extdata", "positions_to_exclude_2.csv", package = "dreams"),
    system.file("extdata", "positions_to_exclude_2.csv", package = "dreams")
  )

  samples_3 <- get_training_data(
    bam_paths = c(read_example_bam_file, read_example_bam_file, read_example_bam_file, read_example_bam_file),
    reference_path = reference_path,
    bed_include_path = NULL,
    factor = factor,
    mm_rate_max = 0.51,
    positions_to_exclude_paths = positions_to_exclude_paths
  )

  common_positions_to_exclude_paths <- c(system.file("extdata", "positions_to_exclude_1.csv", package = "dreams"))

  positions_to_exclude_paths <- c(
    system.file("extdata", "positions_to_exclude_2.csv", package = "dreams"),
    system.file("extdata", "positions_to_exclude_2.csv", package = "dreams")
  )

  samples_3 <- get_training_data(
    bam_paths = c(read_example_bam_file, read_example_bam_file),
    reference_path = reference_path,
    bed_include_path = NULL,
    factor = factor,
    mm_rate_max = 0.51,
    positions_to_exclude_paths = positions_to_exclude_paths,
    common_positions_to_exclude_paths = common_positions_to_exclude_paths
  )

  common_positions_to_exclude_paths <- c(system.file("extdata", "positions_to_exclude_2.csv", package = "dreams"))

  positions_to_exclude_paths <- c(
    system.file("extdata", "positions_to_exclude_2.csv", package = "dreams"),
    system.file("extdata", "positions_to_exclude_1.csv", package = "dreams")
  )

  bed_include_path <- system.file("extdata", "bed_cov3.bed", package = "dreams")

  samples_4 <- get_training_data(
    bam_paths = c(read_example_bam_file, read_example_bam_file),
    reference_path = reference_path,
    bed_include_path = bed_include_path,
    factor = factor,
    mm_rate_max = 0.51,
    positions_to_exclude_paths = positions_to_exclude_paths,
    common_positions_to_exclude_paths = common_positions_to_exclude_paths
  )
})
