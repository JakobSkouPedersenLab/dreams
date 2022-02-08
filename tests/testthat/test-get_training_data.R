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
#       mm_rate_max = 0.49,
#       bed_include_path = NULL
#     )
#
#   n_samples <- nrow(filtered_mismatch_positions_df) * 1
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
#
#
#
#
# test_that("full example", {
#   read_example_bam_file <- system.file("extdata", "mini_example.bam", package = "dreams")
#   reference_path <- system.file("extdata", "ref.fasta", package = "dreams")
#
#   samples <- generate_training_samples(
#     bam_path = read_example_bam_file,
#     reference_path = reference_path,
#     bed_include_path = NULL,
#     factor = 1,
#     mm_rate_max = 0.51
#   )
# })
