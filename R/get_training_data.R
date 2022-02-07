# generate_training_samples <- function(bam_paths, reference_path, coverage_data_path, output_path, bed_include_path = NULL, factor = 1) {
#
#   bam_df = load_BAM(bam_path)
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
#   # Filter data
#   filtered_mismatch_positions_df <-
#     filter_read_positions(
#       read_positions = mismatch_positions_df,
#       coverage_data_path = coverage_data_path,
#       bed_include_path = bed_include_path
#     )
#
#   n_samples = nrow(filtered_mismatch_positions_df) * factor
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
#
#   write.csv()
# }
#
