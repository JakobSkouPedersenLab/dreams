test_that("simple examples 2", {
  # Example 1

  read_example_bam_file <- system.file("extdata", "mini_example.bam", package = "dreams")
  reference_path <- system.file("extdata", "ref.fasta", package = "dreams")

  bam_df = load_BAM(read_example_bam_file, chr = "ref", pos = 3)

  extract_features_from_bam(bam_df, reference_path = reference_path)

})
