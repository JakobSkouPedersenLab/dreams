


test_that("simple examples 2", {
  # Example 1

  read_example_bam_file <- system.file("extdata", "mini_example.bam", package = "dreams")
  reference_path <- system.file("extdata", "ref.fasta", package = "dreams")

  bam_df <- load_BAM(read_example_bam_file, chr = "chr1", pos = 11)
  extract_features_from_bam(bam_df, reference_path = reference_path)

  bam_df <- load_BAM(read_example_bam_file, chr = "chr1", pos = 12)
  extract_features_from_bam(bam_df, reference_path = reference_path)

  bam_df <- load_BAM(read_example_bam_file, chr = "chr1", pos = 13)
  extract_features_from_bam(bam_df, reference_path = reference_path)


  bam_df <- load_BAM(read_example_bam_file, chr = c("chr1","chr1","chr1"), pos = c(11,12,13))
  extract_features_from_bam(bam_df, reference_path = reference_path)


})



test_that("get_reference_seq", {
  # Example 1

  reference_path <- system.file("extdata", "ref.fasta", package = "dreams")
  seq <- get_reference_seq(chr = "chr1", genomic_pos = 6, buffer = 5, reference_path = reference_path)

  expect_equal(seq, c("chr1" = "AAAAAAAAAAA"))

})
