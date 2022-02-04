test_that("calc_string_entropy_k_mer", {
  # No entropy examples
  expect_equal(calc_string_entropy_k_mer(s = "AAA"), 0)
  expect_equal(calc_string_entropy_k_mer(s = c("AAA", "TT")), c(0,0))
  expect_equal(calc_string_entropy_k_mer(s = "AAA", k = 1), 0)
  expect_equal(calc_string_entropy_k_mer(s = "AAA", k = 1, alphabet = c("A", "B")), 0)
  expect_equal(calc_string_entropy_k_mer(s = "AAA", k = 2, alphabet = "A"), 0)
  expect_equal(calc_string_entropy_k_mer(s = "ATG", k = 3), 0)

  # Simple examples
  expect_equal(calc_string_entropy_k_mer(s = "AT", k = 1), 0.30103, tolerance = 0.0001)
  expect_equal(calc_string_entropy_k_mer(s = "ATG", k = 1), 0.47712, tolerance = 0.0001)
  expect_equal(calc_string_entropy_k_mer(s = "ATGCGTC", k = 2), 0.77815, tolerance = 0.0001)
})

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
  reference_path <- system.file("extdata", "ref.fasta", package = "dreams")

  # Example 1
  seq <- get_reference_seq(chr = "chr1", genomic_pos = 6, buffer = 5, reference_path = reference_path)
  expect_equal(seq, "AAAAAAAAAAA")

  # Example 2
  seq <- get_reference_seq(chr = "chr1", genomic_pos = 6, buffer = 3, reference_path = reference_path)
  expect_equal(seq, "AAAAAAA")

  # Example 3
  seq <- get_reference_seq(chr = rep("chr1",2), genomic_pos = c(5,6), buffer = 3, reference_path = reference_path)
  expect_equal(seq, c("AAAAAAA", "AAAAAAA"))

})
