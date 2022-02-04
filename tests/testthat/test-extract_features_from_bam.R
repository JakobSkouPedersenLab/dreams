test_that("simple examples 2", {
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

  bam_df <- load_BAM(read_example_bam_file, chr = "ref", pos = 3)

  extract_features_from_bam(bam_df, reference_path = reference_path)
})

test_that("Choose features", {
  read_example_bam_file <- system.file("extdata", "mini_example.bam", package = "dreams")
  reference_path <- system.file("extdata", "ref.fasta", package = "dreams")

  bam_df <- load_BAM(read_example_bam_file, chr = "chr1", pos = 6)

  extract_features_from_bam(bam_df, reference_path = reference_path)
})
