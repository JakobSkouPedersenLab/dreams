test_that("Load mini bam", {
  mini_bam_file <- system.file("extdata", "mini_example.bam", package = "dreams")
  bam_df <- load_BAM(mini_bam_file)

  # Check number of reads
  expect_equal(nrow(bam_df), 4)

  # Check columns
  expected_standard_columns <- c("qname", "rname", "flag", "strand", "pos", "qwidth", "mapq", "cigar", "mpos", "isize", "seq", "qual", "MD", "chr")
  expected_umi_columns <- c("ce", "cd", "cE", "cD", "chr")

  expect_true(all(expected_standard_columns %in% colnames(bam_df)))
  expect_true(all(expected_umi_columns %in% colnames(bam_df)))

  # Check correct read input
  # Read start positions
  expect_equal(bam_df$pos, c(6, 10, 11, 14))
  # Sequence
  expect_equal(bam_df$seq, c(
    "AAAA", "AAAAAA",
    "AAT", "AAA"
  ))
  # MD Tag
  expect_equal(bam_df$MD, c("4", "6", "2A", "3"))
})

test_that("Load empty bam", {
  empty_bam_file <- system.file("extdata", "empty_example.bam", package = "dreams")
  bam_df <- load_BAM(empty_bam_file)

  # Check number of reads
  expect_equal(nrow(bam_df), 0)

  # Check columns
  expected_standard_columns <- c("qname", "rname", "flag", "strand", "pos", "qwidth", "mapq", "cigar", "mpos", "isize", "seq", "qual", "MD", "chr")
  expect_true(all(expected_standard_columns %in% colnames(bam_df)))
})


test_that("Load example bam", {
  read_example_bam_file <- system.file("extdata", "read_example.bam", package = "dreams")
  bam_df <- load_BAM(read_example_bam_file)

  # Check number of reads
  expect_equal(nrow(bam_df), 2)

  # Check columns
  expected_standard_columns <- c("qname", "rname", "flag", "strand", "pos", "qwidth", "mapq", "cigar", "mpos", "isize", "seq", "qual", "MD", "chr")
  expected_umi_columns <- c("ce", "cd", "cE", "cD", "chr")

  expect_true(all(expected_standard_columns %in% colnames(bam_df)))
  expect_true(all(expected_umi_columns %in% colnames(bam_df)))
})
