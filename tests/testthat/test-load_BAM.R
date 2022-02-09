test_that("Load empty bam", {
  empty_bam_file <- system.file("extdata", "empty_example.bam", package = "dreams")
  bam_df <- load_BAM(empty_bam_file)

  # Check number of reads
  expect_equal(nrow(bam_df), 0)

  # Check columns
  expected_standard_columns <- c("qname", "rname", "flag", "strand", "pos", "qwidth", "mapq", "cigar", "mpos", "isize", "seq", "qual", "MD", "chr")
  expect_true(all(expected_standard_columns %in% colnames(bam_df)))
})

test_that("Mini example - Load", {
  mini_bam_file <- system.file("extdata", "mini_example.bam", package = "dreams")
  bam_df <- load_BAM(mini_bam_file)

  # Check number of reads
  expect_equal(nrow(bam_df), 4)

  # Check columns
  expected_standard_columns <- c("qname", "rname", "flag", "strand", "pos", "qwidth", "mapq", "cigar", "mpos", "isize", "seq", "qual", "MD", "chr")
  expected_umi_columns <- c("ce", "cd", "cE", "cD")

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

test_that("Mini example (no UMI) - Load", {
  mini_bam_no_umifile <- system.file("extdata", "mini_example_no_umi.bam", package = "dreams")
  bam_df <- load_BAM(mini_bam_no_umifile)

  # Check number of reads
  expect_equal(nrow(bam_df), 4)

  # Check UMI not part of columns
  expected_umi_columns <- c("ce", "cd", "cE", "cD")

  expect_false(any(expected_umi_columns %in% colnames(bam_df)))
})

test_that("Mini example - Load positions", {
  mini_bam_file <- system.file("extdata", "mini_example.bam", package = "dreams")

  # Position with no coverage
  empty_position_bam_df <- load_BAM(mini_bam_file, chr = "chr1", pos = 1)
  expect_equal(nrow(empty_position_bam_df), 0)

  # Position with coverage
  position_w_reads_bam_df <- load_BAM(mini_bam_file, chr = "chr1", pos = 13)
  expect_equal(nrow(position_w_reads_bam_df), 2)

  expect_true(all(position_w_reads_bam_df$chr == "chr1"))
  expect_true(all(position_w_reads_bam_df$genomic_pos == 13))

  # Multiple positions with coverage
  positions_w_reads_bam_df <- load_BAM(mini_bam_file, chr = "chr1", pos = c(13, 14))

  expect_equal(positions_w_reads_bam_df$chr, c("chr1", "chr1", "chr1", "chr1"))
  expect_equal(positions_w_reads_bam_df$genomic_pos, c(13, 13, 14, 14))
})

test_that("Load read (field test)", {
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
