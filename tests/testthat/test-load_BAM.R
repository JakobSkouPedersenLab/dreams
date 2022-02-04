test_that("Load simple bam", {
  small_sam_file <- system.file("extdata", "small_sam_file.bam", package = "dreams")
  bam_df = load_BAM(small_sam_file)

  # Check number of reads
  expect_equal(nrow(bam_df), 2)

  # Check columns
  expected_standard_columns = c("qname","rname", "flag","strand","pos","qwidth","mapq","cigar","mpos","isize","seq","qual","MD","chr")
  expected_umi_columns = c("ce","cd","cE","cD","chr")

  expect_true(all(expected_standard_columns %in% colnames(bam_df)))
  expect_true(all(expected_umi_columns %in% colnames(bam_df)))
})
