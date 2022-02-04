test_that("Load simple bam", {
  small_sam_file <- system.file("extdata", "small_sam_file.bam", package="dreams")
  load_BAM(small_sam_file)
})
