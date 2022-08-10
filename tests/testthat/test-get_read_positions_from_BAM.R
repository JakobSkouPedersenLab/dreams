test_that("Empty call", {
  read_example_bam_file <- system.file("extdata", "mini_example.bam", package = "dreams")
  reference_path <- system.file("extdata", "ref.fasta", package = "dreams")

  no_coverage_df <- get_read_positions_from_BAM(
    bam_file_path = read_example_bam_file,
    chr = "chr1",
    genomic_pos = 1,
    reference_path = reference_path
  )

  expect_equal(nrow(no_coverage_df), 0)
})

test_that("Simple example", {
  read_example_bam_file <- system.file("extdata", "mini_example.bam", package = "dreams")
  reference_path <- system.file("extdata", "ref.fasta", package = "dreams")

  two_read_positions_df <- get_read_positions_from_BAM(
    bam_file_path = read_example_bam_file,
    chr = "chr1",
    genomic_pos = 13,
    reference_path = reference_path
  )

  expect_equal(nrow(two_read_positions_df), 2)
  expect_equal(two_read_positions_df$obs, c("A", "T"))
})


test_that("Simple example", {
  read_example_bam_file <- system.file("extdata", "mini_example.bam", package = "dreams")
  reference_path <- system.file("extdata", "ref.fasta", package = "dreams")

  two_read_positions_df <- get_read_positions_from_BAM(
    bam_file_path = read_example_bam_file,
    chr = "chr1",
    genomic_pos = 13,
    reference_path = reference_path
  )
  get_training_data_from_bam(bam_path = read_example_bam_file,
                                 reference_path = reference_path)$read_position_mm_rate %>%
    group_by(chr, genomic_pos) %>%
    mutate(count = n_mismatches,
          total_coverage = sum(coverage),
           mm_rate = n_mismatches/total_coverage) %>%
    pivot_wider(id_cols = c(chr,genomic_pos,total_coverage), names_from = "nucleotide", values_from = mm_rate, values_fill = 0)
})
