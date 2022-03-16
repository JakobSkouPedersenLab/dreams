test_that("Small example", {
  three_mutations_df <-
    data.frame(
      CHROM = c("chr1", "chr2", "chr3"),
      POS = c(10, 13, 17),
      REF = "A",
      ALT = "T"
    )

  # Reads on Mut 1,2 - No reads on Mut 3 - And random read
  some_reads_df <-
    data.frame(
      chr =
        c(
          "chr1", "chr1",
          "chr2", "chr2", "chr2", "chr2",
          "chr5"
        ),
      genomic_pos =
        c(
          10, 10,
          13, 13, 13, 13,
          19
        ),
      ref = "A",
      obs =
        c(
          "A", "T",
          "A", "A", "A", "G",
          "A"
        )
    )

  model_path <- system.file("extdata", "model_test.h5", package = "dreams")
  model <- keras::load_model_hdf5(model_path)

  mutation_calls_df <- call_mutations(mutations_df = three_mutations_df, reads_df = some_reads_df, model = model, beta = 0.01, alpha = 0.05)

  # Mutations info
  expect_equal(nrow(mutation_calls_df), 3)
  expect_equal(mutation_calls_df$count, c(1, 0, 0))
  expect_equal(mutation_calls_df$coverage, c(2, 3, 0))
  expect_equal(mutation_calls_df$full_coverage, c(2, 4, NA))
})


test_that("Snapshot", {
  three_mutations_df <-
    data.frame(
      CHROM = c("chr1", "chr2", "chr3"),
      POS = c(10, 13, 17),
      REF = "A",
      ALT = "T"
    )

  # Reads on Mut 1,2 - No reads on Mut 3 - And random read
  cov_1 <- 17
  count_1 <- 1
  cov_2 <- 43
  count_2 <- 3

  some_reads_df <-
    data.frame(
      chr =
        c(
          rep("chr1", cov_1),
          rep("chr2", cov_2)
        ),
      genomic_pos =
        c(
          rep(10, cov_1),
          rep(13, cov_2)
        ),
      ref = "A",
      obs =
        c(
          rep("A", cov_1 - count_1), rep("T", count_1),
          rep("A", cov_2 - count_2), rep("T", count_2)
        )
    )

  model_path <- system.file("extdata", "model_test.h5", package = "dreams")
  model <- keras::load_model_hdf5(model_path)

  res <- call_mutations(mutations_df = three_mutations_df, reads_df = some_reads_df, model = model, beta = 0.01, alpha = 0.05)
  expect_snapshot(res)

  slow_res <- call_mutations(mutations_df = three_mutations_df, reads_df = some_reads_df, model = model, beta = 0.01, alpha = 0.05, use_warp_speed = FALSE)
  expect_snapshot(slow_res)

  expect_equal(slow_res %>% select(tf_est, p_val), res %>% select(tf_est, p_val))
})
