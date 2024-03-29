test_that("Simple example dreams_vc", {
  read_example_bam_file <- system.file("extdata", "mini_example.bam", package = "dreams")
  reference_path <- system.file("extdata", "ref.fasta", package = "dreams")

  model_path <- system.file("extdata", "model_test.h5", package = "dreams")
  model <- keras::load_model_hdf5(model_path)

  mutations_df <-
    data.frame(
      CHROM = "chr1",
      POS = 10,
      REF = "T",
      alt = "A"
    )


  test <- dreams_vc(
    mutations_df = mutations_df, bam_file_path = read_example_bam_file,
    reference_path = reference_path,
    model = model, alpha = 0.05,
    calculate_confidence_intervals = FALSE
  )
})

test_that("Simple example dreams_vc parallel", {
  read_example_bam_file <- system.file("extdata", "mini_example.bam", package = "dreams")
  reference_path <- system.file("extdata", "ref.fasta", package = "dreams")

  model_path <- system.file("extdata", "model_test.h5", package = "dreams")
  model <- keras::load_model_hdf5(model_path)

  mutations_df <-
    data.frame(
      CHROM = "chr1",
      POS = 10,
      REF = "T",
      ALT = "A"
    )

  test <- dreams_vc_parallel(
    mutations_df = mutations_df,
    bam_file_path = read_example_bam_file,
    reference_path = reference_path,
    model = model,
    alpha = 0.05,
    batch_size = 10,
    calculate_confidence_intervals = FALSE,
    ncores = 1
  )

  test
})


test_that("Simple example dreams_vc - 3 pos", {
  read_example_bam_file <- system.file("extdata", "mini_example.bam", package = "dreams")
  reference_path <- system.file("extdata", "ref.fasta", package = "dreams")

  model_path <- system.file("extdata", "model_test.h5", package = "dreams")
  model <- keras::load_model_hdf5(model_path)

  mutations_df <-
    data.frame(
      CHROM = c("chr1", "chr1", "chr1"),
      POS = c(10, 10, 10),
      REF = c("T", "T", "T"),
      alt = c("A", "G", "C")
    )

  test <- dreams_vc(
    mutations_df = mutations_df,
    bam_file_path = read_example_bam_file,
    reference_path = reference_path, model = model, alpha = 0.05, calculate_confidence_intervals = FALSE
  )
})





# test_that("Simple example dreams_vc - 3 pos - parallel", {
#   read_example_bam_file <- system.file("extdata", "mini_example.bam", package = "dreams")
#   reference_path <- system.file("extdata", "ref.fasta", package = "dreams")
#
#   model_path <- system.file("extdata", "model_test.h5", package = "dreams")
#   model <- keras::load_model_hdf5(model_path)
#
#   mutations_df <-
#     data.frame(
#       CHROM = c("chr1", "chr1", "chr1"),
#       POS = c(10, 10, 10),
#       REF = c("T", "T", "T"),
#       alt = c("A", "G", "C")
#     )
#
#
#   test <- dreams_vc_parallel(
#     mutations_df = mutations_df,
#     bam_file_path = read_example_bam_file,
#     reference_path = reference_path, model = model,
#     beta = 0.01, alpha = 0.05,
#     calculate_confidence_intervals = FALSE, ncores = 8
#   )
#
#   test
#
# })

#
# test_that("Invalid mutations_df", {
#   mutations_no_alt_df <-
#     data.frame(
#       CHROM = "chr1",
#       POS = 10,
#       REF = "T"
#     )
#
#   expect_error(
#     call_mutations(mutations_df = mutations_no_alt_df, read_positions_df = data.frame(), model = NULL, beta = NULL),
#     "mutations_df should have the columns .*"
#   )
# })

test_that("Invalid read_positions_df", {
  mutations_df <-
    data.frame(
      CHROM = "chr1",
      POS = 10,
      REF = "T",
      alt = "A"
    )

  reads_no_pos_df <- data.frame(
    chr = "chr1",
    ref = c("A", "A"),
    obs = c("A", "G")
  )

  expect_error(
    call_mutations(mutations_df = mutations_df, read_positions_df = reads_no_pos_df, model = NULL, beta = NULL),
    "read_positions_df should have the columns \\['chr', genomic_pos', 'ref, 'obs'\\]"
  )
})

test_that("Empty call", {
  empty_mutations_df <-
    data.frame(
      CHROM = numeric(),
      POS = numeric(),
      REF = numeric(),
      ALT = numeric()
    )

  empty_read_positions_df <-
    data.frame(
      chr = numeric(),
      genomic_pos = numeric(),
      ref = numeric(),
      obs = numeric()
    )

  mutation_calls_df <- call_mutations(mutations_df = empty_mutations_df, read_positions_df = empty_read_positions_df, model = NULL, beta = NULL)

  expect_equal(mutation_calls_df, data.frame())
})

test_that("Confidence interval", {
  one_mutation_df <-
    data.frame(
      CHROM = "chr1",
      POS = 10,
      REF = "A",
      ALT = "T"
    )

  # Reads
  some_read_positions_df <-
    data.frame(
      chr = c("chr1", "chr1"),
      genomic_pos =
        c(10, 10),
      ref = "A",
      obs =
        c("A", "T")
    )

  model_path <- system.file("extdata", "model_test.h5", package = "dreams")
  model <- keras::load_model_hdf5(model_path)

  without_ci_df <- call_mutations(mutations_df = one_mutation_df, read_positions_df = some_read_positions_df, model = model, beta = 0.01, alpha = 0.05, calculate_confidence_intervals = FALSE)
  with_ci_df <- call_mutations(mutations_df = one_mutation_df, read_positions_df = some_read_positions_df, model = model, beta = 0.01, alpha = 0.05, calculate_confidence_intervals = TRUE)

  # CI present in result
  expect_true(all(c("tf_min", "tf_max") %in% colnames(without_ci_df)))
  expect_true(all(c("tf_min", "tf_max") %in% colnames(with_ci_df)))

  # CI not calculated
  expect_true(all(is.na(c(without_ci_df$tf_min, without_ci_df$tf_max))))

  # CI is calculated
  expect_true(all(is.numeric(c(with_ci_df$tf_min, with_ci_df$tf_max))))
})

test_that("Only mutants reads", {
  one_mutations_df <-
    data.frame(
      CHROM = "chr1",
      POS = 10,
      REF = "A",
      ALT = "T"
    )

  observed_mutation_read_df <-
    data.frame(
      chr = "chr1",
      genomic_pos = 10,
      ref = "A",
      obs = "T"
    )

  model_path <- system.file("extdata", "model_test.h5", package = "dreams")
  model <- keras::load_model_hdf5(model_path)

  # Expect no error
  expect_error(call_mutations(mutations_df = one_mutations_df, read_positions_df = observed_mutation_read_df, model = model, beta = 0.01), NA)
})

test_that("Small example", {
  three_mutations_df <-
    data.frame(
      CHROM = c("chr1", "chr2", "chr3"),
      POS = c(10, 13, 17),
      REF = "A",
      ALT = "T"
    )

  # Reads on Mut 1,2 - No reads on Mut 3 - And random read
  some_read_positions_df <-
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

  mutation_calls_df <- call_mutations(mutations_df = three_mutations_df, read_positions_df = some_read_positions_df, model = model, beta = 0.01, alpha = 0.05)

  # Mutations info
  expect_equal(nrow(mutation_calls_df), 3)
  expect_equal(mutation_calls_df$count, c(1, 0, 0))
  expect_equal(mutation_calls_df$coverage, c(2, 3, 0))
  expect_equal(mutation_calls_df$full_coverage, c(2, 4, NA))
})


test_that("Bigger example", {
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

  some_read_positions_df <-
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

  res <- call_mutations(
    mutations_df = three_mutations_df, read_positions_df = some_read_positions_df, model = model, beta = 0.01, alpha = 0.05,
    calculate_confidence_intervals = TRUE
  )

  slow_res <- call_mutations(
    mutations_df = three_mutations_df, read_positions_df = some_read_positions_df, model = model, beta = 0.01, alpha = 0.05,
    use_turboem = FALSE,
    calculate_confidence_intervals = TRUE
  )

  expect_equal(slow_res %>% select(tf_est, tf_min, tf_min, p_val), res %>% select(tf_est, tf_min, tf_min, p_val))

  # Snap shots
  # expect_snapshot(res)
  # expect_snapshot(slow_res)
})


test_that("Bigger example with chr/pos wise", {
  three_mutations_df <-
    data.frame(
      CHROM = c("chr1", "chr2", "chr2", "chr3"),
      POS = c(10, 13, 14, 17),
      REF = "A",
      ALT = "T"
    )

  # Reads on Mut 1,2 - No reads on Mut 3 - And random read
  cov_1 <- 17
  count_1 <- 1
  cov_2 <- 43
  count_2 <- 3
  cov_3 <- 30
  count_3 <- 2

  some_read_positions_df <-
    data.frame(
      chr =
        c(
          rep("chr1", cov_1),
          rep("chr2", cov_2),
          rep("chr2", cov_3)
        ),
      genomic_pos =
        c(
          rep(10, cov_1),
          rep(13, cov_2),
          rep(14, cov_3)
        ),
      ref = "A",
      obs =
        c(
          rep("A", cov_1 - count_1), rep("T", count_1),
          rep("A", cov_2 - count_2), rep("T", count_2),
          rep("A", cov_3 - count_3), rep("T", count_3)
        )
    )

  model_path <- system.file("extdata", "model_test.h5", package = "dreams")
  model <- keras::load_model_hdf5(model_path)

  res <- call_mutations(
    mutations_df = three_mutations_df, read_positions_df = some_read_positions_df, model = model, beta = 0.01, alpha = 0.05,
    calculate_confidence_intervals = TRUE
  )

  res_1 <- call_mutations(
    mutations_df = three_mutations_df, read_positions_df = some_read_positions_df, model = model, beta = 0.01, alpha = 0.05,
    calculate_confidence_intervals = TRUE
  )

  res_2 <- call_mutations(
    mutations_df = three_mutations_df, read_positions_df = some_read_positions_df, model = model, beta = 0.01, alpha = 0.05,
    calculate_confidence_intervals = TRUE
  )

  expect_equal(res %>% select(tf_est, tf_min, tf_min, p_val), res_1 %>% select(tf_est, tf_min, tf_min, p_val))
  expect_equal(res %>% select(tf_est, tf_min, tf_min, p_val), res_2 %>% select(tf_est, tf_min, tf_min, p_val))


  # Snap shots
  # expect_snapshot(res)
  # expect_snapshot(slow_res)
})
