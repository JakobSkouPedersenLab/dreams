test_that("Invalid mutations_df", {
  mutations_no_alt_df <-
    data.frame(
      CHROM = "chr1",
      POS = 10,
      REF = "T"
    )

  expect_error(
    call_cancer(mutations_df = mutations_no_alt_df, reads_df = data.frame(), model = NULL, beta = NULL),
    "mutations_df should have the columns .*"
  )
})

test_that("Invalid reads_df", {
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
    call_cancer(mutations_df = mutations_df, reads_df = reads_no_pos_df, model = NULL, beta = NULL),
    "reads_df should have the columns \\['chr', genomic_pos', 'ref, 'obs'\\]"
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

  empty_reads_df <-
    data.frame(
      chr = numeric(),
      genomic_pos = numeric(),
      ref = numeric(),
      obs = numeric()
    )

  res <- call_cancer(mutations_df = empty_mutations_df, reads_df = empty_reads_df, model = NULL, beta = NULL)

  expect_equal(res$mutation_info, data.frame())
  expect_equal(res$cancer_info$mutations_tested, 0)
  expect_equal(res$cancer_info$total_coverage, 0)
  expect_equal(res$cancer_info$p_val, 1)
  expect_equal(res$cancer_info$EM_converged, NA)
})

test_that("Equal columns in simple and empty call", {
  # Empty call
  empty_mutations_df <-
    data.frame(
      CHROM = numeric(),
      POS = numeric(),
      REF = numeric(),
      ALT = numeric()
    )

  empty_reads_df <-
    data.frame(
      chr = numeric(),
      genomic_pos = numeric(),
      ref = numeric(),
      obs = numeric()
    )

  empty_res <- call_cancer(mutations_df = empty_mutations_df, reads_df = empty_reads_df, model = NULL, beta = NULL)

  # One mutation with no count
  one_mutations_df <-
    data.frame(
      CHROM = "chr1",
      POS = 10,
      REF = "A",
      ALT = "T"
    )

  no_count_reads_df <-
    data.frame(
      chr = "chr1",
      genomic_pos = 10,
      ref = "A",
      obs = "A"
    )


  model_path <- system.file("extdata", "model_test.h5", package = "dreams")
  model <- keras::load_model_hdf5(model_path)

  no_count_res <- call_cancer(mutations_df = one_mutations_df, reads_df = no_count_reads_df, model = model, beta = 0.01)

  # One mutation with one count
  one_count_reads_df <-
    data.frame(
      chr = "chr1",
      genomic_pos = 10,
      ref = "A",
      obs = c("A", "T")
    )

  one_count_res <- call_cancer(mutations_df = one_mutations_df, reads_df = one_count_reads_df, model = model, beta = 0.01)


  expect_equal(colnames(empty_res$cancer_info), colnames(no_count_res$cancer_info))
  expect_equal(colnames(empty_res$cancer_info), colnames(one_count_res$cancer_info))
  expect_equal(colnames(no_count_res$cancer_info), colnames(one_count_res$cancer_info))

  expect_equal(colnames(no_count_res$mutation_info), colnames(one_count_res$mutation_info))
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
  expect_error(call_cancer(mutations_df = one_mutations_df, reads_df = observed_mutation_read_df, model = model, beta = 1), NA)

})

test_that("One mutation - no reads", {
  one_mutations_df <-
    data.frame(
      CHROM = "chr1",
      POS = 10,
      REF = "A",
      ALT = "T"
    )

  empty_reads_df <-
    data.frame(
      chr = numeric(),
      genomic_pos = numeric(),
      ref = numeric(),
      obs = numeric()
    )

  model_path <- system.file("extdata", "model_test.h5", package = "dreams")
  model <- keras::load_model_hdf5(model_path)

  res <- call_cancer(mutations_df = one_mutations_df, reads_df = empty_reads_df, model = model, beta = 1)

  # Cancer info
  expect_equal(res$cancer_info$mutations_tested, 1)
  expect_equal(res$cancer_info$total_coverage, 0)
  expect_equal(res$cancer_info$p_val, 1)
  expect_equal(res$cancer_info$EM_converged, NA)

  # Mutations info
  expect_equal(nrow(res$mutation_info), 1)
  expect_equal(res$mutation_info$chr, "chr1")
  expect_equal(res$mutation_info$obs_freq, 0)
  expect_equal(res$mutation_info$P_mut_is_present, NA)
  expect_equal(res$mutation_info$count, 0)
})

test_that("Two mutations - one read in anptother position", {
  two_mutations_df <-
    data.frame(
      CHROM = c("chr1", "chr2"),
      POS = c(10, 13),
      REF = "A",
      ALT = "T"
    )

  some_reads_df <-
    data.frame(
      chr = "chr3",
      genomic_pos = 10,
      ref = "A",
      obs = "T"
    )

  model_path <- system.file("extdata", "model_test.h5", package = "dreams")
  model <- keras::load_model_hdf5(model_path)

  res <- call_cancer(mutations_df = two_mutations_df, reads_df = some_reads_df, model = model, beta = 1)

  # Cancer info
  expect_equal(res$cancer_info$mutations_tested, 2)
  expect_equal(res$cancer_info$total_coverage, 0)
  expect_equal(res$cancer_info$p_val, 1)
  expect_equal(res$cancer_info$EM_converged, NA)

  # Mutations info
  expect_equal(nrow(res$mutation_info), 2)
  expect_equal(res$mutation_info$chr, c("chr1", "chr2"))
  expect_equal(res$mutation_info$obs_freq, c(0, 0))
  expect_equal(res$mutation_info$P_mut_is_present, c(NA, NA))
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
  some_reads_df <-
    data.frame(
      chr =
        c(
          "chr1", "chr1",
          "chr2", "chr2", "chr2",
          "chr5"
        ),
      genomic_pos =
        c(
          10, 10,
          13, 13, 13,
          19
        ),
      ref = "A",
      obs =
        c(
          "A", "T",
          "A", "A", "A",
          "A"
        )
    )

  model_path <- system.file("extdata", "model_test.h5", package = "dreams")
  model <- keras::load_model_hdf5(model_path)

  res <- call_cancer(mutations_df = three_mutations_df, reads_df = some_reads_df, model = model, beta = 0.01)

  # Cancer info
  expect_equal(res$cancer_info$mutations_tested, 3)
  expect_equal(res$cancer_info$total_coverage, 5)
  expect_true(res$cancer_info$EM_converged)

  # Mutations info
  expect_equal(nrow(res$mutation_info), 3)
  expect_equal(res$mutation_info$count, c(1, 0, 0))
  expect_equal(res$mutation_info$coverage, c(2, 3, 0))
})

test_that("Mutation with no data should not change output", {
  one_mutations_df <-
    data.frame(
      CHROM = "chr1",
      POS = 10,
      REF = "A",
      ALT = "T"
    )

  two_mutations_df <-
    bind_rows(
      one_mutations_df,
      data.frame(
        CHROM = "chr2",
        POS = 13,
        REF = "A",
        ALT = "T"
      )
    )

  # Reads on Mut 1,2 - No reads on Mut 3 - And random read
  some_reads_df <-
    data.frame(
      chr = "chr1",
      genomic_pos = 10,
      ref = "A",
      obs =
        c(
          "T", "A", "A", "A", "A", "A"
        )
    )

  model_path <- system.file("extdata", "model_test.h5", package = "dreams")
  model <- keras::load_model_hdf5(model_path)

  res_one_mutations <- call_cancer(mutations_df = one_mutations_df, reads_df = some_reads_df, model = model, beta = 0.01)
  res_two_mutations <- call_cancer(mutations_df = two_mutations_df, reads_df = some_reads_df, model = model, beta = 0.01)

  # Is cancer results the same?
  expect_equal(res_one_mutations$cancer_info$tf_est, res_two_mutations$cancer_info$tf_est)
  expect_equal(res_one_mutations$cancer_info$r_est, res_two_mutations$cancer_info$r_est)
  expect_equal(res_one_mutations$cancer_info$p_val, res_two_mutations$cancer_info$p_val)
})
