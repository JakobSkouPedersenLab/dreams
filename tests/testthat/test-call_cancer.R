test_that("Invalid mutations_df", {
  mutations_no_alt_df <-
    data.frame(
      CHROM = "chr1",
      POS = 10,
      REF = "T"
    )

  expect_error(
    call_cancer(mutations_df = mutations_no_alt_df, reads = data.frame(), model = NULL, beta = NULL),
    "mutations_df should have the columns .*"
  )
})

test_that("Empty mutations_df", {
  empty_mutations_df <-
    data.frame(
      CHROM = numeric(),
      POS = numeric(),
      REF = numeric(),
      ALT = numeric()
    )

  res = call_cancer(mutations_df = empty_mutations_df, reads = data.frame(), model = NULL, beta = NULL)

  expect_null(res$mutation_info)
  expect_equal(res$cancer_info$mutations_tested, 0)
  expect_equal(res$cancer_info$total_coverage, 0)
  expect_equal(res$cancer_info$p_val, 1)
  expect_equal(res$cancer_info$EM_converged, NA)
})

test_that("test ...", {
  mutations_df <-
    data.frame(
      CHROM = "chr1",
      POS = 10,
      REF = "T",
      alt = "A"
    )

  reads_df = data.frame(
    ref = c("A","A"),
    obs = c("A","G")
  )

  model_path = system.file("extdata", "model_test.h5", package = "dreams")

  model = keras::load_model_hdf5(model_path)

  expect_equal(names(call_cancer(mutations_df = mutations_df, reads = reads_df, model = model, beta = 1)), c("cancer_info", "mutation_info"))

})
