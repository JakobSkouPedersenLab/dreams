test_that("test ...", {
  mutations_df <-
    data.frame(
      CHROM = "chr1",
      POS = 10,
      REF = "T",
      alt = "A"
    )

  call_cancer(mutations_df = mutations_df, reads = data.frame(), model = NULL, beta = NULL)
})

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
