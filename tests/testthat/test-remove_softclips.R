
test_that("Simple example", {

  input_df = tibble::tibble(seq = "ATCG",
                  qual = "FFFF",
                  cigar = "1S3M",
                  cd = list(c(3,3,3,3)),
                  ce = list(c(1,0,0,0)))


  output_df = tibble::tibble(seq = "TCG",
                        qual = "FFF",
                        cigar = "3M",
                        cd = list(c(3,3,3)),
                        ce = list(c(0,0,0)))

  clipped_df = remove_softclips(input_df)

  expect_equal(clipped_df,expected = output_df)
})


test_that("Simple example - no UMI", {

  input_df = tibble::tibble(seq = "ATCG",
                            qual = "FFFF",
                            cigar = "1S3M")


  output_df = tibble::tibble(seq = "TCG",
                             qual = "FFF",
                             cigar = "3M")

  clipped_df = remove_softclips(input_df)

  expect_equal(clipped_df,expected = output_df)
})
