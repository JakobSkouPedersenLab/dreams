test_that("expand_cigar() expands CIGAR string into sequence of operations", {
  expect_equal(expand_cigar("2M1I4M"), "MMIMMMM" )
  expect_equal(expand_cigar("10M"), "MMMMMMMMMM" )
  expect_equal(expand_cigar("5M"), "MMMMM" )
  expect_equal(expand_cigar("4M"), "MMMM" )
  expect_equal(expand_cigar("5M1D4M"), "MMMMMDMMMM" )
  expect_equal(expand_cigar("5M1D4I"), "MMMMMDIIII" )
})





