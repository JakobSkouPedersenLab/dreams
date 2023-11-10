test_that("expand_cigar() expands CIGAR string into sequence of operations", {
  expect_equal(expand_cigar("2M1I4M2D"), "MMIMMMM" )
  expect_equal(expand_cigar("10M"), "MMMMMMMMMM" )
  expect_equal(expand_cigar("5M"), "MMMMM" )
  expect_equal(expand_cigar("2H4M2H"), "MMMM" )
  expect_equal(expand_cigar("5M1D4M"), "MMMMMDMMMM" )
  expect_equal(expand_cigar("5M1D4I"), "MMMMMDIIII" )
})





test_that("clean_insertions() adjusts the sequence of CIGAR operations by cleaning up insertion operations", {
  expect_equal(clean_insertions("MMIMMMM"), "MIMMMM" )
  expect_equal(clean_insertions("MMMMIIIMMMMMM"), "MMMIMMMMMM" )
  expect_equal(clean_insertions("MMMMM"), "MMMMM" )
  expect_equal(clean_insertions("MMMM"), "MMMM" )
  expect_equal(clean_insertions("MMMMMDMMMM"), "MMMMMDMMMM" )
  expect_equal(clean_insertions("MMMMMDIIII"), "MMMMMD" )
})



