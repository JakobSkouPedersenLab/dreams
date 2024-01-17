test_that("expand_cigar() expands CIGAR string into sequence of operations", {
  expect_equal(expand_cigar("2M1I4M"), "MMIMMMM" )
  expect_equal(expand_cigar("10M"), "MMMMMMMMMM" )
  expect_equal(expand_cigar("5M"), "MMMMM" )
  expect_equal(expand_cigar("4M"), "MMMM" )
  expect_equal(expand_cigar("5M1D4M"), "MMMMMDMMMM" )
  expect_equal(expand_cigar("5M1D4I"), "MMMMMDIIII" )
})

test_that("Tests clean_insertions()", {
  expect_equal(clean_insertions("MMIMMMM"), "MIMMMM" )
  expect_equal(clean_insertions("MMMMM"), "MMMMM" )
  expect_equal(clean_insertions("MMMMMDMMMM"), "MMMMMDMMMM" )
  expect_equal(clean_insertions("MMMMMDIIII"), "MMMMMD" )
  expect_equal(clean_insertions("MMMMMDIIIID"), "MMMMMDD" )
  expect_equal(clean_insertions("MMMMMIIIIDD"), "MMMMIDD" )
})

test_that("Tests convert_to_cigar()", {
  expect_equal(convert_to_cigar("MIMMMM"), "1M1I4M" )
  expect_equal(convert_to_cigar("MMMMM"), "5M" )
  expect_equal(convert_to_cigar("MMMMMDMMMM"), "5M1D4M" )
  expect_equal(convert_to_cigar("MMMMMD"), "5M1D" )
  expect_equal(convert_to_cigar("MMMMMDDI"), "5M2D1I" )
  expect_equal(convert_to_cigar("MMMMIDD"), "4M1I2D" )
})

test_that("Tests get_indels_start_positions()", {
  expect_equal(get_indels_start_positions("2M1I4M"), c(3))
  expect_equal(get_indels_start_positions("10M"),  numeric(0))
  expect_equal(get_indels_start_positions("5M"), numeric(0) )
  expect_equal(get_indels_start_positions("4M"), numeric(0))
  expect_equal(get_indels_start_positions("5M1D4M2I"), c(6,11) )
  expect_equal(get_indels_start_positions("5M1D4I4D"), c(6,7,11) )
})


