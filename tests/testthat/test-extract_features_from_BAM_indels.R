test_that("Tests get_indels_start_positions()", {
  expect_equal(correct_seq("2M1I1M2D1I1M4I", "AAACACTGGC"), "AICDDI")
  expect_equal(correct_seq("2M1I1M2D1M3I",  "AAACACTG"), "AICDDI")

})
