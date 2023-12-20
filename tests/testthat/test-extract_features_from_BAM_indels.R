test_that("Tests get_indels_start_positions()", {
  expect_equal(correct_seq("MMIMDDIMIII", "AAACACTGG"), "AICDDI")
  expect_equal(correct_seq("MMIMDDMIII",  "AAACACTG"), "AICDDI")

})
