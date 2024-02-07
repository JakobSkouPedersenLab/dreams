test_that("Tests correct_seq()", {
  expect_equal(correct_seq("2M1I1M2D1I1M4I", "AAACACTGGC"), "AICDDI")
  expect_equal(correct_seq("2M1I1M2D1M3I",  "AAACACTG"), "AICDDI")

})



test_that("Tests correct_num()", {
  expect_equal(correct_num("2M1I1M2D1M3I", c(0,0,2,0,0,1,0,2)), c(0,0,0,0,0,0))
  expect_equal(correct_num("2M1I1M2D1M3I",  c(2,3,2,4,5,1,2,3)), c(2,3,4,0,0,5))
  expect_equal(correct_num("2M1I1M2D1I1M4I", c(0,0,2,0,0,1,0,2,1,1)), c(0,0,0,0,0,1))
})
