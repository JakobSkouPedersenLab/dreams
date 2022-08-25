
test_that("simple examples", {
  no_correction_test <- tibble(
    seq = c("OX", "XO"),
    cigar = c("2M", "2M"),
    pos = 1,
    pos_idx = c(1, 2),
    exp_obs = "O"
  )

  simple_insert_test <- tibble(
    seq = c("OIIX", "XIIO"),
    cigar = c("1M2I1M", "1M2I1M"),
    pos = 1,
    pos_idx = c(1, 2),
    exp_obs = "O"
  )

  adv_insert_test <- tibble(
    seq = c("XXOIXXIIX", "XXXIXOIIX", "XXXIXXIIO"),
    cigar = c("3M1I2M2I1M", "3M1I2M2I1M", "3M1I2M2I1M"),
    pos = 1,
    pos_idx = c(3, 5, 6),
    exp_obs = "O"
  )

  simple_del_test <- tibble(
    seq = c("OX", "XO"),
    cigar = c("2H1M2D1M", "2H1M2D1M"),
    pos = 1,
    pos_idx = c(1, 4),
    exp_obs = "O"
  )

  adv_del_test <- tibble(
    seq = c("OXX", "XOX", "XXO"),
    cigar = c("1M2D1M3D1M"),
    pos = 1,
    pos_idx = c(1, 4, 8),
    exp_obs = "O"
  )

  simple_indel_test <- tibble(
    seq = c("OXIIX", "XOIIX", "XXIIO"),
    cigar = c("1M2D1M2I1M"),
    pos = 1,
    pos_idx = c(1, 4, 5),
    exp_obs = "O"
  )

  adv_indel_test <- tibble(
    seq = c("OXIIXIXXX", "XOIIXIXXX", "XXIIOIXXX", "XXIIXIOXX", "XXIIXIXOX", "XXIIXIXXO"),
    cigar = c("1M2D1M2I1M1I2M3D1M"),
    pos = 1,
    pos_idx = c(1, 4, 5, 6, 7, 11),
    exp_obs = "O"
  )

  # Observed case
  obs_test <- tibble(
    seq = c("NGGGACATCCCTTAGGTGACCTCAGCNTGGGNGTGAGCATTAGGTAACTCTCTCCCTTGGCATCTGTGCCCTCTATTCACAGATAACTCTCTCCCCTGTTTCTAGGAGAAAAAAAAGTGCGTTCGCTACATACAAGG"),
    cigar = c("108M1D29M"),
    pos = 113165450,
    genomic_pos = 113165574,
    pos_idx = genomic_pos - pos + 1,
    exp_obs = "C"
  )


  aa <- bind_rows(
    no_correction_test,
    simple_insert_test,
    adv_insert_test,
    simple_del_test,
    adv_del_test,
    simple_indel_test,
    adv_indel_test,
    obs_test
  )


  tests <- aa %>%
    correct_pos_idx_w_cigar() %>%
    mutate(
      obs = str_sub(seq, pos_idx, pos_idx),
      test_accepted = obs == exp_obs
    )

  # All accepted
  expect_true(all(tests$test_accepted))
})

test_that("Getting mismatch positions", {

  # Example 1
  # POS+:1-2345-678
  # REF: T-TTTT-TTT
  # SEQ: TAAA-TAA--
  pos <- 1
  cigar <- "1M1I2M1D1M1I1M2D"
  MD <- "1T0T0^T1T0^TT"
  isize <- 10

  expect_equal(get_mismatch_genomic_pos_list(pos = pos, MDtag = MD), list(c(2, 3, 6)))
  expect_equal(get_match_genomic_pos_list(pos = pos, cigar = cigar, MDtag = MD), list(c(1, 5)))
})


test_that("Getting mismatch positions", {

  # Example 1
  # POS+:1-2345-678
  # REF: T-TTTT-TTT
  # SEQ: TAAA-TAA--
  pos <- 1
  cigar <- "1M1I2M1D1M1I1M2D"
  MD <- "1T0T0^T1T0^TT"
  isize <- 10

  expect_equal(get_mismatch_genomic_pos_list(pos = pos, MDtag = MD), list(c(2, 3, 6)))
  expect_equal(get_match_genomic_pos_list(pos = pos, cigar = cigar, MDtag = MD), list(c(1, 5)))
})


test_that("simple examples 2", {
  # Example 1
  # POS+:1-2345-678
  # REF: T-TTTT-TTT
  # SEQ: TAAA-TAA--
  ex1 <- data.frame(
    pos = 1,
    cigar = "1M1I2M1D1M1I1M2D",
    MD = "1T0T0^T1T0^TT",
    isize = 10
  )

  expect_true(sample_negative_read_positions(bam_df = ex1, n_samples = 1)$genomic_pos %in% c(1, 5))
  expect_equal(extract_mismatch_positions(bam_df = ex1)$genomic_pos, c(2, 3, 6))

  # Example 2
  # POS+:1-234-5
  # REF: T-TTT-T
  # SEQ: TCA-TAA

  ex2 <- data.frame(
    pos = 1,
    cigar = "1M1I1M1D1M1I1M",
    MD = "1T0^T1T",
    isize = 5
  )

  expect_true(sample_negative_read_positions(bam_df = ex2, n_samples = 1)$genomic_pos %in% c(1, 4))
  expect_equal(extract_mismatch_positions(bam_df = ex2)$genomic_pos, c(2, 5))

  # Example 3
  # POS+:123456789
  # REF: T-TTT-T
  # SEQ: TCA-TAA

  ex3 <- data.frame(
    pos = 1,
    cigar = "9M",
    MD = "9",
    isize = 9
  )

  expect_true(sample_negative_read_positions(bam_df = ex3, n_samples = 1)$genomic_pos %in% 1:9)
  expect_equal(extract_mismatch_positions(bam_df = ex3) %>% nrow(), 0)
})


test_that("simple examples 2", {
  # Example 1
  # POS+:123
  # REF: TTT
  # SEQ: T-T
  ex1 <- data.frame(
    pos_idx = 1:4,
    cigar = "1H1M1D1M"
  )

  res1 <- ex1 %>% correct_pos_idx_w_cigar()
  print (res1)
  expect_equal(res1$is_in_deletion, c(T, F, T, F))

  # Example 2
  # POS+:1234567
  # REF: TTTTTTTT
  # SEQ: -T---TTT
  ex2 <- data.frame(
    pos_idx = 1:8,
    cigar = "1H1M3D3M"
  )

  res2 <- ex2 %>% correct_pos_idx_w_cigar()
  expect_equal(res2$is_in_deletion, c(F, T, T, T, F, F, F))



  # Example 3
  # POS+:1--2345
  # REF: T--TTTT
  # SEQ: TAA-TTT
  ex3 <- data.frame(
    pos_idx = 1:6,
    cigar = "1D1M2I1D3M"
  )

  res3 <- ex3 %>% correct_pos_idx_w_cigar()
  expect_equal(res3$is_in_deletion, c(F, T, F, F, F))
})




test_that("Getting mismatch positions with hardclips", {

  # Example 1
  # POS+:12345678
  # REF: TTTTTTTT
  # SEQ: TTTTT
  pos <- 1
  cigar <- "3H5M"
  MD <- "5M"

  get_mismatch_genomic_pos_list(pos = pos, MDtag = MD)
  get_match_genomic_pos_list(pos = pos, cigar = cigar, MDtag = MD)

  expect_equal(get_mismatch_genomic_pos_list(pos = pos, MDtag = MD), list(c(2, 3, 6)))
  expect_equal(get_match_genomic_pos_list(pos = pos, cigar = cigar, MDtag = MD), list(c(1, 5)))
})


test_that("Getting mismatch positions with hardclips", {

  # Example 1
  # POS+:12345678
  # REF: TTTTTTTT
  # SEQ: TTATT
  pos <- 1
  cigar <- "1M3D5M"
  MD <- "3A2"

  mismatches = get_mismatch_genomic_pos_list(pos = pos, MDtag = MD)
  print (mismatches)

  matches = get_match_genomic_pos_list(pos = pos, cigar = cigar, MDtag = MD)
  print (matches)

  expect_equal(get_mismatch_genomic_pos_list(pos = pos, MDtag = MD), list(c(2, 3, 6)))
  expect_equal(get_match_genomic_pos_list(pos = pos, cigar = cigar, MDtag = MD), list(c(1, 5)))
})


test_that("Getting mismatch positions deletions", {

  # Example 1
  # POS+:1-2345-678
  # REF: TTTTTTTTTT
  # SEQ: --TTTATT--
  pos <- 1
  cigar <- "20D20H10M100D"


  get_mismatch_genomic_pos_list(pos = pos, MDtag = MD)
  get_match_genomic_pos_list(pos = pos, cigar = cigar, MDtag = MD)
})


