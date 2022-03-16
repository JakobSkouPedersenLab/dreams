# Small example

    Code
      res
    Output
         chr pos ref alt   tf_est    Q_val      ll_A      ll_0  exp_count count
      1 chr1  10   A   T 1.003677 7.078386 -1.386294 -4.925487 0.01462533     1
      2 chr2  13   A   T 0.000000 0.000000  0.000000  0.000000         NA     0
      3 chr3  17   A   T 0.000000 0.000000  0.000000  0.000000         NA     0
        coverage full_coverage obs_freq EM_converged EM_steps fpeval objfeval
      1        2             2      0.5         TRUE        2      3        3
      2        3             3      0.0        FALSE        0      0        0
      3        0             0      NaN        FALSE        0      0        0
              p_val mutation_detected
      1 0.007801926              TRUE
      2 1.000000000             FALSE
      3 1.000000000             FALSE

# Snapshot

    Code
      res
    Output
         chr pos ref alt    tf_est    Q_val       ll_A      ll_0 exp_count count
      1 chr1  10   A   T 0.1049353 2.464746  -3.803207  -5.03558 0.1243153     1
      2 chr2  13   A   T 0.1272296 8.334869 -10.880590 -15.04802 0.3144445     3
      3 chr3  17   A   T 0.0000000 0.000000   0.000000   0.00000        NA     0
        coverage full_coverage   obs_freq EM_converged EM_steps fpeval objfeval
      1       17            17 0.05882353         TRUE        3      5        5
      2       43            43 0.06976744         TRUE        3      5        5
      3        0             0        NaN        FALSE        0      0        0
              p_val mutation_detected
      1 0.116426529             FALSE
      2 0.003889128              TRUE
      3 1.000000000             FALSE

