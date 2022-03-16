# Snapshot

    Code
      res
    Output
         chr genomic_pos ref alt full_coverage    tf_est EM_converged EM_steps fpeval
      1 chr1          10   A   T            17 0.1049353         TRUE        3      5
      2 chr2          13   A   T            43 0.1272296         TRUE        3      5
      3 chr3          17   A   T            NA 0.0000000        FALSE        0      0
        objfeval exp_count count coverage   obs_freq      ll_0       ll_A    Q_val
      1        4 0.1243153     1       17 0.05882353  -5.03558  -3.803207 2.464746
      2        4 0.3144445     3       43 0.06976744 -15.04802 -10.880590 8.334869
      3        0        NA     0        0        NaN   0.00000   0.000000 0.000000
              p_val mutation_detected
      1 0.116426529             FALSE
      2 0.003889128              TRUE
      3 1.000000000             FALSE

---

    Code
      slow_res
    Output
         chr genomic_pos ref alt full_coverage    tf_est EM_converged EM_steps fpeval
      1 chr1          10   A   T            17 0.1049353         TRUE        9      9
      2 chr2          13   A   T            43 0.1272296         TRUE        8      8
      3 chr3          17   A   T            NA 0.0000000        FALSE        0      0
        objfeval exp_count count coverage   obs_freq      ll_0       ll_A    Q_val
      1        9 0.1243153     1       17 0.05882353  -5.03558  -3.803207 2.464746
      2        8 0.3144445     3       43 0.06976744 -15.04802 -10.880590 8.334869
      3        0        NA     0        0        NaN   0.00000   0.000000 0.000000
              p_val mutation_detected
      1 0.116426529             FALSE
      2 0.003889128              TRUE
      3 1.000000000             FALSE

