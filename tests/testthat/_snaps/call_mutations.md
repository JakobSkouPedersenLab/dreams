# Snapshot of bigger example

    Code
      res
    Output
         chr genomic_pos ref alt full_coverage    tf_est EM_converged EM_steps fpeval
      1 chr1          10   A   T            17 0.1049353         TRUE        3      5
      2 chr2          13   A   T            43 0.1272296         TRUE        3      5
      3 chr3          17   A   T            NA 0.0000000         TRUE        1      0
        objfeval     tf_min    tf_max exp_count count coverage   obs_freq      ll_0
      1        4 0.02778288 0.2885132 0.1243153     1       17 0.05882353  -5.03558
      2        4 0.02819616 0.2870679 0.3144445     3       43 0.06976744 -15.04802
      3        1 0.00000000 0.5063943 0.0000000     0        0        NaN   0.00000
              ll_A    Q_val df       p_val mutation_detected
      1  -3.803207 2.464746  1 0.116426529             FALSE
      2 -10.880590 8.334869  1 0.003889128              TRUE
      3   0.000000 0.000000  1 1.000000000             FALSE

---

    Code
      slow_res
    Output
         chr genomic_pos ref alt full_coverage    tf_est EM_converged EM_steps fpeval
      1 chr1          10   A   T            17 0.1049353         TRUE        9      9
      2 chr2          13   A   T            43 0.1272296         TRUE        8      8
      3 chr3          17   A   T            NA 0.0000000         TRUE        1      0
        objfeval     tf_min    tf_max exp_count count coverage   obs_freq      ll_0
      1        9 0.02778288 0.2885132 0.1243153     1       17 0.05882353  -5.03558
      2        8 0.02819616 0.2870679 0.3144445     3       43 0.06976744 -15.04802
      3        1 0.00000000 0.5063943 0.0000000     0        0        NaN   0.00000
              ll_A    Q_val df       p_val mutation_detected
      1  -3.803207 2.464746  1 0.116426529             FALSE
      2 -10.880590 8.334869  1 0.003889128              TRUE
      3   0.000000 0.000000  1 1.000000000             FALSE

