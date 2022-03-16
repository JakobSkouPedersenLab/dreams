
log_lik <- function(tf, sample_mutations, err_ref_to_mut, err_mut_to_ref) {
  log_lik_mut <- log((1 - err_mut_to_ref) * tf / 2 + err_ref_to_mut * (1 - tf / 2))
  log_lik_ref <- log((1 - err_ref_to_mut) * (1 - tf / 2) + err_mut_to_ref * tf / 2)

  return(sum(log_lik_mut[sample_mutations == 1]) + sum(log_lik_ref[sample_mutations == 0]))
}

update_P_Y1 <- function(x, tf, err_ref_to_mut, err_mut_to_ref) {
  # Get scaled probabilities for fragment state(Y)
  Q_Y0 <- ifelse(x == 0,
    (1 - tf / 2) * (1 - err_ref_to_mut),
    (1 - tf / 2) * err_ref_to_mut
  )

  Q_Y1 <- ifelse(x == 0,
    tf / 2 * err_mut_to_ref,
    tf / 2 * (1 - err_mut_to_ref)
  )

  # Normalize probabilities
  P_Y1 <- Q_Y1 / (Q_Y0 + Q_Y1)

  return(P_Y1)
}

get_tumor_freq <- function(x, err_ref_to_mut, err_mut_to_ref, use_warp_speed, max_it = 5000) {

  # Get starting guess for tf
  expected_signal <- sum(err_ref_to_mut) / length(x)
  excess_signal <- mean(x) - expected_signal

  tf <- max(1e-8, excess_signal * 2)

  # Improve guess if value close to 0 is better
  ll_0 <- log_lik(0, x, err_ref_to_mut = err_ref_to_mut, err_mut_to_ref = err_mut_to_ref)
  ll_eps <- log_lik(1e-8, x, err_ref_to_mut = err_ref_to_mut, err_mut_to_ref = err_mut_to_ref)
  ll_start <- log_lik(tf, x, err_ref_to_mut = err_ref_to_mut, err_mut_to_ref = err_mut_to_ref)

  if (ll_eps >= ll_start) {
    tf <- 1e-8

    if (ll_0 >= ll_eps) {
      tf <- 0
    }
  }

  if (use_warp_speed) {
    em_update <- function(par, sample_mutations, err_ref_to_mut, err_mut_to_ref) {
      tf <- par[1]

      # Update P_Y1
      P_Y1 <- update_P_Y1(x = sample_mutations, tf = tf, err_ref_to_mut = err_ref_to_mut, err_mut_to_ref = err_mut_to_ref)

      # Update tf
      tf_new <- mean(P_Y1) * 2

      return(c(tf_new, P_Y1))
    }

    em_objective <- function(par, sample_mutations, err_ref_to_mut, err_mut_to_ref) {
      tf <- par[1]

      ll <- log_lik(tf, sample_mutations, err_ref_to_mut = err_ref_to_mut, err_mut_to_ref = err_mut_to_ref)

      return(-ll)
    }

    # Start speed up version of EM algorithm
    start_values <- c(tf, rep(0, length(err_ref_to_mut)))
    turboem_res <- turboEM::turboem(
      par = start_values,
      fixptfn = em_update,
      objfn = em_objective,
      pconstr = function(par) {
        par <- par[1]
        lower <- c(0, 0)
        upper <- c(2, 1)
        is_in_parameter_space <- all(lower < par & par < upper)
        return(is_in_parameter_space)
      },
      project = function(par) {
        par <- par[1]
        par_project <- pmax(1e-8, pmin(par, c(2 - 1e-8, 1 - 1e-8)))
        return(par_project)
      },
      method = "squarem",
      sample_mutations = x,
      err_ref_to_mut = err_ref_to_mut,
      err_mut_to_ref = err_mut_to_ref,
      control.run =
        list(
          convtype = "objfn",
          tol = 1e-8
        )
    )

    EM_steps <- turboem_res$itr
    fpeval <- turboem_res$fpeval
    objfeval <- turboem_res$objfeval
    EM_converged <- turboem_res$convergence

    tf <- turboem_res$par[1]
  } else {
    # Start regular EM algorithm
    EM_converged <- FALSE
    EM_steps <- 0

    for (i in 1:max_it) {
      EM_steps <- EM_steps + 1
      tf_old <- tf

      # Update z
      z <- update_P_Y1(x = x, tf = tf_old, err_ref_to_mut = err_ref_to_mut, err_mut_to_ref = err_mut_to_ref)

      # Update tf
      tf <- mean(z) * 2

      if (abs(tf - tf_old) / (tf_old + 1e-8) < 1e-8 & EM_steps >= 5) {
        EM_converged <- TRUE
        break
      }
    }

    fpeval <- EM_steps
    objfeval <- EM_steps
  }


  # Final check of tf value
  ll_final <- log_lik(tf, x, err_ref_to_mut = err_ref_to_mut, err_mut_to_ref = err_mut_to_ref)
  if (ll_0 >= ll_final) {
    tf <- 0
  }

  res <- list(
    tf = tf,
    EM_converged = EM_converged,
    EM_steps = EM_steps,
    fpeval = fpeval,
    objfeval = objfeval
  )

  return(res)
}

em_test_mutation <- function(sample_mutations, err_ref_to_mut, err_mut_to_ref, use_warp_speed) {
  get_tumor_freq_res <-
    get_tumor_freq(
      sample_mutations,
      err_ref_to_mut = err_ref_to_mut,
      err_mut_to_ref = err_mut_to_ref,
      use_warp_speed = use_warp_speed
    )

  tf_est <- get_tumor_freq_res$tf
  EM_converged <- get_tumor_freq_res$EM_converged
  EM_steps <- get_tumor_freq_res$EM_steps
  fpeval <- get_tumor_freq_res$fpeval
  objfeval <- get_tumor_freq_res$fpeval

  ll_0 <- log_lik(0, sample_mutations, err_ref_to_mut = err_ref_to_mut, err_mut_to_ref = err_mut_to_ref)
  ll_A <- log_lik(tf_est, sample_mutations, err_ref_to_mut = err_ref_to_mut, err_mut_to_ref = err_mut_to_ref)
  Q_val <- -2 * (ll_0 - ll_A)
  p_val <- pchisq(Q_val, 1, lower.tail = FALSE)

  return(
    list(
      tf_est = tf_est,
      Q_val = Q_val,
      p_val = p_val,
      ll_A = ll_A,
      ll_0 = ll_0,
      EM_converged = EM_converged,
      EM_steps = EM_steps,
      fpeval = fpeval,
      objfeval = objfeval
    )
  )
}

call_mutations <- function(mutations_df, all_reads, model, beta, alpha = 0.05, use_warp_speed = TRUE) {

  # Clean up mutations
  mutations_df <- mutations_df %>%
    select(
      "chr" = matches("chr|CHR|CHROM"),
      "genomic_pos" = matches("pos|POS"),
      "ref" = matches("ref|REF"),
      "alt" = matches("alt|ALT|obs|OBS")
    )

  res_df <- NULL

  if (nrow(mutations_df) == 0) {
    return(data.frame())
  }

  for (i in 1:nrow(mutations_df)) {
    chr <- mutations_df$chr[i]
    genomic_pos <- as.numeric(as.character(mutations_df$genomic_pos[i]))
    ref <- mutations_df$ref[i]
    alt <- mutations_df$alt[i]

    reads <- all_reads %>% filter(chr == !!chr, genomic_pos == !!genomic_pos)
    full_coverage <- nrow(reads)

    # Prepare EM input
    em_input <- prepare_em_input(mutations_df = mutations_df[i, ], reads_df = reads, model = model, beta = beta)

    X_list <- em_input$X_list
    error_ref_to_mut_list <- em_input$error_ref_to_mut_list
    error_mut_to_ref_list <- em_input$error_mut_to_ref_list


    # if no reads cover position, skip to next mutation
    # TODO
    if (nrow(reads) == 0) {
      new_row <- data.frame(
        chr = chr, pos = genomic_pos, ref = ref, alt = alt,
        tf_est = 0, p_val = 1,
        EM_converged = FALSE, EM_steps = 0, fpeval = 0, objfeval = 0,
        count = n_mut, exp_count = NA,
        coverage = nrow(reads_ref_alt), full_coverage = full_coverage,
        obs_freq = n_mut / nrow(reads_ref_alt),
        Q_val = 0, ll_A = 0, ll_0 = 0,
        mutation_detected = FALSE
      )
      res_df <- rbind(res_df, new_row)
      next
    }

    # No reads reads supporting reference or mutation
    # TODO
    if (n_ref == 0 | n_mut == 0) {
      new_row <- data.frame(
        chr = chr, pos = genomic_pos, ref = ref, alt = alt,
        tf_est = 0, p_val = 1,
        EM_converged = FALSE, EM_steps = 0, fpeval = 0, objfeval = 0,
        count = n_mut, exp_count = NA,
        coverage = nrow(reads_ref_alt), full_coverage = full_coverage,
        obs_freq = n_mut / nrow(reads_ref_alt),
        Q_val = 0, ll_A = 0, ll_0 = 0,
        mutation_detected = FALSE
      )
      res_df <- rbind(res_df, new_row)
      next
    }

    sample_mutations <- X_list[[1]]
    error_ref_to_mut_norm <- error_ref_to_mut_list[[1]]
    error_mut_to_ref_norm <- error_mut_to_ref_list[[1]]


    ##### EM Log-lik

    em_res <- em_test_mutation(sample_mutations, error_ref_to_mut_norm, error_mut_to_ref_norm, use_warp_speed)

    new_row <- data.frame(
      chr = chr, pos = genomic_pos, ref = ref, alt = alt,
      tf_est = em_res$tf_est,
      # tf_min = tf_CI$tf_min, tf_max = tf_CI$tf_max,
      Q_val = em_res$Q_val,
      ll_A = em_res$ll_A,
      ll_0 = em_res$ll_0,
      exp_count = sum(error_ref_to_mut_norm),
      count = n_mut, # TODO: sapply(X_list, sum),
      coverage = nrow(reads_ref_alt), # TODO: sapply(X_list, length),
      full_coverage = full_coverage,
      obs_freq = sapply(X_list, mean),
      EM_converged = em_res$EM_converged,
      EM_steps = em_res$EM_steps,
      fpeval = em_res$fpeval,
      objfeval = em_res$objfeval,
      p_val = em_res$p_val,
      mutation_detected = em_res$p_val <= alpha
    )

    res_df <- rbind(res_df, new_row)
  }
  return(res_df)
}








# -------------------------------------------------------------------------

call_mutations_new <- function(mutations_df, reads_df, model, beta, alpha = 0.05, calculate_confidence_intervals = FALSE, use_warp_speed = TRUE) {
  # Clean up mutations
  mutations_df <- mutations_df %>%
    select(
      "chr" = matches("chr|CHR|CHROM"),
      "genomic_pos" = matches("pos|POS"),
      "ref" = matches("ref|REF"),
      "alt" = matches("alt|ALT|obs|OBS")
    )

  # Stop if mutations do not have the expected columns
  mutations_expected_columns <- c("chr", "genomic_pos", "ref", "alt")
  if (!all(mutations_expected_columns %in% colnames(mutations_df))) {
    stop("mutations_df should have the columns ['chr', genomic_pos', 'ref, 'alt']")
  }

  # Stop if reads do not have the expected columns
  reads_expected_columns <- c("chr", "genomic_pos", "ref", "obs")
  if (!all(reads_expected_columns %in% colnames(reads_df))) {
    stop("reads_df should have the columns ['chr', genomic_pos', 'ref, 'obs']")
  }

  # Prepare inputs for algorithm
  em_input <- prepare_em_input(mutations_df = mutations_df, reads_df = reads_df, model = model, beta = beta)
  X_list <- em_input$X_list
  error_ref_to_mut_list <- em_input$error_ref_to_mut_list
  error_mut_to_ref_list <- em_input$error_mut_to_ref_list

  # If no mutations return empty result
  if (nrow(mutations_df) == 0) {
    return(
      list(
        cancer_info = empty_cancer_info(mutations_df, em_input),
        mutation_info = data.frame()
      )
    )
  }

  # EM algorithm
  em_res <- get_em_parameter_estimates(
    X_list = X_list,
    error_ref_to_mut_list = error_ref_to_mut_list,
    error_mut_to_ref_list = error_mut_to_ref_list,
    use_warp_speed = use_warp_speed
  )

  # Confidence intervals
  if (calculate_confidence_intervals) {
    tf_CI <- get_tf_CI(
      X_list = X_list,
      error_mut_to_ref_list = error_mut_to_ref_list,
      error_ref_to_mut_list = error_ref_to_mut_list,
      r_est = em_res$r_est,
      tf_est = em_res$tf_est,
      alpha = alpha
    )
    r_CI <- get_r_CI(
      X_list = X_list,
      error_mut_to_ref_list = error_mut_to_ref_list,
      error_ref_to_mut_list = error_ref_to_mut_list,
      r_est = em_res$r_est,
      tf_est = em_res$tf_est,
      alpha = alpha
    )
  } else {
    tf_CI <- list(tf_min = NA, tf_max = NA)
    r_CI <- list(r_min = NA, r_max = NA)
  }

  # Test significance
  ll_0 <- log_likelihood(
    X_list = X_list,
    error_mut_to_ref_list = error_mut_to_ref_list,
    error_ref_to_mut_list = error_ref_to_mut_list,
    r = 0,
    tf = 0
  )

  ll_A <- log_likelihood(
    X_list = X_list,
    error_mut_to_ref_list = error_mut_to_ref_list,
    error_ref_to_mut_list = error_ref_to_mut_list,
    r = em_res$r_est,
    tf = em_res$tf_est
  )

  Q_val <- -2 * (ll_0 - ll_A)
  p_val <- stats::pchisq(Q_val, 2, lower.tail = FALSE)

  # Collect cancer information
  mutation_info <-
    data.frame(
      tf_est = em_res$tf_est,
      tf_min = tf_CI$tf_min, tf_max = tf_CI$tf_max,
      Q_val = Q_val,
      ll_A = ll_A,
      ll_0 = ll_0,
      est_mutations_present = sum(em_res$P_mut_is_present),
      exp_count = sapply(error_ref_to_mut_list, sum),
      count = sapply(X_list, sum),
      coverage = sapply(X_list, length),
      obs_freq = sapply(X_list, mean),
      EM_converged = em_res$EM_converged,
      EM_steps = em_res$EM_steps,
      fpeval = em_res$fpeval,
      objfeval = em_res$objfeval,
      p_val = p_val,
      mutation_detected = p_val <= alpha
    )

  return(mutation_info)
}
