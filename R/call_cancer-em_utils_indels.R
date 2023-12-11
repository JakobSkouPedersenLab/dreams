prepare_em_input_indels <- function(mutations_df, read_positions_df, model, beta) {
  obs_is_mut_list <- list()
  error_ref_to_mut_list <- list()
  error_mut_to_ref_list <- list()

  if (nrow(mutations_df) >= 1) {
    for (i in 1:nrow(mutations_df)) {
      chr <- mutations_df[i, "chr"]
      genomic_pos <- mutations_df[i, "genomic_pos"]
      ref <- mutations_df[i, "ref"]
      alt <- mutations_df[i, "alt"]

      mut_reads_ref_alt <-
        read_positions_df %>%
        filter(
          # Filter reads positions to mutation
          .data$chr == !!chr, .data$genomic_pos == !!genomic_pos,
          # Filter mutation type (remove N and "other" than ref/alt alleles)
          .data$obs == !!ref | .data$obs == !!alt
        )

      # Add to list
      X <- mut_reads_ref_alt$obs == alt
      obs_is_mut_list <- append(obs_is_mut_list, list(X))

      # Error rates
      error_ref_df <- predict_error_rates_indels(mut_reads_ref_alt, model, beta)
      error_mut_df <- predict_error_rates_indels(mut_reads_ref_alt %>% mutate(ref = !!alt), model, beta)

      # Pick relevant error rates for ref and alt
      error_ref_to_ref_raw <- error_ref_df[[paste0(ref, "_corrected")]]
      error_ref_to_mut_raw <- error_ref_df[[paste0(alt, "_corrected")]]

      error_mut_to_ref_raw <- error_mut_df[[paste0(ref, "_corrected")]]
      error_mut_to_mut_raw <- error_mut_df[[paste0(alt, "_corrected")]]

      # Normalize
      error_ref_to_mut <- error_ref_to_mut_raw / (error_ref_to_ref_raw + error_ref_to_mut_raw)
      error_mut_to_ref <- error_mut_to_ref_raw / (error_mut_to_ref_raw + error_mut_to_mut_raw)

      # Add to list
      error_ref_to_mut_list <- append(error_ref_to_mut_list, list(error_ref_to_mut))
      error_mut_to_ref_list <- append(error_mut_to_ref_list, list(error_mut_to_ref))
    }
  }

  return(
    list(
      obs_is_mut_list = obs_is_mut_list,
      error_ref_to_mut_list = error_ref_to_mut_list,
      error_mut_to_ref_list = error_mut_to_ref_list
    )
  )
}

em_update <- function(par, obs_is_mut_list, error_ref_to_mut_list, error_mut_to_ref_list) {
  # Unpack par
  tf <- par[1]
  r <- par[2]
  n_frag_vec <- sapply(obs_is_mut_list, length)

  # E-step
  P_list <- update_P(
    tf = tf,
    r = r,
    error_mut_to_ref_list = error_mut_to_ref_list,
    error_ref_to_mut_list = error_ref_to_mut_list,
    obs_is_mut_list = obs_is_mut_list
  )

  P_Z1_vec <- P_list$P_Z1_vec
  P_Y_given_Z1_list <- P_list$P_Y_given_Z1_list

  # M-step
  tf_new <- update_tf(P_Y_given_Z1_list, P_Z1_vec, n_frag_vec)
  r_new <- update_r(P_Z1_vec)

  return(c(tf_new, r_new))
}

em_objective <- function(par, obs_is_mut_list, error_ref_to_mut_list, error_mut_to_ref_list) {
  # Unpack par
  tf <- par[1]
  r <- par[2]

  ll <- log_likelihood(
    obs_is_mut_list = obs_is_mut_list,
    error_mut_to_ref_list = error_mut_to_ref_list,
    error_ref_to_mut_list = error_ref_to_mut_list,
    tf = tf,
    r = r
  )

  # Return negative log-likelihood
  return(-ll)
}

get_starting_values <- function(obs_is_mut_list, error_mut_to_ref_list, error_ref_to_mut_list) {
  raw_observed_signal <- obs_is_mut_list %>% sapply(mean)
  observed_signal <- is.nan(raw_observed_signal) %>% ifelse(0, raw_observed_signal)

  # Find good starting values and (heuristic) limits for parameters
  # tf
  expected_count <- error_ref_to_mut_list %>%
    sapply(sum) %>%
    sum()
  count <- obs_is_mut_list %>%
    sapply(sum) %>%
    sum()
  coverage <- error_ref_to_mut_list %>%
    sapply(length) %>%
    sum()

  # Heuristic max value
  tf_max <- min(2 * max(observed_signal), 1.999)

  # TODO: Make limit function
  # Initial guess of tf
  tf_guess <- (2 * (count - expected_count) / coverage) %>% limit(1e-5, tf_max)

  # Set min to fraction of tf_guess
  tf_min <- tf_guess / 100

  # r
  r_min <- 1 / length(obs_is_mut_list)
  r_max <- mean(observed_signal > 0)

  r_guess <- mean(count > expected_count) %>% limit(r_min, r_max)


  # Grid size
  n_start_guess_low <- 4
  n_start_guess_high <- 4

  # Make grid
  tf_low_seq <- seq(tf_min, tf_guess, length.out = n_start_guess_low)
  tf_high_seq <- seq(tf_guess, tf_max, length.out = n_start_guess_high)

  tf_seq <- c(tf_low_seq, tf_high_seq) %>%
    sort() %>%
    unique()

  r_low_seq <- seq(r_min, r_guess, length.out = n_start_guess_low)
  r_high_seq <- seq(r_guess, r_max, length.out = n_start_guess_high)

  r_seq <- c(r_low_seq, r_high_seq) %>%
    sort() %>%
    unique()

  tf_r_grid <- expand.grid(tf = tf_seq, r = r_seq)

  # Get values of ll function
  ll_seq <- pmap_dbl(
    tf_r_grid,
    function(tf, r) {
      log_likelihood(
        obs_is_mut_list = obs_is_mut_list,
        error_mut_to_ref_list = error_mut_to_ref_list,
        error_ref_to_mut_list = error_ref_to_mut_list,
        r = r,
        tf = tf
      )
    }
  )

  ll_seq[is.nan(ll_seq)] <- min(ll_seq, na.rm = TRUE)

  # Choose best starting values
  tf_start <- tf_r_grid[which.max(ll_seq), 1]
  r_start <- tf_r_grid[which.max(ll_seq), 2]

  return(c(tf_start = tf_start, r_start = r_start))
}

run_full_em <- function(start_values, obs_is_mut_list, error_mut_to_ref_list, error_ref_to_mut_list) {
  tf_t <- start_values["tf_start"]
  r_t <- start_values["r_start"]

  tf_t_hist <- tf_t
  r_t_hist <- r_t

  ll_t_hist <- log_likelihood(
    obs_is_mut_list = obs_is_mut_list,
    error_mut_to_ref_list = error_mut_to_ref_list,
    error_ref_to_mut_list = error_ref_to_mut_list,
    tf = tf_t,
    r = r_t
  )

  P_Y_t_hist <- NULL

  ll_max <- -Inf

  max_iterations <- 10000
  EM_steps <- 0
  EM_converged <- F

  for (i in 1:max_iterations) {
    EM_steps <- EM_steps + 1

    # E-step
    P_t_list <- update_P(
      tf = tf_t,
      r = r_t,
      error_mut_to_ref_list = error_mut_to_ref_list,
      error_ref_to_mut_list = error_ref_to_mut_list,
      obs_is_mut_list = obs_is_mut_list
    )

    P_Y_t_vec <- P_t_list$P_Z1_vec
    P_Y_given_Z1_list <- P_t_list$P_Y_given_Z1_list

    # M-step
    n_frag_vec <- sapply(obs_is_mut_list, length)
    tf_t <- update_tf(P_Y_given_Z1_list, P_Y_t_vec, n_frag_vec)
    r_t <- update_r(P_Y_t_vec)

    ll_t <- log_likelihood(
      obs_is_mut_list = obs_is_mut_list,
      error_mut_to_ref_list = error_mut_to_ref_list,
      error_ref_to_mut_list = error_ref_to_mut_list,
      r = r_t,
      tf = tf_t
    )

    ll_t_hist <- c(ll_t_hist, ll_t)
    r_t_hist <- c(r_t_hist, r_t)
    tf_t_hist <- c(tf_t_hist, tf_t)
    P_Y_t_hist <- rbind(P_Y_t_hist, P_Y_t_vec)

    if (ll_t - ll_max < 1e-08 & EM_steps >= 5) {
      EM_converged <- TRUE
      break
    } else {
      ll_max <- ll_t
    }
  }

  # Gather results
  res <- list(
    tf_est = tf_t,
    r_est = r_t,
    EM_steps = EM_steps,
    fpeval = EM_steps,
    objfeval = EM_steps,
    EM_converged = EM_converged,
    hist = list(
      r_t_hist = r_t_hist,
      tf_t_hist = tf_t_hist,
      P_Y_t_hist = P_Y_t_hist,
      ll_t_hist = ll_t_hist
    )
  )

  return(res)
}

run_turbo_em <- function(start_values, obs_is_mut_list, error_mut_to_ref_list, error_ref_to_mut_list,
                         turboem_method = "squarem") {
  if (!requireNamespace("turboEM", quietly = TRUE)) {
    stop(
      "Package \"turboEM\" must be installed to use this function.",
      call. = FALSE
    )
  }
  turboem_res <- turboEM::turboem(
    par = start_values,
    fixptfn = em_update,
    objfn = em_objective,
    pconstr = function(par) {
      lower <- c(0, 0)
      upper <- c(2, 1)
      is_in_parameter_space <- all(lower <= par & par <= upper)
      return(is_in_parameter_space)
    },
    project = function(par) {
      par_project <- pmax(1e-8, pmin(par, c(2, 1)))
      return(par_project)
    },
    method = turboem_method,
    obs_is_mut_list = obs_is_mut_list,
    error_ref_to_mut_list = error_ref_to_mut_list,
    error_mut_to_ref_list = error_mut_to_ref_list,
    control.run =
      list(
        convtype = "objfn",
        tol = 1e-8
      )
  )

  res <- list(
    tf_est = turboem_res$par[1],
    r_est = turboem_res$par[2],
    EM_steps = turboem_res$itr,
    fpeval = turboem_res$fpeval,
    objfeval = turboem_res$objfeval,
    EM_converged = turboem_res$convergence
  )

  return(res)
}

get_em_parameter_estimates <- function(obs_is_mut_list, error_ref_to_mut_list, error_mut_to_ref_list, use_turboem) {
  observed_mut <- obs_is_mut_list %>%
    sapply(sum) %>%
    sum()
  total_observed_positions <- obs_is_mut_list %>%
    sapply(length) %>%
    sum()
  total_mut_ratio <- observed_mut / total_observed_positions

  # If no or only mut signal, return simple result:
  if (observed_mut == 0) {
    res <- list(
      tf_est = 0,
      r_est = 0,
      P_mut_is_present = rep(0, length(error_mut_to_ref_list)),
      EM_steps = 1,
      fpeval = 0,
      objfeval = 1,
      EM_converged = TRUE
    )

    return(res)
  } else if (total_mut_ratio == 1) {
    res <- list(
      tf_est = 2,
      r_est = 1,
      P_mut_is_present = rep(1, length(error_mut_to_ref_list)),
      EM_steps = 1,
      fpeval = 0,
      objfeval = 1,
      EM_converged = TRUE
    )

    return(res)
  }

  # Get starting values
  start_values <- get_starting_values(obs_is_mut_list, error_mut_to_ref_list, error_ref_to_mut_list)

  # Run EM algorithm
  if (use_turboem) {
    em_res <-
      run_turbo_em(
        start_values = start_values, obs_is_mut_list = obs_is_mut_list,
        error_mut_to_ref_list = error_mut_to_ref_list,
        error_ref_to_mut_list = error_ref_to_mut_list,
        turboem_method = "squarem"
      )
    # If "squarem" fails, rerun with regular EM
    squarem_failed <- is.na(em_res$tf_est)
    if (squarem_failed) {
      em_res <- run_turbo_em(
        start_values = start_values, obs_is_mut_list = obs_is_mut_list,
        error_mut_to_ref_list = error_mut_to_ref_list,
        error_ref_to_mut_list = error_ref_to_mut_list,
        turboem_method = "em"
      )
    }
  } else {
    em_res <- run_full_em(start_values, obs_is_mut_list, error_mut_to_ref_list, error_ref_to_mut_list)
  }

  # Check if 0 is better fit
  ll_0 <- log_likelihood(
    obs_is_mut_list = obs_is_mut_list,
    error_mut_to_ref_list = error_mut_to_ref_list,
    error_ref_to_mut_list = error_ref_to_mut_list,
    r = 0,
    tf = 0
  )

  ll_em <- log_likelihood(
    obs_is_mut_list = obs_is_mut_list,
    error_mut_to_ref_list = error_mut_to_ref_list,
    error_ref_to_mut_list = error_ref_to_mut_list,
    r = em_res$r_est,
    tf = em_res$tf_est
  )

  if (ll_0 >= ll_em) {
    em_res$tf_est <- 0
    em_res$r_est <- 0
  }

  # Get guess for hidden state
  P_Z1_vec <- update_P_Z(
    tf = em_res$tf_est,
    r = em_res$r_est,
    error_mut_to_ref_list = error_mut_to_ref_list,
    error_ref_to_mut_list = error_ref_to_mut_list,
    obs_is_mut_list = obs_is_mut_list
  )

  em_res$P_mut_is_present <- P_Z1_vec

  return(em_res)
}
