em_update <- function(par, X_list, error_ref_to_mut_list, error_mut_to_ref_list) {
  # Unpack par
  tf <- par[1]
  r <- par[2]
  n_frag_vec <- sapply(X_list, length)

  # E-step
  P_list <- update_P(
    tf = tf,
    r = r,
    error_mut_to_ref_list = error_mut_to_ref_list,
    error_ref_to_mut_list = error_ref_to_mut_list,
    X_list = X_list
  )

  P_Z1_vec <- P_list$P_Z1_vec
  P_Y_given_Z1_list <- P_list$P_Y_given_Z1_list

  # M-step
  tf_new <- update_tf(P_Y_given_Z1_list, P_Z1_vec, n_frag_vec)
  r_new <- update_r(P_Z1_vec)

  return(c(tf_new, r_new))
}

em_objective <- function(par, X_list, error_ref_to_mut_list, error_mut_to_ref_list) {
  # Unpack par
  tf <- par[1]
  r <- par[2]

  ll <- log_likelihood(
    X_list = X_list,
    error_mut_to_ref_list = error_mut_to_ref_list,
    error_ref_to_mut_list = error_ref_to_mut_list,
    tf = tf,
    r = r
  )

  # Return negative log-likelihood
  return(-ll)
}

get_starting_values <- function(X_list, error_mut_to_ref_list, error_ref_to_mut_list) {
  raw_observed_signal <- X_list %>% sapply(mean)
  observed_signal <- is.nan(raw_observed_signal) %>% ifelse(0, raw_observed_signal)

  # Find good starting values
  # (Heuristic) Limits for parameters
  tf_max <- min(max(observed_signal) * 2, 2)

  r_max <- max(1 / length(X_list), mean(observed_signal > 0))

  # Grid size
  tf_n_start_guess <- 7
  r_n_start_guess <- 7

  # Make grid
  tf_seq <- tf_max * 10^seq(-3, 0, length.out = tf_n_start_guess)
  r_seq <- seq(1e-8, r_max, length.out = r_n_start_guess)
  tf_r_grid <- expand.grid(tf = tf_seq, r = r_seq)

  # Get values of ll function
  ll_seq <- pmap_dbl(
    tf_r_grid,
    function(tf, r) {
      log_likelihood(
        X_list = X_list,
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

run_full_em <- function(X_list, error_mut_to_ref_list, error_ref_to_mut_list) {
  start_values <- get_starting_values(X_list, error_mut_to_ref_list, error_ref_to_mut_list)

  tf_t <- start_values["tf_start"]
  r_t <- start_values["r_start"]

  tf_t_hist <- tf_t
  r_t_hist <- r_t

  ll_t_hist <- log_likelihood(
    X_list = X_list,
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
      X_list = X_list
    )

    P_Y_t_vec <- P_t_list$P_Z1_vec
    P_Y_given_Z1_list <- P_t_list$P_Y_given_Z1_list

    # M-step
    n_frag_vec <- sapply(X_list, length)
    tf_t <- update_tf(P_Y_given_Z1_list, P_Y_t_vec, n_frag_vec)
    r_t <- update_r(P_Y_t_vec)

    ll_t <- log_likelihood(
      X_list = X_list,
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
    P_mut_is_present = P_Y_t_vec,
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

run_turbo_em <- function(X_list, error_mut_to_ref_list, error_ref_to_mut_list) {
  # Run EM algorithm
  start_values <- get_starting_values(X_list, error_mut_to_ref_list, error_ref_to_mut_list)

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
    method = "squarem",
    X_list = X_list,
    error_ref_to_mut_list = error_ref_to_mut_list,
    error_mut_to_ref_list = error_mut_to_ref_list,
    control.run =
      list(
        convtype = "objfn",
        tol = 1e-8
      )
  )

  tf_est <- turboem_res$par[1]
  r_est <- turboem_res$par[2]

  # Check if 0 is better fit
  ll_0 <- log_likelihood(
    X_list = X_list,
    error_mut_to_ref_list = error_mut_to_ref_list,
    error_ref_to_mut_list = error_ref_to_mut_list,
    r = 0,
    tf = 0
  )

  ll_em <- log_likelihood(
    X_list = X_list,
    error_mut_to_ref_list = error_mut_to_ref_list,
    error_ref_to_mut_list = error_ref_to_mut_list,
    r = r_est,
    tf = tf_est
  )

  if (ll_0 >= ll_em) {
    tf_est <- 0
    r_est <- 0
  }

  # Get guess for hidden state
  P_Y_t_vec <- update_P_Z(
    tf = tf_est,
    r = r_est,
    error_mut_to_ref_list = error_mut_to_ref_list,
    error_ref_to_mut_list = error_ref_to_mut_list,
    X_list = X_list
  )

  res <- list(
    tf_est = tf_est,
    r_est = r_est,
    P_mut_is_present = P_Y_t_vec,
    EM_steps = turboem_res$itr,
    fpeval = turboem_res$fpeval,
    objfeval = turboem_res$objfeval,
    EM_converged = turboem_res$convergence
  )

  return(res)
}

get_em_parameter_estimates <- function(X_list, error_ref_to_mut_list, error_mut_to_ref_list, use_warp_speed) {
  observed_mut <- X_list %>%
    sapply(sum) %>%
    sum()
  total_observed_positions <- X_list %>%
    sapply(length) %>%
    sum()
  mut_ratio <- observed_mut / total_observed_positions

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
  } else if (mut_ratio == 1) {
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

  if (use_warp_speed) {
    return(run_turbo_em(X_list, error_mut_to_ref_list, error_ref_to_mut_list))
  } else {
    return(run_full_em(X_list, error_mut_to_ref_list, error_ref_to_mut_list))
  }
}
