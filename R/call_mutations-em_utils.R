get_starting_value_vc <- function(obs_is_mut, error_ref_to_mut, error_mut_to_ref) {
  expected_signal <- sum(error_ref_to_mut) / length(obs_is_mut)
  excess_signal <- mean(obs_is_mut) - expected_signal

  tf_guess <- max(1e-8, excess_signal * 2)

  # Improve guess if value close to 0 is better
  ll_guess <- log_likelihood_vc(tf_guess, obs_is_mut, error_ref_to_mut = error_ref_to_mut, error_mut_to_ref = error_mut_to_ref)
  ll_eps <- log_likelihood_vc(1e-8, obs_is_mut, error_ref_to_mut = error_ref_to_mut, error_mut_to_ref = error_mut_to_ref)
  ll_0 <- log_likelihood_vc(0, obs_is_mut, error_ref_to_mut = error_ref_to_mut, error_mut_to_ref = error_mut_to_ref)

  # Choose best starting value
  tf_start <- c(tf_guess, 1e-8, 0)[which.max(c(ll_guess, ll_eps, ll_0))]

  return(tf_start)
}

em_update_vc <- function(par, obs_is_mut, error_ref_to_mut, error_mut_to_ref) {
  fixed_r <- 1
  new_par <- em_update(c(par, fixed_r), list(obs_is_mut), list(error_ref_to_mut), list(error_mut_to_ref))
  return(new_par[1])
}

em_objective_vc <- function(par, obs_is_mut, error_ref_to_mut, error_mut_to_ref) {
  fixed_r <- 1
  em_objective(c(par, fixed_r), list(obs_is_mut), list(error_ref_to_mut), list(error_mut_to_ref))
}

log_likelihood_vc <- function(tf, obs_is_mut, error_ref_to_mut, error_mut_to_ref) {
  fixed_r <- 1
  return(log_likelihood(
    obs_is_mut_list = list(obs_is_mut),
    error_mut_to_ref_list = list(error_mut_to_ref),
    error_ref_to_mut_list = list(error_ref_to_mut),
    r = fixed_r,
    tf = tf
  ))
}

get_tf_estimate_vc <- function(obs_is_mut, error_ref_to_mut, error_mut_to_ref, use_warp_speed, max_it = 5000) {

  # If no or only mut signal, return simple result:
  observed_mut <- sum(obs_is_mut)
  total_observed_positions <- length(obs_is_mut)
  total_mut_ratio <- observed_mut / total_observed_positions
  if (observed_mut == 0) {
    res <- list(
      tf_est = 0,
      EM_steps = 1,
      fpeval = 0,
      objfeval = 1,
      EM_converged = TRUE
    )

    return(res)
  } else if (total_mut_ratio == 1) {
    res <- list(
      tf_est = 2,
      EM_steps = 1,
      fpeval = 0,
      objfeval = 1,
      EM_converged = TRUE
    )

    return(res)
  }

  # Get starting guess for tf
  tf_start <- get_starting_value_vc(obs_is_mut, error_ref_to_mut, error_mut_to_ref)

  if (use_warp_speed) {

    # Start speed up version of EM algorithm
    turboem_res <- turboEM::turboem(
      par = tf_start,
      fixptfn = em_update_vc,
      objfn = em_objective_vc,
      pconstr = function(par) {
        lower <- 0
        upper <- 2
        is_in_parameter_space <- all(lower < par & par < upper)
        return(is_in_parameter_space)
      },
      project = function(par) {
        par_project <- pmax(1e-8, pmin(par, 2))
        return(par_project)
      },
      method = "squarem",
      obs_is_mut = obs_is_mut,
      error_ref_to_mut = error_ref_to_mut,
      error_mut_to_ref = error_mut_to_ref,
      control.run =
        list(
          convtype = "objfn",
          tol = 1e-8
        )
    )

    tf_est <- turboem_res$par[1]
    EM_steps <- turboem_res$itr
    fpeval <- turboem_res$fpeval
    objfeval <- turboem_res$objfeval
    EM_converged <- turboem_res$convergence
  } else {
    # Start regular EM algorithm
    EM_converged <- FALSE
    EM_steps <- 0
    tf <- tf_start

    for (i in 1:max_it) {
      EM_steps <- EM_steps + 1
      tf_old <- tf

      # Update tf
      tf <- em_update_vc(par = tf_old, obs_is_mut = obs_is_mut, error_ref_to_mut = error_ref_to_mut, error_mut_to_ref = error_mut_to_ref)

      if (abs(tf - tf_old) / (tf_old + 1e-8) < 1e-8 & EM_steps >= 5) {
        EM_converged <- TRUE
        break
      }
    }

    fpeval <- EM_steps
    objfeval <- EM_steps
    tf_est <- tf
  }


  # Check if 0 is better fit
  ll_0 <- log_likelihood_vc(
    obs_is_mut = obs_is_mut,
    error_mut_to_ref = error_mut_to_ref,
    error_ref_to_mut = error_ref_to_mut,
    tf = 0
  )

  ll_em <- log_likelihood_vc(
    obs_is_mut = obs_is_mut,
    error_mut_to_ref = error_mut_to_ref,
    error_ref_to_mut = error_ref_to_mut,
    tf = tf_est
  )

  if (ll_0 >= ll_em) {
    tf_est <- 0
  }

  # Return results
  res <- list(
    tf_est = tf_est,
    EM_converged = EM_converged,
    EM_steps = EM_steps,
    fpeval = fpeval,
    objfeval = objfeval
  )

  return(res)
}
