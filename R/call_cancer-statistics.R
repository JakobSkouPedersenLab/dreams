update_P_Z <- function(tf, r, error_mut_to_ref_list, error_ref_to_mut_list, obs_is_mut_list) {
  log_Q_Z0 <- pmap_dbl(
    tibble(
      x = obs_is_mut_list,
      err_RM = error_ref_to_mut_list,
    ),
    function(x, err_RM) {
      log(1 - r) +
        sum(log(1 - err_RM[x == 0])) +
        sum(log(err_RM[x == 1]))
    }
  )

  log_Q_Z1 <- pmap_dbl(
    tibble(
      x = obs_is_mut_list,
      err_RM = error_ref_to_mut_list,
      err_MR = error_mut_to_ref_list
    ),
    function(x, err_RM, err_MR) {
      log(r) +
        sum(log(tf / 2 * err_MR[x == 0] +
          (1 - tf / 2) * (1 - err_RM[x == 0]))) +
        sum(log(tf / 2 * (1 - err_MR[x == 1]) +
          (1 - tf / 2) * err_RM[x == 1]))
    }
  )

  # Normalize to probabilities
  max_log_P <- pmax(log_Q_Z0, log_Q_Z1)
  Q_Z0 <- exp(log_Q_Z0 - max_log_P)
  Q_Z1 <- exp(log_Q_Z1 - max_log_P)

  P_Z1 <- Q_Z1 / (Q_Z0 + Q_Z1)

  return(P_Z1)
}

update_P <- function(tf, r, error_mut_to_ref_list, error_ref_to_mut_list, obs_is_mut_list) {

  # Get probabilities for site mutation state (Z)
  P_Z1_vec <- update_P_Z(
    tf = tf,
    r = r,
    error_mut_to_ref_list = error_mut_to_ref_list,
    error_ref_to_mut_list = error_ref_to_mut_list,
    obs_is_mut_list = obs_is_mut_list
  )

  # Get probabilities for fragment state given site state (Y given Z)
  # Get scaled probabilities (Q1 + Q2 = c)
  Q_Y0_given_Z1_list <- pmap(
    tibble(
      x = obs_is_mut_list,
      err_RM = error_ref_to_mut_list,
    ),
    function(x, err_RM) {
      ifelse(x == 0,
        (1 - tf / 2) * (1 - err_RM),
        (1 - tf / 2) * err_RM
      )
    }
  )

  Q_Y1_given_Z1_list <- pmap(
    tibble(
      x = obs_is_mut_list,
      err_MR = error_mut_to_ref_list,
    ),
    function(x, err_MR) {
      ifelse(x == 0,
        tf / 2 * err_MR,
        tf / 2 * (1 - err_MR)
      )
    }
  )

  # Get probabilities
  P_Y_given_Z1_list <- pmap(
    tibble(
      Q_Y0_given_Z1 = Q_Y0_given_Z1_list,
      Q_Y1_given_Z1 = Q_Y1_given_Z1_list
    ),
    function(Q_Y0_given_Z1, Q_Y1_given_Z1) {
      # Normalize probabilities
      Q_sum <- Q_Y0_given_Z1 + Q_Y1_given_Z1

      P_Y0_given_Z1 <- Q_Y0_given_Z1 / Q_sum
      P_Y1_given_Z1 <- Q_Y1_given_Z1 / Q_sum

      # Gather probabilities
      return(cbind(P_Y0_given_Z1, P_Y1_given_Z1))
    }
  )

  # Return probabilities
  return(list(
    P_Z1_vec = P_Z1_vec,
    P_Y_given_Z1_list = P_Y_given_Z1_list
  ))
}

update_r <- function(P_Z1_vec) {

  # r is the mean expected prob. of mutation
  r_new <- max(mean(P_Z1_vec), 1 / length(P_Z1_vec))

  return(r_new)
}

update_tf <- function(P_Y_given_Z1_list, P_Z1_vec, n_frag_vec) {

  # Expected read-positions with mutation
  M <- pmap_dbl(
    tibble(
      P_Y_given_Z1 = P_Y_given_Z1_list,
      P_Z1 = P_Z1_vec,
    ),
    function(P_Y_given_Z1, P_Z1) {
      P_Z1 * sum(P_Y_given_Z1[, 2])
    }
  ) %>% sum()

  # Expected number of read-positions covering positions with mutation present
  N <- sum(P_Z1_vec * n_frag_vec)

  if (N == 0) {
    tf_new <- 0
  } else {
    tf_new <- 2 * M / N # 2* is for heterozygozity i.e. half the tumor is "normal"
  }

  return(tf_new)
}

log_likelihood <- function(obs_is_mut_list, error_mut_to_ref_list, error_ref_to_mut_list, tf, r = 1) {
  # Convert to lists if necessary
  if (!is.list(obs_is_mut_list)) {
    obs_is_mut_list <- list(obs_is_mut_list)
  }
  if (!is.list(error_mut_to_ref_list)) {
    error_mut_to_ref_list <- list(error_mut_to_ref_list)
  }
  if (!is.list(error_ref_to_mut_list)) {
    error_ref_to_mut_list <- list(error_ref_to_mut_list)
  }

  # Check input dimensions
  if (!all(sapply(obs_is_mut_list, length) == sapply(error_mut_to_ref_list, length))) {
    stop("Dimensions of list inputs are not equal")
  }

  # Check model parameters
  if (is.nan(tf) || is.na(tf) || tf < 0 || 2 < tf) {
    stop(paste("Illegal value of tf:", tf))
  }

  ll_m_vec <-
    pmap_dbl(
      tibble(
        X = obs_is_mut_list,
        err_MR = error_mut_to_ref_list,
        err_RM = error_ref_to_mut_list
      ),
      function(X, err_MR, err_RM) {
        err_RM_XM <- err_RM[X == 1]
        err_RM_XR <- err_RM[X == 0]
        err_MR_XM <- err_MR[X == 1]
        err_MR_XR <- err_MR[X == 0]

        # Calculate the vector products in log space
        log_p1_prod <- sum(log((1 - tf / 2) * err_RM_XM + tf / 2 * (1 - err_MR_XM)))
        log_p2_prod <- sum(log((1 - tf / 2) * (1 - err_RM_XR) + tf / 2 * err_MR_XR))
        log_p3_prod <- sum(log(err_RM_XM))
        log_p4_prod <- sum(log(1 - err_RM_XR))

        # Note the max value for exp log trick
        max_val <- max(log_p1_prod + log_p2_prod, log_p3_prod + log_p4_prod)

        ll_m <-
          log(
            r * exp(log_p1_prod + log_p2_prod - max_val) +
              (1 - r) * exp(log_p3_prod + log_p4_prod - max_val)
          ) + max_val

        return(sum(ll_m))
      }
    )

  ll <- sum(ll_m_vec)

  return(ll)
}

#' @importFrom stats qchisq uniroot
get_tf_CI <- function(obs_is_mut_list, error_mut_to_ref_list, error_ref_to_mut_list, r_est, tf_est, alpha) {
  ll_max <- log_likelihood(
    obs_is_mut_list = obs_is_mut_list,
    error_mut_to_ref_list = error_mut_to_ref_list,
    error_ref_to_mut_list = error_ref_to_mut_list,
    r = r_est,
    tf = tf_est
  )

  # Find minimum
  ll_0 <- log_likelihood(
    obs_is_mut_list = obs_is_mut_list,
    error_mut_to_ref_list = error_mut_to_ref_list,
    error_ref_to_mut_list = error_ref_to_mut_list,
    r = r_est,
    tf = 0
  )
  Q_0 <- -2 * (ll_0 - ll_max)

  if (Q_0 <= qchisq(1 - alpha, 1) | tf_est == 0) {
    tf_min <- 0
  } else {
    tf_min <- uniroot(function(x) {
      ll_A <- log_likelihood(
        obs_is_mut_list = obs_is_mut_list,
        error_mut_to_ref_list = error_mut_to_ref_list,
        error_ref_to_mut_list = error_ref_to_mut_list,
        r = r_est,
        tf = x
      )
      Q_A <- -2 * (ll_A - ll_max)
      return(Q_A - qchisq(1 - alpha, 1))
    },
    interval = c(0, tf_est)
    )$root
  }

  # Find maximum
  ll_tf_2 <- log_likelihood(
    obs_is_mut_list = obs_is_mut_list,
    error_mut_to_ref_list = error_mut_to_ref_list,
    error_ref_to_mut_list = error_ref_to_mut_list,
    r = r_est,
    tf = 2
  )
  Q_tf_2 <- -2 * (ll_tf_2 - ll_max)

  if (Q_tf_2 <= qchisq(1 - alpha, 1) | tf_est == 2) {
    tf_max <- 2
  } else {
    upper_start <- min(tf_est * 2 + 1e-8, 2)
    tf_max <- uniroot(function(x) {
      ll_A <- log_likelihood(
        obs_is_mut_list = obs_is_mut_list,
        error_mut_to_ref_list = error_mut_to_ref_list,
        error_ref_to_mut_list = error_ref_to_mut_list,
        r = r_est,
        tf = x
      )
      Q_A <- -2 * (ll_A - ll_max)
      return(Q_A - qchisq(1 - alpha, 1))
    },
    interval = c(tf_est, upper_start),
    extendInt = "upX"
    )$root
  }

  return(list(
    tf_min = tf_min,
    tf_max = tf_max
  ))
}

#' @importFrom stats qchisq uniroot
get_r_CI <- function(obs_is_mut_list, error_mut_to_ref_list, error_ref_to_mut_list, r_est, tf_est, alpha) {
  ll_max <- log_likelihood(
    obs_is_mut_list = obs_is_mut_list,
    error_mut_to_ref_list = error_mut_to_ref_list,
    error_ref_to_mut_list = error_ref_to_mut_list,
    r = r_est,
    tf = tf_est
  )

  # Find minimum
  ll_0 <- log_likelihood(
    obs_is_mut_list = obs_is_mut_list,
    error_mut_to_ref_list = error_mut_to_ref_list,
    error_ref_to_mut_list = error_ref_to_mut_list,
    r = 0,
    tf = tf_est
  )
  Q_0 <- -2 * (ll_0 - ll_max)

  if (Q_0 <= qchisq(1 - alpha, 1)) {
    r_min <- 0
  } else {
    r_min <- uniroot(function(x) {
      ll_A <- log_likelihood(
        obs_is_mut_list = obs_is_mut_list,
        error_mut_to_ref_list = error_mut_to_ref_list,
        error_ref_to_mut_list = error_ref_to_mut_list,
        r = x,
        tf = tf_est
      )
      Q_A <- -2 * (ll_A - ll_max)
      return(Q_A - qchisq(1 - alpha, 1))
    },
    interval = c(0, r_est)
    )$root
  }

  # Find maximum
  ll_1 <- log_likelihood(
    obs_is_mut_list = obs_is_mut_list,
    error_mut_to_ref_list = error_mut_to_ref_list,
    error_ref_to_mut_list = error_ref_to_mut_list,
    r = 1,
    tf = tf_est
  )
  Q_1 <- -2 * (ll_1 - ll_max)

  if (Q_1 <= qchisq(1 - alpha, 1)) {
    r_max <- 1
  } else {
    r_max <- uniroot(function(x) {
      ll_A <- log_likelihood(
        obs_is_mut_list = obs_is_mut_list,
        error_mut_to_ref_list = error_mut_to_ref_list,
        error_ref_to_mut_list = error_ref_to_mut_list,
        r = x,
        tf = tf_est
      )
      Q_A <- -2 * (ll_A - ll_max)
      return(Q_A - qchisq(1 - alpha, 1))
    },
    interval = c(r_est, 1)
    )$root
  }

  return(list(
    r_min = r_min,
    r_max = r_max
  ))
}
