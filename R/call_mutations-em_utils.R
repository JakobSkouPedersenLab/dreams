em_update_vc <- function(par, X_list, error_ref_to_mut_list, error_mut_to_ref_list) {
  fixed_r <- 1
  em_update(c(par, fixed_r), X_list, error_ref_to_mut_list, error_mut_to_ref_list)
}

em_objective_vc <- function(par, X_list, error_ref_to_mut_list, error_mut_to_ref_list) {
  fixed_r <- 1
  em_objective(c(par, fixed_r), X_list, error_ref_to_mut_list, error_mut_to_ref_list)
}
