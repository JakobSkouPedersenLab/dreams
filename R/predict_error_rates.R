limit_0_1 <- function(x) {
  return(pmax(pmin(x, 1), 0))
}

correct_errors_predictions <- function(error_df, beta) {
  error_df %>%
    mutate(
      error_prob_sample = case_when(
        ref == "A" ~ 1 - A,
        ref == "C" ~ 1 - C,
        ref == "T" ~ 1 - T,
        ref == "G" ~ 1 - G
      ),
      error_prob_corrected = beta * error_prob_sample / (beta * error_prob_sample - error_prob_sample + 1),
      A_corrected = ifelse(ref == "A",
        1 - error_prob_corrected,
        error_prob_corrected * A / (error_prob_sample)
      ),
      C_corrected = ifelse(ref == "C",
        1 - error_prob_corrected,
        error_prob_corrected * C / (error_prob_sample)
      ),
      G_corrected = ifelse(ref == "G",
        1 - error_prob_corrected,
        error_prob_corrected * G / (error_prob_sample)
      ),
      T_corrected = ifelse(ref == "T",
        1 - error_prob_corrected,
        error_prob_corrected * T / (error_prob_sample)
      ),
      A_corrected = limit_0_1(ifelse(is.nan(A_corrected), 0, A_corrected)),
      C_corrected = limit_0_1(ifelse(is.nan(C_corrected), 0, C_corrected)),
      G_corrected = limit_0_1(ifelse(is.nan(G_corrected), 0, G_corrected)),
      T_corrected = limit_0_1(ifelse(is.nan(T_corrected), 0, T_corrected))
    )
}

#' Title
#'
#' @param read_positions_df
#' @param model
#' @param beta
#'
#' @return
#' @export
#' @importFrom modelr add_predictions
predict_error_rates <- function(read_positions_df, model, beta) {
  # Predict error rates for read positions from trained DREAM model

  # TODO: add modelr

  prediction <- model %>%
    predict(read_positions_df) %>%
    data.frame() %>%
    rename(
      A = X1,
      T = X2,
      C = X3,
      G = X4
    )

  prediction$ref <- read_positions_df$ref

  corrected_errors <- prediction %>%
    correct_errors_predictions(beta = beta) %>%
    select(-ref)

  predicted_errors <- bind_cols(read_positions_df, corrected_errors)

  return(predicted_errors)
}
