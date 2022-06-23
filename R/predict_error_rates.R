limit <- function(x, min = 0, max = 1) {
  return(pmax(pmin(x, max), min))
}

#' @import dplyr
correct_errors_predictions <- function(error_df, beta) {
  error_df %>%
    mutate(
      error_prob_sample = case_when(
        .data$ref == "A" ~ 1 - A,
        .data$ref == "C" ~ 1 - C,
        .data$ref == "T" ~ 1 - T,
        .data$ref == "G" ~ 1 - G
      ),
      error_prob_corrected = beta * .data$error_prob_sample / (beta * .data$error_prob_sample - .data$error_prob_sample + 1),
      A_corrected = ifelse(.data$ref == "A",
        1 - .data$error_prob_corrected,
        .data$error_prob_corrected * .data$A / (.data$error_prob_sample)
      ),
      C_corrected = ifelse(.data$ref == "C",
        1 - .data$error_prob_corrected,
        .data$error_prob_corrected * .data$C / (.data$error_prob_sample)
      ),
      G_corrected = ifelse(.data$ref == "G",
        1 - .data$error_prob_corrected,
        .data$error_prob_corrected * .data$G / (.data$error_prob_sample)
      ),
      T_corrected = ifelse(.data$ref == "T",
        1 - .data$error_prob_corrected,
        .data$error_prob_corrected * .data$T / (.data$error_prob_sample)
      ),
      A_corrected = ifelse(is.nan(.data$A_corrected), 0, .data$A_corrected) %>% limit(),
      C_corrected = ifelse(is.nan(.data$C_corrected), 0, .data$C_corrected) %>% limit(),
      G_corrected = ifelse(is.nan(.data$G_corrected), 0, .data$G_corrected) %>% limit(),
      T_corrected = ifelse(is.nan(.data$T_corrected), 0, .data$T_corrected) %>% limit()
    )
}

#' @importFrom stats predict
predict_error_rates <- function(read_positions_df, model, beta, batch_size) {

  # Predict error rates for read positions from trained DREAM model
  if (nrow(read_positions_df) == 0) {
    predictions_df <-
      data.frame(
        A = numeric(),
        T = numeric(),
        C = numeric(),
        G = numeric()
      )
  } else {
    predictions_df <- model %>%
      predict(read_positions_df, batch_size = batch_size) %>%
      data.frame() %>%
      rename(
        A = .data$X1,
        T = .data$X2,
        C = .data$X3,
        G = .data$X4
      )
  }

  # Add predictions to read positions and correct error rates with bea
  predicted_errors <-
    read_positions_df %>%
    bind_cols(predictions_df) %>%
    correct_errors_predictions(beta = beta)

  return(predicted_errors)
}
