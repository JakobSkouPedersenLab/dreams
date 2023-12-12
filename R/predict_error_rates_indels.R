limit <- function(x, min = 10^-13, max = 0.999999999999) {
  return(pmax(pmin(x, max), min))
}

#' @import dplyr
correct_errors_predictions_indels <- function(error_df, beta) {
  error_df %>%
    mutate(
      error_prob_sample = case_when(
        .data$ref %in% c("A", "T", "C", "G") ~ 1 - ATCG,
        .data$ref == "D" ~ 1 - D,
        .data$ref == "I" ~ 1 - I
        ),
      error_prob_corrected = beta * .data$error_prob_sample / (beta * .data$error_prob_sample - .data$error_prob_sample + 1),
      ATCG_corrected = ifelse(.data$ref %in% c("A", "T", "C", "G"),
                              1.00 - .data$error_prob_corrected,
                              .data$error_prob_corrected * .data$ATCG / (.data$error_prob_sample)
      ),
      D_corrected = ifelse(.data$ref == "D",
                           1.00 - .data$error_prob_corrected,
                           .data$error_prob_corrected * .data$D / (.data$error_prob_sample)
      ),
      I_corrected = ifelse(.data$ref == "I",
                           1.00 - .data$error_prob_corrected,
                           .data$error_prob_corrected * .data$I / (.data$error_prob_sample)
      ),
      ATCG_corrected = ifelse(is.nan(.data$ATCG_corrected), 0, .data$ATCG_corrected) %>% limit(),
      D_corrected = ifelse(is.nan(.data$D_corrected), 0, .data$D_corrected) %>% limit(),
      I_corrected = ifelse(is.nan(.data$I_corrected), 0, .data$I_corrected) %>% limit()
    )
}

#' @importFrom stats predict
predict_error_rates_indels <- function(read_positions_df,
                                       model, beta) {

  # Predict error rates for read positions from trained DREAM model
  if (nrow(read_positions_df) == 0) {
    predictions_df <-
      data.frame(
        ATCG = numeric(),
        D = numeric(),
        I = numeric()
      )
  } else {
    predictions_df <- model %>%
      predict(read_positions_df) %>%
      data.frame() %>%
      rename(
        ATCG = .data$X1,  # Assuming the first column represents the combined probabilities for A, T, C, G
        D = .data$X2,
        I = .data$X3
      )
  }

  # Add predictions to read positions and correct error rates with bea
  predicted_errors <-
    read_positions_df %>%
    bind_cols(predictions_df) %>%
    correct_errors_predictions_indels(beta = beta)

  return(predicted_errors)
}
