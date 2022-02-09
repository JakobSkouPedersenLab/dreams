predict_error_rates <- function(reads, model, beta) {
  prediction <- model %>%
    predict(reads) %>%
    data.frame() %>%
    rename(
      A = X1,
      T = X2,
      C = X3,
      G = X4
    )

  prediction$ref <- reads$ref

  corrected_errors <- prediction %>%
    correct_errors_predictions(beta = beta) %>%
    select(-ref)

  predicted_errors <- bind_cols(reads, corrected_errors)

  return(predicted_errors)
}
