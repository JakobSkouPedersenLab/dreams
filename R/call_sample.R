
make_sample_calls <- function(bam_file, reference_path, layers,
                              model_features, lr, batch_size, epochs,
                              mutations_df) {
  training_data <- get_training_data(
    bam_paths = bam_file,
    reference_path = reference_path,
    positions_to_exclude_paths = mutations_df
    )

  model <- train_dreams_model(
    training_data,
    layers = layers,
    model_features = model_features,
    lr = lr,
    batch_size = batch_size,
    epochs = epochs,
    ...
  )

  # Call variants using DREAMS-vc
  variant_calls <- dreams_vc(
    bam_path = bam_file,
    mutations_df = mutations_df,
    model = model
  )

  # Call cancer using DREAMS-cc
  cancer_calls <- dreams_cc(
    bam_path = bam_file,
    positions = mutations_df,
    model = model
  )

  return(list(
    "variant_calls" = variant_calls,
    "cancer_calls" = cancer_calls
  ))
}
