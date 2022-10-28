
#' Train error model
#'
#' @param training_data \code{data.frame} Input training data (Generated from [get_training_data()])
#' @param layers Numeric vector. Number of nodes in each layer.
#' @param model_features Vector of feature names. Selected features for model training.
#' @param lr Numeric value between 0 and 1. Learning rate.
#' @param batch_size Integer. Batch size.
#' @param epochs Integer. Number of training epochs.
#' @param model_file_path String. Model output file path. Default is \code{NULL}.
#' @param log_file_path String. Path to model output log file. Default is \code{NULL}.
#' @param decay Numeric value between 0 and 1. Decay rate. Default is 0.
#' @param l2_reg Numeric value between 0 and 1. Level of L2 regularization per layer. Default is 0.
#' @param min_delta Numeric value between 0 and 1. Minimum delta for early stopping. Default is 0.
#' @param patience Integer. Patience when reaching minimum delta. Default is 0.
#' @param validation_split Numeric value between 0 and 1. Validation split ratio. Default is 0.
#' @param ctx3_embed_dim Integer. Number of dimensions to embed trinucleotide context to. Default is 3.
#'
#' @return Trained model in hdf5 format.
#' @keywords model training
#' @importFrom keras save_model_hdf5
#' @export
#' @family Train model
#' @seealso [get_training_data()] Function for getting training data
train_dreams_model <- function(training_data, layers,
                        model_features, lr, batch_size, epochs,
                        model_file_path = NULL, log_file_path = NULL,
                        decay = 0, min_delta = 0, patience = 0, l2_reg = 0,
                        validation_split = 0, ctx3_embed_dim = 3) {
  training_data <- prepare_training_data(
    training_data = training_data,
    model_features = model_features
  )

  prepared_input <- prepare_input_layer(training_data$features,
    ctx3_embed_dim = ctx3_embed_dim
  )

  model_structure <- generate_NN_structure(
    inputs = prepared_input$inputs,
    input_layer = prepared_input$input_layer,
    layers = layers,
    reg = l2_reg
  )

  model <- fit_model(
    features = training_data$features,
    labels = training_data$labels,
    input_structure = model_structure,
    lr = lr,
    decay = decay,
    batch_size = batch_size,
    epochs = epochs,
    min_delta = min_delta,
    patience = patience,
    validation_split = validation_split,
    model_file_path = model_file_path,
    log_file_path = log_file_path
  )

  if (!is.null(model_file_path)) {
    # Save final model
    keras::save_model_hdf5(
      object = model,
      filepath = model_file_path
    )
  }

  print("DONE TRAINING")


  return(model)
}




#' Prepare training data
#'
#' @param training_data training data
#' @param model_features selected features
#'
#' @return prepared training data
#' @keywords internal
#' @importFrom keras to_categorical

prepare_training_data <- function(training_data, model_features) {

  # Load training data
  training_data <- training_data %>%
    select(all_of(model_features), .data$obs) %>%
    sample_frac(1)

  # Split data in features and labels
  training_data_features <-
    training_data %>%
    select(-.data$obs)

  training_data_labels <-
    training_data %>%
    select(.data$obs) %>%
    mutate(obs = as.numeric(factor(.data$obs, levels = c("A", "T", "C", "G"))) - 1) %>%
    as.matrix() %>%
    keras::to_categorical()

  return(list(
    features = training_data_features,
    labels = training_data_labels
  ))
}


#' Title
#'
#' @param training_data_features training data features
#' @param ctx3_embed_dim trinucleotide context embedding dimensions
#'
#' @return input layer
#' @keywords internal
#'
#' @importFrom keras fit layer_dense_features layer_batch_normalization layer_concatenate
#' @importFrom tfdatasets feature_spec step_numeric_column step_categorical_column_with_vocabulary_list step_embedding_column step_indicator_column
prepare_input_layer <- function(training_data_features, ctx3_embed_dim) {


  # Subset selected features into categorical and numeric
  all_numerical_variables <-
    c(
      "read_index",
      "fragment_size",
      "local_GC",
      "umi_count",
      "umi_errors",
      "local_complexity_1",
      "local_complexity_2",
      "n_other_errors",
      "prior_error"
    )
  all_categorical_variables <-
    c(
      "ref",
      "strand",
      "first_in_pair",
      "ctx_minus1",
      "ctx_plus1",
      "chr",
      "genomic_pos"
    )
  all_embedded_variables_ctx3 <-
    c(
      "trinucleotide_ctx"
    )


  model_features <- names(training_data_features)

  numerical_variables <- intersect(model_features, all_numerical_variables)
  categorical_variables <- intersect(model_features, all_categorical_variables)
  embedded_variables_ctx3 <- intersect(model_features, all_embedded_variables_ctx3)

  # Create feature spec
  ft_spec_numeric <- training_data_features %>%
    tfdatasets::feature_spec(all_of(model_features)) %>%
    tfdatasets::step_numeric_column(
      all_of(numerical_variables)
    ) %>%
    keras::fit()

  ctx3_vocabulary <-
    expand.grid(
      replicate(3, c("A", "T", "C", "G"), simplify = FALSE)
    ) %>% apply(1, paste0, collapse = "")

  ft_spec_embed_ctx3 <- training_data_features %>%
    tfdatasets::feature_spec(all_of(model_features)) %>%
    tfdatasets::step_categorical_column_with_vocabulary_list(
      all_of(embedded_variables_ctx3),
      vocabulary_list = ctx3_vocabulary
    ) %>%
    tfdatasets::step_embedding_column(
      all_of(embedded_variables_ctx3),
      dimension = ctx3_embed_dim
    ) %>%
    keras::fit()

  ft_spec_categorical <- training_data_features %>%
    tfdatasets::feature_spec(all_of(model_features)) %>%
    tfdatasets::step_categorical_column_with_vocabulary_list(
      all_of(categorical_variables)
    ) %>%
    tfdatasets::step_indicator_column(
      all_of(categorical_variables)
    ) %>%
    keras::fit()


  # Build NN model
  inputs <- tfdatasets::layer_input_from_dataset(training_data_features)

  input_layer_list <- list()

  if (!is.null(ft_spec_numeric$dense_features())) {
    numerical_layer <- inputs %>%
      keras::layer_dense_features(
        name = "input_numeric",
        ft_spec_numeric$dense_features()
      ) %>%
      keras::layer_batch_normalization()

    input_layer_list <- append(input_layer_list, numerical_layer)
  }


  if (!is.null(ft_spec_embed_ctx3$dense_features())) {
    embed_layer_ctx3 <- inputs %>%
      keras::layer_dense_features(
        name = "input_embed_ctx3",
        ft_spec_embed_ctx3$dense_features(),
        weights = list(matrix(0.01, nrow = length(ctx3_vocabulary), ncol = ctx3_embed_dim))
      )

    input_layer_list <- append(input_layer_list, embed_layer_ctx3)
  }


  if (!is.null(ft_spec_categorical$dense_features())) {
    categorical_layer <- inputs %>%
      keras::layer_dense_features(
        name = "input_categorical",
        ft_spec_categorical$dense_features()
      )

    input_layer_list <- append(input_layer_list, categorical_layer)
  }

  # Concatenate inputs if necessary
  if (length(input_layer_list) > 1) {
    input_layer <- keras::layer_concatenate(input_layer_list)
  } else {
    input_layer <- input_layer_list[[1]]
  }

  return(list(
    inputs = inputs,
    input_layer = input_layer
  ))
}


#' Title
#'
#' @param inputs input structure
#' @param input_layer input layer
#' @param layers layer sizes
#' @param reg regularization
#'
#' @return NN structure
#' @keywords internal
#' @importFrom keras layer_dense

generate_NN_structure <- function(inputs, input_layer, layers, reg = 0) {
  hidden_layers <- input_layer

  for (layer_idx in 1:length(layers)) {
    layer_size <- layers[[layer_idx]]

    hidden_layers <- hidden_layers %>%
      keras::layer_dense(
        name = paste0("hidden_", layer_idx),
        units = layer_size,
        activation = "relu",
        kernel_regularizer = keras::regularizer_l2(reg),
      )
  }

  outputs <- keras::layer_dense(
    object = hidden_layers,
    name = "output",
    units = 4,
    activation = "softmax"
  )

  model <- keras::keras_model(inputs = inputs, outputs = outputs)


  return(model)
}

#' Title
#'
#' @param features input features
#' @param labels input labels
#' @param input_structure input structure
#' @param lr learning rate
#' @param decay decay rate
#' @param batch_size batch size
#' @param epochs number of training epochs
#' @param min_delta minimum delta for early stopping
#' @param patience patience for early stopping
#' @param validation_split validation split
#' @param model_file_path output model file path
#' @param log_file_path output log file path
#'
#' @return fitted model
#' @keywords internal
#' @importFrom keras KerasCallback
#' @importFrom readr write_csv
#' @importFrom R6 R6Class


fit_model <- function(features, labels, input_structure,
                      lr,
                      decay,
                      batch_size,
                      epochs,
                      min_delta,
                      patience,
                      validation_split,
                      model_file_path,
                      log_file_path = NULL) {
  model <- input_structure

  LogMetrics <- R6::R6Class(
    "LogMetrics",
    inherit = keras::KerasCallback,
    public =
      list(
        log_path = NULL,
        init_time = Sys.time(),
        initialize = function(log_path) {
          self$log_path <- log_path
        },
        on_epoch_end = function(epoch, logs = list()) {

          # Write log
          log_df <- data.frame(
            epoch = epoch,
            time = difftime(Sys.time(), self$init_time, units = "secs"),
            data.frame(logs),
            lr = keras::k_get_value(model$optimizer$lr)
          )

          if (epoch == 0) {
            readr::write_csv(x = log_df, file = self$log_path)
          } else {
            readr::write_csv(x = log_df, file = self$log_path, append = TRUE)
          }
        }
      )
  )




  # Callbacks
  # Create a training checkpoint for every 10 epochs
  checkpoint_callback <-
    keras::callback_model_checkpoint(
      filepath = paste0(model_file_path, ".{epoch:02d}_checkpoint.hdf5"),
      monitor = "val_loss",
      save_freq = 1
    )


  callback_list <- list(
    checkpoint_callback
  )

  if (!is.null(log_file_path)) {


    # Log metrics and time
    logger_callback <- LogMetrics$new(log_path = log_file_path)

    # Collect callback in list
    callback_list <- list(
      callback_list,
      logger_callback
    )
  }

  # Add early stopping callback
  if (min_delta > 0) {
    early_stopping <-
      keras::callback_early_stopping(
        monitor = "val_loss",
        min_delta = min_delta,
        patience = patience,
        verbose = 0,
        restore_best_weights = TRUE
      )
    callback_list <- append(callback_list, early_stopping)
  }


  print("COMPILE AND TRAIN MODEL")

  # Optimizer
  opt <- keras::optimizer_adam(
    learning_rate = lr,
    decay = decay
  )

  # Compile model
  model %>%
    keras::compile(
      loss = "categorical_crossentropy",
      optimizer = opt,
      metrics = "accuracy"
    )

  # Train model
  model %>% keras::fit(
    x = features,
    y = labels,
    validation_split = validation_split,
    epochs = epochs,
    batch_size = batch_size,
    verbose = 1,
    callbacks = callback_list
  )

  return(model)
}
