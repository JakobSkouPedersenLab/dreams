
#' Train error model
#'
#' @param training_data \code{data.frame} Input training data (Generated from
#'   [get_training_data()])
#' @param layers Numeric vector. Number of nodes in each layer.
#' @param model_features Vector of feature names. Selected features for model
#'   training.
#' @param lr Numeric value between 0 and 1. Learning rate.
#' @param batch_size Integer. Batch size.
#' @param epochs Integer. Number of training epochs.
#' @param model_file_path String. Model output file path. Default is
#'   \code{NULL}.
#' @param log_file_path String. Path to model output log file. Default is
#'   \code{NULL}.
#' @param l2_reg Numeric value between 0 and 1. Level of L2 regularization per
#'   layer. Default is 0.
#' @param min_delta Numeric value between 0 and 1. Minimum delta for early
#'   stopping. Default is 0.
#' @param patience Integer. Patience when reaching minimum delta. Default is 0.
#' @param validation_split Numeric value between 0 and 1. Validation split
#'   ratio. Default is 0.
#' @param ctx3_embed_dim Integer. Number of dimensions to embed trinucleotide
#'   context to. Default is 3.
#'
#' @return Trained model in hdf5 format.
#' @keywords model training
#' @importFrom keras save_model_hdf5
#' @export
#' @family Train model
#' @seealso [get_training_data()] Function for getting training data
train_dreams_model_indels <- function(training_data, layers,
                        model_features, lr, batch_size, epochs,
                        model_file_path = NULL, log_file_path = NULL,
                        min_delta = 0, patience = 0, l2_reg = 0,
                        validation_split = 0, ctx3_embed_dim = 3) {
  training_data <- training_data$data %>%
    filter(.data$obs %in% c("A", "T", "C", "G", "D", "I"))

  training_data <- prepare_training_data_indels(
    training_data = training_data,
    model_features = model_features
  )

  prepared_input <- prepare_input_layer_indels(training_data$features,
    ctx3_embed_dim = ctx3_embed_dim
  )

  model_structure <- generate_NN_structure_indels(
    inputs = prepared_input$inputs,
    input_layer = prepared_input$input_layer,
    layers = layers,
    reg = l2_reg
  )

  model <- fit_model_indels(
    features = training_data$features,
    labels = training_data$labels,
    input_structure = model_structure,
    lr = lr,
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




#' Prepare Training Data with Indel Information
#'
#' This function processes training data by selecting specific features and the
#' observed values (obs) for genomic data, including nucleotide bases and indels
#' (insertions 'I' and deletions 'D'). It then prepares this data for use in
#' machine learning models, particularly those that require categorical data to
#' be converted into a one-hot encoded format.
#'
#' @param training_data A data frame containing the training data.
#' @param model_features A vector of selected feature names to be used in the
#'   model.
#'
#' @return A list with two elements: `features` containing the selected features
#'   of the training data, and `labels` containing the one-hot encoded labels.
#' @keywords internal
#' @importFrom keras to_categorical
#'
prepare_training_data_indels <- function(training_data, model_features) {

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
    mutate(obs = case_when(
      .data$obs %in% c("A", "T", "C", "G") ~ 0,  # ATCG as one level
      .data$obs == "D" ~ 1,                       # D as another level
      .data$obs == "I" ~ 2                        # I as another level
    )) %>%
    as.matrix() %>%
    keras::to_categorical()

  return(list(
    features = training_data_features,
    labels = training_data_labels
  ))
}


#' Prepare Input Layer for Indel Analysis in Neural Network
#'
#' This function processes training data features and context embedding
#' dimensions to prepare an input layer suitable for a neural network model in
#' indel analysis.
#'
#' @param training_data_features A data frame or similar structure containing
#'   features of the training data.
#' @param ctx3_embed_dim Integer specifying the dimension of the embedding for
#'   trinucleotide context sequences.
#'
#' @return A list containing two elements: `inputs`, the raw input from the
#'   dataset, and `input_layer`, the prepared input layer for the neural network
#'   model.
#' @keywords internal
#'
#' @importFrom keras fit layer_dense_features layer_batch_normalization
#'   layer_concatenate
#' @importFrom tfdatasets feature_spec step_numeric_column
#'   step_categorical_column_with_vocabulary_list step_embedding_column
#'   step_indicator_column
#'
prepare_input_layer_indels <- function(training_data_features, ctx3_embed_dim) {

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
      "prior_error",
      "seq_length"
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


#' Generate Neural Network Structure for Indel Analysis
#'
#' This function builds a neural network structure for indel analysis. It takes
#' an input layer and iteratively adds specified hidden layers, followed by an
#' output layer for classification.
#'
#' @param inputs The input structure for the neural network, usually from the data input layer.
#' @param input_layer The initial layer of the neural network to which additional layers are added.
#' @param layers A list or vector indicating the number of units in each hidden layer of the network.
#' @param reg Regularization parameter applied to each dense layer for preventing overfitting.
#'        Defaults to 0.
#'
#' @return A keras model object representing the constructed neural network.
#' @keywords internal
#' @importFrom keras layer_dense
#'
generate_NN_structure_indels <- function(inputs, input_layer, layers, reg = 0) {
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
    units = 3,
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


fit_model_indels <- function(features, labels, input_structure,
                      lr,
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
      filepath = paste0(model_file_path, "_checkpoint"),
      monitor = "val_loss",
      save_freq = 10
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
    learning_rate = lr
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
