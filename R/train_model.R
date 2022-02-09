#
# train_model <- function(training_data,
#                         model_file_path, log_file_path,
#                         model_features, lr, decay, dropout_rate, batch_size,
#                         epochs, l2_reg, L1_nodes, L2_nodes, L3_nodes,
#                         min_delta, patience, validation_split, ctx3_embed_dim) {
#   LogMetrics <- R6::R6Class(
#     "LogMetrics",
#     inherit = KerasCallback,
#     public =
#       list(
#         log_path = NULL,
#         init_time = Sys.time(),
#         initialize = function(log_path) {
#           self$log_path <- log_path
#         },
#         on_epoch_end = function(epoch, logs = list()) {
#
#           # Write log
#           log_df <- data.frame(
#             epoch = epoch,
#             time = difftime(Sys.time(), self$init_time, units = "secs"),
#             data.frame(logs),
#             lr = k_get_value(model$optimizer$lr)
#           )
#
#           if (epoch == 0) {
#             write_csv(x = log_df, file = self$log_path)
#           } else {
#             write_csv(x = log_df, file = self$log_path, append = TRUE)
#           }
#         }
#       )
#   )
#
#   # Load training data
#   training_data <- training_data %>%
#     select(all_of(model_features), obs) %>%
#     sample_frac(1)
#
#   # Split data in features and labels
#   training_data_features <-
#     training_data %>%
#     select(-obs)
#
#   training_data_labels <-
#     training_data %>%
#     select(obs) %>%
#     mutate(obs = as.numeric(factor(obs, levels = c("A", "T", "C", "G"))) - 1) %>%
#     as.matrix() %>%
#     to_categorical()
#
#
#   # Subset selected features into categorical and numeric
#   all_numerical_variables <-
#     c(
#       "read_index",
#       "fragment_size",
#       "local_GC",
#       "umi_count",
#       "umi_errors",
#       "local_complexity_1",
#       "local_complexity_2",
#       "n_other_errors"
#     )
#   all_categorical_variables <-
#     c(
#       "ref",
#       "strand",
#       "first_in_pair",
#       "ctx_minus1",
#       "ctx_plus1"
#     )
#   all_embedded_variables_ctx3 <-
#     c(
#       "trinucleotide_ctx"
#     )
#
#   numerical_variables <- intersect(model_features, all_numerical_variables)
#   categorical_variables <- intersect(model_features, all_categorical_variables)
#   embedded_variables_ctx3 <- intersect(model_features, all_embedded_variables_ctx3)
#
#   print("CREATE SPEC")
#   # Create feature spec
#   ft_spec_numeric <- training_data_features %>%
#     feature_spec(all_of(model_features)) %>%
#     step_numeric_column(
#       all_of(numerical_variables)
#     ) %>%
#     fit()
#
#   ctx3_vocabulary <-
#     expand.grid(
#       replicate(3, c("A", "T", "C", "G"), simplify = FALSE)
#     ) %>% apply(1, paste0, collapse = "")
#
#   ft_spec_embed_ctx3 <- training_data_features %>%
#     feature_spec(all_of(model_features)) %>%
#     step_categorical_column_with_vocabulary_list(
#       all_of(embedded_variables_ctx3),
#       vocabulary_list = ctx3_vocabulary
#     ) %>%
#     step_embedding_column(
#       all_of(embedded_variables_ctx3),
#       dimension = ctx3_embed_dim
#     ) %>%
#     fit()
#
#   ft_spec_categorical <- training_data_features %>%
#     feature_spec(all_of(model_features)) %>%
#     step_categorical_column_with_vocabulary_list(
#       all_of(categorical_variables)
#     ) %>%
#     step_indicator_column(
#       all_of(categorical_variables)
#     ) %>%
#     fit()
#
#
#   # Build NN model
#   inputs <- layer_input_from_dataset(training_data_features)
#
#   input_layer_list <- list()
#
#   if (!is.null(ft_spec_numeric$dense_features())) {
#     numerical_layer <- inputs %>%
#       layer_dense_features(
#         name = "input_numeric",
#         ft_spec_numeric$dense_features()
#       ) %>%
#       layer_batch_normalization()
#
#     input_layer_list <- append(input_layer_list, numerical_layer)
#   }
#
#
#   if (!is.null(ft_spec_embed_ctx3$dense_features())) {
#     embed_layer_ctx3 <- inputs %>%
#       layer_dense_features(
#         name = "input_embed_ctx3",
#         ft_spec_embed_ctx3$dense_features(),
#         weights = list(matrix(0.01, nrow = length(ctx3_vocabulary), ncol = ctx3_embed_dim))
#       )
#
#     input_layer_list <- append(input_layer_list, embed_layer_ctx3)
#   }
#
#
#   if (!is.null(ft_spec_categorical$dense_features())) {
#     categorical_layer <- inputs %>%
#       layer_dense_features(
#         name = "input_categorical",
#         ft_spec_categorical$dense_features()
#       )
#
#     input_layer_list <- append(input_layer_list, categorical_layer)
#   }
#
#   # Concatenate inputs if necessary
#   if (length(input_layer_list) > 1) {
#     input_layer <- layer_concatenate(input_layer_list)
#   } else {
#     input_layer <- input_layer_list[[1]]
#   }
#
#   outputs <- input_layer %>%
#     layer_dense(
#       name = "hidden_1",
#       units = L1_nodes,
#       activation = "relu",
#       kernel_regularizer = regularizer_l2(l2_reg),
#     ) %>%
#     layer_dropout(rate = dropout_rate) %>%
#     layer_dense(
#       name = "hidden_2",
#       units = L2_nodes,
#       activation = "relu",
#       kernel_regularizer = regularizer_l2(l2_reg)
#     ) %>%
#     layer_dense(
#       name = "hidden_3",
#       units = L3_nodes,
#       activation = "relu",
#       kernel_regularizer = regularizer_l2(l2_reg)
#     ) %>%
#     layer_dropout(rate = dropout_rate) %>%
#     layer_dense(
#       name = "output",
#       units = 4,
#       activation = "softmax"
#     )
#
#   model <- keras_model(inputs = inputs, outputs = outputs)
#
#   # Callbacks
#   # Create a training checkpoint for every 10 epochs
#   checkpoint_callback <-
#     callback_model_checkpoint(
#       filepath = paste0(model_file_path, "_checkpoint"),
#       monitor = "val_loss",
#       save_freq = 10
#     )
#
#   # Log metrics and time
#   logger_callback <- LogMetrics$new(log_path = log_file_path)
#
#   # Collect callback in list
#   callback_list <- list(
#     logger_callback,
#     checkpoint_callback
#   )
#
#   # Add early stopping callback
#   if (min_delta > 0) {
#     early_stopping <-
#       callback_early_stopping(
#         monitor = "val_loss",
#         min_delta = min_delta,
#         patience = patience,
#         verbose = 0,
#         restore_best_weights = TRUE
#       )
#     callback_list <- append(callback_list, early_stopping)
#   }
#
#   print("COMPILE AND TRAIN MODEL")
#
#   # Optimizer
#   opt <- optimizer_adam(
#     learning_rate = lr,
#     decay = 0
#   )
#
#   # Compile model
#   model %>%
#     compile(
#       loss = "categorical_crossentropy",
#       optimizer = opt,
#       metrics = "accuracy"
#     )
#
#   # Train model
#   model %>% fit(
#     x = training_data_features,
#     y = training_data_labels,
#     validation_split = validation_split,
#     epochs = epochs,
#     batch_size = batch_size,
#     verbose = 0,
#     callbacks = callback_list
#   )
#
#   # Save final model
#   save_model_hdf5(
#     object = model,
#     filepath = model_file_path
#   )
#
#   print("DONE TRAINING")
#
#
#   return (model)
# }
