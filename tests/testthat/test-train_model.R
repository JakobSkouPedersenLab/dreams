

########## TEST NEW SETUP #######

test_that("separate train model example", {
  read_example_training_data <- system.file("extdata", "training_data.csv", package = "dreams")

  training_data <- read.csv(read_example_training_data)

  prepared_training_data <- prepare_training_data(training_data,
    model_features = c("ref", "strand")
  )

  input <- prepare_input_layer(prepared_training_data$features,
    ctx3_embed_dim = 3
  )

  NN_model <- generate_NN_structure(input$inputs,
    input_layer = input$input_layer,
    layers = c(8, 4, 2),
    reg = 0
  )



  model <- fit_model(
    features = prepared_training_data$features,
    labels = prepared_training_data$labels,
    input_structure = NN_model,
    lr = 0.0005,
    batch_size = 10,
    epochs = 5,
    model_file_path = "~/Desktop/model.test",
    #log_file_path = "~/Desktop/log.test",
    decay = 0, min_delta = 0, patience = 0, validation_split = 0.1
  )

  predict(
    model,
    prepared_training_data$features
  )

})

test_that("Train model example", {
  read_example_training_data <- system.file("extdata", "training_data.csv", package = "dreams")

  training_data <- read.csv(read_example_training_data) %>% select(ref, obs)

  model <- train_model(
    training_data = training_data,
    model_features = c("ref"),
    layers = c(8, 4, 2),
    lr = 0.0005,
    batch_size = 10,
    epochs = 5,
    model_file_path = "~/Desktop/model.test",
    log_file_path = "~/Desktop/log.test",
    validation_split = 0.1
  )

  predict(
    model,
    training_data
  )
})


test_that("Train model example with early stopping", {
  read_example_training_data <- system.file("extdata", "training_data.csv", package = "dreams")

  training_data <- read.csv(read_example_training_data) %>% select(ref, obs)

  model <- train_model(
    training_data = training_data,
    model_features = c("ref"),
    layers = c(8, 4, 2),
    lr = 0.0005,
    batch_size = 10,
    epochs = 5,
    model_file_path = "~/Desktop/model.test",
    log_file_path = "~/Desktop/log.test",
    validation_split = 0.1,
    min_delta = 0.001
  )

  predict(
    model,
    training_data
  )
})

test_that("Train model example with several feature types", {
  read_example_training_data <- system.file("extdata", "training_data.csv", package = "dreams")

  training_data <- read.csv(read_example_training_data) %>% select(ref, read_index, trinucleotide_ctx, obs)

  model <- train_model(
    training_data = training_data,
    model_features = c("ref", "read_index", "trinucleotide_ctx"),
    layers = c(8, 4, 2),
    lr = 0.0005,
    batch_size = 10,
    epochs = 5,
    model_file_path = "~/Desktop/model.test",
    log_file_path = "~/Desktop/log.test",
    validation_split = 0.1
  )

  predict(
    model,
    training_data
  )
})

test_that("Train model example with one layer", {
  read_example_training_data <- system.file("extdata", "training_data.csv", package = "dreams")

  training_data <- read.csv(read_example_training_data) %>% select(ref, read_index, trinucleotide_ctx, obs)

  model <- train_model(
    training_data = training_data,
    model_features = c("ref", "read_index", "trinucleotide_ctx"),
    layers = c(8),
    lr = 0.0005,
    batch_size = 10,
    epochs = 5,
    model_file_path = "~/Desktop/model.test",
    log_file_path = "~/Desktop/log.test",
    validation_split = 0.1
  )

  predict(
    model,
    training_data
  )
})


test_that("simple example", {

  training_data = data.frame(ref = c("A", "A", "A", "A"),
                             obs = c("A", "C", "G", "T"))

  model <- train_model(
    training_data = training_data,
    model_features = c("ref"),
    layers = c(8, 4, 2),
    lr = 0.0005,
    batch_size = 4,
    epochs = 5,
    model_file_path = "~/Desktop/model.test",
    log_file_path = "~/Desktop/log.test",
    validation_split = 0
  )

  predict(
    model,
    training_data
  )
})


#
# test_that("separate train simple example", {
#   training_data = data.frame(ref = c("A", "A", "A", "A"),
#                              obs = c("A", "C", "G", "T"))
#
#   prepared_training_data <- prepare_training_data(training_data,
#                                                   model_features = c("ref")
#   )
#
#   input <- prepare_input_layer(prepared_training_data$features,
#                                ctx3_embed_dim = 3
#   )
#
#   NN_model <- generate_NN_structure(input$inputs,
#                                     input_layer = input$input_layer,
#                                     layers = c(8, 4, 2),
#                                     reg = 0
#   )
#
#
#
#   model <- fit_model(
#     features = prepared_training_data$features,
#     labels = prepared_training_data$labels,
#     input_structure = NN_model,
#     lr = 0.01,
#     batch_size = 4,
#     epochs = 20,
#     model_file_path = "~/Desktop/model.test",
#     log_file_path = "~/Desktop/log.test",
#     decay = 0, min_delta = 0, patience = 0, validation_split = 0.0
#   )
#
#   prediction = predict(
#     model,
#     prepared_training_data$features
#   )
#
#
#
#   expect_true(between(prediction[1,1], 0.24, 0.26))
#
#
# })
