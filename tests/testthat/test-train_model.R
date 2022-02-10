
library(keras)
library(tfdatasets)
library(readr)
library(tidyverse)

########## TEST OLD SETUP #######
read_example_training_data <- system.file("extdata", "training_data.csv", package = "dreams")

training_data <- read.csv(read_example_training_data)


model <- train_model_old(
  training_data = training_data,
  model_file_path = "~/ctDNA/model",
  log_file_path = "~/ctDNA/test.log",
  model_features = c("ref"),
  lr = 0.0005,
  decay = 0,
  batch_size = 10,
  epochs = 5,
  l2_reg = 0,
  L1_nodes = 16,
  L2_nodes = 8,
  L3_nodes = 4,
  min_delta = 0,
  patience = 0,
  validation_split = 0.1,
  ctx3_embed_dim = 3
)

predict(model, training_data %>% select("ref"))



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

  training_data$features
  training_data$labels



  model <- fit_model(
    features = training_data$features,
    labels = training_data$labels,
    input_structure = NN_model,
    lr = 0.0005,
    batch_size = 10,
    epochs = 5,
    model_file_path = "~/Desktop/model.test",
    log_file_path = "~/Desktop/log.test",
    decay = 0, min_delta = 0, patience = 0, validation_split = 0.1
  )

  predict(
    model,
    training_data$features
  )
})

test_that("Train model example", {
  read_example_training_data <- system.file("extdata", "training_data.csv", package = "dreams")

  training_data <- read.csv(read_example_training_data)

  model <- train_model(
    training_data = training_data,
    model_features = c("ref", "strand"),
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
