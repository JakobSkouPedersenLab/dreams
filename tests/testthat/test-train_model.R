
library(keras)
library(tfdatasets)
library(readr)

negative_data <- read_csv("~/GenomeDK/ctDNA/read_level_prediction/combined_data/wf0.4/duplex/post_op_5_2_repeated_patientCV_1to10_patient_filtered/repeat_1/6301/negative.csv") %>% sample_n(100)
positive_data <- read_csv("~/GenomeDK/ctDNA/read_level_prediction/combined_data/wf0.4/duplex/post_op_5_2_repeated_patientCV_1to10_patient_filtered/repeat_1/6301/positive.csv") %>% sample_n(100)

training_data <- rbind(negative_data, positive_data)

model <- train_model(
  training_data = training_data,
  model_file_path = "~/ctDNA/model",
  log_file_path = "~/ctDNA/test.log",
  model_features = c("ref"),
  lr = 0.0005,
  decay = 0,
  dropout_rate = 0,
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
