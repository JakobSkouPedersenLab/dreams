
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DREAMS

<!-- badges: start -->

[![R-CMD-check](https://github.com/JakobSkouPedersenLab/dreams/workflows/R-CMD-check/badge.svg)](https://github.com/JakobSkouPedersenLab/dreams/actions)
[![DOI](https://zenodo.org/badge/455089263.svg)](https://zenodo.org/badge/latestdoi/455089263)
[![Codecov test
coverage](https://codecov.io/gh/JakobSkouPedersenLab/dreams/branch/main/graph/badge.svg)](https://app.codecov.io/gh/JakobSkouPedersenLab/dreams?branch=main)
[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

<!-- badges: end -->

DREAMS is an analysis pipeline that applies a neural network to
low-frequency variant calling and circulating tumor DNA detection from
next-generation DNA sequencing data.

## Installation

You can install the development version of dreams from
[GitHub](https://github.com/JakobSkouPedersenLab/dreams) with:

``` r
# install.packages("devtools")
devtools::install_github("JakobSkouPedersenLab/dreams")
```

### Additional Setup (If needed)

If you encounter any issues related to TensorFlow integrations within R,
install Keras within the correct python environment to ensure a proper
setup:

``` r
keras::install_keras(envname = "<ENVIRONMENT_NAME>")
```

### Usage

After installation, set the environment at the start of each R session:

``` r
reticulate::use_condaenv("<ENVIRONMENT_NAME>", required = TRUE)
```

## Basic Functions

This section provides an overview of the basic functions available in
the dreams library, including data preparation, model training, variant
calling, and cancer detection.

``` r
library(dreams)

# For training, DREAMS requires one or more BAM files and a reference genome.
training_data = get_training_data(
  bam_paths = "/path/bam_file",
  reference_path = "/path/hg38.fa",
  ...)

beta = training_data$info$beta

# Training the DREAMS Model using a Neural Network
# Basic settings for Keras are required.
model = train_dreams_model(
  training_data,
  layers = c(128, 64, 32),
  model_features = c("read_index", "strand", "trinucleotide_ctx", "first_in_pair", "umi_count"),
  lr = 0.01,
  batch_size = 10000,
  epochs = 100,
  model_file_path = NULL,
  ...)
```

### Feature Selection

The DREAMS model supports a variety of features categorized into
numeric, categorical, and embedded types:

#### Numeric Features

- `read_index`, `fragment_size`, `local_GC`, `umi_count`, `umi_errors`,
  `local_complexity_1`, `local_complexity_2`, `n_other_errors`,
  `prior_error`, `seq_length`

#### Categorical Features

- `ref`, `strand`, `first_in_pair`, `ctx_minus1`, `ctx_plus1`, `chr`,
  `genomic_pos`

#### Embedded Feature

- `trinucleotide_ctx`

Ensure that the dataset used aligns with the selected features and
adjust the parameters such as `layers`, `lr`, `batch_size`, and `epochs`
as needed.

``` r
# Call variants using DREAMS-vc

variant_calls = dreams_vc(
  mutations_df = mutations_df,
  bam_file_path = "/path/test_bam_file",
  reference_path = "/path/hg38.fa",
  beta = beta,
  model = model,
  ...)

# Call cancer using DREAMS-cc

cancer_calls = dreams_cc(
  mutations_df = mutations_df,
  bam_file_path = "/path/test_bam_file",
  reference_path = "/path/hg38.fa",
  beta = beta,
  model = model,
  ...)
```

## Parallel Variant Calling with dreams_vc_parallel

In addition to `dreams_vc`, the DREAMS package offers
`dreams_vc_parallel` for parallel variant calling. This function is
particularly useful when working with large datasets and is designed to
leverage computational resources such as multi-core CPUs for parallel
processing.

### When to use dreams_vc_parallel

- **Large Datasets**: Efficiently handles larger datasets by
  distributing the workload across multiple cores.
- **Parallel Processing Capability**: Ideal when the computational
  environment supports parallel processing (e.g., multi-core systems).

### Example Usage

``` r
# Parallel variant calling
parallel_variant_calls = dreams_vc_parallel(
  mutations_df = mutations_df,
  bam_file_path = "/path/test_bam_file",
  reference_path = "/path/hg38.fa",
  beta = beta,
  model = model,
  ...)
```

### Saving and Loading Models

You can save your trained models for later use and load them as needed.
To save a trained model, use the `save_model_hdf5` function from the
`keras` package. Specify the file path where you want to save the model.
Similarly, when training a model using the `train_dreams_model`
function, you can directly specify a file path where to save the model
using the `model_file_path` argument. This allows for automatic saving
of the model upon training completion.

As default `model_file_path = NULL`, and the model won’t be saved
automatically. You can then manually save the model using
`save_model_hdf5`.

To load a previously saved model, use the `load_model_hdf5` function.

``` r
library(keras)

# Save the model
save_model_hdf5(model, filepath = "path/to/your_model.h5")

# Load the model
loaded_model <- load_model_hdf5(filepath = "path/to/your_model.h5")
```

## About DREAMS

For technical details describing how DREAMS works please see our
[article](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-02920-1).
