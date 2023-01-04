
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dreams

<!-- badges: start -->

[![R-CMD-check](https://github.com/JakobSkouPedersenLab/dreams/workflows/R-CMD-check/badge.svg)](https://github.com/JakobSkouPedersenLab/dreams/actions)
[![DOI](https://zenodo.org/badge/455089263.svg)](https://zenodo.org/badge/latestdoi/455089263)
[![Codecov test
coverage](https://codecov.io/gh/JakobSkouPedersenLab/dreams/branch/main/graph/badge.svg)](https://app.codecov.io/gh/JakobSkouPedersenLab/dreams?branch=main)
[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

<!-- badges: end -->

## Installation

You can install the development version of dreams from
[GitHub](https://github.com/JakobSkouPedersenLab/dreams) with:

``` r
# install.packages("devtools")
devtools::install_github("JakobSkouPedersenLab/dreams")
```

## Basic functions

``` r
library(dreams)

# For training data DREAMS requires one or more bam-files
# and a reference genome

training_data = get_training_data(
  bam_paths = "/path/bam_file",
  reference_path = "/path/hg38.fa",
  ...)

# The model can be trained using a neural network
# - and requires basic settings for keras

model = train_dreams_model(
  training_data,
  layers = c(16,8),
  model_features = c("trinucleotide_ctx", "strand", "read_index"),
  lr = 0.01,
  batch_size = 10000,
  epochs = 100,
  ...)
  
# Call variants using DREAMS-vc

variant_calls = dreams_vc(
  bam_path = "/path/test_bam_file",
  positions = "positions_file",
  model = model,
  ...
)

# Call cancer using DREAMS-cc

cancer_calls = dreams_cc(
  bam_path = "/path/test_bam_file",
  positions = "positions_file",
  model = model,
  ...)
```
