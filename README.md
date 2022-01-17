
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Welcome to the CASTLE

<!-- badges: start -->
<!-- badges: end -->

This package is a user friendly implementation of the **CASTLE**
(**C**oncentration and **A**ssay **S**pecific **T**umor **L**oad
**E**stimator) algorithm, as it is described in the accompanying
[scientific paper](https://doi.org/10.1093/clinchem/hvab274). Using the
right analysis tool for your ddPCR analysis results when detecting low
frequency variants can mean night and day for the performance of your
pipeline. In this package functions for training a model and testing
samples for low frequency variants are provided, as well as
functionality to handle output .csv-files generated by QuantaSoft from
Bio-Rad Laboratories Inc.

## Installation

You can install the R -package from
[GitHub](https://github.com/simondrue/castle/) with:

``` r
# install.packages("devtools")
devtools::install_github("simondrue/castle")
```

## Usage

Basic usage of `castle` consists of the following classic workflow:

``` r
# Load package
library(castle)

# Load training samples
training_data_path <- c("path/to/training/csv_files/")
training_samples <- import_QS_files(training_data_path)

# Train model
trained_model <- train_integrated_ddpcr_model(
  background_samples = training_samples
)

# Test some samples
test_samples <- c("path/to/test/csv_files/")
test_res <- test_tumor_sample_integrated(
  test_samples = test_samples,
  integrated_model = trained_model
)
```

To see a full example of how to do this with real world data, please see
`vignette("castle")` which goes through a full analysis of a real world
dataset.
