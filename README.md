
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Welcome to the **CASTLE**

<!-- badges: start -->
<!-- badges: end -->

Using the right analysis tool for your ddPCR analysis results when
detecting low frequency variants can mean (k)night and day for the
performance. In this package functions for training a model and testing
samples for low frequency variants are provided, and includes
functionality to handle output files generated by QuantaSoft from
Bio-Rad Laboratories Inc.  This package is a user friendly
implementation of the CASTLE algorithm, as it is described in \[1\].

\[1\]: \[article\].

## Installation

You can install the R -package from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("simondrue/castle")
```

## Quickstart guide

In this section we will go through how to load a dataset from
QuantaSoft, train a model for the CASTLE algorithm and test some
samples. The first thing we need to do is load the packages we need for
the analysis:

``` r
library(castle)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
```

### Loading QuantaSoft data

We start by loading a QuantaSoft dataset in -format using the included
-function. Here we will use the example data included in the package,
where we have both negative samples of different concentrations
(training samples) and some test samples available for a KRAS G12D
assay. The training sampled are loaded like this:

``` r
path_to_my_data = system.file("extdata/example_data/", "training_data_KRAS_G12D.csv", package = "castle")
training_samples = import_QS_files(path_to_my_data)
```

The first 5 samples/wells in looks like this:

| FileName                       | Sample     | Well | Ch1TargetType | Ch2TargetType | Target            | DoubleNegativeDroplets | WildtypeOnlyDroplets | MutantOnlyDroplets | DoublePositiveDroplets | TotalDroplets | NumberOfMergedWells | MergedWells |
|:-------------------------------|:-----------|:-----|:--------------|:--------------|:------------------|-----------------------:|---------------------:|-------------------:|-----------------------:|--------------:|--------------------:|:------------|
| training\_data\_KRAS\_G12D.csv | NC20\_150  | A01  | Unknown       | Unknown       | KRAS\_G12D\_(FAM) |                   2415 |                11005 |                  0 |                      0 |         13420 |                   1 | NA          |
| training\_data\_KRAS\_G12D.csv | NC20\_24   | F05  | Unknown       | Unknown       | KRAS\_G12D\_(FAM) |                  12458 |                 4040 |                  0 |                      1 |         16499 |                   1 | NA          |
| training\_data\_KRAS\_G12D.csv | NC20\_3.84 | D10  | Unknown       | Unknown       | KRAS\_G12D\_(FAM) |                  17966 |                  463 |                  0 |                      0 |         18429 |                   1 | NA          |
| training\_data\_KRAS\_G12D.csv | NC20\_60   | C03  | Unknown       | Unknown       | KRAS\_G12D\_(FAM) |                   6762 |                 7283 |                  0 |                      1 |         14046 |                   1 | NA          |
| training\_data\_KRAS\_G12D.csv | NC20\_9.6  | A08  | Unknown       | Unknown       | KRAS\_G12D\_(FAM) |                  14554 |                 1547 |                  0 |                      1 |         16102 |                   1 | NA          |

Note that besides some metadata (filename, target etc.) for each sample,
the datasat contains the columns “WildtypeOnlyDroplets”,
“MutantOnlyDroplets”, “DoubleNegativeDroplets” and
“DoublePositiveDroplets”. These are the counts of the four different
kinds of droplets for each sample based on the two signals in Ch1 and
Ch2. The following model training and statistical tests are based solely
on these counts.

### Training a model for CASTLE

To train the model we are only interested in the negative plasma
samples. Ideally these should cover the range of different
concentrations of DNA that the trained model is expected to encounter in
future test samples. Before training we start by removing the empty
samples (“NTC”) and positive controls from the training dataset we
loaded above:

``` r
# Remove control samples ("NTC" and "Positive Control")
clean_training_samples = 
  training_samples %>%
  filter(
    !(Ch1TargetType %in% c("NTC", "Positive Control")),
    !(Ch2TargetType %in% c("NTC", "Positive Control"))
  )
```

The cleaned data is then used for model training using the function :

``` r
trained_model <- train_integrated_ddpcr_model(
  training_samples = clean_training_samples
)
```

The object is simply a list of parameters used in the statistical model
in CASTLE.

### Testing some samples with CASTLE

For testing we will use a variety of samples from different 2 patients
who have undergone curative surgery of colorectal cancer. We start by
loading the data from these patients. To load the data we will again use
the function . But this time data from each patient is distributed
across multiple wells. These can be merged together using the flag .
Furthermore samples across multiple files can also be merged using the
flag if necessary.

``` r
# Path to example data in the package 
path_to_my_test_samples = system.file(
  "extdata/example_data/test_samples", 
  c(
    "patient_1_plasma_KRAS_G12D.csv",
    "patient_2_plasma_KRAS_G12D.csv"
  ), 
  package = "castle")

# Import test examples 
test_samples = 
  import_QS_files(
    path_to_my_test_samples, 
    merge_wells = TRUE)
```

To test the samples we use the function and use the model we trained
above:

``` r
# Make calls using CASTLE
test_res = test_tumor_sample_integrated(
  test_samples = test_samples, 
  integrated_model = trained_model
)
```

Well start by looking at some of the results for the tumor and PBL
(white blood cells) sample from each patient:

``` r
test_res %>% 
  dplyr::filter(Sample %in% c("Tumor","PBL")) %>%
  dplyr::select(FileName, Sample, p_val, mutation_detected) %>%
  knitr::kable()
```

| FileName                           | Sample |    p\_val | mutation\_detected |
|:-----------------------------------|:-------|----------:|:-------------------|
| patient\_1\_plasma\_KRAS\_G12D.csv | PBL    | 1.0000000 | FALSE              |
| patient\_1\_plasma\_KRAS\_G12D.csv | Tumor  | 0.0000000 | TRUE               |
| patient\_2\_plasma\_KRAS\_G12D.csv | PBL    | 0.4550231 | FALSE              |
| patient\_2\_plasma\_KRAS\_G12D.csv | Tumor  | 0.0000000 | TRUE               |

From this we clearly see that the mutation is present in the tumor of
both of the patients and not in PBL. We can therefore use it as an
indicator of the specific cancer, and thus to check if the patient gets
a relapse after the curative surgery. As a positive test we start by
looking at the “Plasma\_A” sample of both patients, which is a blood
sample taken before the surgery:

``` r
test_res %>% 
  dplyr::filter(Sample %in% c("Plasma_A")) %>%
  dplyr::select(FileName, Sample, p_val, mutation_detected) %>%
  knitr::kable()
```

| FileName                           | Sample    | p\_val | mutation\_detected |
|:-----------------------------------|:----------|-------:|:-------------------|
| patient\_1\_plasma\_KRAS\_G12D.csv | Plasma\_A |      0 | TRUE               |
| patient\_2\_plasma\_KRAS\_G12D.csv | Plasma\_A |      0 | TRUE               |

We see that these samples indicate that the mutation, and thus the
tumor, is present in both patient, as expected. From the result we can
also see the following estimates of properties of the ddPCR analysis and
and the ctDNA in the blood:

``` r
test_res %>% 
  dplyr::filter(Sample %in% c("Plasma_A")) %>%
  dplyr::select(FileName, Sample, 
                mutant_molecules_per_droplet, 
                wildtype_molecules_per_droplet,
                total_mutant_molecules,
                mutant_molecules_per_droplet_CI_lower, mutant_molecules_per_droplet_CI_upper,
                wildtype_molecules_per_droplet_CI_lower, wildtype_molecules_per_droplet_CI_upper) %>%
  knitr::kable(caption = "ddPCR statstics")
```

| FileName                           | Sample    | mutant\_molecules\_per\_droplet | wildtype\_molecules\_per\_droplet | total\_mutant\_molecules | mutant\_molecules\_per\_droplet\_CI\_lower | mutant\_molecules\_per\_droplet\_CI\_upper | wildtype\_molecules\_per\_droplet\_CI\_lower | wildtype\_molecules\_per\_droplet\_CI\_upper |
|:-----------------------------------|:----------|--------------------------------:|----------------------------------:|-------------------------:|-------------------------------------------:|-------------------------------------------:|---------------------------------------------:|---------------------------------------------:|
| patient\_1\_plasma\_KRAS\_G12D.csv | Plasma\_A |                       0.0002420 |                         0.2235039 |                 21.32176 |                                  0.0001289 |                                  0.0004051 |                                    0.2191840 |                                    0.2278834 |
| patient\_2\_plasma\_KRAS\_G12D.csv | Plasma\_A |                       0.0038629 |                         0.0959883 |                376.09432 |                                  0.0033714 |                                  0.0043999 |                                    0.0933878 |                                    0.0986397 |

ddPCR statstics

``` r
test_res %>% 
  dplyr::filter(Sample %in% c("Plasma_A")) %>%
  dplyr::select(FileName, Sample, 
                allele_frequency, 
                allele_frequency_CI_lower, allele_frequency_CI_upper) %>%
  knitr::kable(caption = "Blood (ctDNA) statstics")
```

| FileName                           | Sample    | allele\_frequency | allele\_frequency\_CI\_lower | allele\_frequency\_CI\_upper |
|:-----------------------------------|:----------|------------------:|-----------------------------:|-----------------------------:|
| patient\_1\_plasma\_KRAS\_G12D.csv | Plasma\_A |         0.0010814 |                    0.0005764 |                    0.0018094 |
| patient\_2\_plasma\_KRAS\_G12D.csv | Plasma\_A |         0.0386864 |                    0.0339311 |                    0.0438292 |

Blood (ctDNA) statstics

We can then look what happens after the surgery. To do this we look at
the blood samples marked as “Plasma\_B”:

``` r
test_res %>% 
  dplyr::filter(Sample %in% c("Plasma_B")) %>%
  dplyr::select(FileName, Sample, p_val, mutation_detected) %>%
  knitr::kable()
```

| FileName                           | Sample    |    p\_val | mutation\_detected |
|:-----------------------------------|:----------|----------:|:-------------------|
| patient\_1\_plasma\_KRAS\_G12D.csv | Plasma\_B | 0.3061083 | FALSE              |
| patient\_2\_plasma\_KRAS\_G12D.csv | Plasma\_B | 0.2424636 | FALSE              |

The full list of information available in the analysis results is:

    #>  [1] "FileName"                               
    #>  [2] "Sample"                                 
    #>  [3] "Well"                                   
    #>  [4] "Ch1TargetType"                          
    #>  [5] "Ch2TargetType"                          
    #>  [6] "Target"                                 
    #>  [7] "DoubleNegativeDroplets"                 
    #>  [8] "WildtypeOnlyDroplets"                   
    #>  [9] "MutantOnlyDroplets"                     
    #> [10] "DoublePositiveDroplets"                 
    #> [11] "TotalDroplets"                          
    #> [12] "NumberOfMergedWells"                    
    #> [13] "MergedWells"                            
    #> [14] "mutant_molecules_per_droplet"           
    #> [15] "wildtype_molecules_per_droplet"         
    #> [16] "p_val"                                  
    #> [17] "test_statistic"                         
    #> [18] "mutation_detected"                      
    #> [19] "allele_frequency"                       
    #> [20] "total_mutant_molecules"                 
    #> [21] "total_wildtype_molecules"               
    #> [22] "mutant_molecules_per_droplet_CI_lower"  
    #> [23] "mutant_molecules_per_droplet_CI_upper"  
    #> [24] "allele_frequency_CI_lower"              
    #> [25] "allele_frequency_CI_upper"              
    #> [26] "total_mutant_molecules_CI_lower"        
    #> [27] "total_mutant_molecules_CI_upper"        
    #> [28] "wildtype_molecules_per_droplet_CI_lower"
    #> [29] "wildtype_molecules_per_droplet_CI_upper"

We now see that the test indicates that the mutation is no longer
present in any of the patients, which is consistent with the tumor being
removed by surgery. By testing more blood samples over time it is
possible to test if one of the patients gets a relapse.

## Note

Alternatively the functions and are used for training and testing
respectively. These are based on a faster (and simpler) version of the
CASTLE algorithm, but have very similar performance. See \[1\] for more
information.

## Simulating data

In the package the function for simulating droplet counts from the
statistical model that CASTLE is based on is also provided. Here we will
simulate a positive and negative sample based on the parameters learned
above, and test them using :

``` r
# DNA concentration (molecules per droplet)
l = 0.5 # WT
r_positive = 0.01 # M
r_negative = 0 # M

# Number of droplets
n_drops = 14000

# Positive sample
test_sample_positive <-
  simulate_droplet_counts(
    l = l, r = r_positive, 
    a = trained_model$a_est, b = trained_model$b_est, c = trained_model$c_est, 
    n_drops = n_drops
  )

# Negative sample
test_sample_negative <-
  simulate_droplet_counts(
    l = l, r = r_negative, 
    a = trained_model$a_est, b = trained_model$b_est, c = trained_model$c_est, 
    n_drops = n_drops
  )


simulated_samples = bind_rows(
  data.frame(Sample = "Positive", test_sample_positive),
  data.frame(Sample = "Negative", test_sample_negative)
)


# Test samples
sim_test_res = test_tumor_sample_simple(
  test_samples = simulated_samples,
  simple_model = trained_model
)
```

| Sample   | p\_val | mutation\_detected | total\_mutant\_molecules | total\_mutant\_molecules\_CI\_lower | total\_mutant\_molecules\_CI\_upper |
|:---------|-------:|:-------------------|-------------------------:|------------------------------------:|------------------------------------:|
| Positive |      0 | TRUE               |                 144.3455 |                            115.4577 |                          177.677227 |
| Negative |      1 | FALSE              |                   0.0000 |                              0.0000 |                            5.096182 |
