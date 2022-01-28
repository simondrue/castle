test_that("training model - snapshot", {
  background_samples <- data.frame(
    WildtypeOnlyDroplets = c(10, 20, 30, 0),
    MutantOnlyDroplets = c(1, 2, 3, 0),
    DoubleNegativeDroplets = c(30, 30, 30, 30),
    DoublePositiveDroplets = c(1, 1, 1, 0)
  )

  model <- train_simple_ddpcr_model(
    background_samples = background_samples
  )

  # Snapshot test of output
  expect_snapshot(train_simple_ddpcr_model(
    background_samples = background_samples
  ))
})

test_that("training model - Catch wrong sample types", {
  background_samples <- data.frame(
    WildtypeOnlyDroplets = c(1, 1, 1),
    MutantOnlyDroplets = c(1, 2, 3),
    DoubleNegativeDroplets = c(1, 2, 3),
    DoublePositiveDroplets = c(1, 2, 3)
  )

  # No error
  background_samples$Ch1TargetType <- c("Unknown")
  background_samples$Ch2TargetType <- c("Unknown")
  expect_warning(train_simple_ddpcr_model(background_samples = background_samples), NA)

  # Check channel 1
  background_samples$Ch2TargetType <- c("Unknown")
  for (wrong_target_type in c("Positive Control", "NTC", "Blank")) {
    background_samples$Ch1TargetType <- c("Unknown", "Unknown", wrong_target_type)
    expect_warning(
      train_simple_ddpcr_model(background_samples = background_samples),
      paste0("Ch1.*", wrong_target_type, ".*remove")
    )
  }

  # Check channel 2
  background_samples$Ch1TargetType <- c("Unknown")
  for (wrong_target_type in c("Positive Control", "NTC", "Blank")) {
    background_samples$Ch2TargetType <- c("Unknown", "Unknown", wrong_target_type)
    expect_warning(
      train_simple_ddpcr_model(background_samples = background_samples),
      paste0("Ch2.*", wrong_target_type, ".*remove")
    )
  }
})

test_that("simulation - training - test", {
  set.seed(42)

  # Simulate data
  # Parameters
  l_low <- 0.1
  l_mid <- 0.5
  l_high <- 1

  a <- 0.01
  b <- 0.01
  c <- 0.01
  n_drops <- 14000
  n_samples <- 30

  # Training data
  l_low_df <- simulate_droplet_counts(
    l = l_low, r = 0, a = a, b = b, c = c, n_drops = n_drops, n_samples = n_samples
  )
  l_mid_df <- simulate_droplet_counts(
    l = l_mid, r = 0, a = a, b = b, c = c, n_drops = n_drops, n_samples = n_samples
  )
  l_high_df <- simulate_droplet_counts(
    l = l_high, r = 0, a = a, b = b, c = c, n_drops = n_drops, n_samples = n_samples
  )

  background_samples <- rbind(l_low_df, l_mid_df, l_high_df)

  # Positive sample
  test_sample_positive <-
    simulate_droplet_counts(
      l = l_mid, r = 0.1, a = a, b = b, c = c, n_drops = n_drops
    )

  # Negative sample
  test_sample_negative <-
    simulate_droplet_counts(
      l = l_mid, r = 0, a = a, b = b, c = c, n_drops = n_drops
    )

  # Train model
  trained_model <- train_simple_ddpcr_model(
    background_samples = background_samples
  )

  # Test on "known" samples
  # Positive
  positive_test_res <- test_tumor_sample_simple(
    test_samples = test_sample_positive,
    simple_model = trained_model,
    alpha = 0.05
  )

  expect_true(positive_test_res$mutation_detected)


  # Negative
  negative_test_res <- test_tumor_sample_simple(
    test_samples = test_sample_negative,
    simple_model = trained_model,
    alpha = 0.05
  )

  expect_false(negative_test_res$mutation_detected)

  # Multiple samples
  multiple_test_res <- test_tumor_sample_simple(
    test_samples = dplyr::bind_rows(test_sample_positive, test_sample_negative),
    simple_model = trained_model,
    alpha = 0.05
  )

  # Dimensions
  expect_equal(nrow(multiple_test_res), 2)

  # No CIs
  no_CIs_res <- test_tumor_sample_simple(
    test_samples = dplyr::bind_rows(test_sample_positive, test_sample_negative),
    simple_model = trained_model,
    include_wildtype_CI = FALSE,
    include_mutant_CI = FALSE,
    alpha = 0.05
  )

  expect_false(
    any(
      c("r_CI_lower", "r_CI_upper", "wildtype_molecules_per_droplet_CI_lower", "wildtype_molecules_per_droplet_CI_upper") %in% colnames(no_CIs_res)
    )
  )

  # Consistent with single tests
  expect_true(all(multiple_test_res[1, ] == positive_test_res))
  expect_true(all(multiple_test_res[2, ] == negative_test_res))
})
