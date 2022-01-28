test_that(
  "training model - snapshot",
  {
    background_samples <- data.frame(
      WildtypeOnlyDroplets = c(10, 20, 30, 0),
      MutantOnlyDroplets = c(1, 2, 3, 0),
      DoubleNegativeDroplets = c(30, 30, 30, 30),
      DoublePositiveDroplets = c(1, 1, 1, 0)
    )

    abc_grid_resolution <- 7

    # Check whats in a model
    expect_snapshot(
      train_integrated_ddpcr_model(
        background_samples = background_samples,
        abc_grid_resolution = abc_grid_resolution
      )
    )
  }
)

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
  expect_warning(train_integrated_ddpcr_model(background_samples = background_samples), NA)

  # Check channel 1
  background_samples$Ch2TargetType <- c("Unknown")
  for (wrong_target_type in c("Positive Control", "NTC", "Blank")) {
    background_samples$Ch1TargetType <- c("Unknown", "Unknown", wrong_target_type)
    expect_warning(
      train_integrated_ddpcr_model(background_samples = background_samples),
      paste0("Ch1.*", wrong_target_type, ".*remove")
    )
  }

  # Check channel 2
  background_samples$Ch1TargetType <- c("Unknown")
  for (wrong_target_type in c("Positive Control", "NTC", "Blank")) {
    background_samples$Ch2TargetType <- c("Unknown", "Unknown", wrong_target_type)
    expect_warning(
      train_integrated_ddpcr_model(background_samples = background_samples),
      paste0("Ch2.*", wrong_target_type, ".*remove")
    )
  }
})

test_that(
  "simulation - training - test",
  {
    set.seed(42)

    # Simulate training data
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
    background_samples <- rbind(
      simulate_droplet_counts(
        l = l_low, r = 0, a = a, b = b, c = c, n_drops = n_drops, n_samples = n_samples
      ),
      simulate_droplet_counts(
        l = l_mid, r = 0, a = a, b = b, c = c, n_drops = n_drops, n_samples = n_samples
      ),
      simulate_droplet_counts(
        l = l_high, r = 0, a = a, b = b, c = c, n_drops = n_drops, n_samples = n_samples
      )
    )

    # Train model
    trained_integrated_model <- train_integrated_ddpcr_model(
      background_samples = background_samples,
      abc_grid_resolution = 10
    )

    # Test on "known" samples
    # Positive sample
    test_sample_positive <-
      simulate_droplet_counts(
        l = l_mid, r = 0.1, a = a, b = b, c = c, n_drops = n_drops
      )

    positive_test_res <- test_tumor_sample_integrated(
      test_samples = test_sample_positive,
      integrated_model = trained_integrated_model
    )

    expect_true(positive_test_res$mutation_detected)

    # Negative sample
    test_sample_negative <-
      simulate_droplet_counts(
        l = l_mid, r = 0, a = a, b = b, c = c, n_drops = n_drops
      )

    negative_test_res <- test_tumor_sample_integrated(
      test_samples = test_sample_negative,
      integrated_model = trained_integrated_model
    )

    expect_false(negative_test_res$mutation_detected)

    # High mutational content
    test_sample_high_ctdna <-
      simulate_droplet_counts(
        l = l_mid, r = 2, a = a, b = b, c = c, n_drops = n_drops
      )

    expect_error(
      test_tumor_sample_integrated(
        test_samples = test_sample_high_ctdna,
        integrated_model = trained_integrated_model
      ),
      NA
    )

    # Multiple samples
    multiple_test_res <- test_tumor_sample_integrated(
      test_samples = dplyr::bind_rows(test_sample_positive, test_sample_negative),
      integrated_model = trained_integrated_model
    )

    # Dimensions
    expect_equal(nrow(multiple_test_res), 2)

    # Consistent with single tests
    expect_true(all(multiple_test_res[1, ] == positive_test_res))
    expect_true(all(multiple_test_res[2, ] == negative_test_res))
  }
)
