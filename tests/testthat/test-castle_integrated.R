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

test_that(
  "simulation - training - test",
  {
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
    trained_integrated_model <- train_integrated_ddpcr_model(
      background_samples = background_samples,
      abc_grid_resolution = 10
    )

    # Test on "known" samples
    # Positive
    positive_test_res <- test_tumor_sample_integrated(
      test_samples = test_sample_positive,
      integrated_model = trained_integrated_model
    )

    expect_true(positive_test_res$mutation_detected)


    # Negative
    negative_test_res <- test_tumor_sample_integrated(
      test_samples = test_sample_negative,
      integrated_model = trained_integrated_model
    )

    expect_false(negative_test_res$mutation_detected)

    # Multiple samples
    multiple_test_res <- test_tumor_sample_integrated(
      test_samples = dplyr::bind_rows(test_sample_positive, test_sample_negative),
      integrated_model = trained_integrated_model
    )

    # Dimensions
    expect_equal(nrow(multiple_test_res), 2)
    expect_snapshot(multiple_test_res)


    # Consistent with single tests
    expect_true(all(multiple_test_res[1, ] == positive_test_res))
    expect_true(all(multiple_test_res[2, ] == negative_test_res))
  }
)
