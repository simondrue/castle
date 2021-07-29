test_that("training model - snapshot", {
  training_samples <- data.frame(
    N_WT_only = c(10, 20, 30, 0),
    N_M_only = c(1, 2, 3, 0),
    N_d_neg = c(30, 30, 30, 30),
    N_d_pos = c(1, 1, 1, 0)
  )

  abc_grid_resolution <- 7

  model <- train_integrated_ddpcr_model(
    training_samples = training_samples,
    abc_grid_resolution = abc_grid_resolution
  )

  # Test dimensions of output
  expect_equal(nrow(model$abc_grid), abc_grid_resolution^3)
  expect_equal(ncol(model$abc_grid), 3)
  expect_equal(length(model$abc_grid_train_log_lik), abc_grid_resolution^3)

  # Check whats in a model
  expect_true(
    all(
      c("a_est", "b_est", "c_est", "l_est_vec", "abc_grid", "abc_grid_train_log_lik") %in% names(model)
    )
  )

  # Snapshot test of output
  expect_equal(signif(model$a_est,3), 5.13e-07)
  expect_equal(signif(model$b_est, 3), 8.09e-07)
  expect_equal(signif(model$c_est, 3), 0.113)
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

  training_samples <- rbind(l_low_df, l_mid_df, l_high_df)

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
    training_samples = training_samples,
    abc_grid_resolution = 10
  )

  # Test on "known" samples
  # Positive
  positive_test_res <- test_tumor_sample_integrated(
    test_samples = test_sample_positive,
    integrated_model = trained_integrated_model
  )

  expect_true(positive_test_res$is_tumor_positive)


  # Negative
  negative_test_res <- test_tumor_sample_integrated(
    test_samples = test_sample_negative,
    integrated_model = trained_integrated_model
  )

  expect_false(negative_test_res$is_tumor_positive)

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
})
