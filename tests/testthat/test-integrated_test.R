test_that("training model - snapshot", {
  training_samples <- data.frame(
    N_WT_only = c(10, 20, 30),
    N_M_only = c(1, 2, 3),
    N_d_neg = c(30, 30, 30),
    N_d_pos = c(1, 1, 1)
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
  expect_equal(model$a_est, 1.005659e-06, tolerance = 1e-3)
  expect_equal(model$b_est, 9.966160e-07, tolerance = 1e-3)
  expect_equal(model$c_est, 0.1132368000, tolerance = 1e-3)
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

  # Test negative sample
  test_sample_negative <-
    simulate_droplet_counts(
      l = l_mid, r = 0, a = a, b = b, c = c, n_drops = n_drops
    )

  # Train model
  trained_model <- train_integrated_ddpcr_model(
    training_samples = training_samples,
    abc_grid_resolution = 10
  )

  # Test on "known" samples
  # Positive
  positive_test_res <- test_tumor_sample_integrated(
    test_samples = test_sample_positive,
    model = trained_model
  )

  expect_true(positive_test_res$is_tumor_positive)


  # Negative
  negative_test_res <- test_tumor_sample_integrated(
    test_samples = test_sample_negative,
    model = trained_model
  )

  expect_false(negative_test_res$is_tumor_positive)
})
