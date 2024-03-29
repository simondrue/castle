test_that("sum_log_p", {
  log_p1 <- -1
  log_p2 <- -2

  # Calculations
  expect_equal(sum_log_p(log_p1), log(exp(log_p1)))
  expect_equal(sum_log_p(c(log_p1, log_p2)), log(exp(log_p1) + exp(log_p2)))

  # Avoid underflow (-Inf)
  expect_equal(sum_log_p(c(-1000,-1000)), -999.30685)

  # Output should be one number
  expect_equal(length(sum_log_p(c(log_p1, log_p2))), 1)
})

test_that("check_input_samples", {
  # Empty samples
  empty_sample <- data.frame(
    WildtypeOnlyDroplets = c(0, 1, 0),
    MutantOnlyDroplets = c(0, 1, 0),
    DoubleNegativeDroplets = c(0, 1, 0),
    DoublePositiveDroplets = c(0, 1, 0)
  )

  expect_error(
    check_input_samples(
      samples = empty_sample
    ),
    regexp = "1,3.*total droplet count is 0"
  )

  # Incomplete data
  incomplete_sample <- data.frame(
    WildtypeOnlyDroplets = 1,
    MutantOnlyDroplets = 1,
    DoubleNegativeDroplets = 1
  )

  expect_error(
    check_input_samples(
      samples = incomplete_sample
    ),
    regexp = "are missing from test_samples"
  )

  # No wild type negative
  no_WT_sample <- data.frame(
    WildtypeOnlyDroplets = 1,
    MutantOnlyDroplets = 0,
    DoubleNegativeDroplets = 0,
    DoublePositiveDroplets = 1
  )

  expect_error(
    check_input_samples(
      samples = no_WT_sample
    ),
    regexp = "have no WT-negative droplets"
  )
})
