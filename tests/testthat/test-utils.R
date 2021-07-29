test_that("sum_log_p", {
  log_p1 <- -1
  log_p2 <- -2

  expect_equal(sum_log_p(log_p1), log(exp(log_p1)))
  expect_equal(sum_log_p(c(log_p1, log_p2)), log(exp(log_p1) + exp(log_p2)))
  expect_equal(length(sum_log_p(c(log_p1, log_p2))), 1)
})

test_that("check_samples", {
  # Empty sample
  empty_sample <- data.frame(
    N_WT_only = c(0, 1, 0),
    N_M_only = c(0, 1, 0),
    N_d_neg = c(0, 1, 0),
    N_d_pos = c(0, 1, 0)
  )
  expect_error(
    check_samples(
      samples = empty_sample
    ),
    regexp = "total droplet count is 0"
  )

  # Incomplete data
  incomplete_sample <- data.frame(
    N_WT_only = 1,
    N_M_only = 1,
    N_d_neg = 1
  )

  expect_error(
    check_samples(
      samples = incomplete_sample
    ),
    regexp = "are missing from test_samples"
  )
})
