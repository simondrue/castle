test_that("calculations", {
  log_p1 = -1
  log_p2 = -2

  expect_equal(sum_log_p(log_p1), log(exp(log_p1)))
  expect_equal(sum_log_p(c(log_p1, log_p2)), log(exp(log_p1) + exp(log_p2)))
  expect_equal(length(sum_log_p(c(log_p1, log_p2))), 1)
})

