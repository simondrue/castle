test_that("calculation example", {
  expect_equal(sum_log_p(c(-1,-1)), log(sum(exp(c(-1,-1)))))
})

