find_abd_confidence_intervals <- function(N_WT_only_vec, N_M_only_vec, N_d_neg_vec, N_d_pos_vec,
                                          l_est_vec, a_est, b_est, c_est,
                                          alpha = 0.05) {
  # Helper functions
  ll_ratio_simple <- function(par_0, par) {
    ll_ratio <- ll_ratio(
      par_0 = par_0, par = par,
      l_est_vec = l_est_vec,
      r_est_vec = rep(0, length(l_est_vec)),
      a_est = a_est,
      b_est = b_est,
      c_est = c_est,
      N_WT_only_vec = N_WT_only_vec,
      N_M_only_vec = N_M_only_vec,
      N_d_neg_vec = N_d_neg_vec,
      N_d_pos_vec = N_d_pos_vec
    )

    return(ll_ratio)
  }

  # Find upper and lower bound of a
  if (ll_ratio_simple(TOL_0, par = "a") > stats::qchisq(1 - alpha, 1)) {
    uni_res <- stats::uniroot(function(x) ll_ratio_simple(x, par = "a") - stats::qchisq(1 - alpha, 1),
      interval = c(TOL_0, a_est),
      tol = TOL_0
    )
    a_CI_lower <- uni_res$root
  } else {
    a_CI_lower <- TOL_0
  }

  a_CI_upper <- stats::uniroot(function(x) ll_ratio_simple(x, par = "a") - stats::qchisq(1 - alpha, 1),
    interval = c(a_est, 5),
    tol = TOL_0
  )$root

  # Find upper and lower bound of b
  if (ll_ratio_simple(TOL_0, par = "b") > stats::qchisq(1 - alpha, 1)) {
    uni_res <- stats::uniroot(function(x) ll_ratio_simple(x, par = "b") - stats::qchisq(1 - alpha, 1),
      interval = c(TOL_0, b_est),
      tol = TOL_0
    )
    b_CI_lower <- uni_res$root
  } else {
    b_CI_lower <- TOL_0
  }

  b_CI_upper <- stats::uniroot(function(x) ll_ratio_simple(x, par = "b") - stats::qchisq(1 - alpha, 1),
    interval = c(b_est, 5),
    tol = TOL_0
  )$root

  # Find upper and lower bound of c
  if (ll_ratio_simple(TOL_0, par = "c") > stats::qchisq(1 - alpha, 1)) {
    uni_res <- stats::uniroot(function(x) ll_ratio_simple(x, par = "c") - stats::qchisq(1 - alpha, 1),
      interval = c(TOL_0, c_est),
      tol = TOL_0
    )
    c_CI_lower <- uni_res$root
  } else {
    c_CI_lower <- TOL_0
  }

  c_CI_upper <- stats::uniroot(function(x) ll_ratio_simple(x, par = "c") - stats::qchisq(1 - alpha, 1),
    interval = c(c_est, 5),
    tol = TOL_0
  )$root

  # Return result
  res <- list(
    a_CI_lower = a_CI_lower,
    a_CI_upper = a_CI_upper,
    b_CI_lower = b_CI_lower,
    b_CI_upper = b_CI_upper,
    c_CI_lower = c_CI_lower,
    c_CI_upper = c_CI_upper
  )
  return(res)
}

find_lr_confidence_intervals <- function(N_WT_only, N_M_only, N_d_neg, N_d_pos,
                                         l_est, r_est,
                                         a_est, b_est, c_est,
                                         alpha = 0.05) {
  # Helper functions
  ll_ratio_simple <- function(par_0, par) {
    ll_ratio <- ll_ratio(
      par_0 = par_0,
      par = par,
      l_est_vec = l_est,
      r_est_vec = r_est,
      a_est = a_est,
      b_est = b_est,
      c_est = c_est,
      N_WT_only_vec = N_WT_only,
      N_M_only_vec = N_M_only,
      N_d_neg_vec = N_d_neg,
      N_d_pos_vec = N_d_pos
    )
    return(ll_ratio)
  }

  # Find bounds on l
  # Model as binomial with p = P(WT=0)=exp(-l)
  WT_neg <- N_M_only + N_d_neg
  n_drops <- N_WT_only + N_M_only + N_d_neg + N_d_pos

  binom_res <- stats::binom.test(x = WT_neg, n = n_drops, conf.level = 1 - alpha)

  l_CI_lower <- -log(binom_res$conf.int[[2]])
  l_CI_upper <- -log(binom_res$conf.int[[1]])

  # Find upper and lower bound of r
  if (ll_ratio_simple(par_0 = 0, par = "r") > stats::qchisq(1 - alpha, 1)) {
    uni_res <- stats::uniroot(function(x) ll_ratio_simple(x, par = "r") - stats::qchisq(1 - alpha, 1),
      interval = c(0, r_est),
      tol = TOL_0
    )
    r_CI_lower <- uni_res$root
  } else {
    r_CI_lower <- 0
  }

  r_CI_upper <- stats::uniroot(function(x) ll_ratio_simple(x, par = "r") - stats::qchisq(1 - alpha, 1),
    interval = c(r_est, 1),
    tol = TOL_0,
    extendInt = "upX"
  )$root

  # P val for r=0
  # r MLE under the hypothesis: H0: r=0
  # Note: If r_est < 0 --> r_MLE_H0 = r_est, since r_est is estimate over all values
  r_MLE_H0 <- 0
  r_LLR_test_stat <- ll_ratio_simple(r_MLE_H0, par = "r")
  r_less_than_0_pval <- 1 - stats::pchisq(r_LLR_test_stat, 1)

  # Return result
  res <- data.frame(
    l_CI_lower = l_CI_lower,
    l_CI_upper = l_CI_upper,
    r_CI_lower = r_CI_lower,
    r_CI_upper = r_CI_upper,
    r_less_than_0_pval = r_less_than_0_pval,
    r_LLR_test_stat = r_LLR_test_stat
  )
  return(res)
}
