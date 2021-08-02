get_start_values_abc <- function(N_WT_only_vec, N_M_only_vec, N_d_neg_vec, N_d_pos_vec,
                                 l_est_vec) {
  # Filter out samples with extreme low concentration
  N_WT_only_vec <- N_WT_only_vec[l_est_vec > TOL_0]
  N_M_only_vec <- N_M_only_vec[l_est_vec > TOL_0]
  N_d_neg_vec <- N_d_neg_vec[l_est_vec > TOL_0]
  N_d_pos_vec <- N_d_pos_vec[l_est_vec > TOL_0]
  l_est_vec <- l_est_vec[l_est_vec > TOL_0]

  # Check if there is more than 1 sample to train on
  if (!(length(l_est_vec) > 1)) {
    stop("Starting values for (a,b,c) can only be found with >1 data points!")
  }

  # Estimate of a+c*l_est from method of moments
  y <- log(N_M_only_vec / N_d_neg_vec + 1)

  # Do linear regression of (l_est, y) to find a and c
  c_start <- sum((y - mean(y)) * (l_est_vec - mean(l_est_vec))) /
    sum((l_est_vec - mean(l_est_vec))^2)
  a_start <- mean(y) - c_start * mean(l_est_vec)

  # Make appropriate corrections if any parameter is less than MIN_START_PAR
  if (c_start < MIN_START_PAR & a_start >= MIN_START_PAR) {
    a_start <- mean(y)
    c_start <- MIN_START_PAR
  } else if (c_start >= MIN_START_PAR & a_start < MIN_START_PAR) {
    a_start <- MIN_START_PAR
    c_start <- max(sum(y * l_est_vec) / sum(l_est_vec^2), MIN_START_PAR)
  } else if (c_start < MIN_START_PAR & a_start < MIN_START_PAR) {
    a_start <- MIN_START_PAR
    c_start <- MIN_START_PAR
  }

  # From a and c we find b
  b_start_vec <- -log(
    log(N_WT_only_vec / (N_WT_only_vec + N_d_pos_vec) *
      exp(a_start + c_start * l_est_vec + l_est_vec) *
      (1 - exp(-l_est_vec)) + 1) /
      l_est_vec
  )

  b_start <- max(mean(b_start_vec), MIN_START_PAR)

  return(list(
    a_start = a_start,
    b_start = b_start,
    c_start = c_start
  ))
}

estimate_l <- function(N_WT_only, N_M_only, N_d_neg, N_d_pos) {

  # Estimate l
  l_est <- log((N_WT_only + N_d_pos) / (N_d_neg + N_M_only) + 1)

  return(l_est)
}

estimate_r <- function(N_WT_only, N_M_only, N_d_neg, N_d_pos,
                       l_est,
                       a_est, b_est, c_est) {

  # Get starting value of r (rough estimates)
  r_start_init <- log((N_M_only + N_d_pos) / (N_d_neg + N_WT_only) + 1)
  # Correct with error estimates
  r_start_cor <- r_start_init - (a_est + b_est * l_est + c_est * l_est)
  # Truncate r at 0, as this is the minimum possible value in the model
  r_start <- max(r_start_cor, MIN_START_PAR)

  # Optimize for r
  optim_res <- stats::optim(
    par = log(r_start),
    method = "BFGS",
    fn = full_log_lik_test,
    gr = grad_full_log_lik_test,
    l_est = l_est,
    N_WT_only = N_WT_only,
    N_M_only = N_M_only,
    N_d_neg = N_d_neg,
    N_d_pos = N_d_pos,
    a_est = a_est,
    b_est = b_est,
    c_est = c_est,
    control = list(
      fnscale = -1,
      reltol = RELTOL
    )
  )

  # Get optimal value of r
  r_est <- exp(optim_res$par)

  # Check if r=0 gives a better likelihood:
  # Calculate LL for simple model (H0)
  LL_0 <- full_log_lik(
    l_vec = l_est,
    r_vec = 0,
    a = a_est,
    b = b_est,
    c = c_est,
    N_WT_only_vec = N_WT_only,
    N_M_only_vec = N_M_only,
    N_d_neg_vec = N_d_neg,
    N_d_pos_vec = N_d_pos
  )

  # Calculate LL for simple model (H0)
  LL_A <- full_log_lik(
    l_vec = l_est,
    r_vec = r_est,
    a = a_est,
    b = b_est,
    c = c_est,
    N_WT_only_vec = N_WT_only,
    N_M_only_vec = N_M_only,
    N_d_neg_vec = N_d_neg,
    N_d_pos_vec = N_d_pos
  )

  if (LL_A < LL_0) {
    r_est <- 0
  }

  return(r_est)
}

find_abd_confidence_intervals <- function(N_WT_only_vec, N_M_only_vec, N_d_neg_vec, N_d_pos_vec,
                                          l_est_vec, a_est, b_est, c_est, alpha) {
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

get_r_CI_simple <- function(l_est, r_est,
                            a_est, b_est, c_est,
                            N_WT_only, N_M_only, N_d_neg, N_d_pos,
                            alpha) {
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

  # Find upper and lower bound of r
  if (ll_ratio_simple(par_0 = 0, par = "r") > stats::qchisq(1 - alpha, 1)) {
    uni_res <- stats::uniroot(function(x) ll_ratio_simple(x, par = "r") - stats::qchisq(1 - alpha, 1),
      interval = c(0, r_est),
      tol = TOL_0
    )
    lower <- uni_res$root
  } else {
    lower <- 0
  }

  uni_res <- stats::uniroot(function(x) ll_ratio_simple(x, par = "r") - stats::qchisq(1 - alpha, 1),
    interval = c(r_est, 1),
    tol = TOL_0,
    extendInt = "upX"
  )

  upper <- uni_res$root

  return(list(
    lower = lower,
    upper = upper
  ))
}

test_r_equal_0_simple <- function(N_WT_only, N_M_only, N_d_neg, N_d_pos,
                                  l_est, r_est,
                                  a_est, b_est, c_est,
                                  include_mutant_CI, alpha) {
  # Calculate LL for simple model (H0)
  LL_0 <- full_log_lik(
    l_vec = l_est,
    r_vec = 0,
    a = a_est,
    b = b_est,
    c = c_est,
    N_WT_only_vec = N_WT_only,
    N_M_only_vec = N_M_only,
    N_d_neg_vec = N_d_neg,
    N_d_pos_vec = N_d_pos
  )

  # Calculate LL for simple model (H0)
  LL_A <- full_log_lik(
    l_vec = l_est,
    r_vec = r_est,
    a = a_est,
    b = b_est,
    c = c_est,
    N_WT_only_vec = N_WT_only,
    N_M_only_vec = N_M_only,
    N_d_neg_vec = N_d_neg,
    N_d_pos_vec = N_d_pos
  )

  # P val for r=0
  # r MLE under the hypothesis: H0: r=0
  # Get test statistic and p-value
  r_LLR_test_stat <- -2 * (LL_0 - LL_A)
  p_val <- 1 - stats::pchisq(r_LLR_test_stat, 1)

  # Gather results
  total_droplets <- N_WT_only + N_M_only + N_d_neg + N_d_pos
  allele_frequency <- r_est / (r_est + l_est)
  estimated_total_mutant_molecules <- r_est * total_droplets
  estimated_total_wildtype_molecules <- l_est * total_droplets
  mutation_detected <- p_val < alpha

  # Collect results
  res <- list(
    # Estimated values for r
    mutant_molecules_per_droplet = r_est,
    # Estimated values for l:
    wildtype_molecules_per_droplet = l_est,
    # Test statistics
    p_val = p_val,
    test_statistic = r_LLR_test_stat,
    # Other results:
    mutation_detected = mutation_detected,
    allele_frequency = allele_frequency,
    total_mutant_molecules = estimated_total_mutant_molecules,
    total_wildtype_molecules = estimated_total_wildtype_molecules
  )

  if (include_mutant_CI) {
    r_CI <- get_r_CI_simple(
      l_est = l_est,
      r_est = r_est,
      a_est = a_est,
      b_est = b_est,
      c_est = c_est,
      N_WT_only = N_WT_only,
      N_M_only = N_M_only,
      N_d_neg = N_d_neg,
      N_d_pos = N_d_pos,
      alpha = alpha
    )

    # Extract bounds
    r_CI_lower <- r_CI$lower
    r_CI_upper <- r_CI$upper

    # Calculate CI on the number of "real" tumor molecules
    total_mutant_molecules_CI_lower <- r_CI_lower * total_droplets
    total_mutant_molecules_CI_upper <- r_CI_upper * total_droplets

    # Add CI's to results
    res <- append(res, list(
      # Estimated values for r
      mutant_molecules_per_droplet_CI_lower = r_CI_lower,
      mutant_molecules_per_droplet_CI_upper = r_CI_upper,
      # Allele frequency
      allele_frequency_CI_lower = r_CI_lower / (r_CI_lower + l_est + TOL_0),
      allele_frequency_CI_upper = r_CI_upper / (r_CI_upper + l_est + TOL_0),
      # Total number of molecules
      total_mutant_molecules_CI_lower = total_mutant_molecules_CI_lower,
      total_mutant_molecules_CI_upper = total_mutant_molecules_CI_upper
    ))
  }

  return(res)
}
