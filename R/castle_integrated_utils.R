# Function for calculating log-likelihood on grid
get_abc_grid_log_lik <- function(N_WT_only_vec,
                                 N_M_only_vec,
                                 N_d_neg_vec,
                                 N_d_pos_vec,
                                 l_vec,
                                 r_vec = rep(0, length(l_vec)),
                                 abc_grid) {

  # Put grid into vectors
  a_vec <- abc_grid[, 1]
  b_vec <- abc_grid[, 2]
  c_vec <- abc_grid[, 3]

  # Calculate log-likelihood by looping over individuals and keeping track of sum
  log_lik_sum_vec <- 0
  for (ind in 1:length(N_WT_only_vec)) {
    # Get l and r for individual
    l <- l_vec[ind]
    r <- r_vec[ind]

    # Get the log_P values (note vectorization in a, b and c)
    log_P_d_neg <- log(P_d_neg(l = l, a = a_vec, c = c_vec, r = r))
    log_P_WT_only <- log(P_WT_only(l = l, a = a_vec, b = b_vec, c = c_vec, r = r))
    log_P_d_pos <- log(P_d_pos(l = l, a = a_vec, b = b_vec, c = c_vec, r = r))
    log_P_M_only <- log(P_M_only(l = l, a = a_vec, c = c_vec, r = r))

    # If N=0 set log_P=0 (to avoid log_P=-Inf)
    if (N_d_neg_vec[ind] == 0) {
      log_P_d_neg <- 0
    }
    if (N_WT_only_vec[ind] == 0) {
      log_P_WT_only <- 0
    }
    if (N_d_pos_vec[ind] == 0) {
      log_P_d_pos <- 0
    }
    if (N_M_only_vec[ind] == 0) {
      log_P_M_only <- 0
    }

    # Calculate individual log-likelihood
    log_lik_indiv <- N_d_neg_vec[ind] * log_P_d_neg +
      N_WT_only_vec[ind] * log_P_WT_only +
      N_d_pos_vec[ind] * log_P_d_pos +
      N_M_only_vec[ind] * log_P_M_only

    # Add individual log_lik to sum
    log_lik_sum_vec <- log_lik_sum_vec + log_lik_indiv
  }

  return(log_lik_sum_vec)
}

integrated_log_lik <- function(l, r,
                               N_WT_only,
                               N_M_only,
                               N_d_neg,
                               N_d_pos,
                               abc_grid,
                               abc_grid_train_log_lik) {

  # Calculate test log likelihood on abc_grid
  abc_grid_test_log_lik <- get_abc_grid_log_lik(
    N_WT_only_vec = N_WT_only,
    N_M_only_vec = N_M_only,
    N_d_neg_vec = N_d_neg,
    N_d_pos_vec = N_d_pos,
    l_vec = l,
    r_vec = r,
    abc_grid = abc_grid
  )

  # Calculate integral (volume) on grid, i.e. integrated log_lik
  # Delta in a, b, and c
  delta_a <- diff(sort(unique(abc_grid[, 1])))[1]
  delta_b <- diff(sort(unique(abc_grid[, 2])))[1]
  delta_c <- diff(sort(unique(abc_grid[, 3])))[1]

  # (Log) Riemann sum to approximate integral
  integrated_log_lik <- sum_log_p(abc_grid_test_log_lik + abc_grid_train_log_lik) + log(delta_a) + log(delta_b) + log(delta_c)

  return(integrated_log_lik)
}

integrated_log_lik_test <- function(par, l,
                                    N_WT_only,
                                    N_M_only,
                                    N_d_neg,
                                    N_d_pos,
                                    abc_grid,
                                    abc_grid_train_log_lik) {

  # Unpack parameters
  r <- exp(par)

  # Calculate test log lik on abc_grid
  integrated_log_lik <- integrated_log_lik(
    l = l, r = r,
    N_WT_only = N_WT_only,
    N_M_only = N_M_only,
    N_d_neg = N_d_neg,
    N_d_pos = N_d_pos,
    abc_grid = abc_grid,
    abc_grid_train_log_lik = abc_grid_train_log_lik
  )

  return(integrated_log_lik)
}

get_r_CI_integrated <- function(r_est, l_est,
                                abc_grid,
                                abc_grid_train_log_lik,
                                N_WT_only,
                                N_M_only,
                                N_d_neg,
                                N_d_pos,
                                alpha) {

  # LLR with integrated likelihood function
  # Calculate LL for full (alternative) model - Note this is constant
  LL_A <- integrated_log_lik(
    l = l_est,
    r = r_est,
    N_WT_only = N_WT_only,
    N_M_only = N_M_only,
    N_d_neg = N_d_neg,
    N_d_pos = N_d_pos,
    abc_grid_train_log_lik = abc_grid_train_log_lik,
    abc_grid = abc_grid
  )

  # Function for calculating integrated LLR
  ll_ratio_int <- function(r_0) {
    LL_0 <- integrated_log_lik(
      l = l_est,
      r = r_0,
      N_WT_only = N_WT_only,
      N_M_only = N_M_only,
      N_d_neg = N_d_neg,
      N_d_pos = N_d_pos,
      abc_grid_train_log_lik = abc_grid_train_log_lik,
      abc_grid = abc_grid
    )
    return(-2 * (LL_0 - LL_A))
  }

  # Find confidence intervals
  # Lower bound
  if (ll_ratio_int(TOL_0) > stats::qchisq(1 - alpha, 1)) {
    uni_res <- stats::uniroot(function(x) ll_ratio_int(x) - stats::qchisq(1 - alpha, 1),
      interval = c(TOL_0, r_est),
      tol = TOL_0
    )
    lower <- uni_res$root
  } else {
    lower <- 0
  }

  # Upper bound
  uni_res <- stats::uniroot(function(x) ll_ratio_int(x) - stats::qchisq(1 - alpha, 1),
    interval = c(r_est, 1),
    tol = TOL_0,
    extendInt = "upX"
  )
  upper <- uni_res$root

  # Return result
  res <- list(
    lower = lower,
    upper = upper
  )
  return(res)
}

estimate_r_integrated <- function(N_WT_only, N_M_only, N_d_neg, N_d_pos,
                                  l_est,
                                  integrated_model) {
  abc_grid <- integrated_model$abc_grid
  abc_grid_train_log_lik <- integrated_model$abc_grid_train_log_lik

  # Get start guess
  r_start <- estimate_r(
    N_WT_only = N_WT_only,
    N_M_only = N_M_only,
    N_d_neg = N_d_neg,
    N_d_pos = N_d_pos,
    l_est = l_est,
    a_train_est = integrated_model$a_est,
    b_train_est = integrated_model$b_est,
    c_train_est = integrated_model$c_est
  ) %>% max(MIN_START_PAR)

  # Now optimize the integrated likelihood numerically
  optim_res_int <- stats::optim(
    par = log(r_start),
    l = l_est,
    fn = integrated_log_lik_test,
    gr = NULL,
    method = "BFGS",
    N_WT_only = N_WT_only,
    N_M_only = N_M_only,
    N_d_neg = N_d_neg,
    N_d_pos = N_d_pos,
    abc_grid = abc_grid,
    abc_grid_train_log_lik = abc_grid_train_log_lik,
    control = list(
      fnscale = -1,
      parscale = log(r_start),
      reltol = RELTOL
    )
  )

  # Unpack result and return
  r_est <- exp(optim_res_int$par)

  # If likelihood is greater at r=0 -> change value
  if (optim_res_int$value < integrated_log_lik(
    l = l_est,
    r = 0,
    N_WT_only = N_WT_only,
    N_M_only = N_M_only,
    N_d_neg = N_d_neg,
    N_d_pos = N_d_pos,
    abc_grid_train_log_lik = abc_grid_train_log_lik,
    abc_grid = abc_grid
  )) {
    r_est <- 0
  }

  return(r_est)
}

test_r_equal_0_integrated <- function(l_est, r_est,
                                      N_WT_only, N_M_only, N_d_neg, N_d_pos,
                                      abc_grid_train_log_lik, abc_grid,
                                      include_mutant_CI, alpha) {
  # Test hypothesis H0: r == 0 vs HA: r > 0
  # LLR with integrated likelihood function

  # Calculate LL for full model (HA)
  LL_A <- integrated_log_lik(
    l = l_est,
    r = r_est,
    N_WT_only = N_WT_only,
    N_M_only = N_M_only,
    N_d_neg = N_d_neg,
    N_d_pos = N_d_pos,
    abc_grid_train_log_lik = abc_grid_train_log_lik,
    abc_grid = abc_grid
  )

  # Calculate LL for simple model (H0)
  LL_0 <- integrated_log_lik(
    l = l_est,
    r = 0,
    N_WT_only = N_WT_only,
    N_M_only = N_M_only,
    N_d_neg = N_d_neg,
    N_d_pos = N_d_pos,
    abc_grid_train_log_lik = abc_grid_train_log_lik,
    abc_grid = abc_grid
  )

  # Get test statistic and p-value
  r_LLR_test_stat <- -2 * (LL_0 - LL_A)
  p_val <- 1 - stats::pchisq(r_LLR_test_stat, 1)

  # Gather results
  total_droplets <- N_WT_only + N_M_only + N_d_neg + N_d_pos
  allele_frequency = r_est / (r_est + l_est)
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
    # Find bounds on r
    r_CI <- get_r_CI_integrated(
      r_est = r_est,
      l_est = l_est,
      N_WT_only = N_WT_only,
      N_M_only = N_M_only,
      N_d_neg = N_d_neg,
      N_d_pos = N_d_pos,
      alpha = alpha,
      abc_grid = abc_grid,
      abc_grid_train_log_lik = abc_grid_train_log_lik
    )

    # Extract bounds
    r_CI_lower <- r_CI$lower
    r_CI_upper <- r_CI$upper

    # Calculate CI on the number of "real" tumor molecules
    total_mutant_molecules_CI_lower <- r_CI_lower * total_droplets
    total_mutant_molecules_CI_upper <- r_CI_upper * total_droplets

    1 / (1 + l_est/r_est)

    # Add CI's to results
    res <- append(res, list(
      # Estimated values for r
      mutant_molecules_per_droplet_CI_lower = r_CI_lower,
      mutant_molecules_per_droplet_CI_upper = r_CI_upper,
      # Allele frequency
      allele_frequency_CI_lower = r_CI_lower / (r_CI_lower + l_est),
      allele_frequency_CI_upper = r_CI_upper / (r_CI_upper + l_est),
      # Total number of molecules
      total_mutant_molecules_CI_lower = total_mutant_molecules_CI_lower,
      total_mutant_molecules_CI_upper = total_mutant_molecules_CI_upper
    ))
  }

  return(res)
}
