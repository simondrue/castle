# Function for calculating log lik on grid
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

  # Calculate log_lik by looping over individuals and keeping track of sum
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

    # Calculate individual log_lik
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

  # Calculave test log lik on abc_grid
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

  # Integral
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

train_integrated_ddpcr_model <- function(N_WT_only_vec,
                                         N_M_only_vec,
                                         N_d_neg_vec,
                                         N_d_pos_vec,
                                         abc_grid_resolution = 25) {

  # Train parameters (get MLE on training data)
  par_est <- train_simple_ddpcr_model(
    N_WT_only_vec = N_WT_only_vec,
    N_M_only_vec = N_M_only_vec,
    N_d_neg_vec = N_d_neg_vec,
    N_d_pos_vec = N_d_pos_vec
  )


  l_est_vec <- par_est$l_est_vec
  a_est <- par_est$a_est
  b_est <- par_est$b_est
  c_est <- par_est$c_est


  # Get bounds of a, b and c to make a grid where training lik function has mass
  bounds <- find_abd_confidence_intervals(
    N_WT_only_vec = N_WT_only_vec,
    N_M_only_vec = N_M_only_vec,
    N_d_neg_vec = N_d_neg_vec,
    N_d_pos_vec = N_d_pos_vec,
    l_est_vec = l_est_vec,
    a_est = par_est$a_est,
    b_est = par_est$b_est,
    c_est = par_est$c_est,
    alpha = 1e-5
  )

  # Get symmetric buffer around estimate
  a_buffer <- max(a_est - bounds$a_CI_lower, bounds$a_CI_upper - a_est)
  b_buffer <- max(b_est - bounds$b_CI_lower, bounds$b_CI_upper - b_est)
  c_buffer <- max(c_est - bounds$c_CI_lower, bounds$c_CI_upper - c_est)

  # Get min and max values for error parameters
  a_min <- max(TOL_0, a_est - a_buffer)
  a_max <- a_est + a_buffer

  b_min <- max(TOL_0, b_est - b_buffer)
  b_max <- b_est + b_buffer

  c_min <- max(TOL_0, c_est - c_buffer)
  c_max <- c_est + c_buffer

  # Make intervals
  a_vec <- seq(from = a_min, to = a_max, length.out = abc_grid_resolution)
  b_vec <- seq(from = b_min, to = b_max, length.out = abc_grid_resolution)
  c_vec <- seq(from = c_min, to = c_max, length.out = abc_grid_resolution)

  # Grid of a, b and c
  abc_grid <- expand.grid(a_vec, b_vec, c_vec)

  # Get training log lik
  abc_grid_train_log_lik <- get_abc_grid_log_lik(
    N_WT_only_vec = N_WT_only_vec,
    N_M_only_vec = N_M_only_vec,
    N_d_neg_vec = N_d_neg_vec,
    N_d_pos_vec = N_d_pos_vec,
    l_vec = l_est_vec,
    abc_grid = abc_grid
  )

  # Gather results
  res <- list(
    a_est = a_est,
    b_est = b_est,
    c_est = c_est,
    l_est_vec = l_est_vec,
    abc_grid = abc_grid,
    abc_grid_train_log_lik = abc_grid_train_log_lik
  )

  return(res)
}


find_r_bounds_test_integrated <- function(r_est, l_est,
                                          abc_grid,
                                          abc_grid_train_log_lik,
                                          N_WT_only,
                                          N_M_only,
                                          N_d_neg,
                                          N_d_pos,
                                          alpha = 0.05) {

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
    r_CI_lower <- uni_res$root
  } else {
    r_CI_lower <- 0
  }

  # Upper bound
  uni_res <- stats::uniroot(function(x) ll_ratio_int(x) - stats::qchisq(1 - alpha, 1),
    interval = c(r_est, 1),
    tol = TOL_0,
    extendInt = "upX"
  )
  r_CI_upper <- uni_res$root

  # Return result
  res <- data.frame(
    r_CI_lower = r_CI_lower,
    r_CI_upper = r_CI_upper
  )
  return(res)
}

test_tumor_sample_integrated <- function(N_WT_only,
                                         N_M_only,
                                         N_d_neg,
                                         N_d_pos,
                                         abc_grid,
                                         abc_grid_train_log_lik,
                                         alpha = 0.05,
                                         include_r_confidence_interval = TRUE,
                                         include_l_confidence_interval = TRUE) {

  # Get MLE of a, b and c
  abc_MLE <- abc_grid[which.max(abc_grid_train_log_lik), ][1, ]
  a_train_est <- abc_MLE$Var1
  b_train_est <- abc_MLE$Var2
  c_train_est <- abc_MLE$Var3

  # Estimate test parameters

  # Get starting value from simple method
  par_est <- estimate_parameters_test(
    N_WT_only = N_WT_only,
    N_M_only = N_M_only,
    N_d_neg = N_d_neg,
    N_d_pos = N_d_pos,
    a_train_est = a_train_est,
    b_train_est = b_train_est,
    c_train_est = c_train_est
  )

  l_est <- par_est$l_est
  r_start <- max(par_est$r_est, MIN_START_PAR)

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
  r_less_than_0_pval <- 1 - stats::pchisq(r_LLR_test_stat, 1)


  # Gather results
  total_droplets <- N_WT_only + N_M_only + N_d_neg + N_d_pos
  total_tumor_molecules_expected <- r_est * total_droplets
  is_tumor_positive <- r_less_than_0_pval < alpha

  # Collect results
  res <- list(
    # Estimated values for r
    r_est = r_est,
    # Estimated values for l:
    l_est = l_est,
    # Test statistics
    pval_r_leq_0 = r_less_than_0_pval,
    r_LLR_test_stat = r_LLR_test_stat,
    # Other results:
    is_tumor_positive = is_tumor_positive,
    total_tumor_molecules_expected = total_tumor_molecules_expected,
    total_droplets = total_droplets
  )

  if (include_r_confidence_interval) {
    # Find bounds on r
    bounds_res <- find_r_bounds_test_integrated(
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
    r_CI_lower <- bounds_res$r_CI_lower
    r_CI_upper <- bounds_res$r_CI_upper

    # Calculate CI on the number of "real" tumor molecules
    total_tumor_molecules_CI_lower <- r_CI_lower * total_droplets
    total_tumor_molecules_CI_upper <- r_CI_upper * total_droplets

    # Add CI's to results
    res <- append(res, list(
      # Estimated values for r
      r_CI_lower = r_CI_lower,
      r_CI_upper = r_CI_upper,
      # Total number of molecules
      total_tumor_molecules_CI_lower = total_tumor_molecules_CI_lower,
      total_tumor_molecules_CI_upper = total_tumor_molecules_CI_upper
    ))
  }

  if (include_l_confidence_interval) {
    # Find bounds on l
    # Model as binomial with p = P(WT=0)=exp(-l)
    WT_neg <- N_M_only + N_d_neg
    n_drops <- N_WT_only + N_M_only + N_d_neg + N_d_pos

    binom_res <- stats::binom.test(x = WT_neg, n = n_drops, conf.level = 1 - alpha)

    l_CI_lower <- -log(binom_res$conf.int[[2]])
    l_CI_upper <- -log(binom_res$conf.int[[1]])

    # Add CI's to results
    res <- append(res, list(
      # Estimated values for l:
      l_CI_lower = l_CI_lower,
      l_CI_upper = l_CI_upper
    ))
  }

  return(res)
}
