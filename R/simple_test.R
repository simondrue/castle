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

  # Make appropiate corrections if any parameter is less than MIN_START_PAR
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

train_simple_ddpcr_model <- function(training_samples) {
  # Check if expected columns are present
  if (!all(c("N_WT_only", "N_M_only", "N_d_neg", "N_d_pos") %in% colnames(training_samples))) {
    missing_cols <-
      setdiff(
        c("N_WT_only", "N_M_only", "N_d_neg", "N_d_pos"),
        colnames(training_samples)
      )
    stop(paste(missing_cols, "are missing from training_samples"))
  }

  # Unpack data
  N_WT_only_vec <- training_samples$N_WT_only
  N_M_only_vec <- training_samples$N_M_only
  N_d_neg_vec <- training_samples$N_d_neg
  N_d_pos_vec <- training_samples$N_d_pos

  # Estimate lambdas
  l_est_vec <- log((N_WT_only_vec + N_d_pos_vec) / (N_d_neg_vec + N_M_only_vec) + 1)

  # Get starting values (rough estimates) of a, b and c
  start_values <- get_start_values_abc(
    N_WT_only_vec = N_WT_only_vec,
    N_M_only_vec = N_M_only_vec,
    N_d_neg_vec = N_d_neg_vec,
    N_d_pos_vec = N_d_pos_vec,
    l_est_vec = l_est_vec
  )

  # Set start at mean of individual samples
  start_vec <- log(c(
    start_values$a_start,
    start_values$b_start,
    start_values$c_start
  ))

  # Estimate alpha, beta and gamma by maximizing log-likelihood
  optim_res <- stats::optim(
    par = start_vec,
    method = "BFGS",
    fn = full_log_lik_train,
    gr = grad_full_log_lik_train,
    l_est_vec = l_est_vec,
    N_WT_only_vec = N_WT_only_vec,
    N_M_only_vec = N_M_only_vec,
    N_d_neg_vec = N_d_neg_vec,
    N_d_pos_vec = N_d_pos_vec,
    control = list(
      fnscale = -1,
      reltol = RELTOL
    )
  )

  par_est <- exp(optim_res$par)

  res <- list(
    a_est = par_est[1],
    b_est = par_est[2],
    c_est = par_est[3],
    l_est_vec = l_est_vec
  )
  return(res)
}

estimate_l <- function(N_WT_only, N_M_only, N_d_neg, N_d_pos) {

  # Estimate l
  l_est <- log((N_WT_only + N_d_pos) / (N_d_neg + N_M_only) + 1)

  return(l_est)
}

estimate_r <- function(N_WT_only, N_M_only, N_d_neg, N_d_pos,
                       l_est,
                       a_train_est, b_train_est, c_train_est) {

  # Get starting value of r (rough estimates)
  r_start_init <- log((N_M_only + N_d_pos) / (N_d_neg + N_WT_only) + 1)
  # Correct with error estimates
  r_start_cor <- r_start_init - (a_train_est + b_train_est * l_est + c_train_est * l_est)
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
    a_train_est = a_train_est,
    b_train_est = b_train_est,
    c_train_est = c_train_est,
    control = list(
      fnscale = -1,
      reltol = RELTOL
    )
  )

  # Get optimal value of r
  r_est <- exp(optim_res$par)

  return(r_est)
}

test_tumor_sample_simple <- function(test_samples,
                                     model,
                                     alpha = 0.05) {

  # Unpack parameters
  N_WT_only <- test_samples$N_WT_only
  N_M_only <- test_samples$N_M_only
  N_d_neg <- test_samples$N_d_neg
  N_d_pos <- test_samples$N_d_pos

  a_train_est <- model$a_est
  b_train_est <- model$b_est
  c_train_est <- model$c_est

  # Estimate parameters
  l_est <- estimate_l(
    N_WT_only = N_WT_only,
    N_M_only = N_M_only,
    N_d_neg = N_d_neg,
    N_d_pos = N_d_pos
  )

  r_est <- estimate_r(
    N_WT_only = N_WT_only,
    N_M_only = N_M_only,
    N_d_neg = N_d_neg,
    N_d_pos = N_d_pos,
    l_est = l_est,
    a_train_est = a_train_est,
    b_train_est = b_train_est,
    c_train_est = c_train_est
  )

  bounds_res <- find_lr_confidence_intervals(
    N_WT_only = N_WT_only,
    N_M_only = N_M_only,
    N_d_neg = N_d_neg,
    N_d_pos = N_d_pos,
    l_est = l_est,
    a_est = a_train_est,
    b_est = b_train_est,
    c_est = c_train_est,
    r_est = r_est,
    alpha = alpha
  )

  # Gather results
  total_droplets <- N_WT_only + N_M_only + N_d_neg + N_d_pos
  total_tumor_molecules_expected <- r_est * total_droplets
  total_tumor_molecules_CI_lower <- bounds_res$r_CI_lower * total_droplets
  total_tumor_molecules_CI_upper <- bounds_res$r_CI_upper * total_droplets

  is_tumor_positive <- bounds_res$r_less_than_0_pval < alpha

  # Return results
  res <- list(
    # Estimated values for r:
    r_est = r_est,
    r_CI_lower = bounds_res$r_CI_lower,
    r_CI_upper = bounds_res$r_CI_upper,
    # Estimated values for l:
    l_est = l_est,
    l_CI_lower = bounds_res$l_CI_lower,
    l_CI_upper = bounds_res$l_CI_upper,
    # Test statistics
    pval_r_leq_0 = bounds_res$r_less_than_0_pval,
    r_LLR_test_stat = bounds_res$r_LLR_test_stat,
    # Other results:
    is_tumor_positive = is_tumor_positive,
    total_tumor_molecules_expected = total_tumor_molecules_expected,
    total_tumor_molecules_CI_lower = total_tumor_molecules_CI_lower,
    total_tumor_molecules_CI_upper = total_tumor_molecules_CI_upper,
    total_droplets = total_droplets
  )

  return(res)
}
