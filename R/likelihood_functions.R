# Likelihood functions ----------------------------------------------------

# Full LogLik over multiple samples with shared a, b and c, and individual l and r
full_log_lik <- function(l_vec, r_vec = rep(0, length(l_vec)), a, b, c,
                         N_WT_only_vec, N_M_only_vec, N_d_neg_vec, N_d_pos_vec) {

  # Get the individual log_P_i values
  log_P_d_neg <- log(P_d_neg(l = l_vec, a = a, c = c, r = r_vec))
  log_P_WT_only <- log(P_WT_only(l = l_vec, a = a, b = b, c = c, r = r_vec))
  log_P_d_pos <- log(P_d_pos(l = l_vec, a = a, b = b, c = c, r = r_vec))
  log_P_M_only <- log(P_M_only(l = l_vec, a = a, c = c, r = r_vec))

  # If N=0 set log_P=0 (to avoid log_P=-Inf)
  log_P_d_neg[N_d_neg_vec == 0] <- 0
  log_P_WT_only[N_WT_only_vec == 0] <- 0
  log_P_d_pos[N_d_pos_vec == 0] <- 0
  log_P_M_only[N_M_only_vec == 0] <- 0

  # Calculate log likelihood value as sum_i N_i*log_P_i
  log_lik <- N_d_neg_vec %*% log_P_d_neg +
    N_d_pos_vec %*% log_P_d_pos +
    N_WT_only_vec %*% log_P_WT_only +
    N_M_only_vec %*% log_P_M_only

  return(log_lik)
}

# LogLik for optimizing a, b and c in multiple train samples
full_log_lik_train <- function(par,
                               l_est_vec,
                               N_WT_only_vec, N_M_only_vec, N_d_neg_vec, N_d_pos_vec) {
  # Unpack parameters
  par <- exp(par)
  a_est <- par[1]
  b_est <- par[2]
  c_est <- par[3]

  log_lik <- full_log_lik(
    l_vec = l_est_vec,
    a = a_est,
    b = b_est,
    c = c_est,
    N_WT_only_vec = N_WT_only_vec,
    N_M_only_vec = N_M_only_vec,
    N_d_neg_vec = N_d_neg_vec,
    N_d_pos_vec = N_d_pos_vec
  )
  return(log_lik)
}

# LogLik for optimizing l and r in test sample
full_log_lik_test <- function(par, l_est,
                              N_WT_only, N_M_only, N_d_neg, N_d_pos,
                              a_est, b_est, c_est) {
  r_est <- exp(par)

  log_lik <- full_log_lik(
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

  return(log_lik)
}

ll_ratio <- function(par_0, par,
                     l_est_vec,
                     r_est_vec,
                     a_est, b_est, c_est,
                     N_WT_only_vec,
                     N_M_only_vec,
                     N_d_neg_vec,
                     N_d_pos_vec) {

  # Want to calculate log( L(x_0)/L(x_MLE) ) [ = ll(x_0) - ll(x_MLE) = ll_H0 - ll_H1 ]
  # Using this we can get the confidence of value x_0 for the relevant parameter

  # Set all parameter equal under H0 and H1
  l_est_vec_0 <- l_est_vec
  r_est_vec_0 <- r_est_vec
  a_est_0 <- a_est
  b_est_0 <- b_est
  c_est_0 <- c_est

  # Change only relevant parameter in H0
  if (par == "l") {
    l_est_vec_0 <- par_0
  } else if (par == "r") {
    r_est_vec_0 <- par_0
  } else if (par == "a") {
    a_est_0 <- par_0
  } else if (par == "b") {
    b_est_0 <- par_0
  } else if (par == "c") {
    c_est_0 <- par_0
  }

  log_lik_0 <- full_log_lik(
    l_vec = l_est_vec_0,
    r_vec = r_est_vec_0,
    a = a_est_0,
    b = b_est_0,
    c = c_est_0,
    N_WT_only_vec = N_WT_only_vec,
    N_M_only_vec = N_M_only_vec,
    N_d_neg_vec = N_d_neg_vec,
    N_d_pos_vec = N_d_pos_vec
  )

  log_lik_1 <- full_log_lik(
    l_vec = l_est_vec,
    r_vec = r_est_vec,
    a = a_est,
    b = b_est,
    c = c_est,
    N_WT_only_vec = N_WT_only_vec,
    N_M_only_vec = N_M_only_vec,
    N_d_neg_vec = N_d_neg_vec,
    N_d_pos_vec = N_d_pos_vec
  )

  return(-2 * (log_lik_0 - log_lik_1))
}


# Derivatives of Likelihood functions -------------------------------------

dr_full_log_lik <- function(l, r, a, b, c,
                            N_WT_only, N_M_only, N_d_neg, N_d_pos) {

  # d/dr
  dr_P_d_neg <- dr_P_d_neg(l, a, c, r) / P_d_neg(l, a, c, r)
  dr_P_d_pos <- dr_P_d_pos(l, a, b, c, r) / P_d_pos(l, a, b, c, r)
  dr_P_WT_only <- dr_P_WT_only(l, a, b, c, r) / P_WT_only(l, a, b, c, r)
  dr_P_M_only <- dr_P_M_only(l, a, c, r) / P_M_only(l, a, c, r)

  # If N=0 set dr_P=0 (to avoid dr_P=-Inf)
  dr_P_d_neg[N_d_neg == 0] <- 0
  dr_P_WT_only[N_WT_only == 0] <- 0
  dr_P_d_pos[N_d_pos == 0] <- 0
  dr_P_M_only[N_M_only == 0] <- 0

  dr_ll <- N_d_neg %*% dr_P_d_neg +
    N_d_pos %*% dr_P_d_pos +
    N_WT_only %*% dr_P_WT_only +
    N_M_only %*% dr_P_M_only

  return(dr_ll)
}

da_full_log_lik <- function(l, r, a, b, c,
                            N_WT_only, N_M_only, N_d_neg, N_d_pos) {

  # d/da
  da_P_d_neg <- da_P_d_neg(l, a, c, r) / P_d_neg(l, a, c, r)
  da_P_d_pos <- da_P_d_pos(l, a, b, c, r) / P_d_pos(l, a, b, c, r)
  da_P_WT_only <- da_P_WT_only(l, a, b, c, r) / P_WT_only(l, a, b, c, r)
  da_P_M_only <- da_P_M_only(l, a, c, r) / P_M_only(l, a, c, r)

  # If N=0 set da_P=0 (to avoid da_P=-Inf)
  da_P_d_neg[N_d_neg == 0] <- 0
  da_P_WT_only[N_WT_only == 0] <- 0
  da_P_d_pos[N_d_pos == 0] <- 0
  da_P_M_only[N_M_only == 0] <- 0

  da_ll <- N_d_neg %*% da_P_d_neg +
    N_d_pos %*% da_P_d_pos +
    N_WT_only %*% da_P_WT_only +
    N_M_only %*% da_P_M_only

  return(da_ll)
}

db_full_log_lik <- function(l, r, a, b, c,
                            N_WT_only, N_M_only, N_d_neg, N_d_pos) {

  # d/db
  db_P_d_neg <- db_P_d_neg(l, a, c, r) / P_d_neg(l, a, c, r)
  db_P_d_pos <- db_P_d_pos(l, a, b, c, r) / P_d_pos(l, a, b, c, r)
  db_P_WT_only <- db_P_WT_only(l, a, b, c, r) / P_WT_only(l, a, b, c, r)
  db_P_M_only <- db_P_M_only(l, a, c, r) / P_M_only(l, a, c, r)

  # If N=0 set db_P=0 (to avoid db_P=-Inf)
  db_P_d_neg[N_d_neg == 0] <- 0
  db_P_WT_only[N_WT_only == 0] <- 0
  db_P_d_pos[N_d_pos == 0] <- 0
  db_P_M_only[N_M_only == 0] <- 0

  db_ll <- N_d_neg %*% db_P_d_neg +
    N_d_pos %*% db_P_d_pos +
    N_WT_only %*% db_P_WT_only +
    N_M_only %*% db_P_M_only

  return(db_ll)
}

dc_full_log_lik <- function(l, r, a, b, c,
                            N_WT_only, N_M_only, N_d_neg, N_d_pos) {

  # d/dc
  dc_P_d_neg <- dc_P_d_neg(l, a, c, r) / P_d_neg(l, a, c, r)
  dc_P_d_pos <- dc_P_d_pos(l, a, b, c, r) / P_d_pos(l, a, b, c, r)
  dc_P_WT_only <- dc_P_WT_only(l, a, b, c, r) / P_WT_only(l, a, b, c, r)
  dc_P_M_only <- dc_P_M_only(l, a, c, r) / P_M_only(l, a, c, r)

  # If N=0 set dc_P=0 (to avoid dc_P=-Inf)
  dc_P_d_neg[N_d_neg == 0] <- 0
  dc_P_WT_only[N_WT_only == 0] <- 0
  dc_P_d_pos[N_d_pos == 0] <- 0
  dc_P_M_only[N_M_only == 0] <- 0

  dc_ll <- N_d_neg %*% dc_P_d_neg +
    N_d_pos %*% dc_P_d_pos +
    N_WT_only %*% dc_P_WT_only +
    N_M_only %*% dc_P_M_only

  return(dc_ll)
}

# Gradient of LogLik for optimizing a, b and c in multiple train samples
grad_full_log_lik_train <- function(par,
                                    l_est_vec,
                                    N_WT_only_vec, N_M_only_vec, N_d_neg_vec, N_d_pos_vec) {
  # Unpack parameters
  a_est <- exp(par[1])
  b_est <- exp(par[2])
  c_est <- exp(par[3])

  da_ll <- da_full_log_lik(
    l = l_est_vec,
    r = 0,
    a = a_est,
    b = b_est,
    c = c_est,
    N_WT_only = N_WT_only_vec,
    N_M_only = N_M_only_vec,
    N_d_neg = N_d_neg_vec,
    N_d_pos = N_d_pos_vec
  )

  db_ll <- db_full_log_lik(
    l = l_est_vec,
    r = 0,
    a = a_est,
    b = b_est,
    c = c_est,
    N_WT_only = N_WT_only_vec,
    N_M_only = N_M_only_vec,
    N_d_neg = N_d_neg_vec,
    N_d_pos = N_d_pos_vec
  )

  dc_ll <- dc_full_log_lik(
    l = l_est_vec,
    r = 0,
    a = a_est,
    b = b_est,
    c = c_est,
    N_WT_only = N_WT_only_vec,
    N_M_only = N_M_only_vec,
    N_d_neg = N_d_neg_vec,
    N_d_pos = N_d_pos_vec
  )

  return(c(da_ll, db_ll, dc_ll) * exp(par))
}


# Gradient of LogLik for optimizing r and l in test sample
grad_full_log_lik_test <- function(par, l_est,
                                   N_WT_only, N_M_only, N_d_neg, N_d_pos,
                                   a_est, b_est, c_est) {
  # Unpack parameter
  r_est <- exp(par)

  # Calculate gradient
  d_dr <- dr_full_log_lik(
    l = l_est,
    r = r_est,
    a = a_est,
    b = b_est,
    c = c_est,
    N_WT_only = N_WT_only,
    N_M_only = N_M_only,
    N_d_neg = N_d_neg,
    N_d_pos = N_d_pos
  )

  return(d_dr * exp(par))
}
