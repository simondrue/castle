# Probability function ----------------------------------------------------

P_d_neg <- function(l, a, c, r) {
  p <- exp(-l * (1 + c) - a - r)
  return(pmin(pmax(0, p), 1))
}

P_WT_only <- function(l, a, b, c, r) {
  p <- exp(-l * (1 + c) - a - r) * (exp(l * exp(-b)) - 1)
  return(pmin(pmax(0, p), 1))
}

P_M_only <- function(l, a, c, r) {
  p <- exp(-l) * (1 - exp(-r - a - l * c))
  return(pmin(pmax(0, p), 1))
}

P_d_pos <- function(l, a, b, c, r) {
  p <- 1 - exp(-l * (1 + c) - a - r) * (exp(l * exp(-b)) - 1) - exp(-l)
  return(pmin(pmax(0, p), 1))
}

# Derivatives of probability functions ------------------------------------

# d/dr

dr_P_d_neg <- function(l, a, c, r) {
  return(-exp(-l * (1 + c) - a - r))
}

dr_P_WT_only <- function(l, a, b, c, r) {
  return(-(exp(-l * (1 + c) - a - r) * (exp(l * exp(-b)) - 1)))
}

dr_P_M_only <- function(l, a, c, r) {
  return(exp(-l) * exp(-r - a - l * c))
}

dr_P_d_pos <- function(l, a, b, c, r) {
  return(exp(-l * (1 + c) - a - r) * (exp(l * exp(-b)) - 1))
}

# d/dl

dl_P_d_neg <- function(l, a, c, r) {
  return(-(exp(-l * (1 + c) - a - r) * (1 + c)))
}

dl_P_WT_only <- function(l, a, b, c, r) {
  return(exp(-l * (1 + c) - a - r) * (exp(l * exp(-b)) * exp(-b)) -
           exp(-l * (1 + c) - a - r) * (1 + c) * (exp(l * exp(-b)) - 1))
}

dl_P_M_only <- function(l, a, c, r) {
  return(exp(-l) * (exp(-r - a - l * c) * c) - exp(-l) * (1 - exp(-r - a - l * c)))
}

dl_P_d_pos <- function(l, a, b, c, r) {
  return(-(exp(-l * (1 + c) - a - r) * (exp(l * exp(-b)) * exp(-b)) -
             exp(-l * (1 + c) - a - r) * (1 + c) * (exp(l * exp(-b)) - 1) - exp(-l)))
}

# d/da

da_P_d_neg <- function(l, a, c, r) {
  return(-exp(-l * (1 + c) - a - r))
}

da_P_WT_only <- function(l, a, b, c, r) {
  return(-(exp(-l * (1 + c) - a - r) * (exp(l * exp(-b)) - 1)))
}

da_P_M_only <- function(l, a, c, r) {
  return(exp(-l) * exp(-r - a - l * c))
}

da_P_d_pos <- function(l, a, b, c, r) {
  return(exp(-l * (1 + c) - a - r) * (exp(l * exp(-b)) - 1))
}


# d/db

db_P_d_neg <- function(l, a, c, r) {
  return(0)
}

db_P_WT_only <- function(l, a, b, c, r) {
  return(-(exp(-l * (1 + c) - a - r) * (exp(l * exp(-b)) * (l * exp(-b)))))
}

db_P_M_only <- function(l, a, c, r) {
  return(0)
}

db_P_d_pos <- function(l, a, b, c, r) {
  return(exp(-l * (1 + c) - a - r) * (exp(l * exp(-b)) * (l * exp(-b))))
}

# d/dc

dc_P_d_neg <- function(l, a, c, r) {
  return(-(exp(-l * (1 + c) - a - r) * l))
}

dc_P_WT_only <- function(l, a, b, c, r) {
  return(-(exp(-l * (1 + c) - a - r) * l * (exp(l * exp(-b)) - 1)))
}

dc_P_M_only <- function(l, a, c, r) {
  return(exp(-l) * (exp(-r - a - l * c) * l))
}

dc_P_d_pos <- function(l, a, b, c, r) {
  return(exp(-l * (1 + c) - a - r) * l * (exp(l * exp(-b)) - 1))
}
