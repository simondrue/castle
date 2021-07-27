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
