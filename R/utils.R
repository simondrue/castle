
#' Calculate sum of probabilities in log-space.
#'
#' @param lp_vec numeric vector
#'
#' @return numeric value
#' @export
#'
#' @examples
sum_log_p <- function(lp_vec) {
  lp_max <- max(lp_vec)
  return(lp_max + log(sum(exp(lp_vec - lp_max))))
}
