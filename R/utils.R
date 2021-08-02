#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL

sum_log_p <- function(lp_vec) {
  lp_max <- max(lp_vec)
  return(lp_max + log(sum(exp(lp_vec - lp_max))))
}

check_input_samples <- function(samples) {
  # Check if expected columns are present
  expected_columns = c("WildtypeOnlyDroplets", "MutantOnlyDroplets", "DoubleNegativeDroplets", "DoublePositiveDroplets")
  if (!all(expected_columns %in% colnames(samples))) {
    missing_cols <-
      setdiff(expected_columns, colnames(samples))
    stop(paste("[", paste(missing_cols, collapse = ", "), "], are missing from test_samples.", collapse = ""))
  }

  # Check if any samples are empty
  totals <-
    samples %>%
    dplyr::mutate(
      total_droplets = .data$WildtypeOnlyDroplets + .data$MutantOnlyDroplets + .data$DoubleNegativeDroplets + .data$DoublePositiveDroplets,
      total_WT_negative = .data$MutantOnlyDroplets + .data$DoubleNegativeDroplets
    )

  if (any(totals$total_droplets == 0)) {
    stop(paste(
      "Sample(s) [",
      paste(which(totals$total_droplets == 0), collapse = ","),
      "] are empty (total droplet count is 0). Remove these and rerun analysis."
    ))
  }

  if (any(totals$total_WT_negative == 0)) {
    stop(paste(
      "Sample(s) [",
      paste(which(totals$total_WT_negative == 0), collapse = ","),
      "] have no WT-negative droplets. Remove these and rerun analysis."
    ))
  }
}

get_l_CI <- function(N_M_only, N_d_neg, N_WT_only, N_d_pos, alpha) {
  # Find bounds on l
  # Model as binomial with p = P(WT=0)=exp(-l)
  WT_neg <- N_M_only + N_d_neg
  n_drops <- N_WT_only + N_M_only + N_d_neg + N_d_pos

  binom_res <- stats::binom.test(x = WT_neg, n = n_drops, conf.level = 1 - alpha)

  lower <- -log(binom_res$conf.int[[2]])
  upper <- -log(binom_res$conf.int[[1]])
  return(list(lower = lower, upper = upper))
}
