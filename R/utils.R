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
  if (!all(c("N_WT_only", "N_M_only", "N_d_neg", "N_d_pos") %in% colnames(samples))) {
    missing_cols <-
      setdiff(
        c("N_WT_only", "N_M_only", "N_d_neg", "N_d_pos"),
        colnames(samples)
      )
    stop(paste(missing_cols, "are missing from test_samples."))
  }

  # Check if any samples are empty
  totals <-
    samples %>%
    dplyr::mutate(
      total_droplets = .data$N_WT_only + .data$N_M_only + .data$N_d_neg + .data$N_d_pos,
      total_WT_negative = .data$N_M_only + .data$N_d_neg
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
