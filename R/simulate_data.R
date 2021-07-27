#' Simulating ddPCR data
#'
#' @description Functions for simulating data from the statistical model underlying the \code{castle} package.
#' @param l Numeric value. Concentration of wild type DNA (molecules/droplet).
#' @param r Numeric value. Concentration of mutant DNA (molecules/droplet).
#' @param a,b,c Numeric values. Error rate parameters for the statistical model.
#' @param n_drops Integer value. The number of droplets to simulate in each sample.
#' @param n_samples Integer value. The number of samples to simulate.
#'
#' @return Simulated droplet/molecule counts as \code{data.frame}.
#' @export
#' @examples
#' # Simulate droplet counts
#'
#' # Single cancer negative sample (r = 0)
#' simulate_droplet_counts(1, 0, 0.01, 0.01, 0.01, 14000)
#'
#' # Multiple cancer positive samples (r = 0.1)
#' simulate_droplet_counts(1, 0.1, 0.01, 0.01, 0.01, 14000, 10)
#'
#'
#' # Simulate molecule counts for a cancer positive samples (r = 0.1)
#' simulate_molecule_counts(1, 0.1, 0.01, 0.01, 0.01, 14000)

simulate_droplet_counts <- function(l, r, a, b, c, n_drops, n_samples = 1) {
  # Get the probability for each type of droplet
  p_neg <- P_d_neg(l = l, a = a, c = c, r = r)
  p_pos <- P_d_pos(l = l, a = a, b = b, c = c, r = r)
  p_m <- P_M_only(l = l, a = a, c = c, r = r)
  p_wt <- P_WT_only(l = l, a = a, b = b, c = c, r = r)

  # Simulate droplet counts using multinomial distribution
  sim <- stats::rmultinom(n_samples, size = n_drops, prob = c(p_neg, p_pos, p_wt, p_m))

  # Collect data and output
  df_sim <- data.frame(
    N_d_neg = sim[1, ],
    N_d_pos = sim[2, ],
    N_WT_only = sim[3, ],
    N_M_only = sim[4, ]
  )

  return(df_sim)
}

#' @rdname simulate_droplet_counts
#' @export
#' @importFrom rlang .data
simulate_molecule_counts <- function(l, r, a, b, c, n_drops) {
  # Simulate individual molecule counts for each droplet
  WT_copies <- stats::rpois(n_drops, l)
  M_copies <- stats::rpois(n_drops, r + a + b * WT_copies + c * l)

  # Put data in df and add TRUE/FALSE signal
  indiv_droplet_data <- data.frame(
    # Simulate number of (start) molecules in droplets
    WT_copies = WT_copies,
    M_copies = M_copies
  ) %>%
    # Convert to signal (channel) data from ddPCR
    dplyr::mutate(
      WT_ch = WT_copies > 0,
      M_ch = M_copies > 0
    )

  # Summarize data
  summarized_data <- indiv_droplet_data %>%
    dplyr::summarise(
      N_d_neg = sum(!.data$WT_ch & !.data$M_ch),
      N_d_pos = sum(.data$WT_ch & .data$M_ch),
      N_WT_only = sum(.data$WT_ch & !.data$M_ch),
      N_M_only = sum(!.data$WT_ch & .data$M_ch)
    )

  return(list(
    indiv_droplet_data = indiv_droplet_data,
    summarized_data = summarized_data
  ))
}
