test_that("simulation works", {
  # Parameters
  a = 0.01
  b = 0.01
  c = 0.01

  l = 1
  r = 0.1

  n_droplets = 19

  # Droplet counts
  droplet_counts = simulate_droplet_counts(l = l, r = r, a = a, b = b, c = c, n_drops = n_droplets)
  expect_equal(sum(droplet_counts), n_droplets)

  # Multiple samples
  n_samples = 7
  multisample_droplet_counts = simulate_droplet_counts(l = l, r = r, a = a, b = b, c = c, n_drops = n_droplets, n_samples = n_samples)
  multisample_droplet_counts

  expect_equal(nrow(multisample_droplet_counts), n_samples)
  expect_equal(ncol(multisample_droplet_counts), 4)
  expect_equal(rowSums(multisample_droplet_counts), rep(n_droplets, n_samples))

  # Molecule counts
  molecule_count_res = simulate_molecule_counts(l = l, r = r, a = a, b = b, c = c, n_drops = n_droplets)
  expect_equal(length(molecule_count_res),2)
  expect_equal(nrow(molecule_count_res$indiv_droplet_data),n_droplets)
  expect_equal(ncol(molecule_count_res$indiv_droplet_data),4)

  expect_equal(sum(molecule_count_res$summarized_data),n_droplets)

})
