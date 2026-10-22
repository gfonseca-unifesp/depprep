test_that("k-means and SOM recover well separated classes; pure noise gives one group", {
  set.seed(1)
  cls <- rep(1:3, each = 60)
  centres <- rbind(c(0, 0), c(8, 0), c(4, 7))                 # three groups, apart in every direction
  X <- matrix(stats::rnorm(180 * 2), 180) + centres[cls, ]
  km <- run_downstream(X, "kmeans", truth = list(classes = cls), seed = 1)
  expect_equal(km$value[km$metric == "k_gap"], 3)
  expect_gt(km$value[km$metric == "ari_gap"], 0.9)
  sm <- run_downstream(X, "som", truth = list(classes = cls), seed = 1)
  expect_gt(sm$value[sm$metric == "ari"], 0.8)
  expect_true(all(sm$value[sm$metric %in% c("quantisation_error", "topographic_error")] >= 0))
  noise <- matrix(stats::rnorm(200 * 2), 200)
  kn <- run_downstream(noise, "kmeans", seed = 1)
  expect_equal(kn$value[kn$metric == "k_gap"], 1)
  expect_identical(run_downstream(X, "som", truth = list(classes = cls), seed = 1), sm)
})

test_that("PCA recovers a block structure and finds no axis in noise", {
  # blocks of 4 with correlation 0.8: eigenvalues 3.4, above the broken-stick for 12 variables
  R_true <- kronecker(diag(3), matrix(0.8, 4, 4)); diag(R_true) <- 1
  set.seed(2)
  X <- matrix(stats::rnorm(500 * 12), 500) %*% chol(R_true)
  pc <- run_downstream(X, "pca", truth = list(R_true = R_true))
  expect_lt(pc$value[pc$metric == "subspace_distance"], 0.1)
  expect_equal(pc$value[pc$metric == "components_above_broken_stick"], 3)
  noise <- matrix(stats::rnorm(500 * 12), 500)
  pn <- run_downstream(noise, "pca", truth = list(R_true = diag(12)))
  expect_equal(pn$value[pn$metric == "first_above_broken_stick"], 0)
  expect_true(is.na(pn$value[pn$metric == "subspace_distance"]))
})

test_that("MDS: Procrustes correlation is 1 for the true signal and stress is small", {
  set.seed(3)
  signal <- matrix(stats::rnorm(100 * 4), 100) %*% diag(c(3, 2, 0.1, 0.1))
  md <- run_downstream(signal, "mds", truth = list(signal = signal, position = matrix(stats::runif(200), 100)))
  expect_equal(md$value[md$metric == "procrustes_truth"], 1, tolerance = 1e-8)
  expect_lt(md$value[md$metric == "procrustes_position"], 0.5)
  expect_equal(depprep:::.procrustes_r(signal[, 1:2], signal[, 1:2] %*% matrix(c(0, 1, -1, 0), 2) * 5), 1, tolerance = 1e-10)
  mg <- run_downstream(signal, "mds", distance = "gower", truth = list(signal = signal))
  expect_identical(unique(mg$method), "mds_gower")
  expect_gt(mg$value[mg$metric == "procrustes_truth"], 0.8)
})
