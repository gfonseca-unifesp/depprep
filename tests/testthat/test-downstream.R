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

test_that("constrained clustering recovers contiguous classes in space, time and space-time", {
  set.seed(4)
  xy <- matrix(stats::runif(150 * 2), 150)
  cls <- cut(xy[, 1], c(0, 1 / 3, 2 / 3, 1), labels = FALSE)  # three vertical bands
  X <- matrix(stats::rnorm(150 * 4), 150) + 2 * cls
  cs <- run_downstream(X, "constrained", truth = list(classes = cls, position = xy), seed = 1)
  expect_identical(cs$metric, c("k", "ari_true_k", "rand_kmeans"))
  expect_equal(cs$value[cs$metric == "k"], 3)
  expect_gt(cs$value[cs$metric == "ari_true_k"], 0.8)
  expect_true(cs$value[cs$metric == "rand_kmeans"] >= 0 && cs$value[cs$metric == "rand_kmeans"] <= 1)
  expect_identical(run_downstream(X, "constrained", truth = list(classes = cls, position = xy), seed = 1), cs)

  tt <- sample(200)                                            # rows out of time order
  cls_t <- cut(tt, c(0, 70, 130, 200), labels = FALSE)         # three periods
  Xt <- matrix(stats::rnorm(200 * 3), 200) + 1.5 * cls_t
  ct <- run_downstream(Xt, "constrained", truth = list(classes = cls_t, position = matrix(tt)), seed = 2)
  expect_gt(ct$value[ct$metric == "ari_true_k"], 0.8)

  sites <- matrix(stats::runif(30 * 2), 30)
  pos <- cbind(sites[rep(1:30, 8), ], rep(1:8, each = 30))
  cls_st <- ifelse(pos[, 1] < 0.5, 1, 2)                       # two regions, fixed in time
  Xst <- matrix(stats::rnorm(240 * 3), 240) + 2 * cls_st
  cst <- run_downstream(Xst, "constrained", truth = list(classes = cls_st, position = pos), seed = 3)
  expect_gt(cst$value[cst$metric == "ari_true_k"], 0.8)

  noclass <- run_downstream(X, "constrained", truth = list(position = xy), seed = 1)
  expect_equal(noclass$value[noclass$metric == "k"], 3)
  expect_true(is.na(noclass$value[noclass$metric == "ari_true_k"]))
  expect_error(run_downstream(X, "constrained", truth = list(classes = cls)), "position")
})

test_that("contiguity links connect the graph and the Rand index is correct", {
  set.seed(5)
  xy <- rbind(matrix(stats::runif(20 * 2), 20), matrix(stats::runif(20 * 2) + 10, 20))  # two far clusters
  e <- depprep:::.contiguity_links(xy)
  expect_equal(depprep:::.n_components(40, e), 1)
  expect_true(all(e[, 1] < e[, 2]))
  et <- depprep:::.contiguity_links(matrix(c(3, 1, 2)))
  expect_equal(et[order(et[, 1]), , drop = FALSE], rbind(c(1, 3), c(2, 3)))
  pos <- cbind(rep(c(0, 1), 3), 0, rep(1:3, each = 2))         # 2 sites x 3 times
  est <- depprep:::.contiguity_links(pos)
  expect_equal(nrow(est), 3 + 4)                                 # 1 space link per time + 2 time links per site
  expect_equal(depprep:::.rand_index(c(1, 1, 2, 2), c(2, 2, 1, 1)), 1)
  expect_equal(depprep:::.rand_index(c(1, 1, 2, 2), c(1, 2, 1, 2)), 2 / 6)
})

test_that("PCA is scored against the true signal and large data sets are subsampled", {
  set.seed(5)
  n <- 400; cls <- rep(1:2, each = n / 2)
  sig <- cbind(matrix(stats::rnorm(n * 4), n) + 2 * cls, matrix(stats::rnorm(n * 4), n))
  X <- sig + matrix(stats::rnorm(n * 8, sd = 0.3), n)
  pc <- run_downstream(X, "pca", truth = list(signal = sig))
  d <- pc$value[pc$metric == "subspace_distance_signal"]
  expect_lt(d, 0.2)                                        # the leading axis is the class axis of the signal
  km <- run_downstream(X, "kmeans", truth = list(classes = cls), max_n = 100, seed = 2)
  expect_equal(km$value[km$metric == "subsample"], 0.25)
  expect_gt(km$value[km$metric == "ari_true_k"], 0.7)
  expect_identical(run_downstream(X, "kmeans", truth = list(classes = cls), max_n = 100, seed = 2), km)
  expect_false("subsample" %in% run_downstream(X, "kmeans", truth = list(classes = cls), seed = 2)$metric)
})


test_that("HDBSCAN finds well separated classes and no groups in noise", {
  set.seed(8)
  cls <- rep(1:3, each = 80)
  X <- matrix(stats::rnorm(240 * 4), 240) + 4 * cbind(cls == 1, cls == 2, cls == 3, 0)
  h <- run_downstream(X, "hdbscan", truth = list(classes = cls))
  v <- stats::setNames(h$value, h$metric)
  expect_equal(unname(v["n_clusters"]), 3)
  expect_gt(v["ari_clustered"], 0.95)                      # the clusters are pure
  expect_true(v["noise_share"] > 0 && v["noise_share"] < 0.5)
  expect_gt(v["ari"], 0.3)                                 # noise counted as a group lowers the ARI
  expect_equal(unname(v["min_points"]), 10)
  W <- matrix(stats::rnorm(400 * 5), 400)
  h0 <- run_downstream(W, "hdbscan", seed = 1)
  v0 <- stats::setNames(h0$value, h0$metric)
  expect_equal(unname(v0["clusters_above_1"]), 0)
  expect_true(is.na(v0["ari"]))
  expect_identical(run_downstream(W, "hdbscan", seed = 1), h0)
})
