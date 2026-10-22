test_that("latent decorrelation removes nothing without spatial structure", {
  set.seed(1)
  xy <- matrix(stats::runif(400), 200)
  X <- matrix(stats::rnorm(200 * 4), 200)
  for (m in c("sbss", "maf")) {
    d <- decorrelate_latent(X, xy, method = m, nsim = 99, seed = 2)
    expect_s3_class(d, "depprep_decorrelated")
    expect_lte(d$removed, 1)                             # at most a false positive at 5%
    expect_equal(dim(d$X), dim(X))
    expect_equal(unname(diag(d$R)), rep(1, 4))
  }
})

test_that("a strong shared field is removed by both methods", {
  g <- sim_groups(sim_group("a", 4, "matern", range = 0.5, large_scale_share = 0.8))
  s <- sim_spatial(300, g, seed = 3)
  L <- s$components$large[, 1]
  for (m in c("sbss", "maf")) {
    d <- decorrelate_latent(s$X, s$coords, method = m, nsim = 99, seed = 4)
    expect_gte(d$removed, 1)
    # the shared field is no longer in the data
    expect_lt(max(abs(stats::cor(d$X, L))), 0.3)
    expect_lt(mean(abs(d$R[upper.tri(d$R)])), mean(abs(stats::cor(s$X)[upper.tri(d$R)])))
  }
  expect_identical(decorrelate_latent(s$X, s$coords, "maf", nsim = 99, seed = 4),
                   decorrelate_latent(s$X, s$coords, "maf", nsim = 99, seed = 4))
})
