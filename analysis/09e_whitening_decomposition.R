# Stage 9 — where the pairs live and which step of the whitening removes them (protocol-v1.3.1,
# section 17.2; descriptive). Deterministic: no random numbers.
#
# On the cube with valid fraction 0.25 and the grouping of the main run, for the three
# predicted pairs and three others:
# * by scale, without a model: raw; site means; basin-mean seasonal cycle; seasonal cycle
#   within cells; anomalies from each cell's monthly climatology; the same anomalies as basin
#   means; the same anomalies with the basin mean of each month removed;
# * by step of the whitening: spatial only, AR only and both, each with each variable's own
#   filter and with the filter of either variable on both. The AR filter alone still carries
#   the spatial gradient.
#
# Requires results/fza/ocean_cube.csv, ocean_cells.csv and ocean_groups.csv.
# Output (not tracked): results/fza/ocean_whitening_decomposition.csv

suppressMessages(devtools::load_all(quiet = TRUE))
out_dir <- file.path("results", "fza")
variables <- c("sst_night", "salinity", "precipitation", "wind_u", "wind_v", "par", "kd490", "adg443", "chlorophyll")
log_variables <- c("chlorophyll", "kd490", "adg443")
pairs <- list(c("salinity", "adg443"), c("par", "precipitation"), c("kd490", "chlorophyll"),
              c("adg443", "chlorophyll"), c("salinity", "kd490"), c("sst_night", "precipitation"))

cube <- read.csv(file.path(out_dir, "ocean_cube.csv"), stringsAsFactors = FALSE)
cells <- read.csv(file.path(out_dir, "ocean_cells.csv"))
months <- sort(unique(cube$month))
ok <- !is.na(cube$value) & cube$valid_pixels > 0 & cube$valid_fraction >= 0.25
v <- ifelse(ok, cube$value, NA_real_)
lg <- cube$variable %in% log_variables
v[lg] <- log10(v[lg])
ns <- nrow(cells); nt <- length(months)
X <- matrix(NA_real_, ns * nt, length(variables), dimnames = list(NULL, variables))
X[cbind((match(cube$month, months) - 1) * ns + match(cube$cell, cells$cell), match(cube$variable, variables))] <- v
index <- data.frame(site = rep(seq_len(ns), nt), time = rep(seq_len(nt), each = ns))
coords <- as.matrix(cells[, c("lon", "lat")])
calendar <- as.integer(substr(months, 6, 7))

r_of <- function(a, b) stats::cor(as.vector(a), as.vector(b), use = "pairwise.complete.obs")
to_cube <- function(z) matrix(z, ns, nt)
climatology <- function(M) {
  C <- M
  for (m in 1:12) C[, calendar == m] <- rowMeans(M[, calendar == m, drop = FALSE], na.rm = TRUE)
  C
}
by_scale <- do.call(rbind, lapply(pairs, function(p) {
  A <- to_cube(X[, p[1]]); B <- to_cube(X[, p[2]])
  aA <- A - climatology(A); aB <- B - climatology(B)
  site <- function(M) rowMeans(M, na.rm = TRUE)
  basin <- function(M) colMeans(M, na.rm = TRUE)
  cycle <- function(M) tapply(basin(M), calendar, mean)
  data.frame(a = p[1], b = p[2], step = c("raw", "site means", "basin seasonal cycle", "seasonal cycle within cells",
                                          "anomalies", "anomalies, basin means", "anomalies, local"),
             r = c(r_of(A, B), r_of(site(A), site(B)), r_of(cycle(A), cycle(B)),
                   r_of(climatology(A) - site(A), climatology(B) - site(B)), r_of(aA, aB), r_of(basin(aA), basin(aB)),
                   r_of(sweep(aA, 2, basin(aA)), sweep(aB, 2, basin(aB)))))
}))

groups <- read.csv(file.path(out_dir, "ocean_groups.csv"))
groups <- groups[groups$min_fraction == 0.25, ]
group_of <- stats::setNames(groups$group, groups$variable)
W <- stats::setNames(lapply(unique(groups$group), function(g) {
  cols <- groups$variable[groups$group == g]
  .whitener(X[, cols, drop = FALSE], g, groups$regime[groups$group == g][1], NULL, "spatiotemporal", coords, NULL, index)
}), unique(groups$group))
spatial <- function(z, g) { e <- environment(W[[g]]); e$spatial_whiten(e$sw, e$to_cube(z)) }
temporal <- function(C, g) t(.ar_filter(t(C), environment(W[[g]])$phi))
last_times <- function(M, m) M[, utils::tail(seq_len(ncol(M)), m), drop = FALSE]
r_aligned <- function(A, B) { m <- min(ncol(A), ncol(B)); r_of(last_times(A, m), last_times(B, m)) }
Z <- scale(X)
by_step <- do.call(rbind, lapply(pairs, function(p) {
  a <- Z[, p[1]]; b <- Z[, p[2]]
  ga <- group_of[[p[1]]]; gb <- group_of[[p[2]]]
  sa <- spatial(a, ga); sb <- spatial(b, gb)
  data.frame(a = p[1], b = p[2], step = c("spatial, own filters", "AR, own filters", "both, own filters",
                                          "spatial, filter of a", "spatial, filter of b", "both, filter of a", "both, filter of b"),
             r = c(r_of(sa, sb), r_aligned(temporal(to_cube(a), ga), temporal(to_cube(b), gb)),
                   r_aligned(temporal(sa, ga), temporal(sb, gb)), r_of(sa, spatial(b, ga)), r_of(spatial(a, gb), sb),
                   r_aligned(temporal(sa, ga), temporal(spatial(b, ga), ga)), r_aligned(temporal(spatial(a, gb), gb), temporal(sb, gb))))
}))
out <- rbind(cbind(part = "scale", by_scale), cbind(part = "whitening", by_step))
utils::write.csv(out, file.path(out_dir, "ocean_whitening_decomposition.csv"), row.names = FALSE)
print(stats::xtabs(r ~ step + paste(a, b, sep = "-"), out), digits = 2)
