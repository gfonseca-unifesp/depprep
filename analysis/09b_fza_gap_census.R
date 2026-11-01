# Stage 9 — census of the gaps of the ocean cube (decision D26: how to treat missing data).
#
# The package stops on missing values (stage 8.1). Before choosing what to do with them,
# this script measures how much of the cube is missing and where:
# * per variable, at several minimum valid fractions of a cell (protocol section 12,
#   PENDING): a cell-month has a value when at least that fraction of its water pixels is
#   valid (0 = at least one valid pixel);
# * by calendar month (cloud cover under the ITCZ) and by region: the shelf, where the
#   Amazon plume spreads (depth above 200 m, the median of the ETOPO1 points of the cell),
#   against the slope and deep water;
# * how concentrated the gaps are: the share of the variance of the gap indicator that the
#   months explain, and the share that the cells explain;
# * which gaps are structural: cells where a variable is never observed (the coastal
#   footprint of the 40-km salinity, the land mask of the wind);
# * what each option of D26 would keep:
#   (a) only complete data, within the cells where every variable is observed in some
#       month: cell-months with all variables, cells complete in every month, months
#       complete in every cell;
#   (b) every observed value (no row is dropped; the gaps stay as gaps);
#   (c) longer composites: a cell has a value in a two- or three-month period when it has
#       one in any month of the period.
# It also checks the cube against the feasibility run (analysis/01b): with every pixel of
# a cell, the share of cells with a value must be at least the share seen with one pixel
# per cell, in the months the two have in common.
# No random numbers are used.
#
# Requires results/fza/ocean_cube.csv and ocean_cells.csv (analysis/09a_fza_ocean_cube.R),
# results/fza/coverage_monthly.csv (analysis/01b) and data/raw/fza/etopo1_search_box.csv
# (analysis/01a).
# Outputs (not tracked):
#   - results/fza/gap_census.csv    missing share per threshold, variable and scope
#   - results/fza/gap_options.csv   what each option of D26 keeps, per threshold
#   - results/fza/gap_check.csv     the cube against the feasibility coverage

out_dir <- file.path("results", "fza")
cube <- read.csv(file.path(out_dir, "ocean_cube.csv"), stringsAsFactors = FALSE)
cells <- read.csv(file.path(out_dir, "ocean_cells.csv"))
thresholds <- c(0, 0.25, 0.5, 0.75)
g <- 0.25

# --- Region of each cell: shelf (depth above 200 m) or slope and deep water ------------
topo <- read.csv(file.path("data", "raw", "fza", "etopo1_search_box.csv"), skip = 2, header = FALSE,
                 col.names = c("lat", "lon", "altitude"))
topo$key <- paste(floor(topo$lon / g), floor(topo$lat / g))
depth <- tapply(topo$altitude, topo$key, stats::median)
cells$altitude <- as.vector(depth[paste(floor(cells$lon / g), floor(cells$lat / g))])
cells$region <- ifelse(cells$altitude > -200, "shelf", "slope and deep")
cube$region <- cells$region[match(cube$cell, cells$cell)]
cube$calendar_month <- as.integer(substr(cube$month, 6, 7))

has_value <- function(f) !is.na(cube$value) & cube$valid_pixels > 0 & cube$valid_fraction >= f
share <- function(x, by) { d <- aggregate(x ~ ., data = cbind(x = x, by), FUN = mean); names(d)[ncol(d)] <- "missing"; d }

census <- list(); options <- list()
for (f in thresholds) {
  miss <- !has_value(f)
  by_var <- share(miss, cube[, "variable", drop = FALSE]); by_var$scope <- "all"; by_var$level <- "all"
  by_cm <- share(miss, cube[, c("variable", "calendar_month")]); by_cm$scope <- "calendar month"
  by_cm$level <- as.character(by_cm$calendar_month); by_cm$calendar_month <- NULL
  by_reg <- share(miss, cube[, c("variable", "region")]); by_reg$scope <- "region"; by_reg$level <- by_reg$region; by_reg$region <- NULL
  # concentration: variance of the gap indicator explained by the month and by the cell
  conc <- do.call(rbind, lapply(split(data.frame(miss, month = cube$month, cell = cube$cell), cube$variable), function(d) {
    v <- stats::var(d$miss)
    data.frame(by_month = if (v > 0) stats::var(ave(as.numeric(d$miss), d$month)) / v else NA,
               by_cell = if (v > 0) stats::var(ave(as.numeric(d$miss), d$cell)) / v else NA)
  }))
  conc_rows <- rbind(data.frame(variable = rownames(conc), missing = conc$by_month, scope = "variance explained", level = "months"),
                     data.frame(variable = rownames(conc), missing = conc$by_cell, scope = "variance explained", level = "cells"))
  census[[length(census) + 1]] <- cbind(threshold = f, rbind(by_var, by_cm, by_reg, conc_rows)[, c("variable", "scope", "level", "missing")])

  # structural gaps: cells where a variable is never observed (the coastal footprint of
  # the 40-km salinity, the land mask of the wind) against gaps that come and go
  never <- tapply(!miss, list(cube$cell, cube$variable), any)
  structural <- data.frame(variable = colnames(never), scope = "cells never observed", level = "all",
                           missing = colSums(!never) / nrow(never))
  census[[length(census)]] <- rbind(census[[length(census)]], cbind(threshold = f, structural))
  seen <- as.integer(rownames(never)[apply(never, 1, all)])   # every variable observed some month

  # the options of D26 ((a) within the cells where every variable is observed some month)
  wide <- tapply(!miss, list(paste(cube$cell, cube$month), cube$variable), all)
  complete_row <- apply(wide, 1, all)
  cm <- do.call(rbind, strsplit(names(complete_row), " "))
  in_seen <- as.integer(cm[, 1]) %in% seen
  cell_ok <- tapply(complete_row[in_seen], cm[in_seen, 1], all); month_ok <- tapply(complete_row[in_seen], cm[in_seen, 2], all)
  plume <- complete_row[in_seen & cells$region[match(as.integer(cm[, 1]), cells$cell)] == "shelf"]
  composite <- function(k) {                    # k-month periods: any month with a value
    period <- (match(cube$month, sort(unique(cube$month))) - 1) %/% k
    m <- tapply(!miss, list(cube$variable, cube$cell, period), any)
    c(kept = sum(m), of = length(m), periods = dim(m)[3])
  }
  c2 <- composite(2); c3 <- composite(3)
  options[[length(options) + 1]] <- data.frame(threshold = f,
    option = c("cells with every variable observed some month",
               "(a) complete cell-months", "(a) complete cells", "(a) complete months", "(a) complete cell-months on the shelf",
               "(b) every observed value",
               sprintf("(c) two-month composites (%d periods)", c2[["periods"]]),
               sprintf("(c) three-month composites (%d periods)", c3[["periods"]])),
    count = c(length(seen), sum(complete_row[in_seen]), sum(cell_ok), sum(month_ok), sum(plume), sum(!miss), c2[["kept"]], c3[["kept"]]),
    of = c(nrow(never), sum(in_seen), length(cell_ok), length(month_ok), length(plume), length(miss), c2[["of"]], c3[["of"]]),
    unit = c("cells", "cell-months", "cells", "months", "cell-months", "values", "values", "values"))
  options[[length(options)]]$kept <- options[[length(options)]]$count / options[[length(options)]]$of
}
census <- do.call(rbind, census); options <- do.call(rbind, options)
write.csv(census, file.path(out_dir, "gap_census.csv"), row.names = FALSE)
write.csv(options, file.path(out_dir, "gap_options.csv"), row.names = FALSE)

# --- The cube against the feasibility coverage (one pixel per cell, ERDDAP) -------------
feas <- read.csv(file.path(out_dir, "coverage_monthly.csv"), stringsAsFactors = FALSE)
pairs <- c("chlorophyll-a" = "chlorophyll", "Kd490" = "kd490", "PAR" = "par")
check <- do.call(rbind, lapply(names(pairs), function(v) {
  a <- feas[feas$variable == v, c("month", "coverage")]
  b <- aggregate(x ~ month, data = data.frame(month = cube$month, x = has_value(0))[cube$variable == pairs[[v]], ], FUN = mean)
  d <- merge(a, b)
  data.frame(variable = pairs[[v]], months = nrow(d), feasibility_median = stats::median(d$coverage),
             cube_median = stats::median(d$x), months_cube_below = sum(d$x < d$coverage - 0.01))
}))
write.csv(check, file.path(out_dir, "gap_check.csv"), row.names = FALSE)

print(census[census$scope %in% c("all", "cells never observed"), ], digits = 2, row.names = FALSE)
print(options, digits = 2, row.names = FALSE)
print(check, digits = 3, row.names = FALSE)
cat("cells:", table(cells$region), "\n")
