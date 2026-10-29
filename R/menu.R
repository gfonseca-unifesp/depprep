#' Menu of large-scale structures for the bases
#'
#' A menu restricts or imposes the structures a basis may be built from (decision D16).
#' The user chooses the menu; [build_group_basis()] searches the parameters within it,
#' with the same limits as the automatic choice: a global permutation test, forward
#' selection with permutation tests and double stopping, the cap on `k / n`, and only the
#' group's own variables (marginal information). Every term is a function of the
#' positions only (coordinates or times), so it can be the shared basis of a group;
#' terms that are functions of the variable itself (lags, the mean of the neighbours) are
#' filters, not basis terms, and are not in the menu.
#'
#' Terms, in order of precedence (earlier terms are kept when later ones repeat them):
#' * `term_trend(degree)`: polynomial of time (default degree 3) or trend surface of
#'   the centred coordinates (default degree 2: x, y, x^2, xy, y^2);
#' * `term_cycle(period, harmonics)`: sine and cosine harmonics of a period. Without
#'   `period`, it is the median of the significant peaks of the least-squares
#'   periodogram of the group's variables;
#' * `term_mem(W, margin)`: Moran's eigenvector maps with positive autocorrelation.
#'   `W = "dbmem"` gives the distance-based maps used by default. In space, a
#'   neighbourhood graph and weights can be named (`"del"`, `"gab"`, `"rel"`, `"mst"`,
#'   `"dnear"`, as in [adespatial::listw.candidates()], with `"binary"` or `"flin"`
#'   weights, for example `"gab_flin"`), or `W = "auto"` tries the Delaunay, Gabriel and
#'   relative-neighbourhood graphs with binary and linear weights and keeps the one whose
#'   selection has the largest adjusted R-squared, with the global test at the Sidak
#'   level for the number of graphs tried (Bauman et al. 2018). `margin = "time"` gives
#'   the temporal maps;
#' * `term_aem(direction)`: asymmetric eigenvector maps (Blanchet et al. 2008) on the
#'   Delaunay graph oriented by `direction` (degrees, counter-clockwise from the x
#'   axis), for directional processes such as a current or a plume.
#'
#' In space-time the spatial terms (`term_trend()`, `term_mem()`, `term_aem()`) and the
#' temporal ones (`term_trend()`, `term_cycle()`, `term_mem(margin = "time")`) form the
#' two margins, selected first on the stacked data, and their products are the candidates
#' (as for the automatic choice). A margin exists when a term other than the trend
#' targets it, or when `term_trend(margin = )` names it.
#'
#' @param ... Terms.
#' @param degree Degree of the polynomial or trend surface; `NULL` for the default.
#' @param period Period of the cycle, in time units.
#' @param harmonics Number of harmonics; `NULL` for `floor(period / 2)`.
#' @param W Spatial weights of the maps (see above).
#' @param margin `"space"` or `"time"`; `NULL` for the margin of the data (space when
#'   there are coordinates).
#' @param direction Direction of the process, in degrees.
#' @return `basis_menu()` returns a list of terms of class `depprep_menu`; each `term_*()`
#'   returns one term (class `depprep_term`).
#' @references Bauman, D., Drouet, T., Fortin, M.-J. and Dray, S. (2018). Optimizing the
#'   choice of a spatial weighting matrix in eigenvector-based methods. Ecology 99,
#'   2159-2166.
#'
#'   Blanchet, F. G., Legendre, P. and Borcard, D. (2008). Modelling directional spatial
#'   processes in ecological data. Ecological Modelling 215, 325-336.
#' @examples
#' basis_menu(term_trend(), term_cycle(period = 12))
#' basis_menu(term_mem(W = "auto"))
#' basis_menu(term_trend(), term_aem(direction = 45))
#' @export
basis_menu <- function(...) {
  terms <- list(...)
  if (!length(terms) || !all(vapply(terms, inherits, TRUE, "depprep_term"))) {
    stop("`basis_menu()` takes one or more terms made with term_trend(), term_cycle(), term_mem() or term_aem().", call. = FALSE)
  }
  structure(terms, class = "depprep_menu")
}

#' @rdname basis_menu
#' @export
term_trend <- function(degree = NULL, margin = NULL) .term("trend", margin, degree = degree)

#' @rdname basis_menu
#' @export
term_cycle <- function(period = NULL, harmonics = NULL) {
  if (!is.null(period) && (!is.numeric(period) || length(period) != 1 || period < 2)) stop("`period` must be one number of at least 2.", call. = FALSE)
  .term("cycle", "time", period = period, harmonics = harmonics)
}

#' @rdname basis_menu
#' @export
term_mem <- function(W = "dbmem", margin = NULL) {
  W <- as.character(W)
  ok <- W %in% c("dbmem", "auto") || grepl("^(del|gab|rel|mst|dnear)(_(binary|flin))?$", W)
  if (length(W) != 1 || !ok) stop("Unknown `W`: ", paste(W, collapse = ", "), call. = FALSE)
  .term("mem", margin, W = W)
}

#' @rdname basis_menu
#' @export
term_aem <- function(direction) {
  if (!is.numeric(direction) || length(direction) != 1) stop("`direction` must be one angle in degrees.", call. = FALSE)
  .term("aem", "space", direction = direction)
}

.term <- function(kind, margin, ...) {
  if (!is.null(margin)) margin <- match.arg(margin, c("space", "time"))
  structure(list(kind = kind, margin = margin, ...), class = "depprep_term")
}

#' @export
print.depprep_menu <- function(x, ...) {
  cat("<depprep_menu>\n")
  for (t in x) {
    par <- t[setdiff(names(t), c("kind", "margin"))]
    par <- par[!vapply(par, is.null, TRUE)]
    cat(sprintf("  %s%s%s\n", t$kind, if (length(par)) paste0("(", paste(names(par), unlist(par), sep = " = ", collapse = ", "), ")") else "",
                if (!is.null(t$margin)) paste0(" [", t$margin, "]") else ""))
  }
  invisible(x)
}

# The menu that reproduces the automatic choice of each regime type (build_group_basis()).
.menu_of_type <- function(type, period = NULL, trend = TRUE) {
  tr <- if (trend) list(term_trend()) else list()
  terms <- switch(type,
    isotropic = , anisotropic = c(tr, list(term_mem())),
    seasonal = c(tr, list(term_cycle(period))),
    trend = c(tr, list(term_mem(margin = "time"))),
    ar = list(term_mem(margin = "time")),
    temporal_common = c(tr, list(term_mem(margin = "time"))),
    spatial_only = c(tr, list(term_mem(margin = "space"))),
    separable_seasonal = c(tr, list(term_mem(margin = "space"), term_cycle(period))),
    moving = , separable_ar1 = c(tr, list(term_mem(margin = "space"), term_mem(margin = "time"))))
  do.call(basis_menu, terms)
}

# The margin of each term for a data family: terms without one take the margin of the data.
.term_margins <- function(menu, family) {
  vapply(menu, function(t) if (!is.null(t$margin)) t$margin else if (family == "temporal") "time" else
    if (t$kind == "trend") "both" else "space", "")
}

# Candidate sets of one margin, as a list of alternatives (more than one only with
# term_mem(W = "auto")): the terms' matrices in menu order, made orthonormal together.
# `pos` are the coordinates (space) or the sorted times (time); `Y` the margin's
# response (columns = variables, or stacked variables in space-time).
.margin_candidates <- function(menu, margins, which, pos, Y) {
  use <- margins %in% c(which, "both")
  parts <- lapply(menu[use], function(t) .term_matrix(t, which, pos, Y))
  alt <- which(vapply(parts, is.list, TRUE))
  if (!length(alt)) return(list(.orthonormal(do.call(cbind, parts))))
  lapply(parts[[alt[1]]], function(m) { p <- parts; p[[alt[1]]] <- m
    .orthonormal(do.call(cbind, lapply(p, function(x) if (is.list(x)) x[[1]] else x))) })
}

# Raw matrix of one term on one margin (a named list of matrices for W = "auto").
.term_matrix <- function(t, which, pos, Y) {
  switch(t$kind,
    trend = if (which == "time") {
      tt <- as.numeric(pos)
      degree <- min(if (is.null(t$degree)) 3 else t$degree, length(unique(tt)) - 1)
      if (degree < 1) return(matrix(0, length(tt), 0))
      p <- stats::poly(tt, degree); colnames(p) <- paste0("trend", seq_len(degree)); p
    } else {
      m <- .trend_surface(pos)
      if (identical(t$degree, 1)) m[, seq_len(ncol(as.matrix(pos))), drop = FALSE] else m
    },
    cycle = {
      tt <- as.numeric(pos)
      period <- if (is.null(t$period)) .menu_period(Y, tt) else t$period
      h <- seq_len(if (is.null(t$harmonics)) floor(period / 2) else min(t$harmonics, floor(period / 2)))
      hm <- do.call(cbind, lapply(h, function(i) cbind(sin(2 * pi * i * tt / period), cos(2 * pi * i * tt / period))))
      colnames(hm) <- paste0(rep(c("sin", "cos"), length(h)), rep(h, each = 2), "_p", signif(period, 6))
      hm
    },
    mem = if (which == "time") {
      m <- as.matrix(adespatial::dbmem(matrix(as.numeric(pos)), MEM.autocor = "positive", silent = TRUE))
      colnames(m) <- paste0("tMEM", seq_len(ncol(m))); m
    } else if (t$W == "dbmem") {
      as.matrix(adespatial::dbmem(as.matrix(pos), MEM.autocor = "positive", silent = TRUE))
    } else {
      .mem_of_graphs(as.matrix(pos), if (t$W == "auto") c("del_binary", "del_flin", "gab_binary", "gab_flin", "rel_binary", "rel_flin") else t$W)
    },
    aem = .aem_vectors(as.matrix(pos), t$direction))
}

# Positive Moran's eigenvector maps of named neighbourhood graphs and weights (a matrix
# for one graph, a named list for several).
.mem_of_graphs <- function(coords, codes) {
  codes <- ifelse(grepl("_", codes), codes, paste0(codes, "_binary"))
  out <- lapply(codes, function(code) {
    p <- strsplit(code, "_")[[1]]
    lw <- adespatial::listw.candidates(coords, nb = p[1], weights = p[2])[[1]]
    m <- as.matrix(adespatial::scores.listw(lw, MEM.autocor = "positive"))
    colnames(m) <- paste0("MEM", seq_len(ncol(m)), "_", code)
    m
  })
  if (length(out) == 1) out[[1]] else stats::setNames(out, codes)
}

# Asymmetric eigenvector maps on the Delaunay graph oriented by `direction` (degrees).
# Only the maps with positive autocorrelation on the graph (Moran's I above its
# expectation) are kept, as for the dbMEM: all n - 1 maps span almost every pattern, and
# the global test of the selection would have no power.
.aem_vectors <- function(coords, direction) {
  # building the oriented graph is slow (about a minute at n = 200), and correct() asks for
  # the same maps for the group basis and for each variable's own basis: keep the last few
  coords <- as.matrix(coords)
  for (e in .aem_cache$entries) if (e$direction == direction && identical(e$coords, coords)) return(e$vectors)
  m <- .aem_build(coords, direction)
  .aem_cache$entries <- c(list(list(coords = coords, direction = direction, vectors = m)), utils::head(.aem_cache$entries, 3))
  m
}
.aem_cache <- new.env(parent = emptyenv())

.aem_build <- function(coords, direction) {
  nb <- spdep::tri2nb(coords)
  b <- adespatial::aem.build.binary(nb.object = nb, coords = cbind(seq_len(nrow(coords)), coords),
                                    rot.angle = direction, plot.connexions = FALSE)
  a <- adespatial::aem(aem.build.binary = b)
  m <- a$vectors[, a$values > 1e-8, drop = FALSE]
  W <- spdep::nb2mat(nb, style = "B")
  mc <- scale(m, scale = FALSE)
  moran <- nrow(m) / sum(W) * colSums(mc * (W %*% mc)) / colSums(mc^2)
  m <- m[, moran > -1 / (nrow(m) - 1), drop = FALSE]
  colnames(m) <- paste0("AEM", seq_len(ncol(m)), "_d", direction)
  m
}

# Period of term_cycle() without `period`: the median of the significant peaks (p < 0.05)
# of the least-squares periodogram of the margin's columns.
.menu_period <- function(Y, tt) {
  Y <- as.matrix(Y)
  pk <- lapply(seq_len(ncol(Y)), function(j) .ls_periodogram(Y[, j] - mean(Y[, j]), tt))
  per <- vapply(pk, function(p) if (is.finite(p$p) && p$p < 0.05) p$period else NA_real_, 0)
  if (!any(is.finite(per))) stop("term_cycle(): no significant cycle in the group; give `period`.", call. = FALSE)
  stats::median(per, na.rm = TRUE)
}

# Forward selection over alternative candidate sets: one set is the usual selection;
# several (W = "auto") are each selected with the global test at the Sidak level for the
# number of sets, and the selection with the largest adjusted R-squared is kept.
.select_alternatives <- function(Y, alts, kn_max, alpha, nperm) {
  if (length(alts) == 1) return(.forward_select(Y, alts[[1]], kn_max, alpha, nperm))
  a <- 1 - (1 - alpha)^(1 / length(alts))
  if (1 / (nperm + 1) > a) warning(sprintf("%d permutations cannot reach the Sidak level %.4f for %d graphs; use nperm >= %d.",
                                           nperm, a, length(alts), ceiling(1 / a)), call. = FALSE)
  sels <- lapply(alts, function(cand) .forward_select(Y, cand, kn_max, a, nperm))
  n <- nrow(Y)
  adj <- vapply(sels, function(s) { k <- length(s$terms); if (!k) -Inf else 1 - (1 - s$r2) * (n - 1) / (n - k - 1) }, 0)
  best <- if (all(!is.finite(adj))) 1 else which.max(adj)
  s <- sels[[best]]
  s$W <- if (all(!is.finite(adj))) NA_character_ else names(alts)[best]
  s
}
