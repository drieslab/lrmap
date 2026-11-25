#' Per-cell ligand exposure and local score for one ligand-receptor pair
#'
#' The per-cell quantities behind one row of an [lrmap()] table, for maps: the ligand of each sender cell, the
#' receptor of each receiver cell, the ligand each receiver is exposed to from sender cells at the chosen scale
#' (X_i, with the kernels of [lrmap()]), and the local score R_i * X_i. Averaged over the receivers, the local
#' score is the table's `score` (before leak removal).
#'
#' @param counts raw counts: a genes x cells matrix with gene rownames, or a reader `function(genes)` (see
#'   [lrmap()]); only the ligand and receptor genes are read.
#' @param coords cells x 2 coordinates (microns).
#' @param types cell type per cell.
#' @param sender,receiver sender and receiver cell types.
#' @param ligand,receptor genes of the ligand and the receptor; several genes for a complex (ligand: minimum over
#'   subunits, receptor: mean).
#' @param scale "contact", "paracrine" or "niche".
#' @param lib library size per cell (all genes); default: column sums of `counts` when it is a matrix.
#' @param k,sigma_p,sigma_n,bin_p,bin_n kernel settings, as in [lrmap()].
#' @return data.frame, one row per cell (in the order of `coords`): `type`, `ligand` (sender cells, else NA),
#'   `receptor` (receiver cells, else NA), `exposure` (receiver cells), `local` (= receptor x exposure).
#' @export
lrmap_local <- function(counts, coords, types, sender, receiver, ligand, receptor,
                        scale = c("contact", "paracrine", "niche"), lib = NULL, k = 8, sigma_p = 25, sigma_n = 150,
                        bin_p = 8, bin_n = 25) {
  scale <- match.arg(scale); coords <- as.matrix(coords)[, 1:2, drop = FALSE]; n <- nrow(coords)
  types <- as.character(types); if (length(types) != n) stop("coords and types must describe the same cells")
  g <- unique(c(ligand, receptor))
  X <- if (is.function(counts)) counts(g) else counts[intersect(g, rownames(counts)), , drop = FALSE]
  if (!all(g %in% rownames(X))) stop("genes not found in counts: ", paste(setdiff(g, rownames(X)), collapse = ", "))
  if (is.null(lib)) {
    if (is.function(counts)) stop("lib is required when counts is a reader function")
    lib <- Matrix::colSums(counts)
  }
  Y <- as.matrix(Matrix::t(X[g, , drop = FALSE])) * (100 / pmax(as.numeric(lib), 1))
  L <- if (length(ligand) == 1) Y[, ligand] else do.call(pmin, as.data.frame(Y[, ligand, drop = FALSE]))
  R <- if (length(receptor) == 1) Y[, receptor] else rowMeans(Y[, receptor, drop = FALSE])
  isA <- types == sender; isB <- types == receiver
  if (scale == "contact") {
    nn <- RANN::nn2(coords, k = k + 1); d0 <- stats::median(as.vector(nn$nn.dists[, -1])); rm(nn)
    e <- .contact_edges(coords, k, d0, 3 * d0)
    xa <- .rsum(e$wn * (L * isA)[e$send], e$recv, n)
  } else {
    sg <- if (scale == "paracrine") sigma_p else sigma_n; bn <- if (scale == "paracrine") bin_p else bin_n
    box <- c(min(coords[, 1]), diff(range(coords[, 1])) + 1e-6, min(coords[, 2]), diff(range(coords[, 2])) + 1e-6)
    gr <- .grid(coords, box, sg, bn); nc <- gr$nx * gr$ny
    Bl <- .blur(cbind(.bin_sum(cbind(L * isA), gr$lin, nc), .bin_sum(matrix(1, n, 1), gr$lin, nc)), gr)
    Tf <- Bl[, 2]; F <- Bl[, 1] / pmax(Tf, EPS); F[Tf < 0.1 * stats::median(Tf[gr$lin])] <- 0
    xa <- F[gr$lin]
  }
  data.frame(type = types, ligand = ifelse(isA, L, NA_real_), receptor = ifelse(isB, R, NA_real_),
             exposure = ifelse(isB, xa, NA_real_), local = ifelse(isB, R * xa, NA_real_))
}
