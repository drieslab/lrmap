# ---- leak (transcript spillover between touching cells) ------------------------------------------------------
# Standalone copy of edgeicg's leak model, so lrmap does not depend on edgeicg. For a cell i of type t and gene g
# (composition scale y = counts / library * 100):
#
#   leaked_ig = sum_{j touching i, type(j) = k != t} wn_ij lambda_tk y_jg + cmu_i a_tg,   cmu_i = sum_j wn_ij mu_tk
#
# wn: touching-edge weights normalised per receiver (kNN, exp(-d / median d), or user `spill_edges`). a_tg is the
# stage-1 intercept: the receiver type's baseline in a weighted regression of y on the cell's neighbourhood shares
# (contact, paracrine and niche share of every other type). lambda_tk, mu_tk come from a Huber regression across
# genes of the stage-1 contact coefficient of sender type k on the touching-sender profile (lambda) and the
# receiver baseline (mu). Same definitions, defaults and numerics as edgeicg (exposures.R, engine.R setup and
# pass 1, niche_stats.R .decompose); the niche share uses lrmap's grid blur (the share is a ratio, so the
# kernel normalisation cancels).

LEAK_DEF <- list(k = 8, sigma_p = 25, sigma_n = 150, bin = 10, max_k = 400, min_cells = 60L, min_det = 0.01,
                 huber = 1.345)

# touching edges: kNN (k, weight exp(-d / median d)) or user edges (recv, send, w)
.lk_edges <- function(coords, k, spill_edges = NULL) {
  if (!is.null(spill_edges)) {
    e <- as.data.frame(spill_edges)
    if (!all(c("recv", "send", "w") %in% names(e))) stop("spill_edges needs columns recv, send and w")
    return(list(recv = as.integer(e$recv), send = as.integer(e$send), w = as.numeric(e$w)))
  }
  nn <- RANN::nn2(coords, k = k + 1)
  d <- as.vector(nn$nn.dists[, -1, drop = FALSE])
  list(recv = rep(seq_len(nrow(coords)), k), send = as.vector(nn$nn.idx[, -1, drop = FALSE]), w = exp(-d / stats::median(d)))
}

.lk_share <- function(recv, send, w, tcode, K, n) {
  S <- matrix(0, n, K); s <- rowsum(w, recv + (tcode[send] - 1L) * n)
  S[as.integer(rownames(s))] <- s[, 1]
  S / pmax(rowSums(S), EPS)
}

# neighbourhood shares: contact (kNN), paracrine (Gaussian within 3 sigma, contact neighbours excluded; receivers
# in chunks to bound memory), niche (blurred per-type count grids)
.lk_exposures <- function(coords, tcode, K, o, rs = NULL) {
  n <- nrow(coords); ct <- .lk_edges(coords, o$k)
  C <- .lk_share(ct$recv, ct$send, ct$w, tcode, K, n)
  ck <- ct$recv * (n + 1) + ct$send
  P <- matrix(0, n, K)
  for (ch in split(seq_len(n), ceiling(seq_len(n) / 20000))) {
    nn <- RANN::nn2(coords, coords[ch, , drop = FALSE], k = min(o$max_k, n), searchtype = "radius", radius = 3 * o$sigma_p)
    recv <- rep(ch, ncol(nn$nn.idx)); send <- as.vector(nn$nn.idx); d <- as.vector(nn$nn.dists)
    ok <- send > 0 & send != recv; recv <- recv[ok]; send <- send[ok]; d <- d[ok]
    keep <- !((recv * (n + 1) + send) %in% ck)
    S <- matrix(0, n, K); s <- rowsum(exp(-0.5 * d[keep]^2 / o$sigma_p^2), recv[keep] + (tcode[send[keep]] - 1L) * n)
    S[as.integer(rownames(s))] <- s[, 1]; P[ch, ] <- S[ch, , drop = FALSE]
  }
  P <- P / pmax(rowSums(P), EPS)
  box <- c(min(coords[, 1]), diff(range(coords[, 1])) + 1e-6, min(coords[, 2]), diff(range(coords[, 2])) + 1e-6)
  gr <- .grid(coords, box, o$sigma_n, o$bin); nc <- gr$nx * gr$ny
  Fk <- .blur(.bin_sum(outer(tcode, seq_len(K), "==") * 1, gr$lin, nc), gr, rs)
  N <- Fk[gr$lin, , drop = FALSE]; N <- N / pmax(rowSums(N), EPS)
  list(C = C, P = P, N = N)
}

# stage-1 designs per receiver type (types with >= min_cells cells): A = [1, C, P, N of the other types] with
# constant columns dropped, weights lib / mean(lib); Ginv, M = Aw Ginv, leverage h; contact column per sender
.lk_designs <- function(ex, tcode, K, lib, min_cells) {
  cw <- lib / mean(lib); dz <- vector("list", K); pairs <- NULL
  for (t in seq_len(K)) {
    idx <- which(tcode == t); if (length(idx) < min_cells) next
    others <- setdiff(seq_len(K), t)
    D <- cbind(ex$C[idx, others, drop = FALSE], ex$P[idx, others, drop = FALSE], ex$N[idx, others, drop = FALSE])
    cols <- data.frame(scale = rep(c("contact", "paracrine", "niche"), each = length(others)), k = rep(others, 3))
    keep <- apply(D, 2, stats::var) > 1e-12; D <- D[, keep, drop = FALSE]; cols <- cols[keep, , drop = FALSE]
    A <- cbind(1, D); Aw <- A * cw[idx]; Pm <- diag(ncol(A)); Pm[1, 1] <- 0
    Ginv <- chol2inv(chol(crossprod(Aw, A) + 1e-6 * Pm)); M <- Aw %*% Ginv
    c1 <- vapply(others, function(kk) { w <- which(cols$scale == "contact" & cols$k == kk); if (length(w)) 1L + w else NA_integer_ }, 1L)
    dz[[t]] <- list(idx = idx, A = A, Aw = Aw, Ginv = Ginv, M = M, h = pmin(rowSums(M * A), 0.99), others = others, c1 = c1)
    pairs <- rbind(pairs, data.frame(receiver = t, sender = others))
  }
  list(dz = dz, pairs = pairs, fitted = which(!vapply(dz, is.null, TRUE)))
}

# stage-1 intercepts a (K x genes) for composition-scale Y (cells x genes, sparse)
.lk_intercepts <- function(ds, Y, K) {
  a <- matrix(NA_real_, K, ncol(Y), dimnames = list(NULL, colnames(Y)))
  for (t in ds$fitted) { z <- ds$dz[[t]]; a[t, ] <- as.vector(crossprod(z$Aw %*% z$Ginv[, 1], as.matrix(Y[z$idx, , drop = FALSE]))) }
  a
}

# Huber regression across genes of a pair's contact coefficients on the touching-sender profile (lambda) and the
# receiver baseline (mu); lambda >= 0 (edgeicg .decompose)
.lk_decompose <- function(beta, se, ak, at, huber = 1.345, iters = 50) {
  w0 <- 1 / (se^2 + stats::median(se^2) + EPS); Z <- cbind(ak, at)
  fit <- function(wh) {
    ww <- w0 * wh
    co <- as.vector(solve(crossprod(Z * ww, Z), crossprod(Z * ww, beta)))
    if (co[1] < 0) co <- c(0, sum(ww * at * beta) / (sum(ww * at^2) + EPS))
    co
  }
  wh <- rep(1, length(beta)); co <- fit(wh)
  for (it in seq_len(iters)) {
    r <- as.vector(beta - Z %*% co) / (se + EPS)
    s <- 1.4826 * stats::median(abs(r - stats::median(r))) + 1e-8
    u <- abs(r) / s; wh <- ifelse(u <= huber, 1, huber / u)
    new <- fit(wh)
    if (max(abs(new - co)) < 1e-8) { co <- new; break }
    co <- new
  }
  list(lam = co[1], mu = co[2])
}

#' Estimate leak (spillover) coefficients between touching cells
#'
#' Fits the leak model used by lrmap (`leak = `) and by edgeICG: for each receiver type t and touching sender type
#' k, `lambda` is the fraction of the sender's expression that appears in the receiver and `mu` a share of the
#' receiver's own baseline, both estimated across genes from a first-stage regression of expression on the
#' cell's neighbourhood composition (contact, paracrine and niche shares). Use many genes (all genes of the
#' panel); the coefficients are per type pair, not per gene.
#'
#' @param counts raw counts: genes x cells matrix with gene rownames, or a reader `function(genes)` (as in
#'   [lrmap()]); with a reader, `genes` and `lib` are required.
#' @param coords cells x 2 coordinates (microns). @param types cell type per cell.
#' @param lib library size per cell (default: column sums of `counts`).
#' @param genes genes to fit on (default: all rows of `counts`). Genes detected in fewer than `min_det` of the
#'   cells of every receiver type are dropped.
#' @param spill_edges optional touching edges (`recv`, `send`, `w`; cell indices in the order of `coords`), e.g.
#'   from cell outlines; default: k nearest neighbours with weight exp(-d / median d).
#' @param k,sigma_p,sigma_n,bin neighbourhood settings of the first-stage regression.
#' @param min_cells receiver types with fewer cells are not fitted. @param min_det gene filter (see `genes`).
#' @param block genes read and fitted at a time. @param verbose progress messages.
#' @return data.frame with `receiver`, `sender`, `lambda`, `mu` (attribute `n_genes`).
#' @export
lrmap_leak <- function(counts, coords, types, lib = NULL, genes = NULL, spill_edges = NULL, k = LEAK_DEF$k,
                       sigma_p = LEAK_DEF$sigma_p, sigma_n = LEAK_DEF$sigma_n, bin = LEAK_DEF$bin,
                       min_cells = LEAK_DEF$min_cells, min_det = LEAK_DEF$min_det, block = 300L, verbose = TRUE) {
  t0 <- proc.time()[[3]]
  say <- function(...) if (verbose) message(sprintf("[lrmap_leak %5.0f s] ", proc.time()[[3]] - t0), sprintf(...))
  coords <- as.matrix(coords)[, 1:2, drop = FALSE]; n <- nrow(coords)
  types <- as.character(types); lev <- sort(unique(types)); tcode <- match(types, lev); K <- length(lev)
  reader <- is.function(counts)
  if (reader && (is.null(genes) || is.null(lib))) stop("genes and lib are required when counts is a reader function")
  if (is.null(genes)) genes <- rownames(counts)
  if (is.null(lib)) lib <- Matrix::colSums(counts)
  lib <- pmax(as.numeric(lib), 1)
  o <- utils::modifyList(LEAK_DEF, list(k = k, sigma_p = sigma_p, sigma_n = sigma_n, bin = bin))
  ex <- .lk_exposures(coords, tcode, K, o); ds <- .lk_designs(ex, tcode, K, lib, min_cells); rm(ex)
  pairs <- ds$pairs; P_n <- nrow(pairs); fitted <- ds$fitted
  ct <- .lk_edges(coords, k, spill_edges)
  wn <- ct$w / pmax(.rsum(ct$w, ct$recv, n)[ct$recv], EPS); tr <- tcode[ct$recv]; ts <- tcode[ct$send]
  ep <- match(paste(tr, ts), paste(pairs$receiver, pairs$sender)); ue <- tr != ts & !is.na(ep)
  sw_p <- .rsum(wn[ue], ep[ue], P_n)
  E <- Matrix::sparseMatrix(i = ep[ue], j = ct$send[ue], x = wn[ue] / sw_p[ep[ue]], dims = c(P_n, n))
  nt <- tabulate(tcode, K); Ty <- Matrix::sparseMatrix(i = tcode, j = seq_len(n), x = 1, dims = c(K, n))
  say("%d cells, %d types (%d fitted), %d genes", n, K, length(fitted), length(genes))
  b1 <- se1 <- mbar <- a <- NULL
  for (gb in split(genes, ceiling(seq_along(genes) / block))) {
    X <- if (reader) counts(gb) else counts[intersect(gb, rownames(counts)), , drop = FALSE]
    X <- methods::as(X, "CsparseMatrix"); if (!nrow(X)) next
    Y <- Matrix::Diagonal(x = 100 / lib) %*% Matrix::t(X)
    colnames(Y) <- rownames(X); rm(X)
    det <- as.matrix(Ty %*% (Y > 0)) / nt
    keep <- colSums(det[fitted, , drop = FALSE] >= min_det) > 0
    if (!any(keep)) next
    Y <- Y[, keep, drop = FALSE]; G <- ncol(Y)
    ab <- matrix(NA_real_, K, G); bb <- sb <- matrix(NA_real_, P_n, G)
    for (t in fitted) {
      z <- ds$dz[[t]]; Yt <- as.matrix(Y[z$idx, , drop = FALSE])
      co <- z$Ginv %*% crossprod(z$Aw, Yt); R <- (Yt - z$A %*% co) / (1 - z$h)
      se2 <- crossprod(z$M^2, R^2); ab[t, ] <- co[1, ]
      pp <- match(paste(t, z$others), paste(pairs$receiver, pairs$sender))
      ok <- !is.na(z$c1); bb[pp[ok], ] <- co[z$c1[ok], ]; sb[pp[ok], ] <- sqrt(se2[z$c1[ok], ])
    }
    a <- cbind(a, ab); b1 <- cbind(b1, bb); se1 <- cbind(se1, sb); mbar <- cbind(mbar, as.matrix(E %*% Y))
  }
  if (is.null(a)) stop("no gene passes min_det")
  mbar[sw_p == 0, ] <- NA
  lam <- mu <- rep(NA_real_, P_n)
  for (p in seq_len(P_n)) {
    t <- pairs$receiver[p]
    if (anyNA(mbar[p, ]) || anyNA(a[t, ]) || anyNA(b1[p, ])) next
    dc <- tryCatch(.lk_decompose(b1[p, ], se1[p, ], mbar[p, ], a[t, ], LEAK_DEF$huber), error = function(e) list(lam = 0, mu = 0))
    lam[p] <- dc$lam; mu[p] <- dc$mu
  }
  say("done: %d genes", ncol(a))
  structure(data.frame(receiver = lev[pairs$receiver], sender = lev[pairs$sender], lambda = lam, mu = mu,
                       stringsAsFactors = FALSE), n_genes = ncol(a))
}

# leak table from a data.frame or an edgeicg fit (read by field names; edgeicg is not needed)
.leak_table <- function(leak) {
  if (inherits(leak, "EdgeICG")) {
    P <- leak$pairs
    return(data.frame(receiver = leak$types[P$receiver], sender = leak$types[P$sender], lambda = P$lambda, mu = P$mu))
  }
  df <- as.data.frame(leak)
  if (!all(c("receiver", "sender", "lambda") %in% names(df))) stop("leak needs receiver, sender and lambda")
  if (is.null(df$mu)) df$mu <- 0
  df
}

.leak_mats <- function(leak, lev) {
  df <- .leak_table(leak); K <- length(lev); lam <- mu <- matrix(0, K, K)
  ij <- cbind(match(df$receiver, lev), match(df$sender, lev)); ok <- stats::complete.cases(ij)
  lam[ij[ok, , drop = FALSE]] <- ifelse(is.na(df$lambda[ok]), 0, df$lambda[ok])
  mu[ij[ok, , drop = FALSE]] <- ifelse(is.na(df$mu[ok]), 0, df$mu[ok])
  list(lambda = lam, mu = mu)
}

# leak operator on the touching graph: W y = sum_j wn_ij lambda[t_i, t_j] y_j over touching cells of other types;
# cmu_i = sum_j wn_ij mu[t_i, t_j]; W2 has the squared weights (variance of W y)
.leak_W <- function(ct, tcode, lm) {
  n <- length(tcode); wn <- ct$w / pmax(.rsum(ct$w, ct$recv, n)[ct$recv], EPS)
  het <- tcode[ct$recv] != tcode[ct$send]; tr <- tcode[ct$recv]; ts <- tcode[ct$send]
  x <- wn * lm$lambda[cbind(tr, ts)] * het
  list(W = Matrix::sparseMatrix(i = ct$recv, j = ct$send, x = x, dims = c(n, n)),
       W2 = Matrix::sparseMatrix(i = ct$recv, j = ct$send, x = x^2, dims = c(n, n)),
       cmu = .rsum(wn * lm$mu[cbind(tr, ts)] * het, ct$recv, n))
}
