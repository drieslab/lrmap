#' Multi-scale co-location of ligands and receptors between cell types
#'
#' For each sender type A, receiver type B (A != B unless `same_type`) and ligand-receptor pair (L, R), the
#' score at scale s is the mean over receiver cells i of type B of R_i * X_i, where X_i is the ligand
#' carried by type-A cells around i: sum over neighbours j of type A of w_s(i, j) L_j, divided by the summed
#' weight of all neighbours. Expression is on the composition scale (counts / library * 100). Scales:
#' contact (k nearest cells, weight exp(-d / median d)), paracrine (Gaussian, `sigma_p`) and niche (Gaussian,
#' `sigma_n`); paracrine and niche are evaluated on grids. The score is high when receptor-expressing B cells
#' sit among ligand-expressing A cells.
#'
#' The null distribution comes from toroidal shifts of the whole sender layer (cells with their types and
#' expression) against the receiver layer: each layer keeps its own spatial pattern, only their alignment is
#' broken. z = (score - mean of shifted scores) / their standard deviation; p is one-sided (co-location
#' above chance); q is Benjamini-Hochberg within each scale.
#'
#' The score measures co-location: where and at which range communication through a pair is possible. Within
#' one co-located pair of cell types, ligand-receptor pairs differ mainly through their expression and its
#' spatial pattern; a call does not show that the pair is used. The type-level co-location (`coloc`) is
#' reported separately so pair-level results can be read against it. Giotto's original cell-cell
#' communication test (`spatCellCellcom`) is the contact-only predecessor of this analysis.
#'
#' @param counts genes x cells count matrix (dense or sparse) with gene rownames. Only the genes of `lr` are read.
#' @param coords cells x 2 coordinates (microns).
#' @param types cell type per cell.
#' @param lr ligand-receptor table: columns `ligand` and `receptor`, optionally `ligand_subunits` and
#'   `receptor_subunits` (subunit genes joined by `"_"`). A complex ligand is the minimum over its subunits,
#'   a complex receptor the mean.
#' @param lib library size per cell (default: column sums of `counts`, which is only right if `counts` holds
#'   all genes).
#' @param scales any of "contact", "paracrine", "niche".
#' @param k contact neighbours. @param sigma_p,sigma_n paracrine and niche kernel widths (microns).
#' @param bin_p,bin_n grid bins for paracrine and niche (microns).
#' @param B number of toroidal shifts. @param min_shift minimum shift length (microns).
#' @param min_det a pair is tested for (A, B) when the ligand is detected in at least `min_det` of A cells and
#'   the receptor (every subunit) in at least `min_det` of B cells.
#' @param min_cells minimum cells of the sender and of the receiver type.
#' @param same_type also test A = B.
#' @param leak optional leak coefficients: an `EdgeICG` fit or a data.frame with `receiver`, `sender`
#'   (type names), `lambda` and optionally `mu`. Transcripts leaked from touching cells of other types are then
#'   removed from every cell before scoring.
#' @param chunk ligands blurred at a time (memory: grid bins x chunk numbers).
#' @param seed seed for the shifts. @param verbose progress messages.
#' @return an `LRmap` object: `table` (one row per tested sender, receiver, pair and scale), `coloc`
#'   (type-level co-location per scale), `n_pairs`, `n_covered`, `settings`, `timing`.
#' @export
lrmap <- function(counts, coords, types, lr, lib = NULL, scales = c("contact", "paracrine", "niche"), k = 8,
                  sigma_p = 25, sigma_n = 150, bin_p = 5, bin_n = 10, B = 200L, min_shift = 450, min_det = 0.05,
                  min_cells = 30L, same_type = FALSE, leak = NULL, chunk = 64L, seed = 1L, verbose = TRUE) {
  t0 <- proc.time()[[3]]; tm <- c()
  say <- function(...) if (verbose) message(sprintf("[lrmap %5.0f s] ", proc.time()[[3]] - t0), sprintf(...))
  scales <- match.arg(scales, c("contact", "paracrine", "niche"), several.ok = TRUE)
  coords <- as.matrix(coords)[, 1:2, drop = FALSE]; n <- nrow(coords)
  if (ncol(counts) != n || length(types) != n) stop("counts, coords and types must describe the same cells")
  types <- as.character(types); lev <- sort(unique(types)); tcode <- match(types, lev); K <- length(lev)
  if (is.null(lib)) lib <- Matrix::colSums(counts)
  tb <- .lr_table(lr, rownames(counts)); tab <- tb$table
  if (!nrow(tab)) stop("no ligand-receptor pair has all its genes in counts")
  genes <- unique(unlist(c(tab$lg, tab$rg)))
  X <- t(as.matrix(counts[genes, , drop = FALSE])); D <- (X > 0) * 1
  Y <- X * (100 / pmax(lib, 1)); rm(X)
  if (!is.null(leak)) Y <- .unleak(Y, tcode, K, coords, k, .leak_mats(leak, lev))
  Le <- .entities(tab$ligand, tab$lg, Y, D, "min"); Re <- .entities(tab$receptor, tab$rg, Y, D, "mean")
  rm(Y, D)
  Lm <- cbind(Le$E, .cells = 1); Rm <- cbind(Re$E, .cells = 1)
  nt <- tabulate(tcode, K)
  detL <- rowsum(cbind(Le$D, .cells = 1), tcode) / nt; detR <- rowsum(cbind(Re$D, .cells = 1), tcode) / nt
  mL <- rowsum(Lm, tcode) / nt; mR <- rowsum(Rm, tcode) / nt

  # tests: every (sender A, receiver B, pair) passing the expression filters, plus the type-level rows (.cells)
  li <- c(match(tab$ligand, colnames(Lm)), ncol(Lm)); ri <- c(match(tab$receptor, colnames(Rm)), ncol(Rm))
  ok_t <- nt >= min_cells
  tests <- do.call(rbind, lapply(seq_along(li), function(p) {
    g <- expand.grid(A = which(ok_t & detL[, li[p]] >= min_det), Bt = which(ok_t & detR[, ri[p]] >= min_det))
    if (!same_type) g <- g[g$A != g$Bt, , drop = FALSE]
    if (nrow(g)) cbind(g, p = p) else NULL
  }))
  tests$li <- li[tests$p]; tests$ri <- ri[tests$p]
  say("%d of %d pairs covered; %d sender-receiver-pair tests per scale", tb$n_covered, tb$n,
      sum(tests$p <= nrow(tab)))

  box <- c(min(coords[, 1]), diff(range(coords[, 1])) + 1e-6, min(coords[, 2]), diff(range(coords[, 2])) + 1e-6)
  sh <- .shifts(B, box, min_shift, seed)
  res <- list()
  for (s in scales) {
    t1 <- proc.time()[[3]]
    S <- if (s == "contact") .scores_contact(coords, tcode, K, Lm, Rm, tests, sh, box, k)
         else .scores_grid(coords, tcode, K, Lm, Rm, tests, sh, box, if (s == "paracrine") sigma_p else sigma_n,
                           if (s == "paracrine") bin_p else bin_n, chunk)
    res[[s]] <- .summarise(S, tests, s)
    tm[s] <- proc.time()[[3]] - t1; say("%s done", s)
  }
  R <- do.call(rbind, res)
  out <- data.frame(sender = lev[R$A], receiver = lev[R$Bt], ligand = c(tab$ligand, ".cells")[R$p],
                    receptor = c(tab$receptor, ".cells")[R$p], scale = R$scale, score = R$score,
                    null_mean = R$null_mean, fold = R$score / pmax(R$null_mean, EPS), z = R$z, p = R$pv,
                    p_emp = R$p_emp, q = NA_real_, mean_L = mL[cbind(R$A, R$li)], mean_R = mR[cbind(R$Bt, R$ri)],
                    det_L = detL[cbind(R$A, R$li)], det_R = detR[cbind(R$Bt, R$ri)], n_sender = nt[R$A],
                    n_receiver = nt[R$Bt], stringsAsFactors = FALSE)
  cells <- out$ligand == ".cells"
  for (s in scales) { w <- which(out$scale == s & !cells); out$q[w] <- stats::p.adjust(out$p[w], "BH") }
  for (s in scales) { w <- which(out$scale == s & cells); out$q[w] <- stats::p.adjust(out$p[w], "BH") }
  coloc <- out[cells, c("sender", "receiver", "scale", "score", "null_mean", "fold", "z", "p", "q", "n_sender",
                        "n_receiver")]
  names(coloc)[names(coloc) == "score"] <- "share"
  tm["total"] <- proc.time()[[3]] - t0
  structure(list(table = out[!cells, ], coloc = coloc, types = lev, n_pairs = tb$n, n_covered = tb$n_covered,
                 settings = list(scales = scales, k = k, sigma_p = sigma_p, sigma_n = sigma_n, bin_p = bin_p,
                                 bin_n = bin_n, B = B, min_shift = min_shift, min_det = min_det,
                                 min_cells = min_cells, same_type = same_type, leak = !is.null(leak)),
                 timing = tm), class = "LRmap")
}

# null summaries per test row: S = tests x (1 + B), first column observed
.summarise <- function(S, tests, s) {
  s0 <- S[, 1]; nul <- S[, -1, drop = FALSE]; B <- ncol(nul)
  mu <- rowMeans(nul); sdv <- sqrt(rowSums((nul - mu)^2) / (B - 1))
  z <- (s0 - mu) / pmax(sdv, EPS); z[sdv < EPS & abs(s0 - mu) < EPS] <- 0
  cbind(tests, scale = s, score = s0, null_mean = mu, z = z, pv = stats::pnorm(z, lower.tail = FALSE),
        p_emp = (1 + rowSums(nul >= s0)) / (B + 1))
}

# contact: k nearest cells of the (shifted) sender layer per receiver
.scores_contact <- function(coords, tcode, K, Lm, Rm, tests, sh, box, k) {
  nn <- RANN::nn2(coords, k = k + 1); d <- as.vector(nn$nn.dists[, -1]); d0 <- stats::median(d)
  cap <- 3 * d0
  blk <- split(seq_len(nrow(tests)), (tests$A - 1L) * K + tests$Bt)
  S <- matrix(0, nrow(tests), nrow(sh))
  for (v in seq_len(nrow(sh))) {
    e <- .contact_edges(coords, k, d0, cap, if (v == 1) NULL else sh[v, ], box)
    nv <- tabulate(tcode[e$valid], K)
    key <- (tcode[e$send] - 1L) * K + tcode[e$recv]
    eb <- split(seq_along(key), key)
    for (b in names(blk)) {
      rows <- blk[[b]]; ee <- eb[[b]]; if (is.null(ee)) next
      Bt <- tests$Bt[rows[1]]
      S[rows, v] <- colSums(e$wn[ee] * Rm[e$recv[ee], tests$ri[rows], drop = FALSE] *
                              Lm[e$send[ee], tests$li[rows], drop = FALSE]) / max(nv[Bt], 1)
    }
  }
  S
}

# paracrine / niche: blurred ligand fields of each sender type, read at receivers (shifted sender layer)
.scores_grid <- function(coords, tcode, K, Lm, Rm, tests, sh, box, sigma, bin, chunk) {
  gr <- .grid(coords, box, sigma, bin); nc <- gr$nx * gr$ny
  Tf <- .blur(.bin_sum(matrix(1, nrow(coords), 1), gr$lin, nc), gr)[, 1]
  tau <- 0.1 * stats::median(Tf[gr$lin])
  LV <- vapply(seq_len(nrow(sh)), function(v) .shift_lin(gr, sh[v, ], coords, box), numeric(nrow(coords)))
  VA <- matrix(Tf[LV] >= tau, nrow(coords))                                   # valid receivers per shift
  nv <- rowsum(VA * 1, tcode)                                                  # K x shifts (types present)
  nvK <- matrix(0, K, nrow(sh)); nvK[as.integer(rownames(nv)), ] <- nv
  S <- matrix(0, nrow(tests), nrow(sh))
  for (A in unique(tests$A)) {
    ta <- which(tests$A == A); la <- unique(tests$li[ta]); sa <- which(tcode == A)
    for (ch in split(la, ceiling(seq_along(la) / chunk))) {
      Fq <- .blur(.bin_sum(Lm[sa, ch, drop = FALSE], gr$lin[sa], nc), gr) / pmax(Tf, EPS)
      Fq[Tf < tau, ] <- 0
      rows_c <- ta[tests$li[ta] %in% ch]
      for (Bt in unique(tests$Bt[rows_c])) {
        rows <- rows_c[tests$Bt[rows_c] == Bt]; cb <- which(tcode == Bt)
        RB <- Rm[cb, tests$ri[rows], drop = FALSE]; lp <- match(tests$li[rows], ch)
        for (v in seq_len(nrow(sh))) {
          Xv <- Fq[LV[cb, v], lp, drop = FALSE]
          S[rows, v] <- colSums(RB * Xv * VA[cb, v]) / max(nvK[Bt, v], 1)
        }
      }
    }
  }
  S
}

#' Significant ligand-receptor co-location
#'
#' @param x an `LRmap` object. @param scale scales to keep.
#' @param q_thr q-value threshold. @param min_fold minimum score over the shifted-null mean.
#' @return rows of `x$table` with q < `q_thr` and fold >= `min_fold`, by q.
#' @export
lrmap_results <- function(x, scale = c("contact", "paracrine", "niche"), q_thr = 0.05, min_fold = 1.25) {
  tb <- x$table; keep <- which(tb$scale %in% scale & !is.na(tb$q) & tb$q < q_thr & tb$fold >= min_fold)
  r <- tb[keep, ]; r[order(r$q, -r$z), ]
}

#' @export
print.LRmap <- function(x, ...) {
  cat(sprintf("LRmap: %d cell types, %d of %d ligand-receptor pairs covered, %d shifts%s\n", length(x$types),
              x$n_covered, x$n_pairs, x$settings$B, if (x$settings$leak) ", leak removed" else ""))
  for (s in x$settings$scales) {
    a <- x$table[x$table$scale == s, ]; cl <- x$coloc[x$coloc$scale == s, ]
    cat(sprintf("  %-9s %6d tests, %5d co-located (q < 0.05, fold >= 1.25); type pairs co-located: %d of %d\n",
                s, nrow(a), nrow(lrmap_results(x, s)), sum(cl$q < 0.05, na.rm = TRUE), nrow(cl)))
  }
  invisible(x)
}

# ---- leak -------------------------------------------------------------------------------------------------

.leak_mats <- function(leak, lev) {
  if (inherits(leak, "EdgeICG")) {
    P <- leak$pairs; df <- data.frame(receiver = leak$types[P$receiver], sender = leak$types[P$sender],
                                      lambda = P$lambda, mu = P$mu)
  } else df <- as.data.frame(leak)
  if (!all(c("receiver", "sender", "lambda") %in% names(df))) stop("leak needs receiver, sender and lambda")
  if (is.null(df$mu)) df$mu <- 0
  K <- length(lev); lam <- mu <- matrix(0, K, K)
  ij <- cbind(match(df$receiver, lev), match(df$sender, lev)); ok <- stats::complete.cases(ij)
  lam[ij[ok, , drop = FALSE]] <- ifelse(is.na(df$lambda[ok]), 0, df$lambda[ok])
  mu[ij[ok, , drop = FALSE]] <- ifelse(is.na(df$mu[ok]), 0, df$mu[ok])
  list(lambda = lam, mu = mu)
}

# remove transcripts leaked from touching cells of other types (edgeicg's leak model, contact kNN graph):
# y_i - sum_j wn_ij lambda[t_i, t_j] y_j - (sum_j wn_ij mu[t_i, t_j]) * mean_{t_i}(y)
.unleak <- function(Y, tcode, K, coords, k, lm) {
  n <- nrow(coords); nn <- RANN::nn2(coords, k = k + 1)
  recv <- rep(seq_len(n), k); send <- as.vector(nn$nn.idx[, -1]); d <- as.vector(nn$nn.dists[, -1])
  w <- exp(-d / stats::median(d)); wn <- w / pmax(.rsum(w, recv, n)[recv], EPS)
  het <- tcode[recv] != tcode[send]; tr <- tcode[recv]; ts <- tcode[send]
  Wl <- Matrix::sparseMatrix(i = recv, j = send, x = wn * lm$lambda[cbind(tr, ts)] * het, dims = c(n, n))
  cmu <- .rsum(wn * lm$mu[cbind(tr, ts)] * het, recv, n)
  a <- rowsum(Y, tcode) / tabulate(tcode, K)[sort(unique(tcode))]
  am <- matrix(0, K, ncol(Y)); am[as.integer(rownames(a)), ] <- a
  Y - as.matrix(Wl %*% Y) - cmu * am[tcode, , drop = FALSE]
}
