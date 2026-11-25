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
#' @param counts raw counts: a genes x cells matrix (dense or sparse) with gene rownames, or a reader
#'   `function(genes)` returning a genes x cells matrix for the requested genes that exist (rows named; cells in
#'   the order of `coords`), e.g. reading from a disk-backed store. Only the genes of `lr` are read, once.
#' @param coords cells x 2 coordinates (microns).
#' @param types cell type per cell.
#' @param lr ligand-receptor table: columns `ligand` and `receptor`, optionally `ligand_subunits` and
#'   `receptor_subunits` (subunit genes joined by `"_"`). A complex ligand is the minimum over its subunits,
#'   a complex receptor the mean.
#' @param lib library size per cell (all genes). Default: column sums of `counts` when it is a matrix holding all
#'   genes; required with a reader.
#' @param scales any of "contact", "paracrine", "niche", and "expression": the expression-only reference used
#'   by single-cell analyses (mean ligand in the sender type x mean receptor in the receiver type against random
#'   relabelling of cells, no spatial information), on the same tests and values; compare with
#'   [lrmap_compare()].
#' @param k contact neighbours. @param sigma_p,sigma_n paracrine and niche kernel widths (microns).
#' @param bin_p,bin_n grid bins for paracrine and niche (microns); keep at most a third of the kernel width.
#' @param B number of toroidal shifts. @param min_shift minimum shift length (microns).
#' @param min_det a pair is tested for (A, B) when the ligand is detected in at least `min_det` of A cells and
#'   the receptor (every subunit) in at least `min_det` of B cells.
#' @param min_cells minimum cells of the sender and of the receiver type.
#' @param same_type also test A = B.
#' @param leak optional leak (spillover) correction: a data.frame with `receiver`, `sender` (type names), `lambda`
#'   and optionally `mu` (e.g. from [lrmap_leak()]), an edgeICG fit (its leak coefficients are read; edgeICG is not
#'   needed), or `TRUE` to estimate them with [lrmap_leak()] on all genes of `counts` (matrix input only).
#'   Transcripts leaked from touching cells of other types are removed from every gene before ligand and receptor
#'   values are formed: y - sum_j wn_ij lambda y_j - cmu_i a_t, with a_t the receiver type's first-stage baseline.
#' @param leak_pool how leak coefficients of pairs with little contact are pooled when `leak` is an edgeICG fit or
#'   a plain table (see [lrmap_leak()], argument `pool`); tables from [lrmap_leak()] are used as given.
#' @param spill_edges optional touching edges for the leak correction (`recv`, `send`, `w`; cell indices in the
#'   order of `coords`), e.g. from cell outlines; default: the k nearest neighbours, weight exp(-d / median d).
#' @param chunk ligands blurred at a time; default from the grid size (about `grid_mem` numbers per chunk).
#' @param grid_mem numbers held per blurred chunk (sets the default `chunk`).
#' @param z_cand the ligand-side part of the sampling variance on grid scales is computed for tests with z above
#'   this before it (adding variance only lowers z).
#' @param backend "auto" uses the compiled kernels of the lrmaprs package when it is installed.
#' @param threads threads for the compiled kernels (1: the calling thread only; 0: all cores).
#' @param seed seed for the shifts. @param verbose progress messages.
#' @details Three tests are reported for every spatial row. `q` (with `z`, `fold`): co-location of ligand-expressing
#'   senders and receptor-expressing receivers against shifted sender layers; this includes the co-location of
#'   the two cell types, so any pair expressed in two neighbouring types is called. `q_type` (with `z_type`,
#'   `fold_type`; the default of [lrmap_results()]): co-location beyond that of the two cell types, i.e. whether
#'   ligand-high senders sit nearer receptor-high receivers than other cells of the two types do, calibrated by
#'   the same shifts. `q_shuffle`: the same question against shuffling ligand and receptor values within their
#'   cell types (layout fixed, exact moments); it ignores smooth spatial trends of expression and is
#'   anti-conservative when they are present.
#' @return an `LRmap` object: `table` (one row per tested sender, receiver, pair and scale), `coloc`
#'   (type-level co-location per scale), `n_pairs`, `n_covered`, `settings`, `timing`.
#' @export
lrmap <- function(counts, coords, types, lr, lib = NULL, scales = c("contact", "paracrine", "niche"), k = 8,
                  sigma_p = 25, sigma_n = 150, bin_p = 8, bin_n = 25, B = 200L, min_shift = 450, min_det = 0.05,
                  min_cells = 30L, same_type = FALSE, leak = NULL, leak_pool = c("median", "q1", "none"), spill_edges = NULL, chunk = NULL, grid_mem = 4e7, z_cand = 1.5,
                  backend = c("auto", "rust", "r"), threads = 1L, seed = 1L, verbose = TRUE) {
  backend <- match.arg(backend); leak_pool <- match.arg(leak_pool); has_rs <- requireNamespace("lrmaprs", quietly = TRUE)
  if (backend == "rust" && !has_rs) stop("backend = 'rust' needs the lrmaprs package")
  rs <- if (backend != "r" && has_rs) as.integer(threads) else NULL
  t0 <- proc.time()[[3]]; tm <- c()
  say <- function(...) if (verbose) message(sprintf("[lrmap %5.0f s] ", proc.time()[[3]] - t0), sprintf(...))
  scales <- match.arg(scales, c("contact", "paracrine", "niche", "expression"), several.ok = TRUE)
  coords <- as.matrix(coords)[, 1:2, drop = FALSE]; n <- nrow(coords)
  if (length(types) != n) stop("coords and types must describe the same cells")
  types <- as.character(types); lev <- sort(unique(types)); K <- length(lev)
  if (isTRUE(leak)) {
    if (is.function(counts)) stop("leak = TRUE needs a counts matrix; with a reader, estimate with lrmap_leak()")
    say("estimating leak coefficients on all genes"); leak <- lrmap_leak(counts, coords, types, lib, pool = leak_pool, verbose = FALSE)
  }
  # cells sorted by type: every type is a contiguous block
  ord <- order(match(types, lev)); tcode <- match(types, lev)[ord]; coords <- coords[ord, , drop = FALSE]
  nt <- tabulate(tcode, K); off <- c(0L, cumsum(nt))[seq_len(K)]
  rng <- lapply(seq_len(K), function(t) off[t] + seq_len(nt[t]))

  # read the table's genes once (sparse genes x cells)
  tg <- .lr_table(lr, NULL)$genes
  Xg <- if (is.function(counts)) counts(tg) else counts[intersect(tg, rownames(counts)), , drop = FALSE]
  if (ncol(Xg) != n) stop("counts must have one column per cell")
  if (is.null(lib)) {
    if (is.function(counts)) stop("lib is required when counts is a reader function")
    lib <- Matrix::colSums(counts)
  }
  lib <- pmax(as.numeric(lib)[ord], 1)
  Xg <- methods::as(Xg, "CsparseMatrix")[, ord, drop = FALSE]
  tb <- .lr_table(lr, rownames(Xg)); tab <- tb$table
  if (!nrow(tab)) stop("no ligand-receptor pair has all its genes in counts")
  Y <- Matrix::Diagonal(x = 100 / lib) %*% Matrix::t(Xg[unique(unlist(c(tab$lg, tab$rg))), , drop = FALSE])
  rm(Xg)
  # entities: ligand = minimum over its subunit genes, receptor = mean. No all-cell entity matrices are built:
  # values, sampling variances and detection are formed per cell type from Y, in chunks of entities.
  Ln <- unique(tab$ligand); Rn <- unique(tab$receptor)
  subL <- tab$lg[match(Ln, tab$ligand)]; subR <- tab$rg[match(Rn, tab$receptor)]; nl <- length(Ln); nr <- length(Rn)
  pv <- 100 / lib                                     # Poisson variance on the composition scale: y * 100 / lib
  # detection per type: share of cells with every subunit detected
  det_ent <- function(sub) {
    out <- matrix(0, K, length(sub)); one <- lengths(sub) == 1
    if (any(one)) {
      g1 <- vapply(sub[one], `[`, "", 1); Yd <- Y[, g1, drop = FALSE]; Yd@x <- as.numeric(Yd@x > 0)
      out[, one] <- as.matrix(Matrix::sparseMatrix(i = tcode, j = seq_len(n), x = 1, dims = c(K, n)) %*% Yd)
    }
    for (j in which(!one)) out[, j] <- tabulate(tcode[Matrix::rowSums(Y[, sub[[j]], drop = FALSE] > 0) == length(sub[[j]])], K)
    out / nt
  }
  detL <- det_ent(subL); detR <- det_ent(subR)
  lc <- NULL
  if (!is.null(leak)) {   # leak removal per gene (y - W y - cmu * a_t), before ligand / receptor values are formed
    ct <- if (is.null(spill_edges)) NULL else {
      e <- as.data.frame(spill_edges); inv <- order(ord); list(recv = inv[e$recv], send = inv[e$send], w = e$w) }
    ct <- .lk_edges(coords, k, ct)
    ex <- .lk_exposures(coords, tcode, K, utils::modifyList(LEAK_DEF, list(k = k)), rs)
    # coefficients from lrmap_leak() are used as given; others (an edgeICG fit, a plain table) are pooled here
    lt <- .leak_table(leak)
    if (is.null(lt$status)) {
      sup <- .lk_support(ct, tcode, K)
      ij <- cbind(match(lt$receiver, lev), match(lt$sender, lev)); ok <- stats::complete.cases(ij)
      lt <- lt[ok, , drop = FALSE]; ij <- ij[ok, , drop = FALSE]
      pl <- .leak_pool(ij[, 1], lt$lambda, lt$mu, sup[ij], leak_pool, LEAK_DEF$min_support)
      say("leak: %d of %d receiver-sender pairs pooled (contact support < %g or no usable estimate)",
          sum(pl$status != "ok"), nrow(lt), LEAK_DEF$min_support)
      lt$lambda <- pl$lambda; lt$mu <- pl$mu
    }
    lm <- .leak_mats(lt, lev)
    lw <- .leak_W(ct, tcode, lm)
    a <- .lk_intercepts(.lk_designs(ex, tcode, K, lib, LEAK_DEF$min_cells), Y, K); a[is.na(a)] <- 0; rm(ex)
    lc <- list(WY = lw$W %*% Y, cmu = lw$cmu, a = a, W2 = lw$W2); rm(lw)
  }
  # entity values for cells r (of type t when leak-corrected), chunks of 100 entities
  ent_vals <- function(sub, fun, r, t, corrected = TRUE) {
    do.call(cbind, lapply(split(seq_along(sub), ceiling(seq_along(sub) / 100)), function(ii) {
      g <- unique(unlist(sub[ii])); G <- as.matrix(Y[r, g, drop = FALSE])
      if (corrected && !is.null(lc)) G <- G - as.matrix(lc$WY[r, g, drop = FALSE]) - outer(lc$cmu[r], lc$a[t, g])
      colnames(G) <- g
      matrix(vapply(sub[ii], function(s) if (length(s) == 1) G[, s] else if (fun == "min")
                    do.call(pmin, as.data.frame(G[, s, drop = FALSE])) else rowMeans(G[, s, drop = FALSE]),
                    numeric(length(r))), nrow = length(r))
    }))
  }
  # sampling variance: Poisson of the raw entity (mean of k subunits: / k), plus that of the removed leak (W2)
  ent_var <- function(sub, fun, r) {
    kd <- if (fun == "mean") lengths(sub) else rep(1, length(sub))
    V <- pv[r] * sweep(ent_vals(sub, fun, r, NA, FALSE), 2, kd, "/")
    if (!is.null(lc)) {
      W2r <- lc$W2[r, , drop = FALSE]; nb <- which(Matrix::colSums(W2r) > 0)
      if (length(nb)) V <- V + as.matrix(W2r[, nb, drop = FALSE] %*% (pv[nb] * sweep(ent_vals(sub, fun, nb, NA, FALSE), 2, kd, "/")))
    }
    V
  }

  # tests: every (sender A, receiver B, pair) passing the expression filters, plus the type-level rows (.cells,
  # entity index nl + 1 / nr + 1)
  li <- c(match(tab$ligand, Ln), nl + 1L); ri <- c(match(tab$receptor, Rn), nr + 1L)
  dL <- cbind(detL, 1); dR <- cbind(detR, 1); ok_t <- nt >= min_cells
  tests <- do.call(rbind, lapply(seq_along(li), function(p) {
    g <- expand.grid(A = which(ok_t & dL[, li[p]] >= min_det), Bt = which(ok_t & dR[, ri[p]] >= min_det))
    if (!same_type) g <- g[g$A != g$Bt, , drop = FALSE]
    if (nrow(g)) cbind(g, p = p) else NULL
  }))
  tests$li <- li[tests$p]; tests$ri <- ri[tests$p]
  say("%d of %d pairs covered; %d sender-receiver-pair tests per scale", tb$n_covered, tb$n,
      sum(tests$p <= nrow(tab)))

  # dense blocks: per type, only the entities it is tested with (values and their variances)
  blk <- function(sub, fun, nent, idx, side) {
    lapply(seq_len(K), function(t) {
      cols <- sort(unique(idx[tests[[side]] == t])); if (!length(cols)) return(NULL)
      ent <- cols[cols <= nent]; r <- rng[[t]]
      M <- matrix(1, length(r), length(cols)); Vv <- matrix(0, length(r), length(cols))
      if (length(ent)) { w <- match(ent, cols); M[, w] <- ent_vals(sub[ent], fun, r, t); Vv[, w] <- ent_var(sub[ent], fun, r) }
      list(M = M, V = Vv, cols = cols)
    })
  }
  Lb <- blk(subL, "min", nl, tests$li, "A"); Rb <- blk(subR, "mean", nr, tests$ri, "Bt")
  # expression-only reference: mean and variance of every tested entity over all cells (label-shuffle null)
  est <- function(sub, fun, nent, ents) {
    s1 <- s2 <- numeric(length(ents))
    for (ii in split(seq_along(ents), ceiling(seq_along(ents) / 100))) for (t in seq_len(K)) {
      M <- ent_vals(sub[ents[ii]], fun, rng[[t]], t); s1[ii] <- s1[ii] + colSums(M); s2[ii] <- s2[ii] + colSums(M^2)
    }
    m <- s1 / n; v <- mm <- numeric(nent + 1L); mm[ents] <- m; v[ents] <- pmax(s2 / n - m^2, 0); list(mean = mm, var = v)
  }
  if ("expression" %in% scales) {
    xt <- tests$p <= nrow(tab)
    stL <- est(subL, "min", nl, sort(unique(tests$li[xt]))); stR <- est(subR, "mean", nr, sort(unique(tests$ri[xt])))
  }
  rm(lc, Y); invisible(gc(FALSE))
  bmean <- function(Bl, ncol_) { m <- matrix(NA_real_, K, ncol_)
    for (t in seq_len(K)) if (!is.null(Bl[[t]])) m[t, Bl[[t]]$cols] <- colMeans(Bl[[t]]$M); m }
  mLx <- bmean(Lb, nl + 1L); mRx <- bmean(Rb, nr + 1L)
  tests$lc <- vapply(seq_len(nrow(tests)), function(i) match(tests$li[i], Lb[[tests$A[i]]]$cols), 1L)
  tests$rc <- vapply(seq_len(nrow(tests)), function(i) match(tests$ri[i], Rb[[tests$Bt[i]]]$cols), 1L)
  say("blocks: %.0f MB", (sum(vapply(Lb, function(b) if (is.null(b)) 0 else length(b$M), 1)) +
                            sum(vapply(Rb, function(b) if (is.null(b)) 0 else length(b$M), 1))) * 16 / 2^20)

  box <- c(min(coords[, 1]), diff(range(coords[, 1])) + 1e-6, min(coords[, 2]), diff(range(coords[, 2])) + 1e-6)
  reach <- max(0, c(paracrine = 3 * sigma_p, niche = 3 * sigma_n)[intersect(scales, c("paracrine", "niche"))])
  if (min_shift < reach)
    warning(sprintf(paste("min_shift (%.0f um) is shorter than the widest kernel's reach (3 sigma = %.0f um): shifted",
                          "layouts stay partly aligned and the null is conservative; set min_shift >= %.0f"),
                    min_shift, reach, reach), call. = FALSE)
  sh <- .shifts(B, box, min_shift, seed)
  ctx <- list(coords = coords, tcode = tcode, K = K, off = off, rng = rng, Lb = Lb, Rb = Rb, tests = tests,
              sh = sh, box = box, rs = rs)
  # per test: mean and sum of squares of the ligand over sender cells and the receptor over receiver cells; the
  # type-level (.cells) row of the same sender / receiver
  bst <- function(Bl, typ, col) {
    m <- ss <- numeric(nrow(tests))
    for (t in unique(typ)) {
      w <- which(typ == t); M <- Bl[[t]]$M; cm <- colMeans(M); cs <- pmax(colSums(M^2) - nrow(M) * cm^2, 0)
      m[w] <- cm[col[w]]; ss[w] <- cs[col[w]]
    }
    list(mean = m, ss = ss)
  }
  Lst <- bst(Lb, tests$A, tests$lc); Rst <- bst(Rb, tests$Bt, tests$rc)
  cr <- which(tests$p == nrow(tab) + 1L); crow <- cr[match(paste(tests$A, tests$Bt), paste(tests$A[cr], tests$Bt[cr]))]
  res <- list()
  for (s in setdiff(scales, "expression")) {
    t1 <- proc.time()[[3]]
    S <- if (s == "contact") .scores_contact(ctx, k)
         else .scores_grid(ctx, if (s == "paracrine") sigma_p else sigma_n, if (s == "paracrine") bin_p else bin_n,
                           chunk, grid_mem, z_cand)
    res[[s]] <- cbind(.summarise(S$S, S$V, tests, s), .beyond_type(S$S, S$V, crow, Lst$mean * Rst$mean),
                      .within_shuffle(S$S[, 1], tests, .within_graph(ctx, s, k, if (s == "paracrine") sigma_p else sigma_n,
                                                                      if (s == "paracrine") bin_p else bin_n), K, nt, Lst, Rst))
    tm[s] <- proc.time()[[3]] - t1; say("%s done", s); rm(S); invisible(gc(FALSE))
  }
  if ("expression" %in% scales) res$expression <- .expression_test(tests[xt, ], mLx, mRx, stL, stR, nt, n)
  if (!is.null(res$expression)) res$expression[, c("expected_type", "z_type", "p_type", "z_shuffle", "p_shuffle")] <- NA_real_
  R <- do.call(rbind, res)
  out <- data.frame(sender = lev[R$A], receiver = lev[R$Bt], ligand = c(tab$ligand, ".cells")[R$p],
                    receptor = c(tab$receptor, ".cells")[R$p], scale = R$scale, score = R$score,
                    null_mean = R$null_mean, fold = .ratio(R$score, R$null_mean), z = R$z, p = R$pv,
                    p_emp = R$p_emp, q = NA_real_, mean_L = mLx[cbind(R$A, R$li)], mean_R = mRx[cbind(R$Bt, R$ri)],
                    det_L = dL[cbind(R$A, R$li)], det_R = dR[cbind(R$Bt, R$ri)], n_sender = nt[R$A],
                    n_receiver = nt[R$Bt], fold_type = .ratio(R$score, R$expected_type), z_type = R$z_type,
                    p_type = R$p_type, q_type = NA_real_, z_shuffle = R$z_shuffle, p_shuffle = R$p_shuffle,
                    q_shuffle = NA_real_, stringsAsFactors = FALSE)
  cells <- out$ligand == ".cells"
  for (s in scales) {
    w <- which(out$scale == s & !cells); out$q[w] <- stats::p.adjust(out$p[w], "BH")
    if (s != "expression") { out$q_type[w] <- stats::p.adjust(out$p_type[w], "BH"); out$q_shuffle[w] <- stats::p.adjust(out$p_shuffle[w], "BH") }
  }
  for (s in scales) { w <- which(out$scale == s & cells); out$q[w] <- stats::p.adjust(out$p[w], "BH") }
  coloc <- out[cells, c("sender", "receiver", "scale", "score", "null_mean", "fold", "z", "p", "q", "n_sender",
                        "n_receiver")]
  names(coloc)[names(coloc) == "score"] <- "share"
  tm["total"] <- proc.time()[[3]] - t0
  structure(list(table = out[!cells, ], coloc = coloc, types = lev, n_pairs = tb$n, n_covered = tb$n_covered,
                 settings = list(scales = scales, k = k, sigma_p = sigma_p, sigma_n = sigma_n, bin_p = bin_p,
                                 bin_n = bin_n, B = B, min_shift = min_shift, min_det = min_det,
                                 min_cells = min_cells, same_type = same_type, leak = !is.null(leak),
                                 backend = if (is.null(rs)) "r" else "rust", threads = threads),
                 timing = tm), class = "LRmap")
}

# Expression-only test (as in single-cell analyses, no spatial information): score = mean ligand in the sender type
# x mean receptor in the receiver type; null = random relabelling of cells (type sizes kept). The shuffled group
# means have exact moments (sampling without replacement); the product is approximated as normal with the
# delta-method variance (the O(1/n) covariance between disjoint groups is ignored).
.expression_test <- function(tests, mLx, mRx, stL, stR, nt, n) {
  A <- tests$A; Bt <- tests$Bt; li <- tests$li; ri <- tests$ri
  obs <- mLx[cbind(A, li)] * mRx[cbind(Bt, ri)]
  muL <- stL$mean[li]; muR <- stR$mean[ri]
  vL <- stL$var[li] / nt[A] * (n - nt[A]) / (n - 1); vR <- stR$var[ri] / nt[Bt] * (n - nt[Bt]) / (n - 1)
  e0 <- muL * muR; sdv <- sqrt(muR^2 * vL + muL^2 * vR + vL * vR)
  z <- (obs - e0) / pmax(sdv, EPS); z[sdv < EPS] <- 0
  cbind(tests[, c("A", "Bt", "p", "li", "ri")], scale = "expression", score = obs, null_mean = e0, z = z,
        pv = stats::pnorm(z, lower.tail = FALSE), p_emp = NA_real_)
}

# Pair co-location beyond the co-location of the two cell types, calibrated by the shifts: D = S - mean L x mean R
# x (type-level score of the same sender / receiver), for the observed layout and every shift. Under random
# alignment D is near 0; a pair whose ligand-high senders sit nearer its receptor-high receivers than other cells
# of the two types gives D > 0. Smooth spatial trends of expression enter the shifted D's spread.
.beyond_type <- function(S, V, crow, mm) {
  D <- S - mm * S[crow, , drop = FALSE]
  s0 <- D[, 1]; nul <- D[, -1, drop = FALSE]; mu <- rowMeans(nul)
  sdv <- sqrt(rowSums((nul - mu)^2) / max(ncol(nul) - 1, 1) + V)
  z <- (s0 - mu) / pmax(sdv, EPS); z[sdv < EPS] <- 0
  data.frame(expected_type = mm * S[crow, 1], z_type = z, p_type = stats::pnorm(z, lower.tail = FALSE))
}

# Graph sums per sender / receiver type for the within-type shuffle (observed layout). With a_ij the weight of
# sender j in receiver i's ligand exposure (X_i = sum_j a_ij L_j): s_i = sum_j a_ij, q_i = sum_j a_ij^2,
# c_j = sum_i a_ij. Returns K^2 x (nv, S1 = sum s, S2 = sum s^2, Q = sum q, SSc = sum (c - mean c)^2).
.within_graph <- function(ctx, s, k, sigma, bin) {
  coords <- ctx$coords; tcode <- ctx$tcode; K <- ctx$K; n <- nrow(coords); rng <- ctx$rng
  keys <- unique(ctx$tests[, c("A", "Bt")])
  out <- matrix(NA_real_, K * K, 5, dimnames = list(NULL, c("nv", "S1", "S2", "Q", "SSc")))
  if (s == "contact") {
    nn <- RANN::nn2(coords, k = k + 1); d0 <- stats::median(as.vector(nn$nn.dists[, -1])); rm(nn)
    e <- .contact_edges(coords, k, d0, 3 * d0); ta <- tcode[e$send]; tb <- tcode[e$recv]
    for (r in seq_len(nrow(keys))) {
      A <- keys$A[r]; B <- keys$Bt[r]; w <- which(ta == A & tb == B)
      sv <- .rsum(e$wn[w], e$recv[w], n)[rng[[B]]]; qv <- .rsum(e$wn[w]^2, e$recv[w], n)[rng[[B]]]
      cv <- .rsum(e$wn[w], e$send[w], n)[rng[[A]]]
      out[(A - 1) * K + B, ] <- c(sum(e$valid[rng[[B]]]), sum(sv), sum(sv^2), sum(qv), sum((cv - mean(cv))^2))
    }
    return(out)
  }
  gr <- .grid(coords, ctx$box, sigma, bin); nc <- gr$nx * gr$ny; rs <- ctx$rs
  Tf <- .blur(.bin_sum(matrix(1, n, 1), gr$lin, nc), gr, rs)[, 1]; Tc <- Tf[gr$lin]
  valid <- Tc >= 0.1 * stats::median(Tc); invT <- ifelse(valid, 1 / pmax(Tc, EPS), 0)
  g2 <- gr; g2$ker <- gr$ker^2; g2$Kx@x <- g2$Kx@x^2; g2$Ky@x <- g2$Ky@x^2
  ind <- outer(tcode, seq_len(K), "==") * 1
  FA <- .blur(.bin_sum(ind, gr$lin, nc), gr, rs); QA <- .blur(.bin_sum(ind, gr$lin, nc), g2, rs)
  GB <- .blur(.bin_sum(ind * invT, gr$lin, nc), gr, rs); rm(ind)
  for (r in seq_len(nrow(keys))) {
    A <- keys$A[r]; B <- keys$Bt[r]; rb <- rng[[B]]; lb <- gr$lin[rb]
    sv <- FA[lb, A] * invT[rb]; qv <- QA[lb, A] * invT[rb]^2; cv <- GB[gr$lin[rng[[A]]], B]
    out[(A - 1) * K + B, ] <- c(sum(valid[rb]), sum(sv), sum(sv^2), sum(qv), sum((cv - mean(cv))^2))
  }
  out
}

# Within-type shuffle: ligand values permuted among sender cells, receptor values among receiver cells, layout
# fixed. Exact mean and variance of S = sum_i R_i X_i / nv under both permutations (normal approximation for p).
.within_shuffle <- function(S0, tests, G, K, nt, Lst, Rst) {
  g <- G[(tests$A - 1L) * K + tests$Bt, , drop = FALSE]; nA <- nt[tests$A]; nB <- nt[tests$Bt]; nv <- pmax(g[, "nv"], 1)
  Lm <- Lst$mean; Lss <- Lst$ss; Rm <- Rst$mean; Rss <- Rst$ss; sL2 <- Lss / nA
  E <- Rm * Lm * g[, "S1"] / nv
  t1 <- Rm^2 * g[, "SSc"] * Lss / ((nA - 1) * nv^2)
  Eu2 <- Lm^2 * g[, "S2"] + sL2 * (nA / (nA - 1) * g[, "Q"] - g[, "S2"] / (nA - 1))
  ESu2 <- (Lm * g[, "S1"])^2 + g[, "SSc"] * Lss / (nA - 1)
  V <- t1 + pmax(Eu2 - ESu2 / nB, 0) * Rss / ((nB - 1) * nv^2)
  z <- (S0 - E) / sqrt(pmax(V, EPS)); z[V < EPS] <- 0
  data.frame(z_shuffle = z, p_shuffle = stats::pnorm(z, lower.tail = FALSE))
}

# fold: NA when the expectation is not positive (leak-corrected values can be <= 0)
.ratio <- function(x, e) ifelse(e > EPS, x / pmax(e, EPS), NA_real_)

# null summaries per test row: S = tests x (1 + B), first column observed; V = excess sampling variance of the
# observed score over shifted alignments (heteroscedastic noise tied to the observed layout), added to the
# variance across shifts
.summarise <- function(S, V, tests, s) {
  s0 <- S[, 1]; nul <- S[, -1, drop = FALSE]; B <- ncol(nul)
  mu <- rowMeans(nul); sdv <- sqrt(rowSums((nul - mu)^2) / (B - 1) + V)
  z <- (s0 - mu) / pmax(sdv, EPS); z[sdv < EPS & abs(s0 - mu) < EPS] <- 0
  cbind(tests[, c("A", "Bt", "p", "li", "ri")], scale = s, score = s0, null_mean = mu, z = z,
        pv = stats::pnorm(z, lower.tail = FALSE), p_emp = (1 + rowSums(nul >= s0)) / (B + 1))
}

# excess of the observed alignment's sampling variance over the mean of the shifted ones
.excess <- function(Vs) pmax(Vs[, 1] - rowMeans(Vs[, -1, drop = FALSE]), 0)
NVAR <- 21L

# contact: k nearest cells of the (shifted) sender layer per receiver; score = sum_i R_i sum_j wn_ij L_j / n_B
.scores_contact <- function(ctx, k) {
  coords <- ctx$coords; tcode <- ctx$tcode; K <- ctx$K; tests <- ctx$tests; sh <- ctx$sh
  nn <- RANN::nn2(coords, k = k + 1); d0 <- stats::median(as.vector(nn$nn.dists[, -1])); cap <- 3 * d0; rm(nn)
  nvar <- min(NVAR, nrow(sh))
  if (!is.null(ctx$rs)) {
    em <- matrix(0, 0, 0); g <- function(L, f) lapply(L, function(b) if (is.null(b)) em else b[[f]])
    o <- lrmaprs::rs_contact_scores(coords[, 1], coords[, 2], as.integer(tcode), as.integer(ctx$off),
                                    g(ctx$Lb, "M"), g(ctx$Lb, "V"), g(ctx$Rb, "M"), g(ctx$Rb, "V"),
                                    as.integer(tests$A), as.integer(tests$Bt), as.integer(tests$lc), as.integer(tests$rc),
                                    sh, as.integer(k), d0, cap, ctx$box, as.integer(nvar), ctx$rs)
    return(list(S = o$S, V = .excess(o$Vs)))
  }
  bl <- split(seq_len(nrow(tests)), (tests$A - 1L) * K + tests$Bt)
  S <- matrix(0, nrow(tests), nrow(sh)); Vs <- matrix(0, nrow(tests), nvar)
  for (v in seq_len(nrow(sh))) {
    e <- .contact_edges(coords, k, d0, cap, if (v == 1) NULL else sh[v, ], ctx$box)
    nv <- tabulate(tcode[e$valid], K)
    key <- (tcode[e$send] - 1L) * K + tcode[e$recv]
    eb <- split(seq_along(key), key)
    for (b in names(bl)) {
      rows <- bl[[b]]; ee <- eb[[b]]; if (is.null(ee)) next
      A <- tests$A[rows[1]]; Bt <- tests$Bt[rows[1]]; LB <- ctx$Lb[[A]]; RB <- ctx$Rb[[Bt]]
      ir <- e$recv[ee] - ctx$off[Bt]; js <- e$send[ee] - ctx$off[A]
      Rr <- RB$M[ir, tests$rc[rows], drop = FALSE]; Ls <- LB$M[js, tests$lc[rows], drop = FALSE]
      S[rows, v] <- colSums(e$wn[ee] * Rr * Ls) / max(nv[Bt], 1)
      if (v <= nvar) {                         # sampling variance: score = sum_i R_i X_i, X_i = sum_j wn_ij L_j
        XR <- rowsum(e$wn[ee] * Ls, ir); CL <- rowsum(e$wn[ee] * Rr, js)
        Vs[rows, v] <- (colSums(XR^2 * RB$V[as.integer(rownames(XR)), tests$rc[rows], drop = FALSE]) +
                          colSums(CL^2 * LB$V[as.integer(rownames(CL)), tests$lc[rows], drop = FALSE])) / max(nv[Bt], 1)^2
      }
    }
  }
  list(S = S, V = .excess(Vs))
}

# paracrine / niche: blurred ligand fields of each sender type, read at receivers (sender layer shifted by whole
# bins, wrapped with the grid's period)
.scores_grid <- function(ctx, sigma, bin, chunk, grid_mem, z_cand) {
  coords <- ctx$coords; tcode <- ctx$tcode; K <- ctx$K; tests <- ctx$tests; sh <- ctx$sh; n <- nrow(coords)
  gr <- .grid(coords, ctx$box, sigma, bin); nc <- gr$nx * gr$ny
  if (is.null(chunk)) chunk <- max(1L, min(64L, floor(grid_mem / nc)))
  sxy <- cbind(round(sh[, 1] / bin), round(sh[, 2] / bin))
  linv <- function(idx, v) ((gr$gx[idx] - sxy[v, 1]) %% gr$nx) + ((gr$gy[idx] - sxy[v, 2]) %% gr$ny) * gr$nx + 1
  rs <- ctx$rs
  Tf <- .blur(.bin_sum(matrix(1, n, 1), gr$lin, nc), gr, rs)[, 1]
  tau <- 0.1 * stats::median(Tf[gr$lin])
  nvK <- vapply(seq_len(nrow(sh)), function(v) tabulate(tcode[Tf[linv(seq_len(n), v)] >= tau], K), numeric(K))
  nvK <- matrix(nvK, K)
  nvar <- min(NVAR, nrow(sh)); S <- matrix(0, nrow(tests), nrow(sh)); Vs <- matrix(0, nrow(tests), nvar)
  for (A in unique(tests$A)) {
    ta <- which(tests$A == A); la <- unique(tests$lc[ta]); sa <- ctx$rng[[A]]; LB <- ctx$Lb[[A]]
    for (ch in split(la, ceiling(seq_along(la) / chunk))) {
      Fq <- .blur(.bin_sum(LB$M[, ch, drop = FALSE], gr$lin[sa], nc), gr, rs) / pmax(Tf, EPS)
      Fq[Tf < tau, ] <- 0
      rows_c <- ta[tests$lc[ta] %in% ch]
      for (Bt in unique(tests$Bt[rows_c])) {
        rows <- rows_c[tests$Bt[rows_c] == Bt]; cb <- ctx$rng[[Bt]]; RB <- ctx$Rb[[Bt]]
        Rm <- RB$M[, tests$rc[rows], drop = FALSE]; VRB <- RB$V[, tests$rc[rows], drop = FALSE]
        lp <- match(tests$lc[rows], ch)
        if (!is.null(rs)) {
          o <- lrmaprs::rs_grid_gather(Fq, as.integer(gr$gx[cb]), as.integer(gr$gy[cb]), as.integer(sxy[, 1]),
                                       as.integer(sxy[, 2]), as.integer(gr$nx), as.integer(gr$ny), Rm, VRB,
                                       as.integer(lp), as.numeric(nvK[Bt, ]), as.integer(nvar), rs)
          S[rows, ] <- o$S; Vs[rows, ] <- Vs[rows, ] + o$Vs
          next
        }
        for (v in seq_len(nrow(sh))) {
          lin <- linv(cb, v)
          Xv <- Fq[lin, lp, drop = FALSE]                                    # 0 where invalid (Tf < tau)
          S[rows, v] <- colSums(Rm * Xv) / max(nvK[Bt, v], 1)
          if (v <= nvar) Vs[rows, v] <- Vs[rows, v] + colSums(Xv^2 * VRB) / max(nvK[Bt, v], 1)^2   # receptor side
        }
      }
    }
  }
  # ligand side: dS/dL_j = sum_i R_i K(x_i - shift, x_j) / T(x_i - shift): a blurred R / T field read at senders.
  # Added variance only lowers z, so it is computed for candidate tests (z before it > z_cand) only.
  nul <- S[, -1, drop = FALSE]; mu <- rowMeans(nul)
  sdv <- sqrt(rowSums((nul - mu)^2) / max(ncol(nul) - 1, 1) + .excess(Vs))
  cand <- which((S[, 1] - mu) / pmax(sdv, EPS) > z_cand)
  for (Bt in unique(tests$Bt[cand])) {
    tb <- cand[tests$Bt[cand] == Bt]; rb <- unique(tests$rc[tb]); cb <- ctx$rng[[Bt]]; RB <- ctx$Rb[[Bt]]
    for (v in seq_len(nvar)) {
      lin <- linv(cb, v); ok <- Tf[lin] >= tau; if (!any(ok)) next
      for (ch in split(rb, ceiling(seq_along(rb) / chunk))) {
        G <- .blur(.bin_sum(RB$M[ok, ch, drop = FALSE] / Tf[lin[ok]], lin[ok], nc), gr, rs)
        rows_c <- tb[tests$rc[tb] %in% ch]
        for (A in unique(tests$A[rows_c])) {
          rows <- rows_c[tests$A[rows_c] == A]; sa <- ctx$rng[[A]]
          Cj <- G[gr$lin[sa], match(tests$rc[rows], ch), drop = FALSE]
          Vs[rows, v] <- Vs[rows, v] + colSums(Cj^2 * ctx$Lb[[A]]$V[, tests$lc[rows], drop = FALSE]) /
            max(nvK[Bt, v], 1)^2
        }
      }
    }
  }
  list(S = S, V = .excess(Vs))
}

#' Significant ligand-receptor co-location
#'
#' @param x an `LRmap` object. @param scale scales to keep.
#' @param test which null: "type" (pair co-location beyond the co-location of the two cell types, calibrated by
#'   the shifts; `q_type`, `fold_type`), "shuffle" (ligand and receptor values shuffled within their cell
#'   types, layout fixed; `q_shuffle`, `fold_type`) or "colocation" (shifted sender layer: co-location of
#'   ligand-expressing senders and receptor-expressing receivers, including the co-location of the two types;
#'   `q`, `fold`). The expression scale always uses `q` and `fold`.
#' @param q_thr q-value threshold. @param min_fold minimum fold (see `test`).
#' @return rows of `x$table` passing the thresholds, by q.
#' @export
lrmap_results <- function(x, scale = c("contact", "paracrine", "niche", "expression"), test = c("type", "shuffle", "colocation"),
                          q_thr = 0.05, min_fold = 1.25) {
  test <- match.arg(test); tb <- x$table
  keep <- which(tb$scale %in% scale & .called(tb, test, q_thr, min_fold))
  r <- tb[keep, ]; r[order(r$q, -r$z), ]
}

.called <- function(tb, test, q_thr, min_fold) {
  ex <- tb$scale == "expression"
  qq <- switch(test, type = tb$q_type, shuffle = tb$q_shuffle, colocation = tb$q)
  ff <- if (test == "colocation") tb$fold else tb$fold_type
  qq[ex] <- tb$q[ex]; ff[ex] <- tb$fold[ex]
  !is.na(qq) & qq < q_thr & !is.na(ff) & ff >= min_fold
}

#' @export
print.LRmap <- function(x, ...) {
  cat(sprintf("LRmap: %d cell types, %d of %d ligand-receptor pairs covered, %d shifts%s\n", length(x$types),
              x$n_covered, x$n_pairs, x$settings$B, if (x$settings$leak) ", leak removed" else ""))
  for (s in x$settings$scales) {
    a <- x$table[x$table$scale == s, ]; cl <- x$coloc[x$coloc$scale == s, ]
    cat(sprintf("  %-10s %6d tests, %5d called (beyond type co-location, q < 0.05); type pairs co-located: %d of %d\n",
                s, nrow(a), nrow(lrmap_results(x, s, "type")), sum(cl$q < 0.05, na.rm = TRUE), nrow(cl)))
  }
  invisible(x)
}

#' Expression-only calls versus spatial co-location
#'
#' For a fit with the "expression" scale and at least one spatial scale: one row per tested sender, receiver and
#' ligand-receptor pair, with whether it is called on expression alone and at which spatial scales it is
#' co-located, and whether the two cell types co-locate at all (`coloc`).
#'
#' @param x an `LRmap` object fitted with `scales` including "expression".
#' @param test,q_thr,min_fold call rule for the spatial scales (as in [lrmap_results()]).
#' @return data.frame: sender, receiver, ligand, receptor, `expr_call`, `expr_z`, `spatial_scales` (called
#'   scales joined by "+", "" if none), `spatial_call`, `types_colocated`, and `class`: "both",
#'   "expression only", "spatial only" or "neither".
#' @export
lrmap_compare <- function(x, test = c("type", "shuffle", "colocation"), q_thr = 0.05, min_fold = 1.25) {
  test <- match.arg(test)
  tb <- x$table; sp <- setdiff(x$settings$scales, "expression")
  if (!"expression" %in% x$settings$scales || !length(sp)) stop("x needs the 'expression' scale and a spatial scale")
  key <- function(d) paste(d$sender, d$receiver, d$ligand, d$receptor, sep = "\r")
  call <- .called(tb, test, q_thr, min_fold)
  ex <- tb[tb$scale == "expression", ]; k <- key(ex)
  out <- data.frame(ex[, c("sender", "receiver", "ligand", "receptor")], expr_call = call[tb$scale == "expression"],
                    expr_z = ex$z, stringsAsFactors = FALSE)
  ss <- vapply(k, function(z) "", ""); names(ss) <- NULL
  for (s in sp) {
    w <- tb$scale == s & call; hit <- k %in% key(tb[w, ])
    ss[hit] <- ifelse(ss[hit] == "", s, paste(ss[hit], s, sep = "+"))
  }
  out$spatial_scales <- ss; out$spatial_call <- ss != ""
  cl <- x$coloc[x$coloc$scale %in% sp & !is.na(x$coloc$q) & x$coloc$q < q_thr, ]
  out$types_colocated <- paste(out$sender, out$receiver) %in% paste(cl$sender, cl$receiver)
  out$class <- ifelse(out$expr_call & out$spatial_call, "both", ifelse(out$expr_call, "expression only",
                      ifelse(out$spatial_call, "spatial only", "neither")))
  out
}
