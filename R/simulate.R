#' Simulate a section with ligand-receptor co-location scenarios
#'
#' Cells on a jittered hexagonal lattice. Type A forms round clusters; the layout places type B:
#' * `touch`: a ring of B directly around each A cluster (A and B touch);
#' * `gap`: a ring of C (`gap` microns) around each A cluster, then a 50 micron ring of B;
#' * `far`: A clusters in the left fifth, B clusters in the right fifth;
#' * `indep`: A, B and C clusters at independent random positions.
#' Other cells are D. Gene `L` is expressed in A and gene `R` in B. With `leak_R`, `R` is expressed in A instead
#' (not at all in B), and `leak` > 0 moves transcripts from touching cells of other types into each cell, so B cells
#' next to A pick up `R`. Genes `G1..` have random type profiles and smooth spatial modulation; `lr` pairs them at
#' random as decoys.
#'
#' @param layout one of "touch", "gap", "far", "indep".
#' @param size side of the square section (microns). @param spacing lattice spacing (microns).
#' @param n_clusters A clusters. @param r_a A cluster radius (microns).
#' @param mu_l,mu_r mean counts of L in A and R in B (at the mean library).
#' @param leak fraction of each touching neighbour's expression leaked into a cell (other types only).
#' @param leak_R express R in A instead of B (leak scenario).
#' @param gap width of the C ring in the `gap` layout (microns).
#' @param n_genes decoy genes. @param n_decoy decoy pairs. @param lib_mean mean library size.
#' @param seed random seed.
#' @return list(counts (genes x cells), coords, types, lib, lr, leak (table of the true coefficients)).
#' @export
lrmap_simulate <- function(layout = c("touch", "gap", "far", "indep"), size = 1500, spacing = 15, n_clusters = 20,
                           r_a = 40, mu_l = 4, mu_r = 4, leak = 0, leak_R = FALSE, n_genes = 30, n_decoy = 40, gap = 35,
                           lib_mean = 300, seed = 1) {
  layout <- match.arg(layout); set.seed(seed)
  ys <- seq(0, size, by = spacing * sqrt(3) / 2)
  xy <- do.call(rbind, lapply(seq_along(ys), function(r) cbind(seq(if (r %% 2) 0 else spacing / 2, size, by = spacing), ys[r])))
  xy <- xy + matrix(stats::runif(length(xy), -0.2, 0.2) * spacing, ncol = 2); n <- nrow(xy)
  centres <- function(m, xlim, ylim, sep) {
    out <- matrix(0, 0, 2); tries <- 0
    while (nrow(out) < m && tries < 5000) {
      tries <- tries + 1; p <- c(stats::runif(1, xlim[1], xlim[2]), stats::runif(1, ylim[1], ylim[2]))
      if (!nrow(out) || min(sqrt(colSums((t(out) - p)^2))) >= sep) out <- rbind(out, p)
    }
    out
  }
  dmin <- function(C) { if (!nrow(C)) return(rep(Inf, n)); nn <- RANN::nn2(C, xy, k = 1); nn$nn.dists[, 1] }
  ty <- rep("D", n); m <- r_a + 100
  if (layout %in% c("touch", "gap")) {
    ca <- centres(n_clusters, c(m, size - m), c(m, size - m), 2 * (r_a + max(gap, 0) + 60)); d <- dmin(ca)
    if (layout == "touch") { ty[d < r_a + 30] <- "B" } else { ty[d < r_a + gap + 50] <- "B"; ty[d < r_a + gap] <- "C" }
    ty[d < r_a] <- "A"
    cc <- centres(4, c(m, size - m), c(m, size - m), 0); if (layout == "touch") ty[dmin(cc) < r_a & ty == "D"] <- "C"
  } else if (layout == "far") {
    ca <- centres(ceiling(n_clusters / 2), c(r_a + 20, 0.2 * size), c(r_a + 20, size - r_a - 20), 2 * r_a + 30)
    cb <- centres(ceiling(n_clusters / 2), c(0.8 * size, size - r_a - 20), c(r_a + 20, size - r_a - 20), 2 * r_a + 30)
    cc <- centres(4, c(0.35 * size, 0.65 * size), c(r_a, size - r_a), 2 * r_a)
    ty[dmin(cc) < r_a] <- "C"; ty[dmin(cb) < r_a] <- "B"; ty[dmin(ca) < r_a] <- "A"
  } else {
    for (tt in c("C", "B", "A")) ty[dmin(centres(n_clusters, c(0, size), c(0, size), 2 * r_a)) < r_a] <- tt
  }
  lev <- c("A", "B", "C", "D"); tc <- match(ty, lev)
  lib <- round(lib_mean * exp(stats::rnorm(n, 0, 0.3)))
  # smooth, stationary spatial modulation per gene: random Gaussian bumps placed beyond the section edges too
  smooth <- function() {
    nb <- 30; b <- cbind(stats::runif(nb, -450, size + 450), stats::runif(nb, -450, size + 450))
    f <- rowSums(sapply(seq_len(nb), function(j) exp(-0.5 * ((xy[, 1] - b[j, 1])^2 + (xy[, 2] - b[j, 2])^2) / 150^2)))
    exp(0.6 * (f - mean(f)) / max(stats::sd(f), 1e-9))
  }
  prof <- rbind(L = c(mu_l, 0.05, 0.05, 0.05), R = if (leak_R) c(4 * mu_r, 0, 0.05, 0.05) else c(0.05, mu_r, 0.05, 0.05))
  G <- matrix(stats::rlnorm(n_genes * 4, 0, 1) * (stats::runif(n_genes * 4) < 0.6), n_genes, 4,
              dimnames = list(paste0("G", seq_len(n_genes)), NULL))
  prof <- rbind(prof, G)
  E <- t(prof[, tc, drop = FALSE]) / lib_mean                                     # cells x genes, fraction of library
  E[, -(1:2)] <- E[, -(1:2)] * vapply(seq_len(n_genes), function(g) smooth(), numeric(n))
  if (leak > 0) {
    nn <- RANN::nn2(xy, k = 9); recv <- rep(seq_len(n), 8); send <- as.vector(nn$nn.idx[, -1])
    d <- as.vector(nn$nn.dists[, -1]); w <- exp(-d / stats::median(d)); wn <- w / .rsum(w, recv, n)[recv]
    het <- tc[recv] != tc[send]
    Wl <- Matrix::sparseMatrix(i = recv, j = send, x = wn * leak * het, dims = c(n, n))
    E <- E + as.matrix(Wl %*% E)
  }
  counts <- matrix(stats::rpois(length(E), E * lib), n)
  counts <- Matrix::Matrix(t(counts), sparse = TRUE, dimnames = list(rownames(prof), NULL))
  dec <- data.frame(ligand = sample(rownames(G), n_decoy, TRUE), receptor = sample(rownames(G), n_decoy, TRUE))
  dec <- dec[dec$ligand != dec$receptor, ]
  lr <- rbind(data.frame(ligand = "L", receptor = "R"), dec)
  lr$true <- c(TRUE, rep(FALSE, nrow(dec)))
  lk <- expand.grid(receiver = lev, sender = lev, stringsAsFactors = FALSE)
  lk <- lk[lk$receiver != lk$sender, ]; lk$lambda <- leak; lk$mu <- 0
  list(counts = counts, coords = xy, types = ty, lib = lib, lr = lr[!duplicated(lr[, 1:2]), ], leak = lk)
}
