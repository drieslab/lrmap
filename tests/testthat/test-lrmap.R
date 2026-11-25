# small simulated sections use short shifts on purpose; the min_shift warning is tested on its own below
lrmap <- function(...) suppressWarnings(lrmap::lrmap(...))

sim <- lrmap_simulate("touch", size = 600, n_clusters = 4, seed = 3)
Y <- t(as.matrix(sim$counts)) * (100 / sim$lib); tc <- match(sim$types, sort(unique(sim$types)))
box <- c(min(sim$coords[, 1]), diff(range(sim$coords[, 1])) + 1e-6, min(sim$coords[, 2]), diff(range(sim$coords[, 2])) + 1e-6)
pairLR <- data.frame(ligand = "L", receptor = "R")
fit1 <- function(...) lrmap(sim$counts, sim$coords, sim$types, pairLR, lib = sim$lib, B = 20, min_shift = 150, verbose = FALSE, ...)

test_that("the lr table: subunits, coverage, duplicates", {
  tb <- lrmap:::.lr_table(data.frame(ligand = c("A", "A", "B", "C"), receptor = c("R", "R", "R1_R2", "X"),
                                     receptor_subunits = c(NA, NA, "R1_R2", NA)), c("A", "B", "R", "R1", "R2"))
  expect_equal(tb$n, 3); expect_equal(tb$n_covered, 2); expect_equal(tb$table$rg[[2]], c("R1", "R2"))
})

test_that("complex entities: minimum for ligands, mean for receptors, detection of all subunits", {
  M <- Matrix::Matrix(cbind(a = c(0, 2, 3), b = c(1, 0, 5)), sparse = TRUE)
  l <- lrmap:::.entities_sparse("x", list(c("a", "b")), M, "min"); r <- lrmap:::.entities_sparse("x", list(c("a", "b")), M, "mean")
  expect_equal(as.vector(l$E), c(0, 0, 3)); expect_equal(as.vector(r$E), c(0.5, 1, 4)); expect_equal(as.vector(l$D), c(0, 0, 1))
})

test_that("contact score equals a direct sum over k nearest neighbours", {
  f <- fit1(scales = "contact")
  nn <- RANN::nn2(sim$coords, k = 9); d0 <- median(nn$nn.dists[, -1]); cap <- 3 * d0
  b <- which(tc == 2); num <- 0; nv <- 0
  for (i in b) {
    j <- nn$nn.idx[i, -1]; d <- nn$nn.dists[i, -1]; k <- d <= cap; if (!any(k)) next
    w <- exp(-d[k] / d0); nv <- nv + 1
    num <- num + Y[i, "R"] * sum(w * (tc[j[k]] == 1) * Y[j[k], "L"]) / sum(w)
  }
  expect_equal(f$table$score[f$table$sender == "A" & f$table$receiver == "B"], unname(num / nv), tolerance = 1e-10)
})

test_that("grid score equals the blurred ligand share read at receivers, and the blur is a Gaussian", {
  f <- fit1(scales = "paracrine")
  gr <- lrmap:::.grid(sim$coords, box, 25, 8); nc <- gr$nx * gr$ny
  Tf <- lrmap:::.blur(lrmap:::.bin_sum(matrix(1, nrow(Y), 1), gr$lin, nc), gr)[, 1]; tau <- 0.1 * median(Tf[gr$lin])
  a <- tc == 1; Fa <- lrmap:::.blur(lrmap:::.bin_sum(matrix(Y[a, "L"]), gr$lin[a], nc), gr)[, 1]
  b <- tc == 2; ok <- Tf[gr$lin[b]] >= tau
  ref <- sum((Y[b, "R"] * Fa[gr$lin[b]] / Tf[gr$lin[b]])[ok]) / sum(ok)
  expect_equal(f$table$score[f$table$sender == "A" & f$table$receiver == "B"], ref, tolerance = 1e-10)
  M <- matrix(0, nc, 1); M[sample(nc, 30), 1] <- 1; Bm <- matrix(lrmap:::.blur(M, gr), gr$nx); M2 <- matrix(M, gr$nx)
  i0 <- 7; j0 <- 9; r <- ceiling(3 * 25 / 8); refb <- 0
  for (x in max(1, i0 - r):min(gr$nx, i0 + r)) for (y in max(1, j0 - r):min(gr$ny, j0 + r))
    refb <- refb + M2[x, y] * exp(-0.5 * ((x - i0)^2 + (y - j0)^2) * 64 / 625)
  expect_equal(Bm[i0, j0], refb, tolerance = 1e-10)
})

test_that("the leak operator removes planted leak to first order", {
  s <- lrmap_simulate("touch", size = 600, n_clusters = 4, leak = 0.05, leak_R = TRUE, seed = 4)
  tcs <- match(s$types, sort(unique(s$types)))
  lw <- lrmap:::.leak_W(lrmap:::.lk_edges(s$coords, 8), tcs, lrmap:::.leak_mats(s$leak, sort(unique(s$types))))
  E <- matrix(ifelse(tcs == 1, 16, 0) / 3, ncol = 1); Yl <- E + as.matrix(lw$W %*% E)
  Yc <- Yl - as.matrix(lw$W %*% Yl); b <- tcs == 2
  expect_lt(max(abs(Yc[b] - E[b])), 0.01 * max(Yl[b]))
})

test_that("cell order, reader input and repeat runs give the same results", {
  f <- fit1(); p <- sample(ncol(sim$counts))
  f2 <- lrmap(sim$counts[, p], sim$coords[p, ], sim$types[p], pairLR, lib = sim$lib[p], B = 20, min_shift = 150, verbose = FALSE)
  rd <- function(genes) sim$counts[intersect(genes, rownames(sim$counts)), , drop = FALSE]
  f3 <- lrmap(rd, sim$coords, sim$types, pairLR, lib = sim$lib, B = 20, min_shift = 150, verbose = FALSE)
  expect_equal(f$table, f2$table); expect_equal(f$table, f3$table); expect_equal(f$coloc, f3$coloc)
  tp <- f$table[f$table$sender == "A" & f$table$receiver == "B" & f$table$scale == "contact", ]
  expect_gt(tp$z, 2); expect_equal(nrow(f$coloc), 3 * 12)
  expect_true(all(lrmap_results(f, "niche")$scale == "niche"))
})

test_that("the Rust backend equals the R backend (scores, null, variance), with and without leak", {
  skip_if_not_installed("lrmaprs")
  s <- lrmap_simulate("touch", size = 600, n_clusters = 4, leak = 0.05, leak_R = TRUE, seed = 5)
  for (lk in list(NULL, s$leak)) {
    a <- lrmap(s$counts, s$coords, s$types, s$lr, lib = s$lib, B = 25, min_shift = 150, verbose = FALSE, leak = lk,
               backend = "r", z_cand = -Inf)
    b <- lrmap(s$counts, s$coords, s$types, s$lr, lib = s$lib, B = 25, min_shift = 150, verbose = FALSE, leak = lk,
               backend = "rust", z_cand = -Inf)
    c2 <- lrmap(s$counts, s$coords, s$types, s$lr, lib = s$lib, B = 25, min_shift = 150, verbose = FALSE, leak = lk,
                backend = "rust", threads = 2, z_cand = -Inf)
    expect_equal(b$table, a$table, tolerance = 1e-9); expect_equal(b$coloc, a$coloc, tolerance = 1e-9)
    expect_equal(c2$table, b$table)
  }
})

test_that("lrmap_leak recovers planted leak and matches edgeicg's estimate", {
  s <- lrmap_simulate("touch", size = 900, n_clusters = 8, leak = 0.08, n_genes = 150, seed = 6)
  lk <- lrmap_leak(s$counts, s$coords, s$types, s$lib, min_cells = 30, verbose = FALSE)
  med <- stats::median(lk$lambda_raw, na.rm = TRUE)               # planted 0.08 between every touching type pair
  expect_true(med > 0.04 && med < 0.12)
  skip_if_not_installed("edgeicg")
  f <- edgeicg::edgeicg(s$counts, s$coords, s$types, lib = s$lib, min_cells = 30L, calibrate = FALSE, attribute = FALSE,
                        backend = "r", verbose = FALSE)
  P <- f$pairs; ref <- data.frame(receiver = f$types[P$receiver], sender = f$types[P$sender], lambda = P$lambda, mu = P$mu)
  m <- merge(lk, ref, by = c("receiver", "sender"))
  expect_equal(nrow(m), nrow(ref))
  expect_equal(m$lambda_raw, m$lambda.y, tolerance = 1e-6); expect_equal(m$mu_raw, m$mu.y, tolerance = 1e-6)
})

test_that("leak removal: spill_edges equal to the default graph give the same result; cell order does not matter", {
  s <- lrmap_simulate("touch", size = 600, n_clusters = 4, leak = 0.05, leak_R = TRUE, seed = 5)
  run <- function(...) lrmap(s$counts, s$coords, s$types, s$lr, lib = s$lib, B = 15, min_shift = 150, leak = s$leak,
                             backend = "r", verbose = FALSE, ...)
  a <- run(); e <- lrmap:::.lk_edges(s$coords, 8)
  b <- run(spill_edges = data.frame(recv = e$recv, send = e$send, w = e$w))
  expect_equal(a$table, b$table)
  p <- sample(ncol(s$counts)); inv <- order(p)
  c2 <- lrmap(s$counts[, p], s$coords[p, ], s$types[p], s$lr, lib = s$lib[p], B = 15, min_shift = 150, leak = s$leak,
              backend = "r", verbose = FALSE, spill_edges = data.frame(recv = inv[e$recv], send = inv[e$send], w = e$w))
  expect_equal(a$table, c2$table)
})

test_that("the expression-only null matches label shuffles; far-apart types are called on expression only", {
  s <- lrmap_simulate("far", size = 900, n_clusters = 10, seed = 8)
  f <- lrmap(s$counts, s$coords, s$types, s$lr, lib = s$lib, B = 20, min_shift = 150, backend = "r", verbose = FALSE,
             scales = c("contact", "niche", "expression"))
  ex <- f$table[f$table$scale == "expression", ]
  # Monte Carlo: shuffle the type labels, recompute mean L in A x mean R in B
  y <- t(as.matrix(s$counts[c("L", "R"), ])) * 100 / s$lib
  r1 <- ex[ex$ligand == "L" & ex$receptor == "R" & ex$sender == "A" & ex$receiver == "B", ]
  r2 <- ex[ex$ligand != "L", ][1, ]
  set.seed(3)
  for (r in list(r1, r2)) {
    yy <- if (r$ligand == "L") y else t(as.matrix(s$counts[c(r$ligand, r$receptor), ])) * 100 / s$lib
    perm <- replicate(2000, { tp <- sample(s$types); mean(yy[tp == r$sender, 1]) * mean(yy[tp == r$receiver, 2]) })
    expect_equal(r$null_mean, mean(perm), tolerance = 0.03)
    zmc <- (r$score - mean(perm)) / stats::sd(perm)
    expect_lt(abs(r$z - zmc), 0.1 * max(1, abs(zmc)))
  }
  cmp <- lrmap_compare(f)
  lr1 <- cmp[cmp$sender == "A" & cmp$receiver == "B" & cmp$ligand == "L" & cmp$receptor == "R", ]
  expect_equal(lr1$class, "expression only")
  expect_false(lr1$types_colocated)
})

test_that("min_shift shorter than the widest kernel reach warns", {
  s <- lrmap_simulate("touch", size = 600, n_clusters = 4, seed = 2)
  expect_warning(lrmap::lrmap(s$counts, s$coords, s$types, s$lr, lib = s$lib, B = 5, min_shift = 150,
                              scales = "niche", verbose = FALSE), "min_shift")
  expect_no_warning(lrmap::lrmap(s$counts, s$coords, s$types, s$lr, lib = s$lib, B = 5, min_shift = 200,
                                 scales = c("contact", "paracrine"), verbose = FALSE))
})

test_that("a cell type with a single cell does not break leak removal or the expression scale", {
  s <- lrmap_simulate("touch", size = 600, n_clusters = 4, leak = 0.05, seed = 9)
  s$types[1] <- "E"; s$lr <- rbind(s$lr[, 1:2], data.frame(ligand = "L", receptor = "G1_G2"))
  s$lr$receptor_subunits <- ifelse(s$lr$receptor == "G1_G2", "G1_G2", s$lr$receptor)
  expect_no_error(lrmap(s$counts, s$coords, s$types, s$lr, lib = s$lib, B = 5, min_shift = 150, leak = s$leak,
                        scales = c("contact", "expression"), verbose = FALSE))
})

test_that("fold is NA, not huge, when the expectation is not positive", {
  expect_equal(lrmap:::.ratio(c(1, -2, 3), c(2, 0, -1)), c(0.5, NA, NA))
})

test_that("leak pairs with little contact or out-of-range estimates are not corrected", {
  s <- lrmap_simulate("touch", size = 900, n_clusters = 8, leak = 0.08, n_genes = 150, seed = 6)
  lk <- lrmap_leak(s$counts, s$coords, s$types, s$lib, min_cells = 30, verbose = FALSE)
  expect_true(all(is.na(lk$lambda[lk$status != "ok"])))
  expect_true(all(lk$support[lk$status == "ok"] >= 30))
  lk2 <- lrmap_leak(s$counts, s$coords, s$types, s$lib, min_cells = 30, min_support = 1e9, verbose = FALSE)
  expect_true(all(lk2$status %in% c("low support", "no contact")))
  # lrmap ignores supplied coefficients that are out of range: same result as without them
  bad <- s$leak; bad$lambda <- 5; bad$mu <- 0
  f0 <- lrmap(s$counts, s$coords, s$types, s$lr, lib = s$lib, B = 5, min_shift = 150, scales = "contact", verbose = FALSE)
  f1 <- lrmap(s$counts, s$coords, s$types, s$lr, lib = s$lib, B = 5, min_shift = 150, scales = "contact", leak = bad, verbose = FALSE)
  expect_equal(f0$table$score, f1$table$score)
})
