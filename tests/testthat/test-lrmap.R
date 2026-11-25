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
  lw <- lrmap:::.leak_W(s$coords, tcs, 8, lrmap:::.leak_mats(s$leak, sort(unique(s$types))))
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
