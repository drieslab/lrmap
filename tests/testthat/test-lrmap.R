sim <- lrmap_simulate("touch", size = 600, n_clusters = 4, seed = 3)
Y <- t(as.matrix(sim$counts)) * (100 / sim$lib); tc <- match(sim$types, sort(unique(sim$types)))
box <- c(min(sim$coords[, 1]), diff(range(sim$coords[, 1])) + 1e-6, min(sim$coords[, 2]), diff(range(sim$coords[, 2])) + 1e-6)

test_that("the lr table: subunits, coverage, duplicates", {
  tb <- lrmap:::.lr_table(data.frame(ligand = c("A", "A", "B", "C"), receptor = c("R", "R", "R1_R2", "X"),
                                     receptor_subunits = c(NA, NA, "R1_R2", NA)), c("A", "B", "R", "R1", "R2"))
  expect_equal(tb$n, 3); expect_equal(tb$n_covered, 2); expect_equal(tb$table$rg[[2]], c("R1", "R2"))
})

test_that("contact score equals a direct sum over k nearest neighbours", {
  tests <- data.frame(A = 1L, Bt = 2L, p = 1L, li = 1L, ri = 2L)
  S <- lrmap:::.scores_contact(sim$coords, tc, 4, Y, Y, tests, rbind(c(0, 0), c(200, 100)), box, 8)
  nn <- RANN::nn2(sim$coords, k = 9); d0 <- median(nn$nn.dists[, -1]); cap <- 3 * d0
  b <- which(tc == 2); num <- 0; nv <- 0
  for (i in b) {
    j <- nn$nn.idx[i, -1]; d <- nn$nn.dists[i, -1]; k <- d <= cap; if (!any(k)) next
    w <- exp(-d[k] / d0); nv <- nv + 1
    num <- num + Y[i, 2] * sum(w * (tc[j[k]] == 1) * Y[j[k], 1]) / sum(w)
  }
  expect_equal(S[1, 1], unname(num / nv), tolerance = 1e-10)
})

test_that("the grid blur is the separable Gaussian convolution, and shifts move the sender layer", {
  gr <- lrmap:::.grid(sim$coords, box, sigma = 20, bin = 10)
  M <- matrix(0, gr$nx * gr$ny, 1); M[sample(nrow(M), 30), 1] <- 1
  Bm <- matrix(lrmap:::.blur(M, gr), gr$nx); M2 <- matrix(M, gr$nx)
  i0 <- 7; j0 <- 9; r <- ceiling(3 * 20 / 10)
  ref <- 0
  for (a in max(1, i0 - r):min(gr$nx, i0 + r)) for (b in max(1, j0 - r):min(gr$ny, j0 + r))
    ref <- ref + M2[a, b] * exp(-0.5 * ((a - i0)^2 + (b - j0)^2) * 100 / 400)
  expect_equal(Bm[i0, j0], ref, tolerance = 1e-10)
  sh <- c(137, 251); a <- which(tc == 1); g5 <- lrmap:::.grid(sim$coords, box, 20, 5); nc <- g5$nx * g5$ny
  lin <- lrmap:::.shift_lin(g5, sh, sim$coords, box)
  F_read <- lrmap:::.blur(lrmap:::.bin_sum(matrix(1, length(a), 1), g5$lin[a], nc), g5)[lin, 1]     # original field, shifted read
  moved <- lrmap:::.shift_coords(sim$coords, sh, box); gm <- lrmap:::.grid(moved, box, 20, 5)
  F_move <- lrmap:::.blur(lrmap:::.bin_sum(matrix(1, length(a), 1), gm$lin[a], nc), g5)[g5$lin, 1]  # moved senders
  back <- lrmap:::.shift_coords(sim$coords, -sh, box)            # away from the wrap seam and the box edges
  inner <- function(m) m[, 1] > box[1] + 60 & m[, 1] < box[1] + box[2] - 60 & m[, 2] > box[3] + 60 & m[, 2] < box[3] + box[4] - 60
  ok <- inner(back) & inner(sim$coords)
  expect_gt(cor(F_read[ok], F_move[ok]), 0.99)
  expect_equal(lrmap:::.shift_lin(gr, c(0, 0), sim$coords, box), gr$lin)
})

test_that("leak removal recovers the expression before leak (first order)", {
  s <- lrmap_simulate("touch", size = 600, n_clusters = 4, leak = 0.05, leak_R = TRUE, seed = 4)
  tcs <- match(s$types, sort(unique(s$types)))
  lm <- lrmap:::.leak_mats(s$leak, sort(unique(s$types)))
  E <- matrix(ifelse(tcs == 1, 16, 0) / 3, ncol = 1)               # R before leak (composition scale)
  nn <- RANN::nn2(s$coords, k = 9); n <- nrow(E); recv <- rep(seq_len(n), 8); send <- as.vector(nn$nn.idx[, -1])
  d <- as.vector(nn$nn.dists[, -1]); w <- exp(-d / median(d)); wn <- w / lrmap:::.rsum(w, recv, n)[recv]
  Wl <- Matrix::sparseMatrix(i = recv, j = send, x = wn * 0.05 * (tcs[recv] != tcs[send]), dims = c(n, n))
  Yl <- E + as.matrix(Wl %*% E)
  Yc <- lrmap:::.unleak(Yl, tcs, 4, s$coords, 8, lm)
  b <- tcs == 2
  expect_lt(max(abs(Yc[b] - E[b])), 0.01 * max(Yl[b]))            # second-order remainder only
})

test_that("results: shifts, z and the type-level table", {
  f <- lrmap(sim$counts, sim$coords, sim$types, sim$lr, lib = sim$lib, B = 40, min_shift = 150, verbose = FALSE)
  tp <- f$table[f$table$ligand == "L" & f$table$sender == "A" & f$table$receiver == "B" & f$table$scale == "contact", ]
  expect_gt(tp$z, 2)
  expect_true(all(lrmap_results(f, "niche")$scale == "niche"))
  expect_equal(nrow(f$coloc), 3 * 12)
  f2 <- lrmap(sim$counts, sim$coords, sim$types, sim$lr, lib = sim$lib, B = 40, min_shift = 150, verbose = FALSE)
  expect_equal(f$table, f2$table)
})
