EPS <- 1e-12

# Contact neighbours: k nearest cells of the sender layer for each receiver, weight exp(-d / d0).
# Observed layout (shift = NULL): sender layer = receivers themselves, self excluded.
# Edges farther than `cap` are dropped (relevant for shifted layers that leave the tissue).
.contact_edges <- function(coords, k, d0, cap, shift = NULL, box = NULL) {
  n <- nrow(coords)
  if (is.null(shift)) {
    nn <- RANN::nn2(coords, k = k + 1); idx <- nn$nn.idx[, -1, drop = FALSE]; d <- nn$nn.dists[, -1, drop = FALSE]
  } else {
    sc <- .shift_coords(coords, shift, box)
    nn <- RANN::nn2(sc, coords, k = k); idx <- nn$nn.idx; d <- nn$nn.dists
  }
  recv <- rep(seq_len(n), ncol(idx)); send <- as.vector(idx); d <- as.vector(d)
  keep <- d <= cap; recv <- recv[keep]; send <- send[keep]; w <- exp(-d[keep] / d0)
  W <- .rsum(w, recv, n)
  list(recv = recv, send = send, wn = w / pmax(W[recv], EPS), valid = W > 0)
}

.rsum <- function(x, g, n) { s <- numeric(n); r <- rowsum(x, g); s[as.integer(rownames(r))] <- r[, 1]; s }

# toroidal shift of points within the bounding box
.shift_coords <- function(coords, shift, box) {
  cbind(box[1] + (coords[, 1] - box[1] + shift[1]) %% box[2], box[3] + (coords[, 2] - box[3] + shift[2]) %% box[4])
}

# shifts (microns) with toroidal length >= min_shift; the first row is the observed layout (0, 0)
.shifts <- function(B, box, min_shift, seed) {
  set.seed(seed); out <- matrix(0, 0, 2); tries <- 0
  while (nrow(out) < B && tries < 100) {
    tries <- tries + 1
    s <- cbind(stats::runif(4 * B, 0, box[2]), stats::runif(4 * B, 0, box[4]))
    dx <- pmin(s[, 1], box[2] - s[, 1]); dy <- pmin(s[, 2], box[4] - s[, 2])
    out <- rbind(out, s[sqrt(dx^2 + dy^2) >= min_shift, , drop = FALSE])
  }
  if (nrow(out) < B) stop("no shifts of at least min_shift fit in the section; lower min_shift")
  rbind(c(0, 0), out[seq_len(B), , drop = FALSE])
}

# Grid for a Gaussian scale: bin index per cell and a banded blur operator (zero padded, no wrap)
.grid <- function(coords, box, sigma, bin) {
  nx <- max(1, ceiling(box[2] / bin)); ny <- max(1, ceiling(box[4] / bin))
  gx <- pmin(floor((coords[, 1] - box[1]) / bin), nx - 1); gy <- pmin(floor((coords[, 2] - box[3]) / bin), ny - 1)
  band <- function(m) {
    r <- ceiling(3 * sigma / bin); off <- (-r):r; ker <- exp(-0.5 * (off * bin)^2 / sigma^2)
    i <- rep(seq_len(m), each = length(off)); j <- i + rep(off, m); ok <- j >= 1 & j <= m
    Matrix::sparseMatrix(i = i[ok], j = j[ok], x = rep(ker, m)[ok], dims = c(m, m))
  }
  r <- ceiling(3 * sigma / bin); ker <- exp(-0.5 * ((-r):r * bin)^2 / sigma^2)
  list(nx = nx, ny = ny, gx = gx, gy = gy, lin = gx + gy * nx + 1, Kx = band(nx), Ky = band(ny), bin = bin, ker = ker)
}

# Blur several grids at once. V: (nx * ny) x m (column-major nx x ny per column) -> same shape. Rows (x) by
# the band operator on the left, columns (y) by a block-diagonal band on the right: no transposes.
.blur <- function(V, gr, rs = NULL) {
  if (!is.null(rs)) return(lrmaprs::rs_blur(V, as.integer(gr$nx), as.integer(gr$ny), gr$ker, as.integer(rs)))
  m <- ncol(V); A <- gr$Kx %*% matrix(V, gr$nx)                            # nx x (ny * m)
  A <- A %*% if (m == 1) gr$Ky else Matrix::bdiag(rep(list(gr$Ky), m))      # Ky symmetric
  matrix(as.matrix(A), gr$nx * gr$ny)
}

# Sum values into grid bins: x is cells x m (cells subset `idx`), returns (nx*ny) x m
.bin_sum <- function(x, lin, ncell) {
  out <- matrix(0, ncell, ncol(x)); r <- rowsum(x, lin)
  out[as.integer(rownames(r)), ] <- r; out
}

