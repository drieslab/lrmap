# Ligand-receptor table: columns ligand, receptor; optional ligand_subunits / receptor_subunits ("A_B").
# Returns the unique pairs with subunit gene lists, and how many have all their genes in `genes`.
.lr_table <- function(lr, genes) {
  lr <- as.data.frame(lr, stringsAsFactors = FALSE)
  if (!all(c("ligand", "receptor") %in% names(lr))) stop("lr needs columns ligand and receptor")
  sub <- function(col, main) {
    v <- if (col %in% names(lr)) as.character(lr[[col]]) else as.character(lr[[main]])
    bad <- is.na(v) | v == ""
    v[bad] <- as.character(lr[[main]])[bad]
    strsplit(v, "_", fixed = TRUE)
  }
  out <- data.frame(ligand = as.character(lr$ligand), receptor = as.character(lr$receptor), stringsAsFactors = FALSE)
  out$lg <- sub("ligand_subunits", "ligand"); out$rg <- sub("receptor_subunits", "receptor")
  out <- out[!duplicated(paste(out$ligand, out$receptor)), , drop = FALSE]
  all_g <- unique(unlist(c(out$lg, out$rg)))
  if (is.null(genes)) return(list(genes = all_g))
  ok <- vapply(seq_len(nrow(out)), function(i) all(c(out$lg[[i]], out$rg[[i]]) %in% genes), TRUE)
  list(table = out[ok, , drop = FALSE], n = nrow(out), n_covered = sum(ok), genes = all_g)
}

# Entity matrices (sparse, cells x entities) from composition-scale Y (sparse cells x genes): ligand = minimum
# over subunits, receptor = mean. D: all subunits detected. k: number of subunits.
.entities_sparse <- function(names, subunits, Y, fun) {
  u <- !duplicated(names); names <- names[u]; subunits <- subunits[u]; n <- nrow(Y)
  cols <- lapply(subunits, function(g) {
    if (length(g) == 1) return(list(e = Y[, g, drop = FALSE], d = methods::as((Y[, g, drop = FALSE] > 0) * 1, "CsparseMatrix")))
    M <- as.matrix(Y[, g, drop = FALSE])
    v <- if (fun == "min") do.call(pmin, as.data.frame(M)) else rowMeans(M)
    list(e = Matrix::Matrix(v, n, 1, sparse = TRUE), d = Matrix::Matrix((rowSums(M > 0) == length(g)) * 1, n, 1, sparse = TRUE))
  })
  E <- do.call(cbind, lapply(cols, `[[`, "e")); D <- do.call(cbind, lapply(cols, `[[`, "d"))
  colnames(E) <- colnames(D) <- names
  list(E = methods::as(E, "CsparseMatrix"), D = methods::as(D, "CsparseMatrix"), k = lengths(subunits))
}
