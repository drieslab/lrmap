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
  ok <- vapply(seq_len(nrow(out)), function(i) all(c(out$lg[[i]], out$rg[[i]]) %in% genes), TRUE)
  list(table = out[ok, , drop = FALSE], n = nrow(out), n_covered = sum(ok))
}

# Entity matrices (cells x entities): ligand = minimum over subunits, receptor = mean over subunits.
# Y: cells x genes (composition scale, possibly leak-corrected); D: cells x genes detection (0/1).
.entities <- function(names, subunits, Y, D, fun) {
  u <- !duplicated(names); names <- names[u]; subunits <- subunits[u]
  E <- vapply(subunits, function(g) {
    M <- Y[, g, drop = FALSE]
    if (ncol(M) == 1) M[, 1] else if (fun == "min") do.call(pmin, as.data.frame(M)) else rowMeans(M)
  }, numeric(nrow(Y)))
  Dt <- vapply(subunits, function(g) do.call(pmin, as.data.frame(D[, g, drop = FALSE])), numeric(nrow(Y)))
  E <- matrix(E, nrow(Y)); Dt <- matrix(Dt, nrow(Y)); colnames(E) <- colnames(Dt) <- names
  list(E = E, D = Dt)
}
