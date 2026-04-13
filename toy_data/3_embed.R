# toy_data/03_embed.R
#
# Build code embeddings from a sparse co-occurrence matrix.
# The script first converts co-occurrence counts into a PPMI matrix,
# then applies truncated SVD to get a lower-dimensional embedding space.

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(RSpectra)
})

# Create a directory if it does not already exist
ensure_dir <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

# Convert a co-occurrence count matrix into a PPMI matrix.
#
# Steps:
#   - remove explicit zeros
#   - force symmetry
#   - clear the diagonal
#   - compute PMI for nonzero pairs
#   - keep only positive PMI values
#
# The output stays sparse, which keeps memory usage reasonable.
build_ppmi <- function(C) {
  C <- drop0(C)
  C <- (C + t(C)) / 2
  diag(C) <- 0

  total <- sum(C)
  if (total <= 0) stop("Co-occurrence matrix has zero total count.")

  # Marginal probability for each code
  p_i <- as.numeric(Matrix::rowSums(C)) / total

  # Work only on nonzero entries of the sparse matrix
  S <- summary(C)
  p_ij <- S$x / total
  denom <- p_i[S$i] * p_i[S$j]

  pmi <- log(p_ij / denom)
  pmi[!is.finite(pmi)] <- NA_real_
  keep <- !is.na(pmi) & (pmi > 0)

  if (!any(keep)) stop("PPMI became empty (no positive PMI entries).")

  i <- S$i[keep]
  j <- S$j[keep]
  v <- pmi[keep]

  # Mirror the values so the result is symmetric
  P <- sparseMatrix(i = c(i, j), j = c(j, i), x = c(v, v), dims = dim(C))
  diag(P) <- 0
  drop0(P)
}

# Compute dense embeddings from the sparse PPMI matrix using truncated SVD.
#
# k is the requested embedding size.
# If k is larger than the matrix allows, it is reduced automatically.
#
# The returned embedding matrix uses U * sqrt(D), which is a common way
# to turn the SVD factors into vector representations.
svd_embeddings <- function(P, k = 10, seed = 42) {
  set.seed(seed)
  max_k <- max(1L, min(dim(P)) - 1L)
  if (k > max_k) {
    message("Requested k = ", k, " > max allowable (", max_k, "). Using k = ", max_k, ".")
    k <- max_k
  }

  s <- RSpectra::svds(P, k = k)
  E <- s$u %*% diag(sqrt(pmax(s$d, 0)))
  list(E = E, svals = s$d, k = k)
}

# Save the embedding matrix as a compressed CSV file
save_embeddings <- function(E, out_dir, base = "embeddings") {
  ensure_dir(out_dir)
  fwrite(as.data.table(E),
         file.path(out_dir, paste0(base, ".csv.gz")))
}

# Script entry point
#
# Example:
#   Rscript toy_data/03_embed.R 100
#
# By default, the script expects the earlier steps in the toy pipeline
# to have already written their outputs into toy_data/data/.
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)

  # Default file locations, relative to the repository root
  data_dir   <- file.path("toy_data", "data")
  cooc_path  <- ifelse(length(args) >= 1 && nzchar(args[1]),
                       args[1],
                       file.path(data_dir, "cooc_matrix.rds"))
  vocab_path <- ifelse(length(args) >= 2 && nzchar(args[2]),
                       args[2],
                       file.path(data_dir, "vocab.csv"))
  out_emb    <- ifelse(length(args) >= 3 && nzchar(args[3]),
                       args[3],
                       data_dir)
  k          <- ifelse(length(args) >= 4 && nzchar(args[4]),
                       as.integer(args[4]),
                       10L)

  message("Loading co-occurrence matrix: ", cooc_path)
  C <- readRDS(cooc_path)
  if (!inherits(C, "dgCMatrix")) stop("Expected a dgCMatrix at ", cooc_path)

  message("Reading vocab: ", vocab_path)
  vocab <- fread(vocab_path)
  V <- nrow(vocab)

  # The co-occurrence matrix and vocab should describe the same code set
  if (!all(dim(C) == c(V, V))) {
    warning("Dimension mismatch: cooc dims = (", paste(dim(C), collapse="x"),
            "), vocab size = ", V, ". Proceeding, but check inputs.")
  }

  message("Computing unshifted PPMI ...")
  P <- build_ppmi(C)

  message("Running truncated SVD (k = ", k, ") ...")
  svd_res <- svd_embeddings(P, k = k, seed = 42)
  E <- svd_res$E

  message("Saving embeddings ...")
  save_embeddings(E, out_emb, base = "embeddings")

  message("Done.\nEmbeddings: ", file.path(out_emb, "embeddings.csv.gz"))
}
