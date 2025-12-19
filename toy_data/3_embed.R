# toy_data/03_embed.R
#
# Compute unshifted PPMI from a sparse co-occurrence matrix, then
# apply truncated SVD to obtain k-dimensional embeddings.

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(RSpectra)
})

# ------------------------------- Utilities ------------------------------------

ensure_dir <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

build_ppmi <- function(C) {
  C <- drop0(C)
  C <- (C + t(C)) / 2
  diag(C) <- 0
  
  total <- sum(C)
  if (total <= 0) stop("Co-occurrence matrix has zero total count.")
  
  p_i <- as.numeric(Matrix::rowSums(C)) / total
  
  S <- summary(C)
  p_ij <- S$x / total
  denom <- p_i[S$i] * p_i[S$j]
  
  pmi <- log(p_ij / denom)
  pmi[!is.finite(pmi)] <- NA_real_
  keep <- !is.na(pmi) & (pmi > 0)
  
  if (!any(keep)) stop("PPMI became empty (no positive PMI entries).")
  
  i <- S$i[keep]; j <- S$j[keep]; v <- pmi[keep]
  
  P <- sparseMatrix(i = c(i, j), j = c(j, i), x = c(v, v), dims = dim(C))
  diag(P) <- 0
  drop0(P)
}

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

save_embeddings <- function(E, out_dir, base = "embeddings") {
  ensure_dir(out_dir)
  fwrite(as.data.table(E),
         file.path(out_dir, paste0(base, ".csv.gz")))
}

# ----------------------------- Script entrypoint -------------------------------

# Usage:
#   Rscript toy_data/03_embed.R 100
#
# Defaults assume you ran:
#   01_dataset_toy_gen.R  -> toy_data/data/sequences.csv.gz
#   02_co_occ_calc.R      -> toy_data/data/cooc_matrix.rds + toy_data/data/vocab.csv
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  
  # Default toy paths (repo-relative; run from repo root)
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
