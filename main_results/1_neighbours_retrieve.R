# main_results/01_retrieve_neighbors.R
# ------------------------------------------------------------------------------
# Nearest-neighbour retrieval for a target medical code using cosine similarity.
#
# What this script does:
#   1) Loads released embeddings + vocab mapping
#   2) Selects one ontology (e.g., ICD, ATC, CCAM, LPP, NABM) by code prefix
#   3) Retrieves the top-k most similar codes to a target code (cosine similarity)
#
# Run from the REPOSITORY ROOT:
#   Rscript main_results/01_retrieve_neighbors.R
#
# What you typically change:
#   - target_code
#   - ontology
#   - k

# ------------------------------------------------------------------------------
# Supported ontology prefixes
#
# You can restrict neighbour search to a specific coding system by specifying
# the corresponding prefix via the `ontology` argument.
#
# Common prefixes include:
#   - "ICD" : ICD-10 diagnosis codes
#   - "ATC" : ATC medication codes
#   - "CAM" : CCAM procedure codes (France)
#   - "LPP" : LPP medical devices (France)
#   - "BIO" : NABM laboratory tests (France)
#
# Prefix matching is based on the first 3 characters of the code (check vocab.csv).
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

# ----------------------------- Repo-relative paths -----------------------------

data_dir <- file.path("main_results", "data")
func_dir <- file.path("main_results", "functions")

stop_if_not_repo_root <- function() {
  if (!dir.exists(data_dir) || !dir.exists(func_dir)) {
    stop(
      "Please run this script from the repository root.\n",
      "Expected to find:\n",
      "  - main_results/data/\n",
      "  - main_results/functions/\n"
    )
  }
}

# ------------------------------ Helper functions ------------------------------

# We keep helper code in main_results/functions/ so the main script stays readable.
load_helpers <- function() {
  source(file.path(func_dir, "similarity_measures.R"))  # cosine_dist()
  source(file.path(func_dir, "get_label.R"))            # get_label() (uses label_table)
}

load_vocab_codes <- function(vocab_path) {
  vocab <- fread(vocab_path)
  
  # Preferred convention: explicit index + Event_code
  if (all(c("index", "Event_code") %in% names(vocab))) {
    vocab <- vocab[order(index)]
    return(vocab$Event_code)
  }
  
  # Fallback: accept older or nonstandard vocab formatting
  warning(
    "vocab.csv does not contain both 'index' and 'Event_code'. ",
    "Falling back to second column as code."
  )
  vocab[[2]]
}

safe_label <- function(code) {
  out <- tryCatch(get_label(code), error = function(e) NA)
  
  # get_label() may return NA (logical) or a vector; normalize to length-1 character
  if (length(out) == 0 || is.null(out) || isTRUE(is.na(out))) {
    return(NA_character_)
  }
  
  out <- out[1]
  if (isTRUE(is.na(out))) return(NA_character_)
  as.character(out)
}

# ------------------------------ Core functionality -----------------------------

#' Find nearest neighbours of a target code within a chosen ontology.
#'
#' @param embeddings_path Path to embeddings CSV.GZ (rows = codes, cols = dims)
#' @param vocab_path      Path to vocab.csv (maps row index -> code)
#' @param target_code     Code to query (must exist in vocab)
#' @param ontology        Prefix filter ("ICD", "ATC", "CCA", "LPP", "NAB", ...)
#' @param k               Number of neighbours to return
#' @param verbose         Print progress + preview table
#'
#' @return data.frame with neighbour codes, cosine similarity, and labels
find_neighbours <- function(embeddings_path,
                            vocab_path,
                            target_code,
                            ontology,
                            k = 30,
                            verbose = TRUE) {
  
  if (verbose) message("Loading embeddings: ", embeddings_path)
  embeddings <- as.matrix(fread(embeddings_path))
  
  if (verbose) message("Loading vocab: ", vocab_path)
  code_vec <- load_vocab_codes(vocab_path)
  
  if (nrow(embeddings) != length(code_vec)) {
    stop(
      "Mismatch between embeddings rows (", nrow(embeddings),
      ") and vocab size (", length(code_vec), ")."
    )
  }
  
  if (!target_code %in% code_vec) {
    stop("Target code '", target_code, "' not found in vocab.csv.")
  }
  
  ontology <- toupper(ontology)
  
  # Keep only codes from the requested ontology (+ target itself)
  idx_target <- which(code_vec == target_code)
  idx_onto   <- which(substr(code_vec, 1, 3) == ontology)
  if (length(idx_onto) == 0) stop("No codes found for ontology prefix: ", ontology)
  
  idx_all <- unique(c(idx_target, idx_onto))
  emb_sub <- embeddings[idx_all, , drop = FALSE]
  cod_sub <- code_vec[idx_all]
  
  if (verbose) {
    message(sprintf(
      "Search space: %d codes (%s) | embedding dim: %d",
      length(cod_sub), ontology, ncol(emb_sub)
    ))
    message("Computing cosine similarities...")
  }
  
  target_vec <- emb_sub[cod_sub == target_code, , drop = FALSE]
  
  sim <- vapply(
    seq_len(nrow(emb_sub)),
    function(j) cosine_dist(target_vec[1, ], emb_sub[j, ], sparse = FALSE),
    numeric(1)
  )
  names(sim) <- cod_sub
  
  # Remove self
  sim <- sim[names(sim) != target_code]
  
  # Top-k
  ord <- order(sim, decreasing = TRUE)
  sim_ord <- sim[ord]
  k_use <- min(k, length(sim_ord))
  
  neighbours <- names(sim_ord)[seq_len(k_use)]
  sims <- as.numeric(sim_ord[seq_len(k_use)])
  
  res <- data.frame(
    target = target_code,
    neighbour = neighbours,
    sim = sims,
    neighbour_label = vapply(neighbours, safe_label, character(1)),
    stringsAsFactors = FALSE
  )
  
  if (verbose) {
    tgt_lab <- safe_label(target_code)
    message("\n=== Nearest neighbours ===")
    message("Target: ", target_code, if (!is.na(tgt_lab)) paste0(" — ", tgt_lab) else "")
    message("Ontology filter: ", ontology)
    message("Top-k: ", k_use, "\n")
    
    print(utils::head(res, 10), row.names = FALSE)
    cat("\nTip: write.csv(res, \"neighbours.csv\", row.names = FALSE) to save results.\n")
  }
  
  res
}

# -------------------------------- Entrypoint ----------------------------------

if (sys.nframe() == 0L) {
  
  stop_if_not_repo_root()
  load_helpers()
  
  # Label table used by get_label() (loaded once here)
  label_table <- fread(file.path(data_dir, "label_table.csv"))
  
  # ---- Files (released data) ----
  embeddings_path <- file.path(data_dir, "embeddings_ESND_2FC.csv.gz")
  vocab_path      <- file.path(data_dir, "vocab.csv")
  
  # ---- User parameters (edit these) ----
  target_code <- "ICD-G401"   # example: epilepsy ICD code
  ontology    <- "BIO"        # search neighbours among ATC codes (check script header for possible ontologies)
  k           <- 20
  
  res <- find_neighbours(
    embeddings_path = embeddings_path,
    vocab_path      = vocab_path,
    target_code     = target_code,
    ontology        = ontology,
    k               = k,
    verbose         = TRUE
  )
  
  res
}
