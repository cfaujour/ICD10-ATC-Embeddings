# main_results/01_retrieve_neighbors.R
# Retrieve nearest neighbours for a target code using cosine similarity.
# Intended to be run from the REPOSITORY ROOT (repo-relative paths).
#
# Outputs: prints a preview table; returns a data.frame of neighbours.

suppressPackageStartupMessages({
  library(data.table)
})

# ---- Paths (repo-relative; run from repo root) ----
data_dir <- file.path("main_results", "data")
func_dir <- file.path("main_results", "functions")

# Guardrail: fail fast if not run from repo root
if (!dir.exists(data_dir) || !dir.exists(func_dir)) {
  stop("Please run this script from the repository root (cannot find 'main_results/data' or 'main_results/functions').")
}

# ---- Helper functions ----
source(file.path(func_dir, "similarity_measures.R"))  # provides cosine_dist()
source(file.path(func_dir, "get_label.R"))            # get_label() may use label_table

# ---- Labels (used by get_label()) ----
# Kept global for minimal change; helper expects it.
label_table <- fread(file.path(data_dir, "label_table.csv"))

# =====================================================================
# find_neighbours()
#   - Loads embeddings + vocab
#   - Filters to a given ontology (prefix of code, e.g., "ICD", "ATC")
#   - Computes cosine similarity for one target code
#   - Returns a data.frame of neighbours + similarity + label
# =====================================================================
find_neighbours <- function(embeddings_path,
                            vocab_path,
                            target_code,
                            ontology,
                            k = 30,
                            verbose = TRUE) {
  
  # ---- Load embeddings ----
  if (verbose) message("Loading embeddings from: ", embeddings_path)
  emb_dt <- fread(embeddings_path)
  embeddings <- as.matrix(emb_dt)
  
  if (verbose) {
    message(sprintf(
      "Loaded embeddings: %d codes x %d dimensions.",
      nrow(embeddings), ncol(embeddings)
    ))
  }
  
  # ---- Load vocabulary ----
  if (verbose) message("Loading vocabulary from: ", vocab_path)
  vocab <- fread(vocab_path)
  
  # Preferred convention: explicit index + code
  if (all(c("index", "Event_code") %in% names(vocab))) {
    vocab <- vocab[order(index)]
    code_vec <- vocab$Event_code
  } else {
    warning(
      "vocab.csv does not contain both 'index' and 'Event_code'. ",
      "Falling back to second column as code."
    )
    code_vec <- vocab[[2]]
  }
  
  if (nrow(embeddings) != length(code_vec)) {
    stop(
      "Mismatch between embeddings (", nrow(embeddings),
      ") and vocab size (", length(code_vec), ")."
    )
  }
  
  # ---- Check target exists ----
  if (!target_code %in% code_vec) {
    stop(
      "Target code '", target_code,
      "' not found in vocabulary. Check spelling or prefix."
    )
  }
  
  target_label <- tryCatch(
    get_label(target_code),
    error = function(e) NA_character_
  )
  
  # ---- Filter by ontology prefix ----
  ontology <- toupper(ontology)
  if (verbose) message("Restricting search to ontology: ", ontology)
  
  idx_target <- which(code_vec == target_code)
  idx_onto   <- which(substr(code_vec, 1, 3) == ontology)
  
  if (length(idx_onto) == 0) {
    stop("No codes found for ontology '", ontology, "'.")
  }
  
  idx_all <- unique(c(idx_target, idx_onto))
  embeddings_sub <- embeddings[idx_all, , drop = FALSE]
  code_sub <- code_vec[idx_all]
  
  # ---- Compute cosine similarities ----
  if (verbose) message("Computing cosine similarities...")
  
  target_vec <- embeddings_sub[code_sub == target_code, , drop = FALSE]
  
  sim <- vapply(
    seq_len(nrow(embeddings_sub)),
    function(j) cosine_dist(
      target_vec[1, ],
      embeddings_sub[j, ],
      sparse = FALSE
    ),
    numeric(1)
  )
  names(sim) <- code_sub
  
  # Remove self
  sim <- sim[names(sim) != target_code]
  
  # ---- Build neighbourhood table ----
  ord <- order(sim, decreasing = TRUE)
  sim_ord <- sim[ord]
  
  k_use <- min(k, length(sim_ord))
  neighbours <- names(sim_ord)[seq_len(k_use)]
  sims <- as.numeric(sim_ord[seq_len(k_use)])
  
  res <- data.frame(
    target = target_code,
    neighbour = neighbours,
    sim = sims,
    neighbour_label = NA_character_,
    stringsAsFactors = FALSE
  )
  
  # Label lookup via helper
  res$neighbour_label <- vapply(res$neighbour, get_label, character(1))
  
  if (verbose) {
    
    header <- paste0(
      "\n=== Nearest neighbours for ",
      target_code,
      if (!is.na(target_label)) paste0(" — ", target_label) else "",
      " (ontology: ", ontology, ") ===\n"
    )
    message(header)
    
    print(utils::head(res, 10), row.names = FALSE)
    cat("\n------------------------------------------------------\n")
    cat("Displayed top 10 of ", nrow(res),
        " neighbours (top-", k_use, ").\n", sep = "")
    cat("Tip: write.csv(res, \"neighbours.csv\", row.names = FALSE) to save results.\n")
  }
  
  res
}

# --------------------------- Script entrypoint ---------------------------
# Example call (runs only when executed as a script)
if (sys.nframe() == 0L) {
  
  embeddings_path <- file.path(data_dir, "embeddings_ESND_2FC.csv.gz")
  vocab_path      <- file.path(data_dir, "vocab.csv")
  
  target_code <- "ICD-G401"   # example: epilepsy
  ontology    <- "ATC"        # search neighbours among ATC codes
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
