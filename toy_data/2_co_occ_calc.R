# toy_data/02_co_occ_calc.R
# Build a sparse code–code co-occurrence matrix from timestamped events.
# Two codes co-occur if they appear for the same patient within ±window_days.

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

# ---- Paths (repo-relative; run from repo root) ----
out_dir <- file.path("toy_data",  "data")   # you requested this pattern
in_path <- file.path(out_dir, "sequences.csv.gz")

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

# -------------------------------------------------------------------
# Helper: build a vocabulary (code -> index), with simple frequencies.
# -------------------------------------------------------------------
build_vocab <- function(events_dt, min_count = 1L) {
  stopifnot("Event_code" %in% names(events_dt))
  
  freq <- events_dt[, .N, by = .(Event_code)][order(-N, Event_code)]
  if (min_count > 1L) freq <- freq[N >= min_count]
  
  freq[, index := .I]
  setcolorder(freq, c("Event_code", "index", "N"))
  setnames(freq, c("Event_code", "index", "freq"))
  freq[]
}

# -------------------------------------------------------------------
# Core: build co-occurrence counts using a ±window_days time window.
# -------------------------------------------------------------------
build_cooc <- function(events_dt,
                       window_days = 2L,
                       vocab_dt = NULL,
                       include_self = FALSE,
                       verbose = TRUE) {
  
  required_cols <- c("ID", "Event_code", "t")
  stopifnot(all(required_cols %in% names(events_dt)))
  
  if (is.null(vocab_dt)) vocab_dt <- build_vocab(events_dt, min_count = 1L)
  vocab <- copy(vocab_dt)
  setkey(vocab, Event_code)
  V <- nrow(vocab)
  
  idx_map <- vocab$index
  names(idx_map) <- vocab$Event_code
  
  DT <- copy(events_dt)[, .(ID, Event_code, t)]
  setkeyv(DT, c("ID", "t"))
  DT[, rid := .I]
  
  D2 <- DT[, .(ID,
               Event_code,          # becomes i.Event_code in join output
               tmin = t - as.integer(window_days),
               tmax = t + as.integer(window_days),
               rid)]
  
  if (verbose) message("Joining within ±", window_days, " days ...")
  t0 <- Sys.time()
  
  J <- DT[D2,
          on = .(ID, t >= tmin, t <= tmax),
          allow.cartesian = TRUE,
          nomatch = 0L]
  
  # Keep each unordered pair once (upper triangle in time-order sense)
  J <- J[i.rid > rid]
  
  if (!include_self) {
    J <- J[Event_code != i.Event_code]
  }
  
  J[, i_idx := idx_map[Event_code]]
  J[, j_idx := idx_map[i.Event_code]]
  
  if (verbose) message("Aggregating pairs ...")
  
  P <- J[, .N, by = .(i_idx, j_idx)]
  setnames(P, "N", "count")
  
  # Make symmetric sparse matrix
  ii <- c(P$i_idx, P$j_idx)
  jj <- c(P$j_idx, P$i_idx)
  xx <- c(P$count, P$count)
  
  M <- sparseMatrix(i = ii, j = jj, x = as.numeric(xx),
                    dims = c(V, V), giveCsparse = TRUE)
  
  if (!include_self) diag(M) <- 0
  
  if (verbose) {
    dt <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2)
    message("Matrix built: ", V, " x ", V, " | nnz=", length(M@x), " | ", dt, "s")
  }
  
  list(cooc = M, vocab = vocab)
  
}

# -------------------------------------------------------------------
# Save outputs (all inside toy_data/data/)
# -------------------------------------------------------------------
save_outputs <- function(result, out_dir) {
  ensure_dir(out_dir)
  fwrite(result$vocab, file.path(out_dir, "vocab.csv"))
  saveRDS(result$cooc, file.path(out_dir, "cooc_matrix.rds"))
  invisible(TRUE)
}

# -------------------------------------------------------------------
# Entry point
# Usage:
#   Rscript toy_data/02_co_occ_calc.R 2
# -------------------------------------------------------------------
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  win_days <- if (length(args) >= 1) as.integer(args[1]) else 30L
  
  message("Reading: ", in_path)
  DT <- fread(in_path, colClasses = list(
    character = c("ID", "Event_code", "onto"),
    integer   = c("t")
  ))
  
  vocab <- build_vocab(DT, min_count = 1L)
  res <- build_cooc(DT, window_days = win_days, vocab_dt = vocab,
                    include_self = FALSE, verbose = TRUE)
  
  save_outputs(res, out_dir = out_dir)
  
  message("Wrote outputs to: ", out_dir)
  message("  - vocab.csv")
  message("  - cooc_matrix.rds")
}
