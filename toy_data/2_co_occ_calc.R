# toy_data/02_co_occ_calc.R
# Create a sparse code-to-code co-occurrence matrix from event sequences.
# Two codes are counted as co-occurring when they happen for the same patient
# within a chosen time window (±window_days).


# Output folder and input file, relative to the repository root
out_dir <- file.path("toy_data", "data")
in_path <- file.path(out_dir, "sequences.csv.gz")

# Create a directory if it does not already exist
ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

# Build the vocabulary table used by the matrix.
# Each unique Event_code gets:
#   - its frequency in the data
#   - a numeric index for matrix row/column lookup
build_vocab <- function(events_dt, min_count = 1L) {
  stopifnot("Event_code" %in% names(events_dt))

  freq <- events_dt[, .N, by = .(Event_code)][order(-N, Event_code)]
  if (min_count > 1L) freq <- freq[N >= min_count]

  freq[, index := .I]
  setcolorder(freq, c("Event_code", "index", "N"))
  setnames(freq, c("Event_code", "index", "freq"))
  freq[]
}

# Build a sparse co-occurrence matrix from timestamped events.
#
# Expected columns in events_dt:
#   - ID: patient identifier
#   - Event_code: code observed at that time
#   - t: integer timestamp in days
#
# For each event, we look for other events from the same patient that fall
# within ±window_days, then count code pairs across all patients.
build_cooc <- function(events_dt,
                       window_days = 2L,
                       vocab_dt = NULL,
                       include_self = FALSE,
                       verbose = TRUE) {

  required_cols <- c("ID", "Event_code", "t")
  stopifnot(all(required_cols %in% names(events_dt)))

  # Reuse a provided vocabulary when available so matrix indices stay stable
  if (is.null(vocab_dt)) vocab_dt <- build_vocab(events_dt, min_count = 1L)
  vocab <- copy(vocab_dt)
  setkey(vocab, Event_code)
  V <- nrow(vocab)

  # Named lookup: Event_code -> matrix index
  idx_map <- vocab$index
  names(idx_map) <- vocab$Event_code

  # Keep only the columns needed for co-occurrence counting
  # rid is a unique row id used later to avoid double-counting pairs
  DT <- copy(events_dt)[, .(ID, Event_code, t)]
  setkeyv(DT, c("ID", "t"))
  DT[, rid := .I]

  # For each event, define the time interval in which we search for neighbors
  D2 <- DT[, .(ID,
               Event_code,
               tmin = t - as.integer(window_days),
               tmax = t + as.integer(window_days),
               rid)]

  if (verbose) message("Joining within ±", window_days, " days ...")
  t0 <- Sys.time()

  # Join each event against events from the same patient whose timestamp
  # falls inside the requested window
  J <- DT[D2,
          on = .(ID, t >= tmin, t <= tmax),
          allow.cartesian = TRUE,
          nomatch = 0L]

  # Keep each pair only once
  # This removes mirrored duplicates such as (A, B) and (B, A) from the join
  J <- J[i.rid > rid]

  # Optionally drop pairs where the same code is matched with itself
  if (!include_self) {
    J <- J[Event_code != i.Event_code]
  }

  # Convert event codes to matrix indices
  J[, i_idx := idx_map[Event_code]]
  J[, j_idx := idx_map[i.Event_code]]

  if (verbose) message("Aggregating pairs ...")

  # Count how many times each code pair appears
  P <- J[, .N, by = .(i_idx, j_idx)]
  setnames(P, "N", "count")

  # Expand to a symmetric matrix:
  # if code A co-occurs with code B, then B co-occurs with A
  ii <- c(P$i_idx, P$j_idx)
  jj <- c(P$j_idx, P$i_idx)
  xx <- c(P$count, P$count)

  M <- sparseMatrix(i = ii, j = jj, x = as.numeric(xx),
                    dims = c(V, V), giveCsparse = TRUE)

  # Clear the diagonal unless self-co-occurrence is explicitly requested
  if (!include_self) diag(M) <- 0

  if (verbose) {
    dt <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2)
    message("Matrix built: ", V, " x ", V, " | nnz=", length(M@x), " | ", dt, "s")
  }

  list(cooc = M, vocab = vocab)
}

# Save the vocabulary and sparse matrix to disk
save_outputs <- function(result, out_dir) {
  ensure_dir(out_dir)
  fwrite(result$vocab, file.path(out_dir, "vocab.csv"))
  saveRDS(result$cooc, file.path(out_dir, "cooc_matrix.rds"))
  invisible(TRUE)
}

# Script entry point
# Example:
#   Rscript toy_data/02_co_occ_calc.R 2
#
# If no argument is passed, the default window is 30 days.
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
