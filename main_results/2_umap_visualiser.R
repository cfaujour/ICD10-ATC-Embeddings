# main_results/02_umap_graph.R
# ------------------------------------------------------------------------------
# Explore the released embeddings with two beginner-friendly visual tools:
#
# (1) UMAP (2D scatter plot) for a subset of codes
#     - each point = one medical code
#     - colour = chapter / class label (a coarse grouping used in the paper figure)
#     - shape  = ontology prefix (ICD / ATC / CAM / LPP / BIO / ...)
#
# (2) Local similarity graph around a target code
#     - builds a small network of the closest neighbours of a target code
#     - neighbours can be retrieved across several ontologies (e.g. ICD + ATC)
#
# Important:
# - This script assumes you cloned the repository and kept the folder structure.
# - Run from the REPOSITORY ROOT so relative paths work:
#     Rscript main_results/02_umap_graph.R
#
# What you typically edit (near the bottom, in the "USAGE SETUP" section):
#   - codes_subset: which codes are displayed in the UMAP, should be a subset from vocab.csv code names
#   - target_code:  which code is used as the center of the local graph
#   - ontologies / k_per_onto / sim_threshold : ontologies to include in the graph / top k neighbours in each ontology / similarity value threshold for edge to be included
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(umap)
  library(ggplot2)
  library(igraph)
})

# ----------------------------- Repo-relative paths -----------------------------
# We keep all data files inside main_results/data/ and helper functions inside
# main_results/functions/. These file.path() calls build OS-independent paths.

data_dir <- file.path("main_results", "data")
func_dir <- file.path("main_results", "functions")

# Small helper: ensure the script is executed from the repository root.
# (This avoids confusing "file not found" errors later.)
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

# ------------------------------ Load helper code ------------------------------
# We keep core utilities in separate files to keep this script readable:
# - cosine_dist(): cosine similarity between two embedding vectors
# - get_label():   map a code to a human-readable label (if available)

load_helpers <- function() {
  source(file.path(func_dir, "similarity_measures.R"))  # cosine_dist()
  source(file.path(func_dir, "get_label.R"))            # get_label() (uses label_table)
}

# get_label() may return NA in different "types" (logical vs character),
# depending on the input code. safe_label() ensures we always output a single
# character string (or NA_character_). This prevents vapply() type errors.
safe_label <- function(code) {
  out <- tryCatch(get_label(code), error = function(e) NA)
  if (length(out) == 0 || is.null(out) || isTRUE(is.na(out))) return(NA_character_)
  out <- out[1]
  if (isTRUE(is.na(out))) return(NA_character_)
  as.character(out)
}

# ------------------------------ Load embeddings -------------------------------
# The released embeddings are stored as a CSV.GZ where:
#   - each row corresponds to one code
#   - each column is one embedding dimension
#
# To interpret rows, we also load vocab.csv, which maps row indices -> code strings.

load_embeddings <- function(embeddings_path = file.path(data_dir, "embeddings_ESND_2FC.csv.gz"),
                            vocab_path      = file.path(data_dir, "vocab.csv"),
                            verbose         = TRUE) {
  
  if (verbose) message("Loading embeddings: ", embeddings_path)
  embeddings <- as.matrix(fread(embeddings_path))
  
  if (verbose) message("Loading vocab: ", vocab_path)
  vocab <- fread(vocab_path)
  
  # Preferred vocab format: columns named "index" and "Event_code"
  if (all(c("index", "Event_code") %in% names(vocab))) {
    vocab <- vocab[order(index)]
    dict <- vocab$Event_code
  } else {
    # Fallback for older formatting: assume 2nd column contains the code string
    warning("vocab.csv missing 'index'/'Event_code'. Using 2nd column as code.")
    dict <- vocab[[2]]
  }
  
  # Basic consistency check: embeddings rows must match vocab length
  if (nrow(embeddings) != length(dict)) {
    stop("Mismatch: embeddings rows (", nrow(embeddings), ") vs vocab size (", length(dict), ").")
  }
  
  # Set rownames so we can subset by code easily later
  rownames(embeddings) <- dict
  
  if (verbose) {
    message("Loaded ", nrow(embeddings), " codes × ", ncol(embeddings), " dimensions.\n")
  }
  
  list(embeddings = embeddings, dict = dict)
}

# ------------------------------ Chapter mapping -------------------------------
# We colour the UMAP by a "chapter" or coarse group label used in the paper figure.
# The file chapter_table.csv contains a mapping:
#   code -> y   (where y is a chapter/group label string)

load_chapter_map <- function(path = file.path(data_dir, "chapter_table.csv")) {
  tab <- fread(path)
  if (!all(c("code", "y") %in% names(tab))) {
    stop("chapter_table.csv must contain columns: 'code' and 'y'.")
  }
  setNames(tab$y, tab$code)
}

get_chapter <- function(code, chapter_map) {
  v <- chapter_map[[code]]
  if (is.null(v)) NA_character_ else v
}

# ------------------------------ UMAP functions --------------------------------
# UMAP reduces the high-dimensional embedding vectors to 2 dimensions so we can plot
# them in a scatter plot. We typically run it on a subset (not all codes) for speed
# and readability.

run_umap_subset <- function(embeddings,
                            dict,
                            codes_subset,
                            chapter_map,
                            n_neighbors = 25,
                            min_dist    = 0.5,
                            seed        = 16,
                            verbose     = TRUE) {
  
  # Keep only the codes present in the subset
  idx <- which(dict %in% codes_subset)
  if (!length(idx)) stop("No codes from codes_subset were found in vocab/embeddings.")
  
  dict_sub <- dict[idx]
  emb_sub  <- embeddings[idx, , drop = FALSE]
  
  if (verbose) message("Running UMAP on ", length(dict_sub), " codes...")
  
  # UMAP settings: match paper defaults (cosine distance in embedding space)
  set.seed(seed)
  cfg <- umap::umap.defaults
  cfg$metric      <- "cosine"
  cfg$n_neighbors <- n_neighbors
  cfg$min_dist    <- min_dist
  
  u <- umap::umap(emb_sub, n_components = 2, config = cfg)
  
  # Build a plotting table:
  # - UMAP1/UMAP2 are 2D coordinates
  # - onto is the ontology prefix = first 3 characters of the code string
  data.frame(
    UMAP1   = u$layout[, 1],
    UMAP2   = u$layout[, 2],
    code    = dict_sub,
    onto    = substr(dict_sub, 1, 3),
    chapter = vapply(dict_sub, get_chapter, character(1), chapter_map = chapter_map),
    stringsAsFactors = FALSE
  )
}

plot_umap <- function(df_umap, palette, point_size = 3) {
  ggplot(df_umap, aes(UMAP1, UMAP2)) +
    geom_point(aes(colour = chapter, shape = onto), size = point_size, alpha = 0.9) +
    scale_color_manual(values = palette) +
    theme_minimal(base_size = 13) +
    labs(
      title  = "UMAP of embedding subset",
      x      = "UMAP 1",
      y      = "UMAP 2",
      colour = "Chapter",
      shape  = "Ontology"
    )
}

# -------------------------- Local neighbourhood graph --------------------------
# This section builds a small graph around a target code.
# Steps (high level):
#   1) For each ontology, find the top-k neighbours of the target by cosine similarity
#   2) Combine them into one small node set (target + neighbours)
#   3) Compute pairwise similarities in this small set
#   4) Keep edges above sim_threshold and plot the graph

build_local_neighbour_graph <- function(embeddings,
                                        dict,
                                        target,
                                        ontologies     = c("ICD", "ATC"),
                                        k_per_onto     = 5,
                                        sim_threshold  = 0.5,
                                        verbose        = TRUE) {
  
  if (!target %in% dict) stop("Target code not found: ", target)
  
  if (verbose) {
    message("\nBuilding local graph around: ", target)
    message("Ontologies: ", paste(ontologies, collapse = ", "),
            " | k_per_onto = ", k_per_onto,
            " | sim_threshold = ", sim_threshold)
  }
  
  neighbours <- character(0)
  
  # Neighbour search per ontology prefix
  for (onto in ontologies) {
    idx_onto <- which(substr(dict, 1, 3) == onto)
    if (!length(idx_onto)) next
    
    idx_keep <- unique(c(which(dict == target), idx_onto))
    emb_sub  <- embeddings[idx_keep, , drop = FALSE]
    dict_sub <- rownames(emb_sub)
    
    target_vec <- emb_sub[dict_sub == target, , drop = FALSE]
    
    sims <- vapply(seq_len(nrow(emb_sub)), function(j) {
      cosine_dist(target_vec[1, ], emb_sub[j, ], sparse = FALSE)
    }, numeric(1))
    names(sims) <- dict_sub
    
    sims <- sims[names(sims) != target]
    sims <- sort(sims, decreasing = TRUE)
    
    k_use <- min(k_per_onto, length(sims))
    if (k_use > 0) neighbours <- c(neighbours, names(sims)[seq_len(k_use)])
  }
  
  # Final node list: target + unique neighbours across ontologies
  nodes <- unique(c(target, neighbours))
  emb_g <- embeddings[match(nodes, dict), , drop = FALSE]
  
  # Pairwise similarity on this small set (usually 10–30 nodes)
  n <- nrow(emb_g)
  sim_mat <- matrix(0, n, n, dimnames = list(nodes, nodes))
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      sim_mat[i, j] <- cosine_dist(emb_g[i, ], emb_g[j, ], sparse = FALSE)
    }
  }
  
  # Keep only "strong" edges
  sim_mat[sim_mat < sim_threshold] <- 0
  diag(sim_mat) <- 0
  
  g <- graph_from_adjacency_matrix(sim_mat, mode = "undirected", weighted = TRUE, diag = FALSE)
  g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  
  # Add node metadata for plotting/labels
  V(g)$code     <- V(g)$name
  V(g)$ontology <- substr(V(g)$code, 1, 3)
  V(g)$ontology[V(g)$code == target] <- "Target"
  V(g)$label    <- vapply(V(g)$code, safe_label, character(1))
  
  g
}

plot_local_graph <- function(g,
                             layout_fun  = layout_with_fr,
                             label_type  = c("code", "label"),
                             vertex_size = 16,
                             seed        = 42,
                             main        = "Local similarity graph") {
  
  label_type <- match.arg(label_type)
  
  # layout = where nodes are placed on the page; set.seed() makes it reproducible
  set.seed(seed)
  lay <- layout_fun(g)
  
  # Simple colour mapping by ontology prefix (plus a Target colour)
  base_cols <- c(
    Target = "#e41a1c",
    ICD    = "#377eb8",
    ATC    = "#4daf4a",
    CAM    = "#984ea3",
    LPP    = "#ff7f00",
    BIO    = "#a65628",
    Other  = "grey50"
  )
  
  ontos <- V(g)$ontology
  v_col <- vapply(
    ontos,
    function(o) if (o %in% names(base_cols)) base_cols[o] else base_cols["Other"],
    character(1)
  )
  v_col <- as.character(v_col)
  
  # Make the target node slightly larger and outlined
  is_target <- ontos == "Target"
  v_sz <- ifelse(is_target, vertex_size * 1.3, vertex_size)
  v_fr <- ifelse(is_target, "black", NA)
  
  # Choose what text to display on nodes
  v_lab <- if (label_type == "code") V(g)$code else V(g)$label
  
  # Edge width scaled by similarity weight (if there are edges)
  w <- E(g)$weight
  e_w <- if (length(w) && max(w, na.rm = TRUE) > min(w, na.rm = TRUE)) {
    1 + 4 * (w - min(w, na.rm = TRUE)) / (max(w, na.rm = TRUE) - min(w, na.rm = TRUE))
  } else {
    rep(2, length(w))
  }
  
  plot(
    g,
    layout             = lay,
    vertex.size        = v_sz,
    vertex.color       = v_col,
    vertex.frame.color = v_fr,
    vertex.label       = v_lab,
    vertex.label.cex   = 0.7,
    vertex.label.color = "black",
    edge.width         = e_w,
    edge.color         = adjustcolor("grey40", alpha.f = 0.8),
    main               = main
  )
}

# ==============================================================================
# USAGE SETUP (runs only when executed as a script)
# ==============================================================================
if (sys.nframe() == 0L) {
  
  stop_if_not_repo_root()
  load_helpers()
  
  # get_label() expects label_table to exist in the global environment
  label_table <- fread(file.path(data_dir, "label_table.csv"))
  
  # Load released embeddings and code list
  emb_data   <- load_embeddings()
  embeddings <- emb_data$embeddings
  dict       <- emb_data$dict
  
  # Load chapter/class mapping used for colouring in the paper figure
  chapter_map <- load_chapter_map()
  
  # ----------------------------------------------------------------------------
  # (1) UMAP EXAMPLE
  # ----------------------------------------------------------------------------
  # codes_subset controls which codes appear in the UMAP plot.
  # For large vocabularies, plotting everything is slow and unreadable, so we
  # recommend selecting a subset (chapters, ontologies, or manual lists).
  
  codes_subset <- dict[substr(dict, 1, 5) %in% c(
    "ICD-C", "ICD-D", "ATC-L", "ICD-G", "ATC-N", "ICD-H", "ATC-S", "ICD-J", "ATC-R"
  )]
  
  umap_df <- run_umap_subset(
    embeddings   = embeddings,
    dict         = dict,
    codes_subset = codes_subset,
    chapter_map  = chapter_map,
    n_neighbors  = 25,
    min_dist     = 0.5,
    seed         = 16
  )
  
  # Optional filtering (matches the paper code subset)
  umap_df <- subset(umap_df, !chapter %in% c("ICD-NA", "NA"))
  umap_df <- umap_df[umap_df$chapter != "ICD-D50-D89", ]
  
  # Fixed palette used in the paper figure (order matters)
  pal <- c(
    "#969ae0", "#0ce846", "#f7d386", "#eda8d7",
    "#0e158a", "#277d28", "#8a2d80", "#e3099f", "#fca40a"
  )
  
  # Ensure we have enough colours for the chapters shown
  n_ch <- length(unique(umap_df$chapter))
  if (n_ch > length(pal)) {
    stop("UMAP subset contains ", n_ch, " chapters but palette has only ", length(pal),
         " colours. Extend 'pal' or reduce the subset.")
  }
  
  p <- plot_umap(umap_df, palette = pal, point_size = 3)
  print(p)
  
  # Save the UMAP figure in main_results/data/ so it is easy to find
  umap_png <- file.path(data_dir, "umap_chapter_subset.png")
  ggsave(umap_png, p, width = 8, height = 6, dpi = 300, bg = "white")
  message("\nSaved UMAP figure: ", umap_png)
  
  # ----------------------------------------------------------------------------
  # (2) LOCAL GRAPH EXAMPLE
  # ----------------------------------------------------------------------------
  # target_code controls which code is at the center of the graph.
  # ontologies controls which coding systems we include in neighbour retrieval.
  
  target_code <- "ICD-B169"  # example target used in the paper
  
  g_local <- build_local_neighbour_graph(
    embeddings    = embeddings,
    dict          = dict,
    target        = target_code,
    ontologies    = c("ICD", "ATC"),   # e.g. add "CAM", "LPP", "BIO" if desired
    k_per_onto    = 5,
    sim_threshold = 0.5
  )
  
  plot_local_graph(
    g_local,
    label_type = "code",
    main = paste("Local similarity graph around", target_code)
  )
}
