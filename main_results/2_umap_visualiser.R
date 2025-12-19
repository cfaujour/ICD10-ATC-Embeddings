# main_results/02_umap_graph.R
# Explore released embeddings:
#  (1) UMAP 2D projection for a selected subset (colour = chapter, shape = ontology)
#  (2) Local similarity graph around a target code (small subgraph)
#
# Intended usage: run from REPOSITORY ROOT (repo-relative paths).

suppressPackageStartupMessages({
  library(data.table)
  library(umap)
  library(ggplot2)
  library(igraph)
})

# ---- Paths (repo-relative; run from repo root) ----
data_dir <- file.path("main_results", "data")
func_dir <- file.path("main_results", "functions")


ensure_dir <- function(path) if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(data_dir) || !dir.exists(func_dir)) {
  stop("Please run this script from the repository root (cannot find 'main_results/data' or 'main_results/functions').")
}


# ---- Helpers ----
source(file.path(func_dir, "similarity_measures.R"))  # provides cosine_dist()
source(file.path(func_dir, "get_label.R"))            # get_label() may use label_table

# NOTE: cosine_dist() is expected to return cosine SIMILARITY (higher = more similar).

# Global label table (kept for minimal changes if get_label.R expects it)
label_table <- fread(file.path(data_dir, "label_table.csv"))

# ---- Load embeddings + vocab ----
load_embeddings <- function(embeddings_path = file.path(data_dir, "embeddings_ESND_2FC.csv.gz"),
                            vocab_path      = file.path(data_dir, "vocab.csv"),
                            verbose         = TRUE) {
  if (verbose) message("Loading embeddings: ", embeddings_path)
  embeddings <- as.matrix(fread(embeddings_path))
  
  if (verbose) message("Loading vocab: ", vocab_path)
  vocab <- fread(vocab_path)
  
  if (all(c("index", "Event_code") %in% names(vocab))) {
    vocab <- vocab[order(index)]
    dict <- vocab$Event_code
  } else {
    warning("vocab.csv missing 'index'/'Event_code'; falling back to 2nd column as code.")
    dict <- vocab[[2]]
  }
  
  if (nrow(embeddings) != length(dict)) {
    stop("Mismatch: embeddings rows (", nrow(embeddings), ") vs vocab size (", length(dict), ").")
  }
  
  rownames(embeddings) <- dict
  if (verbose) message("Loaded ", nrow(embeddings), " codes x ", ncol(embeddings), " dims.")
  list(embeddings = embeddings, dict = dict)
}

# ---- Chapter mapping (colouring) ----
chapter_table <- fread(file.path(data_dir, "chapter_table.csv"))
if (!all(c("code", "y") %in% names(chapter_table))) {
  stop("chapter_table.csv must contain columns: 'code' and 'y' (chapter label).")
}

chapter_map <- setNames(chapter_table$y, chapter_table$code)

get_chapter <- function(code) {
  val <- chapter_map[[code]]
  if (is.null(val)) NA_character_ else val
}


# ---- UMAP on subset ----
run_umap_subset <- function(embeddings,
                            dict,
                            codes_subset,
                            n_neighbors = 25,
                            min_dist    = 0.5,
                            seed        = 16,
                            verbose     = TRUE) {
  idx <- which(dict %in% codes_subset)
  if (length(idx) == 0) stop("No codes from codes_subset found in vocabulary.")
  
  dict_sub <- dict[idx]
  emb_sub  <- embeddings[idx, , drop = FALSE]
  
  if (verbose) message("UMAP on ", length(dict_sub), " codes.")
  set.seed(seed)
  
  cfg <- umap::umap.defaults
  cfg$metric      <- "cosine"
  cfg$n_neighbors <- n_neighbors
  cfg$min_dist    <- min_dist
  
  u <- umap::umap(emb_sub, n_components = 2, config = cfg)
  
  df <- data.frame(
    UMAP1   = u$layout[, 1],
    UMAP2   = u$layout[, 2],
    code    = dict_sub,
    onto    = substr(dict_sub, 1, 3),
    chapter = vapply(dict_sub, get_chapter, character(1)),
    stringsAsFactors = FALSE
  )
  
  df
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

# ---- Local graph around a target code ----
build_local_neighbour_graph <- function(embeddings,
                                        dict,
                                        target,
                                        ontologies     = c("ICD", "ATC"),
                                        k_per_onto     = 5,
                                        sim_threshold  = 0.5,
                                        verbose        = TRUE) {
  if (!target %in% dict) stop("Target code '", target, "' not found in embeddings.")
  
  if (verbose) message("Local graph around: ", target)
  
  neighbours <- character(0)
  
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
    neighbours <- c(neighbours, names(sims)[seq_len(min(k_per_onto, length(sims)))])
  }
  
  nodes <- unique(c(target, neighbours))
  emb_g <- embeddings[match(nodes, dict), , drop = FALSE]
  
  n <- nrow(emb_g)
  sim_mat <- matrix(0, n, n, dimnames = list(nodes, nodes))
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      sim_mat[i, j] <- cosine_dist(emb_g[i, ], emb_g[j, ], sparse = FALSE)
    }
  }
  sim_mat[sim_mat < sim_threshold] <- 0
  diag(sim_mat) <- 0
  
  g <- graph_from_adjacency_matrix(sim_mat, mode = "undirected", weighted = TRUE, diag = FALSE)
  g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  
  V(g)$code     <- V(g)$name
  V(g)$ontology <- substr(V(g)$code, 1, 3)
  V(g)$ontology[V(g)$code == target] <- "Target"
  V(g)$label    <- vapply(V(g)$code, get_label, character(1))
  
  g
}

plot_local_graph <- function(g,
                             layout_fun  = layout_with_fr,
                             label_type  = c("code", "label"),
                             vertex_size = 16,
                             seed        = 42,
                             main        = "Local similarity graph") {
  label_type <- match.arg(label_type)
  set.seed(seed)
  lay <- layout_fun(g)
  
  base_cols <- c(
    Target = "#e41a1c",
    ICD    = "#377eb8",
    ATC    = "#4daf4a",
    Other  = "grey50"
  )
  
  ontos <- V(g)$ontology
  v_col <- vapply(ontos, function(o) if (!is.null(base_cols[[o]])) base_cols[[o]] else base_cols[["Other"]], character(1))
  is_target <- ontos == "Target"
  
  v_lab <- if (label_type == "code") V(g)$code else V(g)$label
  v_sz  <- ifelse(is_target, vertex_size * 1.3, vertex_size)
  v_fr  <- ifelse(is_target, "black", NA)
  
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

# ================================================================
# USAGE (runs only when executed as a script)
# ================================================================
if (sys.nframe() == 0L) {
  
  # --- Load once ---
  emb_data   <- load_embeddings()
  embeddings <- emb_data$embeddings
  dict       <- emb_data$dict
  
  # --- Subset selection  ---
  # Current subset matches the mains paper figure
  set.seed(13)
  codes_subset <- dict[substr(dict, 1, 5) %in% c(
    "ICD-C", "ICD-D", "ATC-L", "ICD-G", "ATC-N", "ICD-H", "ATC-S", "ICD-J", "ATC-R"
  )]
  
  # --- UMAP parameters (match paper) ---
  umap_df <- run_umap_subset(
    embeddings   = embeddings,
    dict         = dict,
    codes_subset = codes_subset,
    n_neighbors  = 25,
    min_dist     = 0.5,
    seed         = 16
  )
  
  # Match paper filtering exactly 
  umap_df <- subset(umap_df, !chapter %in% c("ICD-NA", "NA"))
  umap_df <- umap_df[umap_df$chapter != "ICD-D50-D89", ]
  
  # Fixed palette (add colors if you have more groups)
  
  pal <- c(
    "#969ae0",
    "#0ce846",  
    "#f7d386",  
    "#eda8d7",  
    "#0e158a",  
    "#277d28",
    "#8a2d80",
    "#e3099f",
    "#fca40a"
  )
  
  
  # Guardrail: ensure palette is long enough for chapters present
  n_ch <- length(unique(umap_df$chapter))
  if (n_ch > length(pal)) {
    stop("Palette has ", length(pal), " colors but UMAP subset contains ", n_ch, " chapters. Extend 'pal' or reduce subset.")
  }
  
  p <- plot_umap(umap_df, palette = pal, point_size = 3)
  print(p)
  
  umap_png <- file.path(data_dir, "umap_chapter_subset.png")
  ggsave(umap_png, p, width = 8, height = 6, dpi = 300, bg = "white")
  message("Saved UMAP figure: ", umap_png)
  
  # --- Local graph example ---
  target_code <- "ICD-B169"  # example target used in the paper
  g_local <- build_local_neighbour_graph(
    embeddings    = embeddings,
    dict          = dict,
    target        = target_code,
    ontologies    = c("ICD", "ATC"),
    k_per_onto    = 5,
    sim_threshold = 0.5
  )
  
  plot_local_graph(
    g_local,
    label_type = "code",
    main = paste("Local similarity graph around", target_code)
  )
}
