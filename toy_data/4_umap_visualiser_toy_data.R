# Visualize code embeddings with UMAP or t-SNE, colored by concept/topic.
#
# Inputs:
#   1) embeddings_path : CSV.GZ produced by 02_ppmi_svd.R (embeddings.csv.gz)
#   2) vocab_path      : vocab.csv (Event_code, index, freq)
#   3) code_topics_path: code_topics.csv (Event_code, onto, topic)
#   4) out_dir         : directory for figures
#   5) method          : "umap" or "tsne" (default: "umap")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(uwot)          # for UMAP
  library(Rtsne)         # for t-SNE
  library(plotly)        # 3D / interactive plots
  library(htmlwidgets)   # to save HTML widgets
})

ensure_dir <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

# -------------------- CLI args & defaults --------------------

args <- commandArgs(trailingOnly = TRUE)

data_dir <- file.path("toy_data", "data")

emb_path         <- ifelse(length(args) >= 1, args[1], file.path(data_dir, "embeddings.csv.gz"))
vocab_path       <- ifelse(length(args) >= 2, args[2], file.path(data_dir, "vocab.csv"))
code_topics_path <- ifelse(length(args) >= 3, args[3], file.path(data_dir, "code_topics.csv"))
out_dir          <- ifelse(length(args) >= 4, args[4], data_dir)
method           <- ifelse(length(args) >= 5, args[5], "umap")
dims             <- ifelse(length(args) >= 6, as.integer(args[6]), 2L)  # 3D or 2D mapping

method <- tolower(method)
if (!method %in% c("umap", "tsne")) {
  stop("Method must be 'umap' or 'tsne'. Got: ", method)
}

if (!dims %in% c(2L, 3L)) {
  stop("dims must be 2 or 3. Got: ", dims)
}

message("Embeddings:   ", emb_path)
message("Vocab:        ", vocab_path)
message("Code topics:  ", code_topics_path)
message("Output dir:   ", out_dir)
message("Method:       ", method)
message("Dims:         ", dims)

ensure_dir(out_dir)

# -------------------- Load data & join labels --------------------

emb_dt <- fread(emb_path)
if (ncol(emb_dt) < 2) stop("Embeddings file seems to have too few columns.")

emb_dt[, index := .I]

vocab <- fread(vocab_path)

code_topics <- fread(code_topics_path)
code_topics <- unique(code_topics[, .(Event_code, topic)])

dt <- merge(emb_dt, vocab[, .(index, Event_code, freq)], by = "index", all.x = TRUE)
dt <- merge(dt, code_topics, by = "Event_code", all.x = TRUE)

dt[is.na(topic), topic := "background"]

embed_cols <- grep("^V[0-9]+$", names(dt), value = TRUE)
X <- as.matrix(dt[, ..embed_cols])

# -------------------- Dimensionality reduction --------------------

set.seed(42)

if (method == "umap") {
  message("Running UMAP ...")
  umap_res <- uwot::umap(
    X,
    n_neighbors = 15,
    min_dist = 0.1,
    n_components = dims,
    metric = "euclidean",
    verbose = TRUE
  )
  for (d in seq_len(dims)) {
    dt[[paste0("dim", d)]] <- umap_res[, d]
  }
  coord_label <- "UMAP"
} else {
  message("Running t-SNE ...")
  tsne_res <- Rtsne(
    X,
    dims = dims,
    perplexity = min(30, floor(nrow(X) / 4)),
    verbose = TRUE,
    check_duplicates = FALSE
  )
  for (d in seq_len(dims)) {
    dt[[paste0("dim", d)]] <- tsne_res$Y[, d]
  }
  coord_label <- "t-SNE"
}

# -------------------- Plot --------------------
if (dims == 2L) {
  p <- ggplot(dt, aes(x = dim1, y = dim2, color = topic)) +
    geom_point(alpha = 0.8, size = 2) +
    labs(
      title = paste(coord_label, "of code embeddings"),
      subtitle = "Points = codes, colored by assigned concept/topic",
      x = paste0(coord_label, " 1"),
      y = paste0(coord_label, " 2"),
      color = "Concept"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )
}

if (dims == 3L) {
  p <- plot_ly(
    dt,
    x = ~dim1, y = ~dim2, z = ~dim3,
    color = ~topic,
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 3, opacity = 0.8)
  ) |>
    plotly::layout(
      title = paste(coord_label, "of code embeddings (3D)"),
      scene = list(
        xaxis = list(title = paste0(coord_label, " 1")),
        yaxis = list(title = paste0(coord_label, " 2")),
        zaxis = list(title = paste0(coord_label, " 3"))
      )
    )
  
}

#Saving figure/widget

base <- paste0(coord_label, "_", method, "_", dims, "d")
png_path  <- file.path(out_dir, paste0(base, ".png"))
html_path <- file.path(out_dir, paste0(base, ".html"))

if (dims == 2L) {
  ggsave(filename = png_path, plot = p, width = 8, height = 6, dpi = 200, bg= "white")

}

if (dims == 3L) {
  htmlwidgets::saveWidget(p, file = html_path, selfcontained = TRUE)
}

message("Saved outputs:")
if (dims == 2L) {
  message("  - ", png_path)
} else {
  message("  - ", html_path)
}