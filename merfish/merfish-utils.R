## ===========================================================================
## merfish-utils.R
## ===========================================================================
## Headless analysis primitives for the MERFISH batch pipeline. Ported from the
## Phenomenalist MERFISH Shiny module (merfish/app.R) and stripped of all Shiny
## reactivity so they can run unattended on a SLURM array task.
##
## Sourced by RunMerfish.R. Nothing here touches the filesystem or prints —
## these are pure compute helpers plus a couple of small IO wrappers.
##
## Dependencies: Matrix, RANN, igraph, ggplot2 (figures), viridis,
##               RColorBrewer, pheatmap; optional RSpectra (UMAP init).
##               SpatialExperiment / SingleCellExperiment for the .rds export.
## ===========================================================================

## ---------------------------------------------------------------------------
## IO ------------------------------------------------------------------------
## ---------------------------------------------------------------------------

# Read a delimited table, transparently handling .gz and tab/comma separators.
# Uses data.table::fread when available (fast, gz-aware), else base read.delim.
.read_table <- function(path, row_names = NULL) {
  ext <- tolower(tools::file_ext(sub("\\.gz$", "", path)))
  sep <- if (ext %in% c("tsv", "txt")) "\t" else ","
  if (requireNamespace("data.table", quietly = TRUE)) {
    df <- data.table::fread(path, sep = sep, header = TRUE, data.table = FALSE,
                            check.names = FALSE)
    if (!is.null(row_names)) {
      rn <- df[[row_names]]
      df[[row_names]] <- NULL
      rownames(df) <- as.character(rn)
    }
    return(df)
  }
  con <- if (grepl("\\.gz$", path)) gzfile(path) else path
  rn_arg <- if (is.null(row_names)) NULL else row_names
  read.delim(con, sep = sep, header = TRUE, row.names = rn_arg,
             check.names = FALSE, stringsAsFactors = FALSE)
}

# Load a MERFISH expression matrix as a numeric cells x genes matrix.
# First column is assumed to hold cell IDs (becomes rownames). If `transpose`
# the input is genes x cells and is flipped to cells x genes.
read_expression <- function(path, transpose = FALSE) {
  df <- .read_table(path, row_names = 1)
  m  <- as.matrix(df)
  storage.mode(m) <- "numeric"
  if (transpose) m <- t(m)
  m
}

# Load cell metadata; first column is the cell ID used to align with expression.
read_metadata <- function(path) {
  .read_table(path)
}

## ---------------------------------------------------------------------------
## QC ------------------------------------------------------------------------
## ---------------------------------------------------------------------------

# Per-cell QC statistics from a cells x genes count matrix.
qc_stats <- function(expr) {
  data.frame(
    cell         = rownames(expr),
    total_counts = rowSums(expr),
    n_genes      = rowSums(expr > 0),
    stringsAsFactors = FALSE
  )
}

# Build a logical keep-mask from QC thresholds. Any threshold left NULL/NA is
# skipped. area / negctrl filters require their metadata column to be present.
qc_filter <- function(expr, meta, params) {
  qs   <- qc_stats(expr)
  keep <- rep(TRUE, nrow(qs))

  if (!is.null(params$qc_min_counts)) keep <- keep & qs$total_counts >= params$qc_min_counts
  if (!is.null(params$qc_max_counts)) keep <- keep & qs$total_counts <= params$qc_max_counts
  if (!is.null(params$qc_min_genes))  keep <- keep & qs$n_genes      >= params$qc_min_genes
  if (!is.null(params$qc_max_genes))  keep <- keep & qs$n_genes      <= params$qc_max_genes

  if (!is.null(params$area_col) && params$area_col %in% colnames(meta)) {
    areas <- as.numeric(meta[[params$area_col]])
    if (!is.null(params$qc_min_area)) keep <- keep & areas >= params$qc_min_area
    if (!is.null(params$qc_max_area)) keep <- keep & areas <= params$qc_max_area
    if (!is.null(params$qc_max_density)) {
      dens <- qs$total_counts / (areas + 1e-10)
      keep <- keep & dens <= params$qc_max_density
    }
    if (!is.null(params$qc_min_density)) {
      dens <- qs$total_counts / (areas + 1e-10)
      keep <- keep & dens >= params$qc_min_density
    }
  }

  if (!is.null(params$negctrl_col) && params$negctrl_col %in% colnames(meta) &&
      !is.null(params$qc_max_negctrl_ratio)) {
    neg_ratio <- as.numeric(meta[[params$negctrl_col]]) / (qs$total_counts + 1e-10)
    keep <- keep & neg_ratio <= params$qc_max_negctrl_ratio
  }

  keep[is.na(keep)] <- FALSE
  keep
}

## ---------------------------------------------------------------------------
## Normalization -------------------------------------------------------------
## ---------------------------------------------------------------------------

# Normalize a cells x genes count matrix. `volumes` (optional) is a per-cell
# vector used by the cell-volume method.
normalize_expr <- function(expr, method = "lognorm", volumes = NULL) {
  n_cells <- nrow(expr); n_genes <- ncol(expr)
  if (method == "lognorm") {
    lib <- rowSums(expr); lib[lib == 0] <- 1
    norm <- log1p(sweep(expr, 1, lib, "/") * 10000)
  } else if (method == "cp10k") {
    lib <- rowSums(expr); lib[lib == 0] <- 1
    norm <- sweep(expr, 1, lib, "/") * 10000
  } else if (method == "cellvol") {
    if (is.null(volumes)) {
      warning("cellvol normalization requested but no volume column; using lognorm")
      lib <- rowSums(expr); lib[lib == 0] <- 1
      norm <- log1p(sweep(expr, 1, lib, "/") * 10000)
    } else {
      volumes[volumes <= 0 | is.na(volumes)] <- 1e-6
      norm <- log1p(sweep(expr, 1, volumes, "/") * median(volumes))
    }
  } else if (method == "none") {
    norm <- expr
  } else { # sct-like Pearson residuals
    lib       <- rowSums(expr)
    gene_mean <- colMeans(expr)
    theta     <- 100
    norm <- matrix(0, nrow = n_cells, ncol = n_genes)
    mean_lib <- mean(lib)
    for (j in seq_len(n_genes)) {
      mu <- (lib * gene_mean[j]) / mean_lib
      norm[, j] <- (expr[, j] - mu) / sqrt(mu + mu^2 / theta)
    }
    clip <- sqrt(n_cells)
    norm[norm >  clip] <-  clip
    norm[norm < -clip] <- -clip
  }
  dimnames(norm) <- dimnames(expr)
  norm
}

## ---------------------------------------------------------------------------
## Feature selection ---------------------------------------------------------
## ---------------------------------------------------------------------------

# Return the names of the selected highly variable genes.
select_hvg <- function(norm, method = "variance", n_hvg = 2000) {
  n_genes <- ncol(norm)
  if (method == "all") return(colnames(norm))
  if (method == "vst") {
    gm <- colMeans(norm); gv <- apply(norm, 2, var)
    score <- tryCatch({
      lo  <- loess(log1p(gv) ~ log1p(gm))
      fit <- exp(predict(lo)) - 1
      gv / (fit + 1e-10)
    }, error = function(e) gv)
  } else {
    score <- apply(norm, 2, var)
  }
  k <- min(n_hvg, n_genes)
  colnames(norm)[order(score, decreasing = TRUE)[seq_len(k)]]
}

## ---------------------------------------------------------------------------
## Dimensionality reduction --------------------------------------------------
## ---------------------------------------------------------------------------

# PCA on HVG-subset, scaled per the chosen method. Returns scores + sdev.
run_pca <- function(norm, hvg, scale_method = "zscore", n_pcs = 20) {
  x <- norm[, hvg, drop = FALSE]
  if (scale_method == "zscore") {
    x <- scale(x, center = TRUE, scale = TRUE)
  } else if (scale_method == "center") {
    x <- scale(x, center = TRUE, scale = FALSE)
  }
  x[is.nan(x)] <- 0
  n_pcs <- min(n_pcs, ncol(x) - 1, nrow(x) - 1)
  pca   <- prcomp(x, center = FALSE, scale. = FALSE, rank. = n_pcs)
  list(scores = pca$x[, seq_len(n_pcs), drop = FALSE],
       sdev   = pca$sdev[seq_len(n_pcs)])
}

# Lightweight UMAP-style 2D embedding: fuzzy KNN graph + normalized-Laplacian
# spectral layout (RSpectra if available, else classical MDS fallback). This is
# the same portable approximation the Shiny module uses.
run_umap <- function(pca_scores, n_neighbors = 15, min_dist = 0.3, seed = 1) {
  set.seed(seed)
  nn       <- RANN::nn2(pca_scores, k = n_neighbors + 1)
  nn_idx   <- nn$nn.idx[, -1, drop = FALSE]
  nn_dists <- nn$nn.dists[, -1, drop = FALSE]
  n        <- nrow(pca_scores)

  sigma <- apply(nn_dists, 1, function(d) d[min(n_neighbors, length(d))])
  sigma[sigma < 1e-8] <- 1e-8

  trip_i <- rep(seq_len(n), each = ncol(nn_idx))
  trip_j <- as.vector(t(nn_idx))
  trip_v <- as.vector(t(exp(-nn_dists / sigma)))
  W <- Matrix::sparseMatrix(i = trip_i, j = trip_j, x = trip_v, dims = c(n, n))
  W <- (W + Matrix::t(W)) / 2

  D_inv_sqrt <- Matrix::Diagonal(x = 1 / sqrt(Matrix::rowSums(W) + 1e-10))
  L_norm     <- D_inv_sqrt %*% W %*% D_inv_sqrt

  emb <- tryCatch({
    if (!requireNamespace("RSpectra", quietly = TRUE)) stop("no RSpectra")
    eig <- RSpectra::eigs_sym(L_norm, k = 3, which = "LM")
    eig$vectors[, 2:3, drop = FALSE]
  }, error = function(e) {
    d <- dist(pca_scores[, seq_len(min(5, ncol(pca_scores))), drop = FALSE])
    cmdscale(d, k = 2)
  })
  emb <- emb + matrix(rnorm(n * 2, 0, min_dist * 0.1), ncol = 2)
  colnames(emb) <- c("UMAP_1", "UMAP_2")
  list(embedding = emb, knn_idx = nn_idx, knn_dists = nn_dists)
}

## ---------------------------------------------------------------------------
## Clustering ----------------------------------------------------------------
## ---------------------------------------------------------------------------

# Shared-nearest-neighbor graph + Leiden community detection.
snn_cluster <- function(knn_idx, resolution = 0.8,
                        objective = c("modularity", "CPM")) {
  objective <- match.arg(objective)
  n <- nrow(knn_idx); k <- ncol(knn_idx)
  edge_list <- vector("list", n)
  for (i in seq_len(n)) {
    for (j in knn_idx[i, ]) {
      if (j > i) {
        shared <- length(intersect(knn_idx[i, ], knn_idx[j, ]))
        if (shared > 0) {
          edge_list[[i]] <- rbind(edge_list[[i]],
                                  data.frame(from = i, to = j, weight = shared / k))
        }
      }
    }
  }
  edges <- do.call(rbind, edge_list)
  if (is.null(edges) || nrow(edges) == 0) return(rep(1L, n))
  g  <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = seq_len(n))
  cl <- igraph::cluster_leiden(g, weights = igraph::E(g)$weight,
                               resolution = resolution,
                               objective_function = objective, n_iterations = 10)
  igraph::membership(cl)
}

## ---------------------------------------------------------------------------
## Spatial statistics --------------------------------------------------------
## ---------------------------------------------------------------------------

# Cluster x cluster neighborhood enrichment Z-scores over a spatial KNN graph.
neighborhood_enrichment <- function(clusters, nn_idx) {
  n_cells   <- length(clusters)
  cl_levels <- sort(unique(clusters))
  n_cl      <- length(cl_levels)
  obs <- matrix(0, n_cl, n_cl, dimnames = list(cl_levels, cl_levels))
  for (i in seq_len(n_cells)) {
    for (cl in clusters[nn_idx[i, ]]) {
      obs[clusters[i], cl] <- obs[clusters[i], cl] + 1
    }
  }
  cl_sizes <- table(clusters)[cl_levels]
  total    <- sum(obs)
  exp_mat  <- outer(cl_sizes, cl_sizes) / sum(cl_sizes)^2 * total
  (obs - exp_mat) / sqrt(exp_mat + 1e-10)
}

# Fast Moran's I (distance-weighted) for one gene over a spatial KNN graph.
morans_i_fast <- function(expr_vec, nn_idx, nn_dists) {
  n  <- length(expr_vec)
  x  <- expr_vec - mean(expr_vec)
  ss <- sum(x^2)
  if (ss < 1e-15) return(list(I = 0, p = 1))
  W_sum <- 0; lag_sum <- 0
  for (i in seq_len(n)) {
    w <- 1 / (nn_dists[i, ] + 1e-6)
    W_sum   <- W_sum + sum(w)
    lag_sum <- lag_sum + sum(w * x[nn_idx[i, ]])
  }
  I  <- (n / W_sum) * (lag_sum / ss)
  EI <- -1 / (n - 1)
  z  <- (I - EI) / (0.5 / sqrt(n))
  list(I = I, p = 2 * pnorm(-abs(z)))
}

# Detect spatially variable genes among the top-variance genes.
detect_svg <- function(norm, coords, k = 10, n_top = 200) {
  nn       <- RANN::nn2(coords, k = k + 1)
  nn_idx   <- nn$nn.idx[, -1, drop = FALSE]
  nn_dists <- nn$nn.dists[, -1, drop = FALSE]
  gv       <- apply(norm, 2, var)
  n_test   <- min(n_top, ncol(norm))
  genes    <- colnames(norm)[order(gv, decreasing = TRUE)[seq_len(n_test)]]
  res <- data.frame(gene = genes, morans_I = NA_real_, p_value = NA_real_,
                    stringsAsFactors = FALSE)
  for (i in seq_along(genes)) {
    mi <- morans_i_fast(norm[, genes[i]], nn_idx, nn_dists)
    res$morans_I[i] <- mi$I
    res$p_value[i]  <- mi$p
  }
  res$p_adj <- p.adjust(res$p_value, method = "BH")
  res[order(res$p_adj, -res$morans_I), ]
}

## ---------------------------------------------------------------------------
## Differential expression ---------------------------------------------------
## ---------------------------------------------------------------------------

# Wilcoxon DE of group1 vs group2 over every gene.
run_de <- function(expr_mat, g1_idx, g2_idx, gene_names) {
  res <- data.frame(gene = gene_names, avg_log2FC = NA_real_,
                    pct_1 = NA_real_, pct_2 = NA_real_,
                    p_val = NA_real_, p_adj = NA_real_, stringsAsFactors = FALSE)
  for (i in seq_along(gene_names)) {
    g1 <- expr_mat[g1_idx, i]; g2 <- expr_mat[g2_idx, i]
    res$avg_log2FC[i] <- log2((mean(g1) + 1e-9) / (mean(g2) + 1e-9))
    res$pct_1[i] <- mean(g1 > 0); res$pct_2[i] <- mean(g2 > 0)
    res$p_val[i] <- if (sd(c(g1, g2)) > 1e-15)
      tryCatch(wilcox.test(g1, g2)$p.value, error = function(e) 1) else 1
  }
  res$p_adj <- p.adjust(res$p_val, method = "BH")
  res[order(res$p_adj, -abs(res$avg_log2FC)), ]
}

# One-vs-rest marker genes for every cluster, stacked into one table. This is
# the headless analog of the app's interactive pairwise DE.
cluster_markers <- function(norm, clusters, min_cells = 3) {
  cl_levels <- sort(unique(clusters))
  out <- list()
  for (cl in cl_levels) {
    g1 <- which(clusters == cl); g2 <- which(clusters != cl)
    if (length(g1) < min_cells || length(g2) < min_cells) next
    de <- run_de(norm, g1, g2, colnames(norm))
    de$cluster <- cl
    out[[cl]] <- de
  }
  if (length(out) == 0) return(NULL)
  do.call(rbind, out)
}

## ---------------------------------------------------------------------------
## Figures (dark Phenomenalist style, headless) ------------------------------
## ---------------------------------------------------------------------------

phen_palette_discrete <- function(n) {
  base <- c("#00c2ff","#e040fb","#00e676","#ffab00","#ff5252","#18ffff","#b388ff",
            "#69f0ae","#ffd740","#ff8a80","#40c4ff","#ea80fc","#b9f6ca","#ffe57f",
            "#ff867c","#80d8ff","#ce93d8","#a5d6a7","#fff176","#ef9a9a")
  if (n <= length(base)) base[seq_len(n)] else grDevices::colorRampPalette(base)(n)
}

theme_phenomenalist <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) %+replace% ggplot2::theme(
    text             = ggplot2::element_text(color = "#e8eaed"),
    plot.title       = ggplot2::element_text(face = "bold", color = "#e8eaed"),
    plot.background  = ggplot2::element_rect(fill = "#12151c", color = NA),
    panel.background = ggplot2::element_rect(fill = "#12151c", color = NA),
    panel.grid.major = ggplot2::element_line(color = "#FFFFFF0A", linewidth = 0.3),
    panel.grid.minor = ggplot2::element_blank(),
    axis.text        = ggplot2::element_text(color = "#8b919e"),
    axis.title       = ggplot2::element_text(color = "#8b919e", face = "bold"),
    legend.text      = ggplot2::element_text(color = "#8b919e"),
    legend.title     = ggplot2::element_text(color = "#e8eaed", face = "bold")
  )
}

# Write the standard figure panel set into <out>/figures.
write_figures <- function(out_dir, qs, keep, coords_filt, umap, clusters,
                          width = 8, height = 6, dpi = 300) {
  fig_dir <- file.path(out_dir, "figures")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  save <- function(name, p) ggplot2::ggsave(file.path(fig_dir, name), p,
                                            width = width, height = height,
                                            dpi = dpi, bg = "#12151c")

  qdf <- data.frame(
    value  = c(qs$total_counts, qs$n_genes),
    metric = rep(c("Total Counts", "Genes Detected"), each = nrow(qs)),
    pass   = rep(keep, 2))
  save("qc_violin.pdf",
       ggplot2::ggplot(qdf, ggplot2::aes(metric, value, fill = pass)) +
         ggplot2::geom_violin(alpha = 0.7, scale = "width", color = NA) +
         ggplot2::scale_fill_manual(values = c("TRUE" = "#00c2ff", "FALSE" = "#ff5252")) +
         ggplot2::labs(title = "QC Distributions", x = NULL, y = "Value") +
         theme_phenomenalist())

  if (!is.null(clusters)) {
    n_cl <- length(unique(clusters))
    sdf  <- data.frame(x = coords_filt[, 1], y = coords_filt[, 2],
                       cluster = factor(clusters))
    save("spatial_clusters.pdf",
         ggplot2::ggplot(sdf, ggplot2::aes(x, y, color = cluster)) +
           ggplot2::geom_point(size = 0.5, alpha = 0.7) +
           ggplot2::scale_color_manual(values = phen_palette_discrete(n_cl)) +
           ggplot2::coord_fixed() + ggplot2::scale_y_reverse() +
           ggplot2::labs(title = "Tissue Space — Clusters", x = "X (µm)", y = "Y (µm)") +
           theme_phenomenalist() +
           ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3))))
    if (!is.null(umap)) {
      udf <- data.frame(UMAP_1 = umap[, 1], UMAP_2 = umap[, 2], cluster = factor(clusters))
      save("umap_clusters.pdf",
           ggplot2::ggplot(udf, ggplot2::aes(UMAP_1, UMAP_2, color = cluster)) +
             ggplot2::geom_point(size = 0.5, alpha = 0.7) +
             ggplot2::scale_color_manual(values = phen_palette_discrete(n_cl)) +
             ggplot2::coord_fixed() +
             ggplot2::labs(title = "UMAP — Clusters") +
             theme_phenomenalist() +
             ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3))))
    }
  }
  invisible(fig_dir)
}
