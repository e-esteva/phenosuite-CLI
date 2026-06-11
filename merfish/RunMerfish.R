## ===========================================================================
## RunMerfish.R
## ===========================================================================
## Orchestrates the headless MERFISH spatial-transcriptomics pipeline for a
## single sample:
##
##   Import -> QC -> Normalize -> HVG -> PCA -> UMAP -> Leiden cluster ->
##   Neighborhood enrichment -> Spatially variable genes -> Cluster markers ->
##   Export (CSVs, SpatialExperiment .rds, figures, run manifest)
##
## Sourced by run-merfish.R after merfish-utils.R. Exposes a single entry point,
## RunMerfish(), driven entirely by an explicit parameter list (no globals).
## ===========================================================================

# `params` is the validated option list assembled by run-merfish.R. Required
# fields: expression_file, metadata_file, out_dir, x_col, y_col. Everything
# else has a default applied upstream.
RunMerfish <- function(params) {

  t0  <- Sys.time()
  log <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
  out <- params$out_dir
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  sample_id <- params$sample_id

  ## ---- 1. Import -----------------------------------------------------------
  log("Reading expression: ", params$expression_file)
  expr <- read_expression(params$expression_file, transpose = params$transpose)
  log("Reading metadata:   ", params$metadata_file)
  meta <- read_metadata(params$metadata_file)

  if (!params$x_col %in% colnames(meta)) stop("x_col '", params$x_col, "' not in metadata")
  if (!params$y_col %in% colnames(meta)) stop("y_col '", params$y_col, "' not in metadata")

  # Align cells: expression rownames against the metadata's first (ID) column.
  meta_ids <- as.character(meta[[1]])
  common   <- intersect(rownames(expr), meta_ids)
  if (length(common) < 10) {
    stop("fewer than 10 cells matched between expression rownames and metadata ",
         "column '", colnames(meta)[1], "'. Check cell IDs / --transpose.")
  }
  expr <- expr[common, , drop = FALSE]
  meta <- meta[match(common, meta_ids), , drop = FALSE]
  coords <- cbind(x = as.numeric(meta[[params$x_col]]),
                  y = as.numeric(meta[[params$y_col]]))
  log(sprintf("Loaded %d cells x %d genes (%d aligned).",
              nrow(expr), ncol(expr), length(common)))

  ## ---- 2. QC ---------------------------------------------------------------
  qs   <- qc_stats(expr)
  keep <- qc_filter(expr, meta, params)
  log(sprintf("QC: %d / %d cells pass (%.1f%%).",
              sum(keep), length(keep), 100 * mean(keep)))
  if (sum(keep) < 10) stop("fewer than 10 cells pass QC; loosen thresholds.")

  expr_f   <- expr[keep, , drop = FALSE]
  meta_f   <- meta[keep, , drop = FALSE]
  coords_f <- coords[keep, , drop = FALSE]

  ## ---- 3. Normalize --------------------------------------------------------
  volumes <- NULL
  if (identical(params$norm_method, "cellvol") &&
      !is.null(params$vol_col) && params$vol_col %in% colnames(meta_f)) {
    volumes <- as.numeric(meta_f[[params$vol_col]])
  }
  log("Normalizing (", params$norm_method, ")")
  norm <- normalize_expr(expr_f, method = params$norm_method, volumes = volumes)

  ## ---- 4. HVG + 5. PCA + 6. UMAP ------------------------------------------
  log("Selecting HVGs (", params$hvg_method, ", n=", params$n_hvg, ")")
  hvg <- select_hvg(norm, method = params$hvg_method, n_hvg = params$n_hvg)
  log("PCA (", params$n_pcs, " PCs, scale=", params$scale_method, ")")
  pca <- run_pca(norm, hvg, scale_method = params$scale_method, n_pcs = params$n_pcs)
  log("UMAP embedding (k=", params$umap_neighbors, ")")
  um  <- run_umap(pca$scores, n_neighbors = params$umap_neighbors,
                  min_dist = params$umap_min_dist, seed = params$seed)

  ## ---- 7. Clustering -------------------------------------------------------
  log("Leiden clustering (k=", params$cluster_k, ", res=", params$cluster_res,
      ", obj=", params$leiden_objective, ")")
  nn_cl    <- RANN::nn2(pca$scores, k = params$cluster_k + 1)
  knn_idx  <- nn_cl$nn.idx[, -1, drop = FALSE]
  clusters <- as.character(snn_cluster(knn_idx, resolution = params$cluster_res,
                                       objective = params$leiden_objective))
  log(sprintf("Found %d clusters.", length(unique(clusters))))

  ## ---- 8. Neighborhood enrichment -----------------------------------------
  log("Neighborhood enrichment (k=", params$nhood_k, ")")
  nn_sp   <- RANN::nn2(coords_f, k = params$nhood_k + 1)
  nhood_z <- neighborhood_enrichment(clusters, nn_sp$nn.idx[, -1, drop = FALSE])

  ## ---- 9. Spatially variable genes ----------------------------------------
  log("Detecting SVGs (k=", params$svg_k, ", top=", params$svg_n_top, ")")
  svg <- detect_svg(norm, coords_f, k = params$svg_k, n_top = params$svg_n_top)
  log(sprintf("  %d SVGs at p_adj < 0.05.", sum(svg$p_adj < 0.05, na.rm = TRUE)))

  ## ---- 10. Cluster markers (one-vs-rest DE) -------------------------------
  log("Cluster markers (one-vs-rest Wilcoxon)")
  markers <- cluster_markers(norm, clusters)

  ## ---- 11. Export ----------------------------------------------------------
  pre <- function(name) file.path(out, paste0(sample_id, "_", name))

  meta_out <- meta_f
  meta_out$cluster   <- clusters
  meta_out$Phenotype <- clusters            # PCF / downstream compatibility
  meta_out$UMAP_1    <- um$embedding[, 1]
  meta_out$UMAP_2    <- um$embedding[, 2]
  write.csv(meta_out,                pre("metadata_clusters.csv"), row.names = FALSE)
  write.csv(norm,                    pre("normalized_expression.csv"))
  write.csv(svg,                     pre("svg.csv"), row.names = FALSE)
  write.csv(as.data.frame(nhood_z),  pre("neighborhood_enrichment.csv"))
  if (!is.null(markers)) write.csv(markers, pre("cluster_markers.csv"), row.names = FALSE)

  if (isTRUE(params$export_spe)) {
    log("Exporting SpatialExperiment (.rds)")
    ok <- requireNamespace("SpatialExperiment", quietly = TRUE) &&
          requireNamespace("S4Vectors", quietly = TRUE)
    if (ok) {
      expr_mat   <- t(norm)
      assay_list <- list(exprs = expr_mat, logcounts = expr_mat,
                         counts = t(expr_f))
      sp_coords  <- coords_f; colnames(sp_coords) <- c("x", "y")
      cd <- meta_f
      cd$cluster <- clusters; cd$cluster_merfish <- clusters; cd$Phenotype <- clusters
      spe <- SpatialExperiment::SpatialExperiment(
        assays        = assay_list,
        colData       = S4Vectors::DataFrame(cd),
        spatialCoords = sp_coords)
      SingleCellExperiment::reducedDim(spe, "UMAP") <- um$embedding
      SingleCellExperiment::reducedDim(spe, "PCA")  <- pca$scores
      saveRDS(spe, pre("spe.rds"))
    } else {
      log("  SpatialExperiment/S4Vectors unavailable — skipping .rds export.")
    }
  }

  if (isTRUE(params$export_figures)) {
    log("Writing figures")
    tryCatch(
      write_figures(out, qs[keep, , drop = FALSE], keep[keep],
                    coords_f, um$embedding, clusters,
                    width = params$fig_width, height = params$fig_height,
                    dpi = params$fig_dpi),
      error = function(e) log("  figure export failed: ", conditionMessage(e)))
  }

  ## ---- Run manifest --------------------------------------------------------
  manifest <- list(
    sample_id        = sample_id,
    expression_file  = normalizePath(params$expression_file),
    metadata_file    = normalizePath(params$metadata_file),
    out_dir          = normalizePath(out),
    n_cells_raw      = nrow(expr),
    n_genes          = ncol(expr),
    n_cells_pass_qc  = sum(keep),
    n_clusters       = length(unique(clusters)),
    n_svg_sig        = sum(svg$p_adj < 0.05, na.rm = TRUE),
    parameters       = params[setdiff(names(params), c())],
    git_sha          = Sys.getenv("PHENOSUITE_GIT_SHA", "unknown"),
    started          = format(t0, "%Y-%m-%d %H:%M:%S"),
    finished         = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    elapsed_seconds  = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1),
    R_version        = R.version.string
  )
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(manifest, pre("run_manifest.json"),
                         auto_unbox = TRUE, pretty = TRUE, null = "null")
  } else {
    saveRDS(manifest, pre("run_manifest.rds"))
  }

  log(sprintf("Done: %s  (%.1fs)  ->  %s",
              sample_id, manifest$elapsed_seconds, out))
  invisible(manifest)
}
