# neighborhood_analysis-utils.R
# ─────────────────────────────────────────────────────────────────────────────
# Core analysis functions for NeighborhoodR (KNN niche matrix + LOO stability
# sweep + MiniBatchKMeans neighborhood assignment).
#
# Sourced by run-neighborhood_analysis.R (CLI).
# All heavy compute runs in Python when sklearn/scipy/numpy are available;
# RANN + base kmeans used as a pure-R fallback.
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(reticulate)
  library(SpatialExperiment)
  library(SingleCellExperiment)
  library(RANN)
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(jsonlite)
  library(glue)
})

`%||%` <- function(a, b) if (!is.null(a) && length(a) && !all(is.na(a))) a else b

# ─── Python backend ──────────────────────────────────────────────────────────

PY_CODE <- r"(
import numpy as np
from scipy.sparse import coo_matrix as _coo
from sklearn.neighbors import NearestNeighbors as _NN
from sklearn.cluster import MiniBatchKMeans as _MBK
import gc as _gc, os as _os

_PSUTIL = False
try:
    import psutil as _ps
    _PSUTIL = True
except ImportError:
    pass

def memory_mb():
    if _PSUTIL:
        return round(_ps.Process(_os.getpid()).memory_info().rss / 1024**2, 1)
    return None

def deep_clean():
    _gc.collect()
    try:
        import ctypes
        ctypes.CDLL("libc.so.6").malloc_trim(0)
    except Exception:
        pass
    return True

def build_pooled_niche_matrix(coords_list, ct_enc_list, n_cts, k1):
    n_cts = int(n_cts); k1 = int(k1)
    sizes = [np.asarray(c, dtype=np.float32).shape[0] for c in coords_list]
    total = int(sum(sizes))
    niche = np.zeros((total, n_cts), dtype=np.float32)
    ptr = 0
    for coords_r, ct_enc_r in zip(coords_list, ct_enc_list):
        coords = np.asarray(coords_r, dtype=np.float32)
        ct_enc = np.asarray(ct_enc_r, dtype=np.int32).ravel()
        n = coords.shape[0]; k = min(k1, n - 1)
        nbrs = _NN(n_neighbors=k + 1, algorithm='ball_tree', n_jobs=-1, metric='euclidean')
        nbrs.fit(coords)
        _, idx = nbrs.kneighbors(coords)
        idx = idx[:, 1:]
        row  = np.repeat(np.arange(n, dtype=np.int32), k)
        col  = ct_enc[idx.ravel()]
        vals = np.ones(len(row), dtype=np.float32)
        block = _coo((vals, (row, col)), shape=(n, n_cts)).toarray()
        block /= k
        niche[ptr:ptr + n] = block; ptr += n
        del coords, ct_enc, idx, row, col, vals, block; _gc.collect()
    return niche

def loo_stability_sweep(niche_mat, ct_enc_r, samp_enc_r, n_samples,
                        k_sweep_r, loo_n, loo_mode, agg, n_cts,
                        sample_group_r=None):
    nm = niche_mat
    ct_enc   = np.asarray(ct_enc_r,  dtype=np.int32).ravel()
    samp_enc = np.asarray(samp_enc_r, dtype=np.int32).ravel()
    k_sweep  = [int(k) for k in k_sweep_r]
    n_s = int(n_samples); n_cts = int(n_cts); loo_n = float(loo_n)
    rng    = np.random.default_rng(0)
    n_iter = min(n_s, 20)
    if loo_mode == 'group':
        # sample_group_r: 0-indexed group-per-sample vector (e.g. timepoint,
        # from --condition-col/--condition-map). Holds out loo_n samples from
        # *every* group each fold, instead of drawing loo_n/pct from the
        # pooled sample list, so a stratified design (N groups x R
        # replicates) actually gets tested as intended.
        sample_group = np.asarray(sample_group_r, dtype=np.int32).ravel()
        groups       = np.unique(sample_group)
        group_idx    = {g: np.where(sample_group == g)[0] for g in groups}
        min_group_sz = min(len(idx) for idx in group_idx.values())
        n_hold_per_group = max(1, min(int(loo_n), min_group_sz - 1))
        if n_hold_per_group < 1 or len(groups) < 2:
            return {'k': k_sweep, 'stability': [0.0] * len(k_sweep)}
        loo_sets = [
            np.concatenate([rng.choice(idx, n_hold_per_group, replace=False)
                             for idx in group_idx.values()])
            for _ in range(n_iter)
        ]
    else:
        n_hold = max(1, min(int(loo_n), n_s - 1)) if loo_mode == 'count' \
                 else max(1, int(np.floor(n_s * loo_n / 100.0)))
        loo_sets = [rng.choice(n_s, n_hold, replace=False) for _ in range(n_iter)]
    agg_fn = np.median if agg == 'median' else np.mean
    def _nh_freq(assign, ct, k2):
        freq = np.zeros((k2, n_cts), dtype=np.float32)
        for j in range(k2):
            m = assign == j
            if not m.any(): continue
            bc = np.bincount(ct[m], minlength=n_cts).astype(np.float32)
            freq[j] = bc / m.sum()
        return freq
    def _nearest(data, centers):
        chunk = 8192; out = np.empty(data.shape[0], dtype=np.int32)
        for s in range(0, data.shape[0], chunk):
            e = min(s + chunk, data.shape[0])
            d = data[s:e, np.newaxis, :] - centers[np.newaxis, :, :]
            out[s:e] = np.argmin(np.einsum('ijk,ijk->ij', d, d), axis=1)
        return out
    scores = []
    for k2 in k_sweep:
        deltas = []
        for held in loo_sets:
            held_set   = set(held.tolist())
            train_mask = np.array([s not in held_set for s in samp_enc.tolist()])
            test_mask  = ~train_mask
            train_data = nm[train_mask]
            if len(train_data) < k2: continue
            km = _MBK(n_clusters=k2, n_init=5, random_state=42,
                      batch_size=min(4096, len(train_data)),
                      max_iter=100, max_no_improvement=10)
            try: tr_assign = km.fit_predict(train_data)
            except Exception: continue
            all_assign = np.empty(len(nm), dtype=np.int32)
            all_assign[train_mask] = tr_assign
            if test_mask.any():
                all_assign[test_mask] = _nearest(nm[test_mask],
                                                  km.cluster_centers_.astype(np.float32))
            ft = _nh_freq(tr_assign,  ct_enc[train_mask], k2)
            ff = _nh_freq(all_assign, ct_enc,             k2)
            deltas.append(float(agg_fn(np.abs(ff - ft))))
        scores.append(float(np.mean(deltas)) if deltas else float('nan'))
    _gc.collect()
    return {'k': k_sweep, 'stability': scores}

def final_kmeans(niche_mat, k2, n_init=10, random_state=42):
    k2 = int(k2); bs = min(4096, max(k2 * 10, 1024))
    km = _MBK(n_clusters=k2, n_init=int(n_init), random_state=int(random_state),
              batch_size=bs, max_iter=300, max_no_improvement=30)
    labels  = (km.fit_predict(niche_mat).astype(np.int32) + 1).tolist()
    centers = km.cluster_centers_.tolist()
    _gc.collect()
    return {'labels': labels, 'centers': centers}
)"

PY_AVAILABLE <- tryCatch({
  py_module_available("sklearn") &&
  py_module_available("scipy")   &&
  py_module_available("numpy")
}, error = function(e) FALSE)

if (PY_AVAILABLE) {
  py_run_string(PY_CODE)
  message("[NeighborhoodR] Python backend active (sklearn ",
          py_eval("__import__('sklearn').__version__"), ")")
} else {
  message("[NeighborhoodR] Python unavailable — using R backend (RANN + kmeans)")
}

# ─── Memory helpers ───────────────────────────────────────────────────────────

mem_mb_r <- function() {
  tryCatch({ m <- gc(verbose=FALSE, reset=FALSE); round(sum(m[, 2L]), 1) },
           error = function(e) NA_real_)
}
mem_log <- function(msg) message(sprintf("[mem] R=%.0f MB | %s", mem_mb_r(), msg))

# ─── Colour palette ──────────────────────────────────────────────────────────

PALETTE_20 <- c(
  "#4e79a7","#f28e2b","#e15759","#76b7b2","#59a14f",
  "#edc948","#b07aa1","#ff9da7","#9c755f","#bab0ac",
  "#d4e157","#26c6da","#ab47bc","#7e57c2","#66bb6a",
  "#ffa726","#ef5350","#42a5f5","#26a69a","#c0ca33"
)
pal <- function(n) {
  if (n <= length(PALETTE_20)) PALETTE_20[seq_len(n)]
  else colorRampPalette(PALETTE_20)(n)
}

# ─── R-fallback analysis functions ───────────────────────────────────────────

.r_build_niche_matrix <- function(spe, ct_col, k1) {
  coords   <- spatialCoords(spe)
  ctypes   <- as.character(colData(spe)[[ct_col]])
  cts_uniq <- sort(unique(ctypes))
  n_cells  <- ncol(spe)
  k1       <- min(k1, n_cells - 1L)
  nn       <- RANN::nn2(coords, k = k1 + 1L)
  nbr_idx  <- nn$nn.idx[, -1L, drop = FALSE]
  mat <- matrix(0.0, nrow = n_cells, ncol = length(cts_uniq),
                dimnames = list(colnames(spe), cts_uniq))
  for (i in seq_len(n_cells)) {
    tab      <- table(factor(ctypes[nbr_idx[i, ]], levels = cts_uniq))
    mat[i, ] <- as.numeric(tab) / k1
  }
  mat
}

.r_assign_to_centers <- function(data, centers) {
  apply(data, 1L, function(row)
    which.min(apply(centers, 1L, function(ctr) sum((row - ctr)^2))))
}

compute_nh_freq <- function(assignments, ctypes, k2, all_cts) {
  freq <- matrix(0.0, nrow = k2, ncol = length(all_cts),
                 dimnames = list(paste0("N", seq_len(k2)), all_cts))
  for (j in seq_len(k2)) {
    idx <- which(assignments == j)
    if (!length(idx)) next
    tab       <- table(factor(ctypes[idx], levels = all_cts))
    freq[j, ] <- as.numeric(tab) / length(idx)
  }
  freq
}

.r_loo_stability_sweep <- function(niche_mat, ctypes_vec, sample_labels,
                                   k_sweep, loo_n, loo_mode, agg_fn,
                                   sample_group = NULL) {
  all_cts        <- colnames(niche_mat)
  unique_samples <- unique(sample_labels)
  n_s            <- length(unique_samples)
  n_iters        <- min(n_s, 20L)

  if (loo_mode == "group") {
    # sample_group: named vector, names = sample names, values = group label
    # (e.g. timepoint, from --condition-col/--condition-map). Holds out
    # loo_n samples from *every* group each fold, instead of drawing
    # loo_n/pct from the pooled sample list.
    groups       <- sample_group[unique_samples]
    group_names  <- unique(groups)
    group_idx    <- lapply(group_names, function(g) unique_samples[groups == g])
    min_group_sz <- min(lengths(group_idx))
    n_hold_per_group <- max(1L, min(as.integer(loo_n), min_group_sz - 1L))
    if (n_hold_per_group < 1L || length(group_names) < 2L) {
      return(data.frame(k = k_sweep, stability = rep(0, length(k_sweep))))
    }
    loo_sets <- lapply(seq_len(n_iters), function(i) {
      set.seed(i)
      unlist(lapply(group_idx, function(idx) sample(idx, n_hold_per_group)))
    })
  } else {
    n_hold <- if (loo_mode == "count") min(as.integer(loo_n), n_s - 1L)
              else max(1L, floor(n_s * loo_n / 100))
    loo_sets <- lapply(seq_len(n_iters), function(i) {
      set.seed(i); sample(unique_samples, n_hold)
    })
  }

  scores <- numeric(length(k_sweep))
  for (ki in seq_along(k_sweep)) {
    k2 <- k_sweep[ki]
    iter_deltas <- numeric(n_iters)
    for (li in seq_along(loo_sets)) {
      held       <- loo_sets[[li]]
      train_mask <- !sample_labels %in% held
      test_mask  <-  sample_labels %in% held
      train_data <- niche_mat[train_mask, , drop = FALSE]
      if (nrow(train_data) < k2) { iter_deltas[li] <- NA_real_; next }
      km <- tryCatch(
        kmeans(train_data, centers = k2, nstart = 5L, iter.max = 100L),
        error = function(e) NULL)
      if (is.null(km)) { iter_deltas[li] <- NA_real_; next }
      te_assign <- if (sum(test_mask))
        .r_assign_to_centers(niche_mat[test_mask, , drop = FALSE], km$centers)
        else integer(0)
      all_assign             <- integer(nrow(niche_mat))
      all_assign[train_mask] <- km$cluster
      all_assign[test_mask]  <- te_assign
      ft <- compute_nh_freq(km$cluster, ctypes_vec[train_mask], k2, all_cts)
      ff <- compute_nh_freq(all_assign, ctypes_vec, k2, all_cts)
      iter_deltas[li] <- agg_fn(abs(ff - ft))
    }
    scores[ki] <- mean(iter_deltas, na.rm = TRUE)
  }
  gc()
  data.frame(k = k_sweep, stability = scores)
}

# ─── Unified niche-matrix builder ─────────────────────────────────────────────

build_niche_data <- function(spe_list, sample_names, ct_cols, k1) {
  all_cts <- sort(unique(unlist(lapply(sample_names, function(s) {
    sort(unique(as.character(colData(spe_list[[s]])[[ct_cols[[s]]]])))
  }))))

  ct_raw_list <- ct_enc_list <- coords_list <- vector("list", length(sample_names))
  n_cells_v   <- integer(length(sample_names))

  for (i in seq_along(sample_names)) {
    s      <- sample_names[i]
    spe    <- spe_list[[s]]
    ctypes <- as.character(colData(spe)[[ct_cols[[s]]]])
    ct_raw_list[[i]] <- ctypes
    ct_enc_list[[i]] <- match(ctypes, all_cts) - 1L
    coords_list[[i]] <- spatialCoords(spe)
    n_cells_v[i]     <- ncol(spe)
  }

  cell_types_v  <- unlist(ct_raw_list)
  ct_encoded_v  <- unlist(ct_enc_list)
  sample_labels <- rep(sample_names, n_cells_v)
  samp_encoded  <- rep(seq_along(sample_names) - 1L, n_cells_v)

  if (PY_AVAILABLE) {
    py_nm <- py$build_pooled_niche_matrix(
      coords_list, ct_enc_list, length(all_cts), as.integer(k1)
    )
    return(list(niche_mat=py_nm, is_python=TRUE,
                cell_types_v=cell_types_v, ct_encoded_v=ct_encoded_v,
                sample_labels=sample_labels, samp_encoded=samp_encoded,
                all_cts=all_cts, n_samples=length(sample_names)))
  }

  mats <- lapply(seq_along(sample_names), function(i) {
    s <- sample_names[i]
    m <- .r_build_niche_matrix(spe_list[[s]], ct_cols[[s]], k1)
    pad_cts <- setdiff(all_cts, colnames(m))
    if (length(pad_cts)) {
      pad <- matrix(0, nrow(m), length(pad_cts), dimnames=list(NULL, pad_cts))
      m   <- cbind(m, pad)
    }
    m[, all_cts, drop=FALSE]
  })
  niche_mat <- do.call(rbind, mats)
  colnames(niche_mat) <- all_cts
  rm(mats); gc()
  list(niche_mat=niche_mat, is_python=FALSE,
       cell_types_v=cell_types_v, ct_encoded_v=ct_encoded_v,
       sample_labels=sample_labels, samp_encoded=samp_encoded,
       all_cts=all_cts, n_samples=length(sample_names))
}

# ─── Celltype-column auto-detection ──────────────────────────────────────────

find_celltype_col <- function(spe) {
  cols <- colnames(colData(spe))
  hits <- grep("annotation|celltype|cell_type|cluster|label",
               cols, ignore.case=TRUE, value=TRUE)
  if (length(hits)) hits[1L] else NULL
}

# ─── Plot helpers ─────────────────────────────────────────────────────────────

.spatial_plot <- function(spe, sname, colour_by="neighborhood",
                          ct_col=NULL, pt_size=1.2, alpha=0.8) {
  df <- as.data.frame(spatialCoords(spe))
  colnames(df) <- c("x", "y")
  df$colour <- if (colour_by == "neighborhood") {
    colData(spe)$neighborhood
  } else {
    if (!is.null(ct_col) && ct_col %in% colnames(colData(spe)))
      as.character(colData(spe)[[ct_col]]) else "unknown"
  }
  lvls    <- sort(unique(df$colour))
  pal_vec <- setNames(pal(length(lvls)), lvls)
  ggplot(df, aes(x, y, colour=colour)) +
    geom_point(size=pt_size, alpha=alpha) +
    scale_colour_manual(values=pal_vec) +
    labs(title=paste0(sname, " — ", colour_by),
         colour=colour_by, x="x", y="y") +
    theme_minimal(base_size=11) +
    theme(legend.position="right",
          plot.title=element_text(face="bold"))
}

.barplot_composition <- function(spe_list, sample_names, ct_cols) {
  rows <- lapply(sample_names, function(s) {
    spe    <- spe_list[[s]]
    ct_col <- ct_cols[[s]]
    if (is.null(ct_col) || !ct_col %in% colnames(colData(spe))) return(NULL)
    data.frame(neighborhood = colData(spe)$neighborhood,
               celltype      = as.character(colData(spe)[[ct_col]]),
               sample        = s)
  })
  df <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(df) || !nrow(df)) return(NULL)

  df_freq <- df |>
    count(neighborhood, celltype) |>
    group_by(neighborhood) |>
    mutate(prop = n / sum(n)) |>
    ungroup()

  pal_vec <- setNames(pal(length(unique(df_freq$celltype))),
                      sort(unique(df_freq$celltype)))
  ggplot(df_freq, aes(neighborhood, prop, fill=celltype)) +
    geom_col(position="stack", width=0.75) +
    scale_fill_manual(values=pal_vec) +
    scale_y_continuous(labels=percent_format()) +
    labs(x="Neighborhood", y="Proportion", fill="Cell type",
         title="Celltype composition per neighborhood") +
    theme_minimal(base_size=11) +
    theme(axis.text.x=element_text(angle=30, hjust=1),
          plot.title=element_text(face="bold"))
}

# ─── Main analysis entry point ────────────────────────────────────────────────

run_neighborhood_analysis <- function(
  spe_paths,
  celltype_cols  = NULL,
  out_dir        = ".",
  label          = format(Sys.Date(), "%Y%m%d"),
  k1             = 10L,
  k2             = NULL,
  k2_min         = 3L,
  k2_max         = NULL,
  loo_mode       = "count",
  loo_n          = 1,
  agg_fn_name    = "median",
  condition_col  = NULL,
  condition_map  = NULL,
  seed           = 42L,
  make_plots     = TRUE,
  verbose        = TRUE
) {
  msg <- function(...) if (verbose) message("[NeighborhoodR] ", ...)

  out_dir <- file.path(out_dir, label)
  dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
  out <- list()

  # 1. Load SPEs
  msg("Loading ", length(spe_paths), " SpatialExperiment object(s)...")
  sample_names <- names(spe_paths)
  spe_list <- lapply(spe_paths, function(p) {
    spe <- readRDS(p)
    stopifnot(is(spe, "SpatialExperiment") || is(spe, "SingleCellExperiment"))
    spe
  })
  names(spe_list) <- sample_names
  gc(); mem_log(sprintf("after load (%d samples)", length(spe_list)))

  # 2. Resolve celltype columns
  if (is.null(celltype_cols)) {
    ct_cols <- setNames(lapply(sample_names, function(s) {
      col <- find_celltype_col(spe_list[[s]])
      if (is.null(col))
        stop("Cannot auto-detect celltype column for '", s,
             "'. Pass --celltype-col.")
      msg("  ", s, ": auto-detected celltype column '", col, "'")
      col
    }), sample_names)
  } else if (length(celltype_cols) == 1 && is.null(names(celltype_cols))) {
    ct_cols <- setNames(rep(celltype_cols, length(sample_names)), sample_names)
  } else {
    ct_cols <- as.list(celltype_cols[sample_names])
  }

  # 3. Build niche matrix
  msg("Building niche matrix (K1=", k1, ")...")
  nd <- build_niche_data(spe_list, sample_names, ct_cols, k1)
  mem_log("after niche matrix build")

  n_cts    <- length(nd$all_cts)
  k2_max_v <- k2_max %||% (n_cts - 1L)
  k2_max_v <- max(k2_min, min(k2_max_v, n_cts - 1L))
  k_sweep  <- seq(k2_min, k2_max_v)

  # 4. LOO stability sweep or fixed K2
  agg_fn   <- if (agg_fn_name == "median") median else mean
  sweep_df <- NULL

  if (!is.null(k2)) {
    msg("Skipping LOO sweep — using fixed K2=", k2)
    optimal_k2 <- as.integer(k2)
  } else {
    msg("Running LOO stability sweep (K2=", k2_min, " to ", k2_max_v, ")...")

    # loo_mode="group": resolve one label per sample from condition_col /
    # condition_map (same source used for the final condition annotation
    # below), so held-out samples come from every group each fold instead
    # of the pooled sample list.
    sample_group_v <- NULL
    if (loo_mode == "group") {
      groups <- vapply(sample_names, function(s) {
        if (!is.null(condition_col) &&
            condition_col %in% colnames(colData(spe_list[[s]]))) {
          v <- as.character(colData(spe_list[[s]])[[condition_col]])
          v <- v[!is.na(v) & nchar(v)]
          if (length(v)) v[1] else NA_character_
        } else if (!is.null(condition_map) && !is.null(condition_map[[s]])) {
          condition_map[[s]]
        } else {
          NA_character_
        }
      }, character(1))
      names(groups) <- sample_names
      if (anyNA(groups))
        stop("--loo-mode=group requires every sample to have a condition label ",
             "— set --condition-col or provide every sample in --condition-map.")
      if (length(unique(groups)) < 2L)
        stop("--loo-mode=group requires at least 2 distinct condition values across samples.")
      sample_group_v <- groups
    }

    if (nd$is_python) {
      sample_group_enc <- if (!is.null(sample_group_v))
        as.integer(match(sample_group_v[sample_names], unique(sample_group_v)) - 1L)
      else NULL
      res <- py$loo_stability_sweep(
        nd$niche_mat, nd$ct_encoded_v, nd$samp_encoded,
        nd$n_samples, k_sweep, loo_n, loo_mode, agg_fn_name, n_cts,
        sample_group_r = sample_group_enc
      )
      sweep_df <- data.frame(k=as.integer(unlist(res$k)),
                             stability=as.numeric(unlist(res$stability)))
    } else {
      sweep_df <- .r_loo_stability_sweep(
        nd$niche_mat, nd$cell_types_v, nd$sample_labels,
        k_sweep, loo_n, loo_mode, agg_fn,
        sample_group = sample_group_v
      )
    }
    optimal_k2 <- sweep_df$k[which.min(sweep_df$stability)]
    msg("Optimal K2 = ", optimal_k2,
        "  (stability = ", round(min(sweep_df$stability, na.rm=TRUE), 4), ")")
  }

  # 5. Final assignment
  msg("Running final MiniBatchKMeans (K2=", optimal_k2, ", seed=", seed, ")...")
  set.seed(seed)
  if (nd$is_python) {
    km_res      <- py$final_kmeans(nd$niche_mat, optimal_k2,
                                    n_init=10L, random_state=as.integer(seed))
    assignments <- as.integer(unlist(km_res$labels))
    py$deep_clean()
  } else {
    km          <- kmeans(nd$niche_mat, centers=optimal_k2,
                          nstart=25L, iter.max=300L)
    assignments <- km$cluster
    rm(km)
  }
  rm(nd); gc(); mem_log("after final kmeans")

  # 6. Write neighborhood + sample into colData, cbind
  msg("Writing neighborhood assignments...")
  ptr <- 1L
  for (s in sample_names) {
    n_c <- ncol(spe_list[[s]])
    colData(spe_list[[s]])$neighborhood <- paste0("N", assignments[ptr:(ptr+n_c-1L)])
    colData(spe_list[[s]])$sample        <- rep(s, n_c)
    if (!is.null(condition_col) &&
        condition_col %in% colnames(colData(spe_list[[s]]))) {
      # condition already in colData — no action needed
    } else if (!is.null(condition_map) && !is.null(condition_map[[s]])) {
      colData(spe_list[[s]])$condition_map <- rep(condition_map[[s]], n_c)
    }
    ptr <- ptr + n_c
  }
  joint_spe <- do.call(cbind, unname(spe_list))
  rm(spe_list, assignments); gc(); mem_log("after cbind")

  # 7. Save outputs
  rds_path <- file.path(out_dir, paste0(label, "_joint_spe.rds"))
  msg("Saving joint SPE → ", rds_path)
  saveRDS(joint_spe, rds_path)
  out$joint_spe <- rds_path

  cd      <- as.data.frame(colData(joint_spe)[, c("sample","neighborhood"), drop=FALSE])
  sum_df  <- do.call(rbind, lapply(sample_names, function(s) {
    tab <- table(cd$neighborhood[cd$sample == s])
    data.frame(sample=s, neighborhood=names(tab), n_cells=as.integer(tab))
  }))
  sum_path <- file.path(out_dir, paste0(label, "_assignment_summary.csv"))
  write.csv(sum_df, sum_path, row.names=FALSE)
  out$assignment_summary <- sum_path

  if (!is.null(sweep_df)) {
    sweep_path <- file.path(out_dir, paste0(label, "_sweep_results.csv"))
    write.csv(sweep_df, sweep_path, row.names=FALSE)
    out$sweep_results <- sweep_path
  }

  prov_path <- file.path(out_dir, paste0(label, "_provenance.json"))
  writeLines(toJSON(list(
    tool           = "NeighborhoodR-CLI",
    label          = label,
    run_date       = as.character(Sys.time()),
    r_version      = paste(R.version$major, R.version$minor, sep="."),
    python_backend = PY_AVAILABLE,
    samples        = sample_names,
    celltype_cols  = ct_cols,
    k1             = k1,
    optimal_k2     = optimal_k2,
    k2_fixed       = !is.null(k2),
    k2_min         = k2_min,
    k2_max         = k2_max_v,
    loo_mode       = loo_mode,
    loo_n          = loo_n,
    agg_fn         = agg_fn_name,
    condition_col  = condition_col %||% NA,
    seed           = seed,
    git_sha        = Sys.getenv("PHENOSUITE_GIT_SHA", "unknown"),
    image_digest   = Sys.getenv("PHENOSUITE_IMAGE_DIGEST", "unknown")
  ), pretty=TRUE, auto_unbox=TRUE), prov_path)
  out$provenance <- prov_path

  # 8. Plots
  if (make_plots) {
    msg("Generating plots...")
    snames_j <- unique(as.character(colData(joint_spe)$sample))
    for (s in snames_j) {
      spe_s    <- joint_spe[, colData(joint_spe)$sample == s]
      p <- .spatial_plot(spe_s, s, colour_by="neighborhood",
                         ct_col=ct_cols[[s]], pt_size=1.2, alpha=0.8)
      png_path <- file.path(out_dir,
                            paste0(label, "_spatial_", make.names(s), ".png"))
      ggsave(png_path, p, width=8, height=7, dpi=150)
      out[[paste0("spatial_", make.names(s))]] <- png_path
    }

    spe_list_j <- setNames(
      lapply(snames_j, function(s) joint_spe[, colData(joint_spe)$sample == s]),
      snames_j)
    bp <- .barplot_composition(spe_list_j, snames_j, ct_cols)
    if (!is.null(bp)) {
      bp_path <- file.path(out_dir, paste0(label, "_composition_barplot.png"))
      ggsave(bp_path, bp, width=10, height=6, dpi=150)
      out$composition_barplot <- bp_path
    }
  }

  msg("Done. Outputs written to: ", out_dir)
  invisible(out)
}
