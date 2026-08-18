# Writes Akoya Vectra cell_seg_data-style CSVs (Sample Name, Cell X Position,
# Cell Y Position, Phenotype, Tissue Category) so automated-phenotyping
# outputs drop directly into the existing PCF toolkit
# (utils/spatial-shiny/vectra_lib_v4.py's extract_data(), which selects these
# columns by name). Mirrors the convention already used by
# modify_clusters/dev/server.R's prepare_pcf_inputs(), except Sample Name is
# derived from a real per-cell sample column instead of row indices.

# Returns the first candidate colData column name that's actually present, or NULL.
detect_sample_column <- function(cd, candidates = c("sample_id", "sample", "Sample", "Sample_Name", "image_id")) {
  hit <- candidates[candidates %in% names(cd)]
  if (length(hit) == 0) return(NULL)
  hit[1]
}

# Per-cell Tissue Category: prefers real tissue/region colData signals, then
# a caller-supplied placeholder (e.g. a free-text "Tissue type" input), and
# only then falls to `sample_id` — SpatialExperiment auto-populates a
# `sample_id` column ("sample01") on every object even when the user never
# sets one, so it must rank below the user's actual input or that input
# would never be reachable in practice.
resolve_tissue_category <- function(spe, fallback = "Unknown") {
  cd <- colData(spe)
  for (col in c("Classifier_Label", "Analysis_Region")) {
    if (col %in% names(cd) && length(unique(cd[[col]])) > 0) {
      return(as.character(cd[[col]]))
    }
  }
  if (!is.null(fallback) && nzchar(trimws(fallback))) {
    return(rep(fallback, ncol(spe)))
  }
  if ("sample_id" %in% names(cd)) return(as.character(cd[["sample_id"]]))
  rep("Unknown", ncol(spe))
}

# Writes one Vectra CSV per detected sample if >1 unique sample is found
# (extract_data() assumes one CSV = one sample, reading Sample Name only from
# the first row), otherwise a single file.
write_vectra_csv <- function(spe, phenotype_col, tissue_fallback = "Unknown",
                              sample_col = NULL, out_dir, file_prefix = "vectra") {
  cd     <- colData(spe)
  coords <- as.data.frame(spatialCoords(spe))

  if (is.null(sample_col)) sample_col <- detect_sample_column(cd)
  sample_vec <- if (!is.null(sample_col)) as.character(cd[[sample_col]]) else rep(file_prefix, ncol(spe))

  out <- data.frame(
    `Sample Name`      = sample_vec,
    `Cell X Position`  = coords[[1]],
    `Cell Y Position`  = coords[[2]],
    Phenotype          = as.character(cd[[phenotype_col]]),
    `Tissue Category`  = resolve_tissue_category(spe, tissue_fallback),
    check.names        = FALSE,
    stringsAsFactors   = FALSE
  )

  safe_name <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

  samples <- unique(sample_vec)
  paths <- character(0)
  if (length(samples) > 1) {
    for (s in samples) {
      path <- file.path(out_dir, glue::glue("{file_prefix}_{safe_name(s)}_{phenotype_col}.csv"))
      write.csv(out[sample_vec == s, ], path, row.names = FALSE)
      paths <- c(paths, path)
    }
  } else {
    path <- file.path(out_dir, glue::glue("{file_prefix}_{phenotype_col}.csv"))
    write.csv(out, path, row.names = FALSE)
    paths <- path
  }
  invisible(paths)
}
