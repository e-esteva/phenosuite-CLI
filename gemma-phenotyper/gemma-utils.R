library(SpatialExperiment)
library(SingleCellExperiment)
library(jsonlite)
library(curl)

DEFAULT_GEMMA_ONTOLOGY <- paste(
  "CD3+CD4+ -> CD4 T cell",
  "CD3+CD8+ -> CD8 T cell",
  "CD3+FOXP3+ -> Treg",
  "CD20+  -> B cell",
  "CD68+  -> Macrophage",
  "EPCAM+ -> Tumor cell",
  "CD31+  -> Endothelial cell",
  "Other  -> Other",
  sep = " | "
)

# Builds a starter ontology scaffold from a marker panel — one placeholder
# row per (already-cleaned) marker name, so every reference in the default
# text is guaranteed to exist in the uploaded object. Used only as the
# fallback when no ontology_manifest.json ships with the checkpoint (see
# read_ontology_manifest()).
generate_default_ontology <- function(marker_names) {
  marker_names <- unique(marker_names)
  rows <- paste0(marker_names, "+ -> ?")
  paste(c(rows, "Other -> Other"), collapse = " | ")
}

# Fine-tuning necessarily uses some marker -> cell-type mapping to build its
# training prompts; that mapping should ship alongside the checkpoint as
# ontology_manifest.json (sibling to config.json) so inference doesn't have
# to guess it:
#   { "markers": ["CD3", "CD4", ...], "ontology": "CD3+CD4+ -> CD4 T cell | ..." }
# where "markers" is the (cleaned) marker vocabulary the model was trained to
# recognize and "ontology" is the exact mapping text used at training time,
# in the same "A+B+ -> Label | ..." format this app's ontology box expects.
# Returns NULL if the checkpoint doesn't provide one (older/external models),
# so callers can fall back to generate_default_ontology().
read_ontology_manifest <- function(model_dir) {
  path <- file.path(model_dir, "ontology_manifest.json")
  if (!file.exists(path)) return(NULL)
  manifest <- tryCatch(jsonlite::fromJSON(path), error = function(e) NULL)
  if (is.null(manifest) || is.null(manifest$markers) || is.null(manifest$ontology))
    return(NULL)
  list(markers = as.character(manifest$markers), ontology = as.character(manifest$ontology))
}

# Assembles the ontology text shown to the user: markers the checkpoint was
# fine-tuned on use its trained mapping verbatim (untouched); markers present
# in the uploaded object but outside that training vocabulary get a "?"
# scaffold row appended, so the user only ever has to write ontology rules
# for genuinely novel markers. With no manifest (trained = NULL), every
# marker is treated as novel.
build_ontology_text <- function(marker_names, trained = NULL) {
  marker_names <- unique(marker_names)
  if (is.null(trained)) return(generate_default_ontology(marker_names))

  uncovered <- setdiff(marker_names, trained$markers)
  if (length(uncovered) == 0) return(trained$ontology)

  paste(c(trained$ontology, paste0(uncovered, "+ -> ?")), collapse = " | ")
}

# Strips assay-specific suffixes (e.g. "_Cytoplasm_Intensity", "_Nucleus_Intensity")
# down to the base marker name, matching the convention used in automated_phenotyping.
clean_marker_name <- function(x) {
  vapply(strsplit(x, "(?i)_Cytoplasm|_Nucleus_", perl = TRUE), `[`, character(1), 1)
}

# Markers that report cell STATE rather than lineage. Ki67 leading a cluster's
# elevated list was being read as lineage evidence; naming its role instead
# lets the model use it as a qualifier ("Proliferating CD4 T cell") and take
# the lineage from the remaining markers.
#
# Deliberately conservative — only markers whose role is unambiguous. Several
# common markers are context-dependent and are intentionally NOT listed,
# because calling them "activation" would suppress real lineage calls:
#   CD25    activation, but also Treg-defining with FOXP3
#   HLA-DR  activation on T cells, lineage on B cells / APCs
#   CD38    activation, but also plasma-cell lineage
#   CD103   residency, but also a cDC1 subset marker
# Override via the `state_markers` argument if your panel needs it.
DEFAULT_STATE_MARKERS <- list(
  proliferation = c("KI67", "MKI67", "PCNA", "MCM2", "TOP2A", "PHH3",
                    "PHOSPHOHISTONEH3", "BRDU", "CCNB1", "AURKB"),
  `activation/checkpoint` = c("PD1", "PDCD1", "TIM3", "HAVCR2", "LAG3", "TIGIT",
                              "CTLA4", "CD69", "ICOS", "OX40", "TNFRSF4",
                              "CD137", "TNFRSF9", "41BB"),
  `cytotoxic effector` = c("GZMB", "GRANZYMEB", "PRF1", "PERFORIN", "GNLY",
                           "GRANULYSIN"),
  `apoptosis/survival` = c("BCL2", "CLEAVEDCASPASE3", "CASP3", "CC3")
)

# Marker names vary in punctuation and case across panels (Ki67 / Ki-67 /
# MKI67, PD_1 / PD-1 / PD1), so compare on an alphanumeric-only uppercase form.
.norm_marker <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))

# Given the marker names being reported as elevated, returns a sentence naming
# any that are state rather than lineage markers — or "" when there are none.
.state_marker_clause <- function(elevated_names, state_markers = DEFAULT_STATE_MARKERS) {
  if (length(elevated_names) == 0) return("")
  norm <- .norm_marker(elevated_names)
  hits <- character(0)
  for (role in names(state_markers)) {
    m <- elevated_names[norm %in% .norm_marker(state_markers[[role]])]
    if (length(m) > 0)
      hits <- c(hits, paste0(paste(m, collapse = " and "), " (", role, ")"))
  }
  if (length(hits) == 0) return("")
  paste0(" Note — the following report cell STATE, not lineage: ",
         paste(hits, collapse = ", "),
         ". Determine the lineage from the other elevated markers and use these",
         " only as a qualifier (e.g. \"Proliferating CD4 T cell\").")
}

# Describes the scale a matrix is actually on, measured rather than assumed.
# The assay name alone doesn't tell you: "exprs" is arcsinh in some pipelines
# and per-marker z-scores in others, and asserting the wrong one in the prompt
# is worse than saying nothing. Columns are sampled so this stays cheap on
# objects with millions of cells.
.describe_scale <- function(mat) {
  n   <- ncol(mat)
  idx <- if (n > 5000) sample.int(n, 5000) else seq_len(n)
  sub <- suppressWarnings(as.matrix(mat[, idx, drop = FALSE]))
  if (!is.numeric(sub) || length(sub) == 0) return("units unknown")

  mu  <- rowMeans(sub, na.rm = TRUE)
  sdv <- apply(sub, 1, stats::sd, na.rm = TRUE)
  rng <- suppressWarnings(range(sub, na.rm = TRUE))
  if (!all(is.finite(rng))) return("units unknown")

  z_like <- all(abs(mu) < 0.15, na.rm = TRUE) &&
            all(abs(sdv - 1) < 0.25, na.rm = TRUE)
  if (z_like)
    return("per-marker z-scores across all cells (0 = that marker's average, +2 = two SD above)")

  if (rng[1] < -0.05)
    return(sprintf("per-marker centred values, negative = below that marker's average (range %.1f to %.1f)",
                   rng[1], rng[2]))
  sprintf("transformed intensities, higher = more expression (range %.1f to %.1f)",
          rng[1], rng[2])
}

# Priority: exprs > scaled > z-scored counts. The chosen assay is used as-is —
# no further transformation — so whatever scaling the object already carries is
# what reaches the prompt.
.resolve_intensity_matrix <- function(spe, marker_cols) {
  avail <- assayNames(spe)
  if ("exprs" %in% avail) {
    mat <- assay(spe, "exprs"); nm <- "exprs"
  } else if ("scaled" %in% avail) {
    mat <- assay(spe, "scaled"); nm <- "scaled"
  } else if ("counts" %in% avail) {
    message("No 'exprs' or 'scaled'; computing feature-wise z-score of counts.")
    raw  <- assay(spe, "counts")
    sds  <- apply(raw, 1, sd); sds[sds == 0] <- 1
    mat  <- (raw - rowMeans(raw)) / sds; nm <- "z-scored counts"
  } else {
    stop("No usable assay. Available: ", paste(assayNames(spe), collapse = ", "))
  }
  missing <- setdiff(marker_cols, rownames(mat))
  if (length(missing) > 0)
    stop("marker_cols not in assay rownames: ", paste(missing, collapse = ", "))

  mat   <- mat[marker_cols, , drop = FALSE]
  scale <- .describe_scale(mat)
  list(mat = mat, assay = nm, scale = scale,
       label = paste0(nm, "; ", scale))
}

# Serialises SPE to JSONL prompts for Gemma inference.
# cluster mode (cluster_col set): one prompt per cluster, fast.
# cell mode (cluster_col = NULL): one prompt per cell, slow.
prepare_inference_input <- function(spe,
                                    marker_cols,
                                    ontology_text,
                                    cluster_col = NULL,
                                    morph_cols  = NULL,
                                    out_path    = tempfile(fileext = ".jsonl")) {
  res     <- .resolve_intensity_matrix(spe, marker_cols)
  int_mat <- res$mat
  lbl     <- res$label
  cd      <- colData(spe)

  cluster_mode <- !is.null(cluster_col) && cluster_col %in% names(cd)
  if (!is.null(cluster_col) && !cluster_mode)
    warning("cluster_col '", cluster_col, "' not in colData; falling back to cell mode.")

  con          <- file(out_path, "w"); on.exit(close(con))
  marker_names <- clean_marker_name(marker_cols)

  if (cluster_mode) {
    clusters <- unique(cd[[cluster_col]])
    for (cl in clusters) {
      idx      <- which(cd[[cluster_col]] == cl)
      cl_means <- rowMeans(int_mat[, idx, drop = FALSE])
      mkr      <- paste(mapply(function(m, v) sprintf("%s=%.3f", m, v),
                               marker_names, cl_means), collapse = ", ")
      prompt <- paste0(
        "Phenotype this cell cluster.\n",
        "Cluster: ", cl, " (n=", length(idx), " cells)\n",
        "Mean markers (", lbl, "): ", mkr, "\n",
        "Ontology: ", ontology_text
      )
      writeLines(toJSON(list(
        cluster_id = as.character(cl),
        n_cells    = length(idx),
        messages   = list(
          list(role = "system",
               content = "You are a cellular phenotyping assistant for mIF data. Return JSON only."),
          list(role = "user", content = prompt)
        )
      ), auto_unbox = TRUE), con)
    }
    message("Cluster mode: wrote ", length(clusters), " prompts to ", out_path)

  } else {
    for (i in seq_len(ncol(spe))) {
      mkr <- paste(mapply(function(m, v) sprintf("%s=%.3f", m, v),
                          marker_names, int_mat[, i]), collapse = ", ")
      mph <- if (!is.null(morph_cols))
        paste(sapply(morph_cols, function(m)
          sprintf("%s=%.2f", m, cd[[m]][[i]])), collapse = ", ") else ""
      prompt <- paste0(
        "Phenotype this cell.\n",
        if ("sample_id" %in% names(cd)) paste0("Sample: ", cd[["sample_id"]][[i]], "\n"),
        "Markers (", lbl, "): ", mkr, "\n",
        if (nchar(mph) > 0) paste0("Morphology: ", mph, "\n") else "",
        "Ontology: ", ontology_text
      )
      writeLines(toJSON(list(
        cell_index = i,
        messages   = list(
          list(role = "system",
               content = "You are a cellular phenotyping assistant for mIF data. Return JSON only."),
          list(role = "user", content = prompt)
        )
      ), auto_unbox = TRUE), con)
    }
    message("Cell mode: wrote ", ncol(spe), " prompts to ", out_path)
  }
  out_path
}

# Per-cluster signal QC.
#
# `exprs` is z-scored per MARKER across cells, which equalises channels but not
# cells: a cluster of large or densely-packed cells can sit above average on
# many channels at once (segmentation spillover, autofluorescence, sheer
# protein content), while a dim cluster sits below on all of them. That shows
# up as a cluster reading "positive" for a biologically impossible share of the
# panel, or for none of it.
#
# This is reported, NOT corrected. Dividing it out would be per-cell
# normalisation, which is defensible in scRNA-seq (library size is an
# independent per-cell technical artifact) but not here: one staining round
# covers every cell, so the markers are the independent axis and a cell's total
# signal carries real biology. Genuine spillover is a segmentation/compensation
# problem, to be fixed upstream rather than scaled away here.
#
# Returns one row per cluster:
#   mean_z_all_markers   mean value across the whole panel (the offset)
#   n_markers_over_1sd   average count of markers above +1 per cell
#   frac_panel_elevated  that count as a fraction of the panel
cluster_signal_qc <- function(spe, marker_cols, cluster_col) {
  res <- .resolve_intensity_matrix(spe, marker_cols)
  m   <- res$mat
  n_panel <- nrow(m)
  cl  <- as.character(colData(spe)[[cluster_col]])
  over <- colSums(m > 1.0)
  data.frame(
    cluster_id          = names(tapply(colMeans(m), cl, mean)),
    mean_z_all_markers  = round(as.numeric(tapply(colMeans(m), cl, mean)), 3),
    n_markers_over_1sd  = round(as.numeric(tapply(over, cl, mean)), 2),
    frac_panel_elevated = round(as.numeric(tapply(over, cl, mean)) / n_panel, 3),
    stringsAsFactors    = FALSE
  )
}

# Reads a marker glossary file: "MARKER = gloss" per line, '#' comments ignored.
# Returns a named character vector keyed on the normalised marker name.
#
# A zero-shot model only knows markers it has read about. Panel-specific or
# ambiguously-named ones it guesses at, confidently: on a 39-plex lymph node
# panel gemma3:12b-it-qat answered "Regulatory T cells" for TCF4 and "Plasma
# cells" for its alias E2-2, both at confidence 0.95, when TCF4/E2-2 is the
# defining pDC transcription factor. A pure 3,554-cell pDC cluster was
# consequently called "Macrophage" in every prompt variant tested. Supplying the
# mapping is the only thing that fixes a wrong prior — no rewording does.
read_marker_glossary <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return(NULL)
  ln <- readLines(path, warn = FALSE)
  ln <- ln[!grepl("^\\s*#", ln) & grepl("=", ln)]
  if (length(ln) == 0) return(NULL)
  key <- .norm_marker(trimws(sub("=.*$", "", ln)))
  val <- trimws(sub("^[^=]*=", "", ln))
  keep <- nzchar(key) & nzchar(val)
  if (!any(keep)) return(NULL)
  setNames(val[keep], key[keep])
}

# Glossary lines for the markers actually being reported as elevated. Scoped
# this way so a long glossary costs no tokens on clusters that don't use it.
.glossary_clause <- function(elevated_names, glossary) {
  if (is.null(glossary) || length(elevated_names) == 0) return("")
  hit <- .norm_marker(elevated_names) %in% names(glossary)
  if (!any(hit)) return("")
  entries <- paste0(elevated_names[hit], " = ",
                    unname(glossary[.norm_marker(elevated_names[hit])]))
  paste0(" Marker notes for this panel: ", paste(entries, collapse = "; "), ".")
}

# Zero-shot prompts — based on automated_phenotyping's "Symmetric single
# choice" GPT prompt: instead of an explicit ontology table, each prompt names
# only the most distinctive markers for that cluster/cell (top/bottom 10% by
# value, widened from that prompt's original 5% so lineage-defining markers
# are less likely to be cut from large panels), plus the tissue type if given.
#
# Markers are split into ELEVATED and REDUCED rather than listed as one signed
# series. A flat "CD68:-0.102, Mac2Gal3:-0.038, ..." list was being read as
# evidence FOR the named lineages regardless of sign — a marker-negative
# cluster came back a confident "Macrophage" purely because CD68 appeared in
# the text. Splitting the two makes the sign impossible to miss and turns the
# low markers into usable exclusion evidence.
#
# `min_z` is an absolute floor: a marker only counts as elevated if it clears
# it. Top-10% alone is purely relative, so a cluster with nothing actually
# expressed still surfaced its least-negative markers formatted identically to
# a genuine CD68-high cluster. When nothing clears the floor the prompt says
# so and asks for "Unclassified" instead of inviting a guess.
prepare_ollama_prompts <- function(spe,
                                   marker_cols,
                                   cluster_col = NULL,
                                   tissue      = NULL,
                                   min_z       = 0.5,
                                   state_markers = DEFAULT_STATE_MARKERS,
                                   spillover_frac = 0.25,
                                   marker_glossary = NULL,
                                   center_profile = FALSE,
                                   out_path    = tempfile(fileext = ".jsonl")) {
  res     <- .resolve_intensity_matrix(spe, marker_cols)
  int_mat <- res$mat
  cd      <- colData(spe)

  cluster_mode <- !is.null(cluster_col) && cluster_col %in% names(cd)
  if (!is.null(cluster_col) && !cluster_mode)
    warning("cluster_col '", cluster_col, "' not in colData; falling back to cell mode.")

  tissue_clause <- if (!is.null(tissue) && nzchar(trimws(tissue)))
    paste0(" This is a ", trimws(tissue), ".") else ""

  # State the scale the numbers are on. Without it the model is reading bare
  # values with no idea whether 3.4 is high — and the units differ by object
  # (arcsinh exprs vs scaled vs z-scored counts).
  scale_clause <- paste0(" Values are ", res$scale, ".")

  json_instr <- paste0(
    'Respond with ONLY JSON in the form {"cell_type": "<answer, at most 5 words>", ',
    '"confidence": <0-1>} — no markdown, no explanation.'
  )

  fmt <- function(nm, v) paste0(nm, " (", sprintf("%+.2f", v), ")", collapse = ", ")

  n_panel <- length(marker_cols)

  # Warns when a profile is elevated across an implausible share of the panel.
  # No cell type expresses a quarter of a lineage panel, so this is a symptom
  # of spillover / unusually bright cells, not co-expression. Reported to the
  # model rather than corrected — see cluster_signal_qc().
  spillover_clause <- function(vals) {
    n_up <- sum(vals > 1.0, na.rm = TRUE)
    if (n_panel == 0 || n_up / n_panel < spillover_frac) return("")
    paste0(" Caution: signal is elevated for ", n_up, " of ", n_panel,
           " markers here, which usually reflects segmentation spillover or",
           " unusually bright cells rather than true co-expression. Weigh the",
           " elevated markers accordingly.")
  }

  # Opt-in, default off. Subtracting a profile's own across-marker mean removes
  # the brightness offset, but it is per-cell normalisation by another name:
  # valid in scRNA-seq where library size is an independent per-cell artifact,
  # not here where one staining round covers every cell and the markers are the
  # independent axis. It also re-inflates dim profiles (a marker at -0.04 in a
  # profile averaging -1.08 becomes +1.04 and reads as "elevated"), so the
  # absolute min_z floor is still applied to the RAW values afterwards.
  maybe_center <- function(vals) if (isTRUE(center_profile)) vals - mean(vals, na.rm = TRUE) else vals

  # Returns the marker section of the prompt: elevated markers (top 10% AND
  # above min_z) and reduced markers (bottom 10%) as separately labelled lists.
  marker_clause <- function(names_, vals) {
    sel  <- maybe_center(vals)          # ranking basis
    hi_q <- sel > quantile(sel, .90, na.rm = TRUE)
    lo_q <- sel < quantile(sel, .10, na.rm = TRUE)
    up   <- which(hi_q & vals >= min_z) # floor always on the RAW value
    dn   <- which(lo_q)

    if (length(up) == 0) {
      # Nothing is actually expressed. Say that plainly and offer the honest
      # answer, rather than handing over a list of near-zero values that reads
      # like positive evidence.
      dn_txt <- if (length(dn) > 0)
        paste0(" The lowest markers are ", fmt(names_[dn][order(vals[dn])], sort(vals[dn])), ".")
        else ""
      return(list(
        body = paste0("No markers are elevated in this cluster (none reach ",
                      sprintf("%+.2f", min_z), ").", dn_txt),
        empty = TRUE))
    }

    up <- up[order(vals[up], decreasing = TRUE)]
    dn <- dn[order(vals[dn])]
    body <- paste0("Elevated markers: ", fmt(names_[up], vals[up]), ".")
    if (length(dn) > 0)
      body <- paste0(body, " Reduced markers: ", fmt(names_[dn], vals[dn]), ".")
    body <- paste0(body, .glossary_clause(names_[up], marker_glossary),
                   .state_marker_clause(names_[up], state_markers),
                   spillover_clause(vals))
    list(body = body, empty = FALSE)
  }

  # Assembles one full user prompt from a marker vector.
  build_prompt <- function(names_, vals, preamble) {
    mc <- marker_clause(names_, vals)
    # Softened escape hatch. The hard form ("answer Unclassified if the markers
    # do not support a specific cell type") pushed ~1.08 of recall mass into
    # Unclassified — most of it clusters that DID have usable markers, notably
    # T cells (0.62 abstained) and Proliferating B cells (0.66). Biasing toward
    # a committed best-guess recovered those (T cells 0.02 -> 0.62,
    # Proliferating B cells 0.24 -> 0.86) and raised mean accuracy 0.64 -> 0.67
    # with zero abstention, which is what a manual-review workflow wants.
    # The option is retained for genuinely contradictory profiles.
    ask <- if (mc$empty)
      paste0(" Give your best call from what evidence there is;",
             " answer \"Unclassified\" only if nothing supports any lineage.")
      else
      paste0(" What cell type is this?",
             " Base the call on the elevated markers; treat reduced markers as",
             " evidence against those lineages. Give your best call.",
             " Answer \"Unclassified\" ONLY if the elevated markers are genuinely",
             " contradictory or too weak to support any lineage.")
    paste0(preamble, mc$body, tissue_clause, scale_clause, ask, " ", json_instr)
  }

  con          <- file(out_path, "w"); on.exit(close(con))
  marker_names <- clean_marker_name(marker_cols)

  if (cluster_mode) {
    clusters <- unique(cd[[cluster_col]])
    for (cl in clusters) {
      idx      <- which(cd[[cluster_col]] == cl)
      cl_means <- rowMeans(int_mat[, idx, drop = FALSE])
      prompt   <- build_prompt(marker_names, cl_means, "")
      writeLines(toJSON(list(
        cluster_id = as.character(cl),
        n_cells    = length(idx),
        messages   = list(
          list(role = "system", content = "You are an expert immunologist."),
          list(role = "user", content = prompt)
        )
      ), auto_unbox = TRUE), con)
    }
    message("Cluster mode (zero-shot prompts): wrote ", length(clusters), " prompts to ", out_path)

  } else {
    for (i in seq_len(ncol(spe))) {
      prompt <- build_prompt(marker_names, int_mat[, i], "")
      writeLines(toJSON(list(
        cell_index = i,
        messages   = list(
          list(role = "system", content = "You are an expert immunologist."),
          list(role = "user", content = prompt)
        )
      ), auto_unbox = TRUE), con)
    }
    message("Cell mode (zero-shot prompts): wrote ", ncol(spe), " prompts to ", out_path)
  }
  out_path
}

# Broadcasts cluster-level Gemma predictions to per-cell colData columns.
broadcast_cluster_labels <- function(spe, preds_df, cluster_col) {
  cluster_ids <- as.character(colData(spe)[[cluster_col]])
  lookup <- setNames(as.character(preds_df$cell_type), as.character(preds_df$cluster_id))
  conf   <- setNames(as.numeric(preds_df$confidence),  as.character(preds_df$cluster_id))
  lab    <- unname(lookup[cluster_ids])

  spe[["gemma_cell_type"]]  <- lab
  spe[["gemma_confidence"]] <- unname(conf[cluster_ids])

  # Cluster-indexed label, e.g. "C01-Treg". Without this the cluster -> cell-type
  # mapping is only recoverable by joining the predictions CSV back onto the
  # object, which downstream tools (manual annotation, plotting) won't do. Same
  # "{cluster}-{label}" convention automated_phenotyping uses for
  # annotated_clusters_marked, so both apps read the same way.
  spe[["gemma_cell_type_marked"]] <- ifelse(is.na(lab), NA_character_,
                                            paste0(cluster_ids, "-", lab))
  # Keep the cluster column that produced the call, so the provenance survives
  # even if the user later re-clusters or renames columns.
  spe[["gemma_cluster_col"]] <- cluster_col
  spe[["gemma_cluster_id"]]  <- cluster_ids
  spe
}

# Returns "merged", "adapter", "unknown", or "not_found"
detect_model_type <- function(dir_path) {
  if (!dir.exists(dir_path)) return("not_found")
  files <- list.files(dir_path, recursive = FALSE)
  if ("adapter_config.json" %in% files) return("adapter")
  if ("config.json" %in% files && any(grepl("safetensors|pytorch_model", files)))
    return("merged")
  return("unknown")
}

# Locates the model weights directory within an extracted zip, tolerating
# extra nesting or junk entries (e.g. __MACOSX/) that a flat single-subdir
# check would miss. Falls back to `root` if no marker file is found.
find_model_dir <- function(root) {
  hits <- list.files(root, pattern = "^(config|adapter_config)\\.json$",
                      recursive = TRUE, full.names = TRUE)
  hits <- hits[!grepl("__MACOSX", hits)]
  if (length(hits) == 0) return(root)
  dirname(hits[[1]])
}

# ---------------------------------------------------------------------------
# Ollama backend — zero-shot alternative to the fine-tuned merged/adapter
# checkpoints above. Calls a running Ollama server's chat API over HTTP
# instead of loading weights into this process (no transformers/torch/PEFT
# needed for this path), so it works with an off-the-shelf tag like
# "gemma3:12b-it-qat" that was never fine-tuned on this ontology.
# ---------------------------------------------------------------------------

# Character vector of tags currently pulled on the server, e.g.
# "gemma3:12b-it-qat", or NULL if the server couldn't be reached.
ollama_list_tags <- function(host, timeout = 10) {
  h <- curl::new_handle(timeout = timeout, connecttimeout = 5)
  resp <- tryCatch(
    curl::curl_fetch_memory(paste0(sub("/+$", "", host), "/api/tags"), handle = h),
    error = function(e) NULL
  )
  if (is.null(resp) || resp$status_code != 200) return(NULL)
  parsed <- tryCatch(jsonlite::fromJSON(rawToChar(resp$content)), error = function(e) NULL)
  if (is.null(parsed) || is.null(parsed$models)) return(character(0))
  as.character(parsed$models$name)
}

# Pings the server and confirms the requested tag has actually been pulled
# there. Returns list(ok = TRUE/FALSE, msg = "...") for the Step 2 status badge.
validate_ollama <- function(host, model_tag) {
  host      <- trimws(host)
  model_tag <- trimws(model_tag)
  if (nchar(host) == 0)
    return(list(ok = FALSE, msg = "Enter an Ollama host URL."))
  if (nchar(model_tag) == 0)
    return(list(ok = FALSE, msg = "Enter a model tag (e.g. gemma3:12b-it-qat)."))

  tags <- ollama_list_tags(host)
  if (is.null(tags))
    return(list(ok = FALSE, msg = paste("Could not reach Ollama at", host)))

  # tolerate the implicit ":latest" Ollama appends to a bare (colon-less) tag
  hit <- model_tag %in% tags || paste0(model_tag, ":latest") %in% tags
  if (!hit)
    return(list(ok = FALSE, msg = sprintf(
      "%s not found on server (%d model(s) pulled) — run `ollama pull %s` first.",
      model_tag, length(tags), model_tag)))

  list(ok = TRUE, msg = sprintf("%s · reachable at %s", model_tag, host))
}

# One /api/chat call. `messages` is a list of list(role=, content=) — the same
# shape prepare_inference_input() already writes into the JSONL. Ollama's
# format = "json" forces valid-JSON output, matching the "Return JSON only"
# system prompt used for both backends. `temperature`, if given, is passed
# through Ollama's generation `options` (default 0.8 server-side if omitted).
ollama_chat <- function(host, model, messages, temperature = NULL, timeout = 600) {
  body <- list(model = model, messages = messages, stream = FALSE, format = "json")
  if (!is.null(temperature)) body$options <- list(temperature = temperature)
  h <- curl::new_handle(timeout = timeout, connecttimeout = 10)
  curl::handle_setheaders(h, "Content-Type" = "application/json")
  curl::handle_setopt(h, postfields = jsonlite::toJSON(body, auto_unbox = TRUE))
  resp <- curl::curl_fetch_memory(paste0(sub("/+$", "", host), "/api/chat"), handle = h)
  if (resp$status_code != 200)
    stop("Ollama request failed (HTTP ", resp$status_code, "): ", rawToChar(resp$content))
  parsed <- jsonlite::fromJSON(rawToChar(resp$content))
  parsed$message$content
}

# Reasoning-trace wrappers to discard before parsing — the R counterpart of
# infer.py's _THINK_RE. Gemma 4 emits "<|channel>thought ... <channel|>" when
# thinking is enabled; <think> tags are the common convention elsewhere.
# Ollama's format="json" usually suppresses these, but a model that leaks one
# shouldn't silently cost us the whole prediction.
.THINK_PATTERNS <- c(
  "<\\|?channel\\|?>\\s*thought\\b.*?<\\|?/?channel\\|?>",
  "<think>.*?</think>"
)

# Best-effort JSON object out of a model reply; returns a named list, or NULL
# if nothing usable is there. Mirrors _extract_json() in infer.py so both
# backends fail (and succeed) on the same inputs.
extract_json_object <- function(text) {
  if (is.null(text) || length(text) != 1 || is.na(text)) return(NULL)
  s <- trimws(text)
  for (p in .THINK_PATTERNS) s <- gsub(p, "", s, perl = TRUE, ignore.case = TRUE)
  s <- trimws(s)
  s <- sub("^```[A-Za-z]*\\n?", "", s)
  s <- sub("\\n?```$", "", s)
  s <- trimws(s)

  as_obj <- function(x) {
    o <- tryCatch(jsonlite::fromJSON(x, simplifyVector = TRUE), error = function(e) NULL)
    if (is.list(o) && !is.null(names(o))) o else NULL
  }

  hit <- as_obj(s)
  if (!is.null(hit)) return(hit)

  # Fall back to the first balanced {...} span, so a "Sure, here is the JSON:"
  # preamble or a trailing remark doesn't cost us the answer. Brace counting is
  # string-aware so braces inside values don't throw off the depth.
  chars  <- strsplit(s, "", fixed = TRUE)[[1]]
  starts <- which(chars == "{")
  for (st in starts) {
    depth <- 0L; in_str <- FALSE; esc <- FALSE
    for (i in seq(st, length(chars))) {
      ch <- chars[[i]]
      if (in_str) {
        if (esc)            esc <- FALSE
        else if (ch == "\\") esc <- TRUE
        else if (ch == '"')  in_str <- FALSE
      } else if (ch == '"') { in_str <- TRUE
      } else if (ch == "{")  { depth <- depth + 1L
      } else if (ch == "}")  {
        depth <- depth - 1L
        if (depth == 0L) {
          cand <- paste(chars[st:i], collapse = "")
          obj  <- as_obj(cand)
          if (!is.null(obj)) return(obj)
          break
        }
      }
    }
  }
  NULL
}

# TRUE for predictions that did not produce a usable label — i.e. the rows
# labelled "unknown" because the request failed or the reply couldn't be
# parsed. Both are transient in principle (a re-request can succeed), which is
# what makes a retry pass worthwhile.
failed_prediction_idx <- function(preds) {
  if (is.null(preds) || nrow(preds) == 0) return(integer(0))
  if (is.null(preds$error)) return(integer(0))
  which(!is.na(preds$error))
}

# Re-runs only the prompts whose predictions failed, merging any successes back
# into output_jsonl. Used by the CLI (automatically) and the Shiny app (behind
# a button), so both recover the same way.
#
# Returns list(before=, after=, recovered=, attempts=).
retry_failed_ollama <- function(input_jsonl, output_jsonl, host, model,
                                temperature = NULL, timeout = 600,
                                max_attempts = 1L, on_progress = NULL) {
  prompt_lines <- readLines(input_jsonl, warn = FALSE)
  prompt_lines <- prompt_lines[nchar(trimws(prompt_lines)) > 0]

  preds <- jsonlite::stream_in(file(output_jsonl), verbose = FALSE)
  failed <- failed_prediction_idx(preds)
  before <- length(failed)
  if (before == 0)
    return(list(before = 0L, after = 0L, recovered = 0L, attempts = 0L))

  attempt <- 0L
  while (attempt < max_attempts && length(failed) > 0) {
    attempt <- attempt + 1L
    still   <- integer(0)

    for (k in seq_along(failed)) {
      i  <- failed[[k]]
      ex <- jsonlite::fromJSON(prompt_lines[[i]], simplifyVector = FALSE)

      err  <- NA_character_
      text <- tryCatch(
        ollama_chat(host, model, ex$messages, temperature = temperature, timeout = timeout),
        error = function(e) { err <<- conditionMessage(e); NA_character_ }
      )
      pred <- if (!is.na(text)) extract_json_object(text) else NULL

      if (!is.null(pred) && !is.null(pred$cell_type)) {
        preds$cell_type[[i]]  <- as.character(pred$cell_type)
        preds$confidence[[i]] <- if (is.null(pred$confidence)) NA_real_
                                 else as.numeric(pred$confidence)
        preds$error[[i]] <- NA_character_
        if (!is.null(preds$raw)) preds$raw[[i]] <- NA_character_
      } else {
        preds$error[[i]] <- if (!is.na(err)) "request_failure" else "parse_failure"
        if (!is.null(preds$raw))
          preds$raw[[i]] <- if (!is.na(err)) err else if (is.na(text)) "" else text
        still <- c(still, i)
      }
      if (!is.null(on_progress)) on_progress(k, length(failed), attempt)
    }
    failed <- still
  }

  # Rewrite the predictions file so downstream consumers see the merged result.
  con <- file(output_jsonl, "w"); on.exit(close(con))
  for (i in seq_len(nrow(preds))) {
    row <- as.list(preds[i, , drop = FALSE])
    row <- row[!vapply(row, function(x) length(x) == 0 || all(is.na(x)), logical(1))]
    writeLines(jsonlite::toJSON(row, auto_unbox = TRUE), con)
  }

  list(before = before, after = length(failed),
       recovered = before - length(failed), attempts = attempt)
}

# N-sample voting for clusters that came back "Unclassified".
#
# The model's own `confidence` field is not informative (it returns ~0.95 for
# almost everything, including clusters with no elevated markers at all), so
# this samples the SAME prompt n_samples times and uses the agreement rate as
# an empirical confidence instead. It needs vote_temperature > 0: the main run
# is greedy (temperature 0), which would return ten identical samples.
#
# Two things can produce "Unclassified", and they mean different things here:
#   * the min_z floor fired  — the prompt literally says no markers are
#     elevated, so a majority label would be the model inventing a lineage from
#     nothing. This is exactly what the floor exists to prevent, so such
#     clusters are re-sampled for DIAGNOSIS but never overridden.
#   * the model declined despite having elevated markers — here a stable
#     majority is meaningful, and is allowed to override.
#
# Adds vote_* columns to the predictions file either way, so the distribution
# is visible even when the label is left alone.
vote_unclassified_ollama <- function(input_jsonl, output_jsonl, host, model,
                                     n_samples = 10L, vote_temperature = 0.7,
                                     min_agreement = 0.6, timeout = 600,
                                     override_floor = FALSE, on_progress = NULL) {
  prompt_lines <- readLines(input_jsonl, warn = FALSE)
  prompt_lines <- prompt_lines[nchar(trimws(prompt_lines)) > 0]
  preds <- jsonlite::stream_in(file(output_jsonl), verbose = FALSE)

  is_uncl <- grepl("unclassif", as.character(preds$cell_type), ignore.case = TRUE)
  idx <- which(is_uncl)
  if (length(idx) == 0)
    return(list(n = 0L, resolved = 0L, preds = preds))

  for (col in c("vote_top","vote_agreement","vote_n","vote_basis","vote_detail")) {
    if (is.null(preds[[col]]))
      preds[[col]] <- if (col %in% c("vote_agreement","vote_n")) NA_real_ else NA_character_
  }

  norm <- function(x) {
    x <- tolower(trimws(x)); x <- gsub("[^a-z0-9 ]", " ", x)
    trimws(gsub("\\s+", " ", x))
  }
  resolved <- 0L

  for (k in seq_along(idx)) {
    i  <- idx[[k]]
    ex <- jsonlite::fromJSON(prompt_lines[[i]], simplifyVector = FALSE)
    user_txt <- ex$messages[[2]]$content
    floor_fired <- grepl("No markers are elevated", user_txt, fixed = TRUE)

    votes <- character(0)
    for (j in seq_len(n_samples)) {
      txt <- tryCatch(ollama_chat(host, model, ex$messages,
                                  temperature = vote_temperature, timeout = timeout),
                      error = function(e) NA_character_)
      o <- if (!is.na(txt)) extract_json_object(txt) else NULL
      if (!is.null(o) && !is.null(o$cell_type)) votes <- c(votes, as.character(o$cell_type))
      if (!is.null(on_progress)) on_progress(k, length(idx), j, n_samples)
    }
    if (length(votes) == 0) next

    tab   <- sort(table(norm(votes)), decreasing = TRUE)
    top   <- names(tab)[1]
    agree <- as.numeric(tab[[1]]) / length(votes)
    # report the winner in the model's own casing
    top_disp <- votes[which(norm(votes) == top)[1]]

    preds$vote_top[[i]]       <- top_disp
    preds$vote_agreement[[i]] <- round(agree, 3)
    preds$vote_n[[i]]         <- length(votes)
    preds$vote_basis[[i]]     <- if (floor_fired) "no_elevated_markers" else "model_declined"
    preds$vote_detail[[i]]    <- paste(sprintf("%s:%d", names(tab), as.integer(tab)), collapse = "; ")

    may_override <- (!floor_fired || isTRUE(override_floor)) &&
                    !grepl("unclassif", top, ignore.case = TRUE) &&
                    agree >= min_agreement
    if (may_override) {
      preds$cell_type[[i]] <- top_disp
      preds$confidence[[i]] <- round(agree, 3)   # empirical, not self-reported
      resolved <- resolved + 1L
    }
  }

  con <- file(output_jsonl, "w"); on.exit(close(con))
  for (i in seq_len(nrow(preds))) {
    row <- as.list(preds[i, , drop = FALSE])
    row <- row[!vapply(row, function(x) length(x) == 0 || all(is.na(x)), logical(1))]
    writeLines(jsonlite::toJSON(row, auto_unbox = TRUE), con)
  }
  list(n = length(idx), resolved = resolved, preds = preds)
}

# Mirrors infer.py's run loop but calls Ollama's HTTP API per prompt instead
# of loading a model. on_progress(i, n), if given, fires after each
# prediction so the caller can drive a Shiny progress bar.
run_ollama_inference <- function(input_jsonl, output_jsonl, host, model,
                                 temperature = NULL, on_progress = NULL,
                                 timeout = 600) {
  lines <- readLines(input_jsonl, warn = FALSE)
  lines <- lines[nchar(trimws(lines)) > 0]
  n     <- length(lines)

  con <- file(output_jsonl, "w"); on.exit(close(con))

  for (i in seq_along(lines)) {
    ex <- jsonlite::fromJSON(lines[i], simplifyVector = FALSE)

    # Keep the failure reason: a swallowed message here is undiagnosable in a
    # batch log, and request failures (timeouts especially) look identical to
    # malformed-JSON failures once the text is gone.
    err  <- NA_character_
    text <- tryCatch(
      ollama_chat(host, model, ex$messages, temperature = temperature, timeout = timeout),
      error = function(e) { err <<- conditionMessage(e); NA_character_ }
    )

    pred <- if (!is.na(text)) extract_json_object(text) else NULL
    if (is.null(pred) || is.null(pred$cell_type)) {
      pred <- list(cell_type = "unknown", confidence = 0,
                   error = if (!is.na(err)) "request_failure" else "parse_failure",
                   raw   = if (!is.na(err)) err else text)
    }
    pred$cluster_id <- ex$cluster_id
    pred$cell_index <- ex$cell_index

    writeLines(jsonlite::toJSON(pred, auto_unbox = TRUE), con)
    if (!is.null(on_progress)) on_progress(i, n)
  }

  output_jsonl
}

# Collapses near-duplicate cell_type labels (e.g. "Treg cell" / "T regulatory"
# / "Regulatory T cell") into a single canonical spelling, via one extra
# Ollama call over the unique label set — the same harmonisation step
# automated_phenotyping runs after its per-cluster GPT calls, adapted from
# OpenAI to Ollama. Returns a named character vector (raw label -> canonical
# label) for every unique input label, or NULL if harmonisation fails (the
# caller should then just keep the raw labels rather than error out).
harmonize_ollama_labels <- function(host, model, labels, temperature = NULL,
                                    timeout = 600) {
  unique_labels <- unique(labels)
  if (length(unique_labels) <= 1) return(NULL)

  prompt <- paste0(
    "Below is a list of cell type annotation labels produced by automated immunophenotyping. ",
    "Many entries are near-duplicates that differ only in pluralisation, spacing, capitalisation, ",
    "or minor wording (e.g. 'CD4 T cell' and 'CD4 T cells', or 'NK cell' and 'Natural Killer cell'). ",
    "Return a JSON object mapping EVERY label in the input list to its canonical standardised form. ",
    "Use Title Case. Do not merge biologically distinct populations. ",
    "Return ONLY the JSON object, no markdown fences, no explanation.\n\n",
    "Labels: ", paste(unique_labels, collapse = ", ")
  )
  messages <- list(
    list(role = "system",
         content = "You are an expert immunologist. Return ONLY valid JSON, no markdown, no explanation."),
    list(role = "user", content = prompt)
  )

  raw <- tryCatch(ollama_chat(host, model, messages, temperature = temperature,
                              timeout = timeout),
                  error = function(e) { message("harmonisation request failed: ", conditionMessage(e)); NA_character_ })
  if (is.na(raw)) return(NULL)

  mapping <- extract_json_object(raw)
  if (is.null(mapping)) return(NULL)
  mapping <- unlist(mapping)
  if (length(mapping) == 0) return(NULL)
  mapping
}
