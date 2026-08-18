# RunGemmaPhenotyper.R
# ====================
# Core pipeline for the Gemma mIF phenotyper, batch/headless counterpart to
# PhenoSuite's gemma_phenotyper Shiny app.
#
# Annotates every cluster (or every cell) in a SpatialExperiment /
# SingleCellExperiment with an LLM-assigned cell type. Two backends:
#
#   hf        — HuggingFace-format weights (merged checkpoint or LoRA adapter,
#               fine-tuned or an off-the-shelf base) run through infer.py
#               (transformers/PEFT) as a subprocess. The only backend usable on
#               an HPC node with no outbound network.
#   ollama    — HTTP call to a running Ollama server. Convenient for local dev;
#               needs a reachable server, so usually not viable on a cluster.
#
# Prompt style is independent of backend:
#   zero-shot — tissue-aware, top/bottom 10% most distinctive markers only.
#               For models that were never fine-tuned on your ontology.
#   ontology  — every marker as marker=value plus an explicit ontology table.
#               For checkpoints fine-tuned against that ontology.
#
# Harmonisation (collapsing near-duplicate labels) runs on both backends: via a
# second HTTP call for ollama, or inside the same infer.py process for hf, so
# the model is loaded only once.
#
# Called by run-gemma-phenotyper.R; sources gemma-utils.R for prompt
# construction, inference, and harmonisation.

RunGemmaPhenotyper <- function(spe_file,
                               out_dir,
                               backend          = c("hf", "ollama"),
                               prompt_style     = c("zero-shot", "ontology"),
                               label            = NULL,
                               cluster_col      = NULL,
                               marker_cols      = NULL,
                               mode             = c("cluster", "cell"),
                               ollama_host      = "http://localhost:11434",
                               ollama_model     = NULL,
                               temperature      = 0,
                               harmonize        = TRUE,
                               timeout          = 600,
                               retry            = 1L,
                               tissue           = NULL,
                               min_z            = 0.5,
                               marker_glossary  = NULL,
                               vote_samples     = 0L,
                               vote_temperature = 0.7,
                               vote_min_agreement = 0.6,
                               vote_override_floor = FALSE,
                               model_dir        = NULL,
                               base_model       = NULL,
                               load_in_8bit     = FALSE,
                               ontology         = NULL,
                               python_bin       = "python3",
                               infer_script     = NULL,
                               export_loom      = FALSE) {

  backend      <- match.arg(backend)
  prompt_style <- match.arg(prompt_style)
  mode         <- match.arg(mode)

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  if (is.null(label)) label <- tools::file_path_sans_ext(basename(spe_file))

  ## -- Load object ----------------------------------------------------------
  message("[gemma] Loading ", spe_file)
  spe <- readRDS(spe_file)
  if (!inherits(spe, "SpatialExperiment") && !inherits(spe, "SingleCellExperiment"))
    stop("input must be a SpatialExperiment or SingleCellExperiment .rds (got ",
         class(spe)[1], ")")

  ## -- Resolve markers / cluster column -------------------------------------
  if (is.null(marker_cols)) {
    marker_cols <- rownames(spe)
    message("[gemma] No --markers given; using all ", length(marker_cols), " features")
  } else {
    missing <- setdiff(marker_cols, rownames(spe))
    if (length(missing) > 0)
      stop("markers not found in object: ", paste(missing, collapse = ", "))
  }

  if (mode == "cluster") {
    if (is.null(cluster_col))
      stop("--cluster-col is required in cluster mode (columns available: ",
           paste(names(colData(spe)), collapse = ", "), ")")
    if (!cluster_col %in% names(colData(spe)))
      stop("cluster column '", cluster_col, "' not in colData. Available: ",
           paste(names(colData(spe)), collapse = ", "))
    n_units <- length(unique(colData(spe)[[cluster_col]]))
  } else {
    cluster_col <- NULL
    n_units     <- ncol(spe)
  }
  message("[gemma] ", mode, " mode: ", n_units, " prompt(s) over ", ncol(spe), " cells")

  ## -- Build prompts --------------------------------------------------------
  tmp_in  <- file.path(out_dir, paste0(label, "_prompts.jsonl"))
  tmp_out <- file.path(out_dir, paste0(label, "_predictions.jsonl"))
  hf_hmap_path <- file.path(out_dir, paste0(label, "_harmonisation.json"))

  if (prompt_style == "zero-shot") {
    gloss <- if (!is.null(marker_glossary)) read_marker_glossary(marker_glossary) else NULL
    if (!is.null(gloss)) message("[gemma] marker glossary: ", length(gloss), " entries")
    prepare_ollama_prompts(spe, marker_cols = marker_cols, cluster_col = cluster_col,
                           tissue = tissue, min_z = min_z, marker_glossary = gloss,
                           out_path = tmp_in)
  } else {
    if (is.null(ontology)) ontology <- DEFAULT_GEMMA_ONTOLOGY
    prepare_inference_input(spe, marker_cols = marker_cols, ontology_text = ontology,
                            cluster_col = cluster_col, out_path = tmp_in)
  }

  ## -- Inference ------------------------------------------------------------
  if (backend == "ollama") {
    message("[gemma] Querying ", ollama_model, " at ", ollama_host,
            " (temperature=", temperature, ")")
    run_ollama_inference(tmp_in, tmp_out, ollama_host, ollama_model,
                         temperature = temperature, timeout = timeout,
                         on_progress = function(i, n) {
                           step <- max(1L, n %/% 10L)
                           if (i %% step == 0 || i == n)
                             message("[gemma]   ", i, " / ", n)
                         })
  } else {
    mtype <- detect_model_type(model_dir)
    if (mtype %in% c("not_found", "unknown"))
      stop("could not resolve model at '", model_dir,
           "' (expected config.json or adapter_config.json); detect_model_type() said: ", mtype)
    message("[gemma] ", mtype, " checkpoint at ", model_dir)

    args <- c(infer_script,
              if (mtype == "merged") c("--model_path", model_dir)
              else c("--adapter_path", model_dir, "--base_model", base_model),
              "--input_jsonl", tmp_in, "--output_jsonl", tmp_out)
    if (isTRUE(load_in_8bit)) args <- c(args, "--load_in_8bit")
    # Harmonise inside the same process so the checkpoint loads once, not twice.
    if (isTRUE(harmonize)) args <- c(args, "--harmonize_out", hf_hmap_path)
    # infer.py retries in-process so the checkpoint isn't reloaded.
    if (retry > 0) args <- c(args, "--retry", as.integer(retry))
    # infer.py votes in-process too, so the checkpoint is loaded once.
    if (vote_samples > 1) {
      args <- c(args, "--vote_samples", as.integer(vote_samples),
                      "--vote_temperature", vote_temperature,
                      "--vote_min_agreement", vote_min_agreement)
      if (isTRUE(vote_override_floor)) args <- c(args, "--vote_override_floor")
    }

    status <- system2(python_bin, args = args, stdout = "", stderr = "", wait = TRUE)
    if (!identical(status, 0L))
      stop("infer.py exited with status ", status)
  }

  # Ollama retries here (the HF backend already retried inside infer.py, while
  # the model was still loaded).
  if (backend == "ollama" && retry > 0) {
    rr <- tryCatch(
      retry_failed_ollama(tmp_in, tmp_out, ollama_host, ollama_model,
                          temperature = temperature, timeout = timeout,
                          max_attempts = as.integer(retry)),
      error = function(e) { message("[gemma] retry pass failed: ", conditionMessage(e)); NULL })
    if (!is.null(rr) && rr$before > 0)
      message("[gemma] Retry: recovered ", rr$recovered, " of ", rr$before,
              " failed prediction(s) in ", rr$attempts, " pass(es)")
  }

  # Ollama votes here; the hf backend already voted inside infer.py.
  if (backend == "ollama" && vote_samples > 1) {
    vres <- tryCatch(
      vote_unclassified_ollama(tmp_in, tmp_out, ollama_host, ollama_model,
                               n_samples = as.integer(vote_samples),
                               vote_temperature = vote_temperature,
                               min_agreement = vote_min_agreement,
                               override_floor = vote_override_floor,
                               timeout = timeout),
      error = function(e) { message("[gemma] voting failed: ", conditionMessage(e)); NULL })
    if (!is.null(vres) && vres$n > 0)
      message("[gemma] Voting: ", vres$n, " Unclassified cluster(s) re-sampled, ",
              vres$resolved, " resolved")
  }

  preds <- jsonlite::stream_in(file(tmp_out), verbose = FALSE)
  if (nrow(preds) == 0) stop("inference produced no predictions")

  n_failed <- if (is.null(preds$error)) 0L else sum(!is.na(preds$error))
  if (n_failed > 0) {
    message("[gemma] WARNING: ", n_failed, " of ", nrow(preds),
            " prediction(s) failed and were labelled 'unknown'")
    for (e in unique(preds$error[!is.na(preds$error)])) {
      idx <- which(preds$error == e)
      message("[gemma]   ", e, ": ", length(idx), " — e.g. ",
              substr(as.character(preds$raw[idx[1]]), 1, 160))
    }
    if (any(preds$error == "request_failure", na.rm = TRUE))
      message("[gemma]   (request_failure often means the ", timeout,
              "s timeout was hit — raise --timeout for large models on slow hardware)")
  }

  ## -- Harmonise ------------------------------------------------------------
  ## Collapses near-duplicate labels ("Treg cell" / "T regulatory" / ...) that
  ## arise because each cluster is labelled in an independent call. For hf the
  ## mapping was already produced in-process by infer.py; for ollama it needs a
  ## second HTTP call here.
  if (isTRUE(harmonize)) {
    message("[gemma] Harmonising ", length(unique(preds$cell_type)), " unique label(s)")
    hmap <- tryCatch({
      if (backend == "ollama") {
        harmonize_ollama_labels(ollama_host, ollama_model, preds$cell_type,
                                temperature = temperature, timeout = timeout)
      } else if (file.exists(hf_hmap_path)) {
        m <- unlist(jsonlite::fromJSON(hf_hmap_path))
        if (length(m) == 0) NULL else m
      } else {
        message("[gemma] no harmonisation sidecar at ", hf_hmap_path)
        NULL
      }
    },
      error = function(e) { message("[gemma] harmonisation failed: ", conditionMessage(e)); NULL })
    if (!is.null(hmap)) {
      hit <- preds$cell_type %in% names(hmap)
      preds$cell_type[hit] <- unname(hmap[preds$cell_type[hit]])
      message("[gemma]   -> ", length(unique(preds$cell_type)), " after harmonisation")
      write.csv(data.frame(Original = names(hmap), Harmonised = unname(hmap)),
                file.path(out_dir, paste0(label, "_harmonisation-map.csv")), row.names = FALSE)
    } else {
      message("[gemma] harmonisation skipped; keeping raw labels")
    }
  }

  ## -- Inject annotations ---------------------------------------------------
  if (mode == "cluster") {
    spe <- broadcast_cluster_labels(spe, preds, cluster_col)
  } else {
    spe[["gemma_cell_type"]]  <- as.character(preds$cell_type)
    spe[["gemma_confidence"]] <- as.numeric(preds$confidence)
  }

  ## -- Outputs --------------------------------------------------------------
  message("[gemma] Writing outputs to ", out_dir)
  saveRDS(spe, file.path(out_dir, paste0(label, "_spe_gemma_annotated.rds")))
  write.csv(as.data.frame(colData(spe)),
            file.path(out_dir, paste0(label, "_gemma_colData.csv")))
  if (!is.null(cluster_col)) {
    qc <- tryCatch(cluster_signal_qc(spe, marker_cols, cluster_col),
                   error = function(e) NULL)
    if (!is.null(qc)) preds <- merge(preds, qc, by = "cluster_id", all.x = TRUE, sort = FALSE)
  }
  write.csv(preds, file.path(out_dir, paste0(label, "_gemma_predictions.csv")),
            row.names = FALSE)

  # Vectra-format CSV needs spatial coords; SCE objects without them are skipped.
  if (inherits(spe, "SpatialExperiment")) {
    tryCatch(
      write_vectra_csv(spe, phenotype_col = "gemma_cell_type",
                       tissue_fallback = if (is.null(tissue)) "Unknown" else tissue,
                       out_dir = out_dir, file_prefix = paste0("vectra_gemma_", label)),
      error = function(e) message("[gemma] vectra export skipped: ", conditionMessage(e)))
  } else {
    message("[gemma] vectra export skipped: input is not a SpatialExperiment")
  }

  if (isTRUE(export_loom)) {
    tryCatch({
      cell_ids <- colnames(spe)
      if (is.null(cell_ids)) cell_ids <- paste0("cell_", seq_len(ncol(spe)))
      loom_mat <- .resolve_intensity_matrix(spe, marker_cols)$mat
      col_attrs <- list(CellID           = cell_ids,
                        gemma_cell_type  = as.character(spe[["gemma_cell_type"]]),
                        gemma_cell_type_marked = as.character(spe[["gemma_cell_type_marked"]]),
                        gemma_cluster_id = as.character(spe[["gemma_cluster_id"]]),
                        gemma_confidence = as.numeric(spe[["gemma_confidence"]]))
      if (inherits(spe, "SpatialExperiment")) {
        col_attrs$x <- spatialCoords(spe)[, 1]
        col_attrs$y <- spatialCoords(spe)[, 2]
      }
      write_loom(matrix = loom_mat,
                 row_attrs = list(Gene = clean_marker_name(rownames(loom_mat))),
                 col_attrs = col_attrs,
                 file_path = file.path(out_dir, paste0(label, "_gemma_annotated.loom")))
    }, error = function(e) message("[gemma] loom export skipped: ", conditionMessage(e)))
  }

  summary_tbl <- as.data.frame(table(spe[["gemma_cell_type"]], useNA = "ifany"),
                               stringsAsFactors = FALSE)
  names(summary_tbl) <- c("cell_type", "n_cells")
  summary_tbl <- summary_tbl[order(-summary_tbl$n_cells), ]
  write.csv(summary_tbl, file.path(out_dir, paste0(label, "_gemma_celltype-summary.csv")),
            row.names = FALSE)

  message("[gemma] Done: ", nrow(summary_tbl), " cell type(s) across ", ncol(spe), " cells")
  invisible(list(spe = spe, predictions = preds, summary = summary_tbl))
}
