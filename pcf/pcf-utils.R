# pcf-utils.R
# ─────────────────────────────────────────────────────────────────────────────
# Core analysis functions for the PCF (pair correlation function) module —
# headless port of the PhenoSuite `pcf-v2` Shiny app.
#
# Sourced by run-pcf.R (CLI).
#
# The GUI ran R → reticulate → vectra_lib_v4.py → rpy2 → spatstat, i.e. it
# bounced back into R for every actual computation. The maths therefore always
# lived in spatstat; this port calls spatstat directly and drops the rpy2 /
# reticulate round-trip entirely, so nothing but R is needed at runtime.
#
# Ported logic and its source:
#   utils/spatial-shiny/vectra_lib_v4.py :: extract_data()      -> extract_vectra_data()
#   utils/spatial-shiny/vectra_lib_v4.py :: pcf()               -> compute_pcf()
#   utils/spatial-shiny/vectra_lib_v4.py :: plot_pcf_curves()   -> pcf_curve_tables()
#   utils/spatial-shiny/plot-pfcs-dev.R  :: plot_pcfs_R()       -> plot_pcf_curve_grid()
#   pcf-v2/production/server.R           :: AUC violin branches -> pcf_auc_table() +
#                                                                  plot_pcf_auc_violins()
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(spatstat.geom)
  library(spatstat.explore)
  library(ggplot2)
  library(ggpubr)
  library(gridExtra)
  library(glue)
  library(jsonlite)
})

`%||%` <- function(a, b) if (!is.null(a) && length(a) && !all(is.na(a))) a else b

# spatstat.explore 3.8-1 changed the pcf defaults: `zerocor` gained a
# non-"none" default and the divisor menu grew. The GUI ran against the older
# behaviour, so the legacy settings are requested explicitly where the
# installed spatstat understands them — this keeps curves comparable with
# results produced by the app, and silences the "defaults have changed" note.
.pcf_legacy_args <- function() {
  fmls <- names(formals(spatstat.explore::pcfinhom))
  args <- list()
  if ("divisor" %in% fmls) args$divisor <- "r"
  if ("zerocor" %in% fmls) args$zerocor <- "none"
  args
}

# Number of radius steps evaluated per curve. The GUI coerced spatstat's
# ~513-point default grid to 500 so every curve is the same length and can be
# stacked into a matrix across samples and interactions.
PCF_N_STEPS <- 500L

# ─── Input parsing ───────────────────────────────────────────────────────────

# Vectra cell_seg_data tables ship as either tab- or comma-delimited .csv.
# Mirrors vectra_lib_v4.py::read_csv_tsv(): try tab first, fall back to comma
# when the expected columns are absent. Legacy exports use underscores in
# column names ('Sample_Name'); those are normalised to spaces.
PCF_REQUIRED_COLS <- c("Cell X Position", "Cell Y Position")

read_vectra_table <- function(path) {
  read_one <- function(sep) {
    tryCatch(
      read.delim(path, sep = sep, check.names = FALSE, stringsAsFactors = FALSE),
      error = function(e) NULL
    )
  }
  normalise <- function(df) {
    if (is.null(df)) return(NULL)
    if ("Sample_Name" %in% names(df)) names(df) <- gsub("_", " ", names(df), fixed = TRUE)
    df
  }
  df <- normalise(read_one("\t"))
  if (is.null(df) || !all(PCF_REQUIRED_COLS %in% names(df))) {
    df2 <- normalise(read_one(","))
    if (!is.null(df2) && all(PCF_REQUIRED_COLS %in% names(df2))) df <- df2
  }
  if (is.null(df))
    stop("could not read '", path, "' as a tab- or comma-delimited table")
  missing_cols <- setdiff(PCF_REQUIRED_COLS, names(df))
  if (length(missing_cols))
    stop("'", basename(path), "' is missing required column(s): ",
         paste(missing_cols, collapse = ", "))
  df
}

# One entry per CSV: PCF treats one file as one sample, exactly as the GUI did
# (extract_data() reads 'Sample Name' from the first row only).
#
# The GUI's score-file branch (thresholding stain intensities against
# *_score_data.txt) is not ported: pcf() only ever reads Cell X/Y Position,
# Tissue Category and Phenotype, so the stain columns it produced never
# reached any PCF computation.
extract_vectra_data <- function(paths, phenotype_col = "Phenotype",
                                drop_na = TRUE, verbose = TRUE) {
  samples <- list()
  for (p in paths) {
    file <- read_vectra_table(p)
    if (!phenotype_col %in% names(file))
      stop("'", basename(p), "' has no '", phenotype_col, "' column")

    tissue <- if ("Tissue Category" %in% names(file)) {
      as.character(file[["Tissue Category"]])
    } else {
      if (verbose)
        message("[PCF]   ", basename(p),
                ": no 'Tissue Category' column — filling 'Unknown' ",
                "(unused by the PCF computation)")
      rep("Unknown", nrow(file))
    }

    df <- data.frame(
      `Cell X Position` = suppressWarnings(as.numeric(file[["Cell X Position"]])),
      `Cell Y Position` = suppressWarnings(as.numeric(file[["Cell Y Position"]])),
      `Tissue Category` = tissue,
      Phenotype         = as.character(file[[phenotype_col]]),
      check.names       = FALSE,
      stringsAsFactors  = FALSE
    )
    if (drop_na)
      df <- df[!is.na(df$Phenotype) & nzchar(trimws(df$Phenotype)), , drop = FALSE]
    df <- df[is.finite(df[["Cell X Position"]]) &
             is.finite(df[["Cell Y Position"]]), , drop = FALSE]
    if (nrow(df) == 0)
      stop("'", basename(p), "' has no usable cells after dropping NA ",
           "phenotypes / coordinates")

    sample_name <- if ("Sample Name" %in% names(file) && nrow(file) > 0) {
      as.character(file[["Sample Name"]])[1]
    } else {
      tools::file_path_sans_ext(basename(p))
    }

    samples[[length(samples) + 1L]] <- list(
      file_path   = p,
      file_name   = tools::file_path_sans_ext(basename(p)),
      sample_name = sample_name,
      data        = df
    )
    if (verbose)
      message(sprintf("[PCF]   %-40s %6d cells, %2d phenotypes",
                      basename(p), nrow(df), length(unique(df$Phenotype))))
  }
  samples
}

# Cell types present in *every* input file — the GUI's checkbox list
# (server.R::mydata() keeps only phenotypes whose per-file indicator column
# sums to the file count).
shared_phenotypes <- function(samples) {
  per_file <- lapply(samples, function(s) unique(s$data$Phenotype))
  all_ph   <- unique(unlist(per_file))
  keep     <- vapply(all_ph, function(ph)
    all(vapply(per_file, function(x) ph %in% x, logical(1))), logical(1))
  all_ph[keep]
}

# ─── PCF computation ─────────────────────────────────────────────────────────

# Port of vectra_lib_v4.py::pcf(). Iterates the upper triangle (including the
# diagonal) of ['All'] + cell_types and calls the matching inhomogeneous
# spatstat estimator for each pair:
#
#   All  vs All   -> pcfinhom()        (regular)
#   All  vs type  -> pcfdot.inhom()    (dot)
#   type vs type  -> pcfcross.inhom()  (cross)
#
# Every estimator uses correction="isotropic" on the same 500-point radius
# grid (0 .. radius/resolution pixels), and each curve is divided by the
# intensity-based normalisation the GUI used, so curves are comparable across
# samples of different size and density.
compute_pcf <- function(samples, cell_types = NULL, phenotype = "Phenotype",
                        count_threshold = 10L, radius = 30, resolution = 0.377,
                        verbose = TRUE) {

  if (count_threshold < 1L)
    stop("count_threshold must be >= 1 (a cell type with zero cells cannot be ",
         "given to spatstat)")

  if (is.null(cell_types)) {
    all_ph     <- unlist(lapply(samples, function(s) s$data[[phenotype]]))
    cell_types <- names(sort(table(all_ph), decreasing = TRUE))
  }
  cell_types <- c("All", setdiff(unique(cell_types), "All"))

  rmax   <- radius / resolution              # microns -> pixels
  step   <- rmax / PCF_N_STEPS
  r_vals <- seq(0, by = step, length.out = PCF_N_STEPS)

  legacy    <- .pcf_legacy_args()
  rows      <- list()
  pcf_list  <- list()
  norm_list <- list()

  for (s in samples) {
    sel <- s$data
    ph  <- as.character(sel[[phenotype]])
    x   <- round(sel[["Cell X Position"]])
    y   <- round(sel[["Cell Y Position"]])

    X <- ppp(x, y,
             c(min(x) - 1, max(x) + 1),
             c(min(y) - 1, max(y) + 1),
             marks = factor(ph))
    area_sq <- ((max(x) - min(x)) * (max(y) - min(y)))^2

    count_of <- function(ct) if (ct == "All") length(ph) else sum(ph == ct)

    if (verbose)
      message("[PCF] ", s$file_name, ": ", npoints(X), " points, ",
              length(cell_types), " cell types (",
              length(cell_types) * (length(cell_types) + 1L) / 2L, " pairs)")

    for (i in seq_along(cell_types)) {
      cell_one <- cell_types[i]
      for (j in seq(i, length(cell_types))) {
        cell_two <- cell_types[j]

        n_one <- count_of(cell_one)
        n_two <- count_of(cell_two)

        below <- (n_one < count_threshold && cell_one != "All") ||
                 (n_two < count_threshold && cell_two != "All")

        if (below) {
          # Skipped pair. The GUI appended a short row here, which pandas
          # padded with NaN and so silently shifted every column after
          # 'Cell_two'; the fields are written to their proper columns below
          # and skipped rows are dropped again before plotting.
          if (verbose)
            message("[PCF]   skip ", cell_one, " vs ", cell_two,
                    " (counts ", n_one, "/", n_two,
                    " < count_threshold=", count_threshold, ")")
          rows[[length(rows) + 1L]] <- data.frame(
            Patient          = s$file_name,
            `Sample Name`    = s$sample_name,
            Cell_one         = cell_one,
            Cell_two         = cell_two,
            PCFsum           = NA_real_,
            normalization    = NA_real_,
            min_count        = min(n_one, n_two),
            count_one        = n_one,
            count_two        = n_two,
            skipped          = TRUE,
            check.names      = FALSE,
            stringsAsFactors = FALSE
          )
          pcf_list[[length(pcf_list) + 1L]]   <- NA
          norm_list[[length(norm_list) + 1L]] <- NA
          next
        }

        if (cell_two == "All") {                      # regular pcf
          lambda1     <- density(X, at = "points", leaveoneout = FALSE)
          lambdasums  <- sum(1 / lambda1)^2
          fv          <- do.call(pcfinhom, c(list(
                                  X, lambda1, correction = "isotropic",
                                  r = r_vals, renormalise = FALSE), legacy))
          cell_count  <- c(npoints(X), npoints(X))
        } else if (cell_one == "All") {               # dot pcf
          c2          <- X[marks(X) == cell_two]
          lambda1     <- density(X,  at = "points", leaveoneout = FALSE)
          lambda2     <- density(c2, at = "points", leaveoneout = FALSE)
          lambdasums  <- sum(1 / lambda1) * sum(1 / lambda2)
          fv          <- do.call(pcfdot.inhom, c(list(
                                  X, i = cell_two,
                                  lambdaI = lambda2, lambdadot = lambda1,
                                  r = r_vals, correction = "isotropic"), legacy))
          cell_count  <- c(npoints(X), npoints(c2))
        } else {                                      # cross pcf
          c1          <- X[marks(X) == cell_one]
          c2          <- X[marks(X) == cell_two]
          lambda1     <- density(c1, at = "points", leaveoneout = FALSE)
          lambda2     <- density(c2, at = "points", leaveoneout = FALSE)
          lambdasums  <- sum(1 / lambda1) * sum(1 / lambda2)
          fv          <- do.call(pcfcross.inhom, c(list(
                                  X, i = cell_one, j = cell_two,
                                  lambdaI = lambda1, lambdaJ = lambda2,
                                  r = r_vals, correction = "isotropic"), legacy))
          cell_count  <- c(npoints(c1), npoints(c2))
        }

        pcf_vals    <- as.numeric(fv$iso)
        pcf_vals[1] <- 0                    # r = 0 is Inf by construction
        normalization <- lambdasums / area_sq
        norm_vals     <- pcf_vals / normalization

        rows[[length(rows) + 1L]] <- data.frame(
          Patient          = s$file_name,
          `Sample Name`    = s$sample_name,
          Cell_one         = cell_one,
          Cell_two         = cell_two,
          PCFsum           = sum(norm_vals),
          normalization    = normalization,
          min_count        = min(cell_count),
          count_one        = cell_count[1],
          count_two        = cell_count[2],
          skipped          = FALSE,
          check.names      = FALSE,
          stringsAsFactors = FALSE
        )
        pcf_list[[length(pcf_list) + 1L]]   <- pcf_vals
        norm_list[[length(norm_list) + 1L]] <- norm_vals
      }
    }
  }

  out          <- do.call(rbind, rows)
  out$PCF      <- pcf_list
  out$normPCF  <- norm_list
  attr(out, "r_pixels")   <- r_vals
  attr(out, "r_um")       <- r_vals * resolution
  attr(out, "cell_types") <- cell_types
  out
}

# ─── Curve tables ────────────────────────────────────────────────────────────

# vectra_lib_v4.py::interaction_subset() — both orderings of a pair, filtered
# on min_count.
interaction_subset <- function(pcf_df, cell1, cell2, min_count = 0) {
  keep <- (((pcf_df$Cell_one == cell1 & pcf_df$Cell_two == cell2) |
            (pcf_df$Cell_one == cell2 & pcf_df$Cell_two == cell1)) &
           pcf_df$min_count > min_count)
  pcf_df[keep, , drop = FALSE]
}

# vectra_lib_v4.py::plot_pcf_curves() — one table per cell type holding every
# interaction that cell type takes part in, named c('All', cell_types) as
# server.R did before writing ppc.rds.
pcf_curve_tables <- function(pcf_df, cell_types, min_count = 0) {
  cts <- c("All", setdiff(unique(cell_types), "All"))
  ppc <- lapply(cts, function(ct) {
    parts <- lapply(cts, function(other) interaction_subset(pcf_df, ct, other,
                                                            min_count = min_count))
    do.call(rbind, parts)
  })
  names(ppc) <- cts
  ppc
}

# ─── Curve grid plot ─────────────────────────────────────────────────────────

# Port of plot-pfcs-dev.R::plot_pcfs_R(): one panel per cell type, one line per
# interaction partner, mean +/- variance ribbon across samples when more than
# one sample contributes.
#
# Two fixes over the original: the partner cell type is read per row rather
# than inferred from row position (which assumed every pair survived the
# count threshold in every sample), and each panel is titled with its own cell
# type instead of indexing the column names by panel number.
plot_pcf_curve_grid <- function(ppc, out_file, resolution, radius, verbose = TRUE) {
  stepsize <- (radius / resolution) / PCF_N_STEPS
  mu       <- "µ"
  cts      <- names(ppc)

  grobs <- list()
  for (i in seq_along(ppc)) {
    this_ct <- cts[i]
    tbl     <- ppc[[i]]
    if (is.null(tbl) || nrow(tbl) == 0) next

    keep <- vapply(tbl$PCF, function(v) is.numeric(v) && length(v) > 0, logical(1))
    if (!any(keep)) next
    if (verbose && any(!keep))
      message("[PCF]   ", this_ct, ": dropping ", sum(!keep),
              " interaction(s) with no curve (below count_threshold)")
    tbl <- tbl[keep, , drop = FALSE]

    partner <- ifelse(tbl$Cell_one == this_ct, tbl$Cell_two, tbl$Cell_one)
    mat     <- do.call(cbind, tbl$PCF)
    n_r     <- nrow(mat)
    radius_um <- seq_len(n_r) * resolution * stepsize

    partners  <- unique(partner)
    multi     <- max(table(partner)) > 1L

    means <- vapply(partners, function(p)
      rowMeans(mat[, partner == p, drop = FALSE]), numeric(n_r))
    vars  <- vapply(partners, function(p) {
      sub <- mat[, partner == p, drop = FALSE]
      if (ncol(sub) < 2L) rep(0, n_r) else apply(sub, 1, stats::var)
    }, numeric(n_r))

    df <- data.frame(
      Radius   = rep(radius_um, times = length(partners)),
      variable = factor(rep(partners, each = n_r), levels = partners),
      mean     = as.vector(means),
      var      = as.vector(vars),
      stringsAsFactors = FALSE
    )

    p <- ggplot(df, aes(x = Radius, y = mean, group = variable)) +
      geom_line(aes(colour = variable), linewidth = 1)
    if (multi)
      p <- p + geom_ribbon(aes(ymin = mean - var, ymax = mean + var,
                               fill = variable), alpha = 0.2, colour = NA)
    p <- p +
      geom_hline(yintercept = 1) +
      xlab(glue("Radius ({mu}m)")) + ylab("PCF") + ggtitle(this_ct) +
      ylim(0, 2) +
      theme_bw() +
      theme(legend.title = element_blank(), legend.key = element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"))

    grobs[[length(grobs) + 1L]] <- p
  }

  if (length(grobs) == 0) {
    if (verbose) message("[PCF] no curves to plot — skipping curve grid")
    return(invisible(NULL))
  }
  # arrangeGrob() converts each ggplot to a grob, which needs a graphics device
  # for text metrics — without one R opens the default pdf device and leaves an
  # Rplots.pdf in whatever directory the job ran from. pdf(NULL) gives it a
  # device that writes nothing; ggsave() opens its own for the real output.
  g <- local({
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    gridExtra::arrangeGrob(grobs = grobs)
  })
  ggsave(out_file, g, width = 20, height = 20, limitsize = FALSE)
  invisible(out_file)
}

# ─── AUC table + violins ─────────────────────────────────────────────────────

# Long-format helper (replaces reshape2::melt). Building the frame directly
# keeps cell-type names intact — melt() mangles 'CD8+ T cells' into
# 'CD8..T.cells', which the GUI then had to patch back up with two gsub()s.
.pcf_melt <- function(df, id = NULL) {
  value_cols <- setdiff(names(df), id)
  out <- data.frame(
    variable = factor(rep(value_cols, each = nrow(df)), levels = value_cols),
    value    = unlist(df[value_cols], use.names = FALSE),
    stringsAsFactors = FALSE
  )
  if (!is.null(id)) out[[id]] <- rep(df[[id]], times = length(value_cols))
  out
}

# server.R's AUC assembly: for reference cell type `ref`, pull the normPCF
# series of every ref-vs-other interaction and stack samples end to end, one
# column per partner cell type plus a Sample column.
#
# Two deviations from the GUI, both deliberate:
#   * values are the raw normPCF series in every case. The GUI's one-sample
#     branch (python pcf_AUC) scaled them by `resolution` while its
#     multi-sample branch did not, so single- and multi-sample runs came out
#     in different units and could not be compared in pcf-builder.
#   * rows are assembled per sample so the Sample column always lines up with
#     its values, including when a pair was skipped in one sample only.
pcf_auc_table <- function(ppc, cell_types, ref, sample_names, label = NULL,
                          verbose = TRUE) {
  if (ref == "All")
    stop("--ref-celltype cannot be 'All' — 'All' is the comparison baseline ",
         "every other interaction is tested against")
  if (!ref %in% names(ppc))
    stop("reference cell type '", ref, "' not found among: ",
         paste(names(ppc), collapse = ", "))

  cell2_list <- c("All", setdiff(cell_types, c(ref, "All")))
  tbl        <- ppc[[ref]]

  n_r <- local({
    lens <- vapply(tbl$normPCF,
                   function(v) if (is.numeric(v)) length(v) else 0L, integer(1))
    if (!any(lens > 0))
      stop("no PCF curves were computed for reference cell type '", ref,
           "' — every interaction it takes part in was skipped. Lower ",
           "--count-threshold or --min-count, or pick a more abundant ",
           "--ref-celltype.")
    max(lens)
  })

  series <- function(sample, other) {
    hit <- tbl$Patient == sample &
           ((tbl$Cell_one == other & tbl$Cell_two == ref) |
            (tbl$Cell_one == ref   & tbl$Cell_two == other))
    v <- tbl$normPCF[hit]
    v <- v[vapply(v, is.numeric, logical(1))]
    if (length(v) == 0) rep(NA_real_, n_r) else as.numeric(v[[1]])
  }

  blocks <- lapply(sample_names, function(s) {
    block <- as.data.frame(
      lapply(cell2_list, function(other) series(s, other)),
      col.names = cell2_list, check.names = FALSE
    )
    names(block) <- cell2_list
    block$Sample <- rep(sub("-vectra_cell_seg_data.*$", "", s), n_r)
    block
  })

  out <- do.call(rbind, blocks)
  # Single-sample runs are labelled by the run, matching the GUI — one run is
  # one group when these tables are pooled in pcf-builder.
  if (length(sample_names) == 1L && !is.null(label) && nzchar(label))
    out$Sample <- rep(label, nrow(out))

  # Drop empty columns and empty rows before this table is written. Both only
  # ever hold NA — the interactions behind them were skipped below
  # count_threshold, and they are still recorded in <label>_pcf_summary.csv
  # with skipped=TRUE. Carrying them into the CSV would put a cell type with
  # no data in pcf-builder's celltype dropdown, and a sample with no data in
  # its reference-group dropdown, where selecting it makes
  # stat_compare_means() fail and leaves the reference line undrawn.
  value_cols <- setdiff(names(out), "Sample")
  empty_cols <- vapply(out[value_cols], function(v) all(is.na(v)), logical(1))
  if (any(empty_cols)) {
    if (verbose)
      message("[PCF]   dropping cell type(s) with no ", ref,
              " interaction in any sample: ",
              paste(value_cols[empty_cols], collapse = ", "))
    value_cols <- value_cols[!empty_cols]
    if (length(value_cols) == 0)
      stop("every '", ref, "' interaction was skipped — no AUC table to write. ",
           "Lower --count-threshold or --min-count, or pick a more abundant ",
           "--ref-celltype.")
    out <- out[, c(value_cols, "Sample"), drop = FALSE]
  }

  keep_rows <- rowSums(!is.na(out[, value_cols, drop = FALSE])) > 0
  if (any(!keep_rows)) {
    dropped <- unique(out$Sample[!keep_rows])
    if (verbose)
      message("[PCF]   dropping sample(s) with no ", ref,
              " interaction at all: ", paste(dropped, collapse = ", "))
    out <- out[keep_rows, , drop = FALSE]
  }

  # Contiguous row names: pcf-builder loads these files with
  # read.csv(..., row.names = 1), which needs the first column unique.
  rownames(out) <- NULL
  out
}

# pcf-builder's "rename samples into groups" step runs
#   gsub(<sample name>, <new name>, Sample)
# once per sample, i.e. it treats the value we write into the Sample column as
# a regular expression and as a substring. Two of our own naming choices can
# therefore mis-group data over there, so they are flagged here rather than
# discovered as a silently wrong plot:
#   * a name containing regex metacharacters never matches itself ('run+2'),
#     so renaming that group quietly does nothing;
#   * a name that is a prefix of another ('s1' inside 's10') rewrites part of
#     the longer one.
# The CSV itself is valid either way — only the app's optional rename is hit.
.check_builder_sample_names <- function(sample_values, verbose = TRUE) {
  if (!verbose) return(invisible(NULL))
  names_ <- unique(as.character(sample_values))

  meta <- names_[grepl("[][(){}.^$*+?|\\\\]", names_)]
  if (length(meta))
    message("[PCF]   NOTE: sample name(s) with regex metacharacters — ",
            "pcf-builder's group rename will not match them: ",
            paste(meta, collapse = ", "))

  prefixes <- names_[vapply(names_, function(n)
    any(n != names_ & startsWith(names_, n)), logical(1))]
  if (length(prefixes))
    message("[PCF]   NOTE: sample name(s) that are a prefix of another — ",
            "pcf-builder's group rename would also rewrite the longer one: ",
            paste(prefixes, collapse = ", "))

  invisible(NULL)
}

# The violin panels from server.R. Every interaction is compared against the
# 'All' baseline (Wilcoxon via stat_compare_means(ref.group='All')), with a
# red mean marker and a horizontal line at the 'All' mean.
.pcf_violin <- function(long_df, x, title, colour, signif_vs_all = TRUE,
                        legend = FALSE) {
  base_mean <- mean(long_df$value[long_df$variable == "All"], na.rm = TRUE)
  p <- ggviolin(long_df, x = x, y = "value", color = colour, add = "boxplot") +
    geom_hline(yintercept = base_mean) +
    xlab("") + ylab("norm PCF") + ggtitle(title)
  # The GUI also added a bare geom_signif(map_signif_level = TRUE, ...) layer
  # here. Without a `comparisons` argument stat_signif has nothing to test and
  # fails at draw time ("Computation failed in `stat_signif()`"), so that layer
  # never rendered a bracket in the app either — it only printed a warning per
  # panel. The p-values that do get drawn come from stat_compare_means() below,
  # which tests every interaction against the 'All' baseline.
  if (signif_vs_all) p <- p + stat_compare_means(ref.group = "All")
  if (!legend)       p <- p + theme(legend.position = "none")
  # Headroom on the flipped value axis, clipping off: stat_compare_means()
  # writes its labels past the data range and they are otherwise cut at the
  # panel edge after coord_flip().
  p + scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
    coord_flip(clip = "off") +
    stat_summary(fun = "mean", geom = "point", color = "red")
}

plot_pcf_auc_violins <- function(auc_df, ref, label, out_dir, verbose = TRUE) {
  out      <- list()
  samples  <- unique(auc_df$Sample)
  values   <- auc_df[, setdiff(names(auc_df), "Sample"), drop = FALSE]
  # Interactions skipped below count_threshold arrive here as all-NA columns
  # (see pcf_auc_table). Dropping them up front keeps ggplot from failing
  # loudly on empty panels.
  long_all <- .pcf_melt(values)
  long_all <- long_all[is.finite(long_all$value), , drop = FALSE]

  if (nrow(long_all) == 0) {
    if (verbose)
      message("[PCF]   no finite AUC values for '", ref,
              "' — skipping violin plots")
    return(out)
  }

  if (length(samples) == 1L) {
    path <- file.path(out_dir, glue("{label}-PCF_AUC_violins.pdf"))
    ggsave(path, .pcf_violin(long_all, "variable",
                             glue("{ref} Interactions"), "variable"),
           width = 10, height = 10)
    out$violins <- path
    return(out)
  }

  path <- file.path(out_dir, glue("{label}-PCF_AUC_violins-global.pdf"))
  ggsave(path, .pcf_violin(long_all, "variable",
                           glue("{ref} Interactions"), "variable"),
         width = 10, height = 10)
  out$violins_global <- path

  long_by <- .pcf_melt(auc_df, id = "Sample")
  long_by <- long_by[is.finite(long_by$value), , drop = FALSE]
  path <- file.path(out_dir, glue("{label}-PCF_AUC_violins-bySample.pdf"))
  ggsave(path, .pcf_violin(long_by, "Sample", glue("{ref} Interactions"),
                           "variable", signif_vs_all = FALSE, legend = TRUE),
         width = 10, height = 10)
  out$violins_by_sample <- path

  ind_dir <- file.path(out_dir, "individual-samples")
  dir.create(ind_dir, showWarnings = FALSE, recursive = TRUE)
  for (s in samples) {
    sub <- long_by[long_by$Sample == s, , drop = FALSE]
    if (nrow(sub) == 0) {
      if (verbose)
        message("[PCF]   ", s, ": every ", ref,
                " interaction was skipped — no per-sample violin written")
      next
    }
    path <- file.path(ind_dir, glue("{make.names(s)}-PCF_AUC_violins.pdf"))
    ggsave(path, .pcf_violin(sub, "variable", glue("{ref} Interactions"),
                             "variable"),
           width = 10, height = 10)
  }
  out$violins_individual <- ind_dir
  if (verbose) message("[PCF]   per-sample violins → ", ind_dir)
  out
}

# ─── Tidy curve export ───────────────────────────────────────────────────────

# One row per (sample, interaction, radius step). The GUI only ever wrote the
# curves as an .rds; this is the same data in a form other tools can read.
pcf_curves_long <- function(pcf_df, resolution) {
  r_um <- attr(pcf_df, "r_um")
  keep <- vapply(pcf_df$PCF, function(v) is.numeric(v) && length(v) > 0, logical(1))
  sub  <- pcf_df[keep, , drop = FALSE]
  if (nrow(sub) == 0) return(NULL)
  n_r  <- length(r_um)
  data.frame(
    Patient       = rep(sub$Patient,       each = n_r),
    `Sample Name` = rep(sub$`Sample Name`, each = n_r),
    Cell_one      = rep(sub$Cell_one,      each = n_r),
    Cell_two      = rep(sub$Cell_two,      each = n_r),
    radius_um     = rep(r_um, times = nrow(sub)),
    PCF           = unlist(sub$PCF,     use.names = FALSE),
    normPCF       = unlist(sub$normPCF, use.names = FALSE),
    check.names   = FALSE,
    stringsAsFactors = FALSE
  )
}

# ─── Orchestrator ────────────────────────────────────────────────────────────

run_pcf <- function(
  vectra_files,
  out_dir         = ".",
  label           = format(Sys.Date(), "%Y%m%d"),
  cell_types      = NULL,
  ref_celltype    = NULL,
  radius          = 30,
  resolution      = 0.377,
  count_threshold = 10L,
  min_count       = 0,
  phenotype_col   = "Phenotype",
  make_plots      = TRUE,
  verbose         = TRUE
) {
  msg <- function(...) if (verbose) message("[PCF] ", ...)

  out_dir <- file.path(out_dir, label)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out <- list()

  # 1. Load
  msg("Loading ", length(vectra_files), " Vectra annotation file(s)...")
  samples      <- extract_vectra_data(vectra_files, phenotype_col = phenotype_col,
                                      verbose = verbose)
  sample_names <- vapply(samples, function(s) s$file_name, character(1))

  # 2. Resolve cell types — default to the phenotypes shared by every file
  shared <- shared_phenotypes(samples)
  if (length(shared) == 0)
    stop("no phenotype is present in all ", length(samples),
         " input files — pass --celltypes explicitly")
  if (is.null(cell_types)) {
    cell_types <- shared
    msg("Cell types (shared across all inputs): ", paste(cell_types, collapse = ", "))
  } else {
    unknown <- setdiff(cell_types, shared)
    if (length(unknown))
      msg("WARNING: cell type(s) not present in every input file: ",
          paste(unknown, collapse = ", "))
    msg("Cell types (from --celltypes): ", paste(cell_types, collapse = ", "))
  }
  cell_types <- setdiff(unique(cell_types), "All")
  if (length(cell_types) < 1)
    stop("at least one cell type is required")

  if (is.null(ref_celltype)) {
    ref_celltype <- cell_types[1]
    msg("No --ref-celltype given; using '", ref_celltype, "'")
  }
  if (!ref_celltype %in% cell_types)
    stop("--ref-celltype '", ref_celltype, "' is not among the analysed cell ",
         "types: ", paste(cell_types, collapse = ", "))

  # 3. PCF
  msg("Computing PCF (radius=", radius, " microns, resolution=", resolution,
      " microns/px, count_threshold=", count_threshold, ")...")
  pcf_df <- compute_pcf(samples, cell_types = cell_types,
                        phenotype = "Phenotype",
                        count_threshold = count_threshold,
                        radius = radius, resolution = resolution,
                        verbose = verbose)
  n_skipped <- sum(pcf_df$skipped)
  if (n_skipped)
    msg(n_skipped, " of ", nrow(pcf_df),
        " interaction(s) skipped below count_threshold=", count_threshold)

  # 4. Tables
  summary_path <- file.path(out_dir, glue("{label}_pcf_summary.csv"))
  write.csv(pcf_df[, c("Patient", "Sample Name", "Cell_one", "Cell_two",
                       "PCFsum", "normalization", "min_count",
                       "count_one", "count_two", "skipped")],
            summary_path, row.names = FALSE)
  out$pcf_summary <- summary_path
  msg("Wrote interaction summary → ", summary_path)

  curves <- pcf_curves_long(pcf_df, resolution)
  if (!is.null(curves)) {
    curves_path <- file.path(out_dir, glue("{label}_pcf_curves.csv"))
    write.csv(curves, curves_path, row.names = FALSE)
    out$pcf_curves <- curves_path
    msg("Wrote curve table → ", curves_path)
  }

  # 5. Per-cell-type curve tables (the GUI's ppc.rds)
  ppc      <- pcf_curve_tables(pcf_df, cell_types, min_count = min_count)
  ppc_path <- file.path(out_dir, glue("{label}_ppc.rds"))
  saveRDS(ppc, ppc_path)
  out$ppc <- ppc_path

  # 6. AUC table — the input format pcf-builder reads (row names included, as
  #    pcf-builder loads these with read.csv(..., row.names = 1))
  msg("Building PCF AUC table (reference cell type: ", ref_celltype, ")...")
  auc_df   <- pcf_auc_table(ppc, cell_types, ref_celltype, sample_names,
                            label = label, verbose = verbose)
  .check_builder_sample_names(auc_df$Sample, verbose = verbose)
  auc_path <- file.path(out_dir, glue("{label}-PCF_AUCs.csv"))
  # row.names are written on purpose: pcf-builder loads these files with
  # read.csv(..., row.names = 1).
  write.csv(auc_df, auc_path)
  out$pcf_aucs <- auc_path
  msg("Wrote AUC table → ", auc_path)

  # 7. Plots
  if (make_plots) {
    msg("Generating PCF curve grid...")
    grid_path <- file.path(out_dir, glue("{label}-PCF-plots.pdf"))
    grid_out  <- plot_pcf_curve_grid(ppc, grid_path, resolution, radius,
                                     verbose = verbose)
    if (!is.null(grid_out)) out$pcf_plots <- grid_path

    msg("Generating PCF AUC violins...")
    out <- c(out, plot_pcf_auc_violins(auc_df, ref_celltype, label, out_dir,
                                       verbose = verbose))
  }

  # 8. Provenance
  prov_path <- file.path(out_dir, glue("{label}_provenance.json"))
  writeLines(toJSON(list(
    tool            = "PCF-CLI",
    label           = label,
    run_date        = as.character(Sys.time()),
    r_version       = paste(R.version$major, R.version$minor, sep = "."),
    spatstat_geom   = as.character(utils::packageVersion("spatstat.geom")),
    spatstat_explore = as.character(utils::packageVersion("spatstat.explore")),
    vectra_files    = unname(vectra_files),
    samples         = unname(sample_names),
    sample_names    = unname(vapply(samples, function(s) s$sample_name, character(1))),
    cell_types      = cell_types,
    ref_celltype    = ref_celltype,
    radius_um       = radius,
    resolution      = resolution,
    count_threshold = count_threshold,
    min_count       = min_count,
    n_radius_steps  = PCF_N_STEPS,
    n_interactions  = nrow(pcf_df),
    n_skipped       = n_skipped,
    make_plots      = make_plots,
    git_sha         = Sys.getenv("PHENOSUITE_GIT_SHA", "unknown"),
    image_digest    = Sys.getenv("PHENOSUITE_IMAGE_DIGEST", "unknown")
  ), pretty = TRUE, auto_unbox = TRUE), prov_path)
  out$provenance <- prov_path

  msg("Done. Outputs written to: ", out_dir)
  invisible(out)
}
