#options(repos = c(CRAN = "https://cloud.r-project.org"))
#remotes::install_version("matrixStats", version="1.1.0")
plot_heatmap.mod=function(x, group_by, assay = "logcounts", out_dir = NULL,size.col=4,size.row=4,auto=F,segment=F){
  #remotes::install_version("matrixStats", version="1.1.0")
  
  if (!is(x, "SpatialExperiment")) {
    stop("input is not a SpatialExperiment object")
  }
  if (!is.character(group_by)) {
    stop("`group_by` is not a character string")
  }
  if (!all(group_by %in% names(colData(x)))) {
    stop("not all `group_by` values are present in the object")
  }
  if (!is.character(assay)) {
    stop("`assay` is not a character string")
  }
  if (!assay %in% assayNames(x)) {
    stop("input SpatialExperiment object does not have a `", 
         assay, "` assay")
  }
  if (!is.null(out_dir)) {
    if (!dir.exists(out_dir)) {
      stop("output directory `", out_dir, "` does not exist")
    }
  }
  gradient_colors <- rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
  for (g in group_by) {
    e <- scuttle::summarizeAssayByGroup(x, ids = colData(x)[[g]], 
                                        assay.type = assay, statistics = "median")
    
    e <- SummarizedExperiment::assay(e, i = "median")
    
    if(sum(rowSums(e)==0) == 0){
      e <- scale(t(e))
    }else{
      e_ <- scuttle::summarizeAssayByGroup(x, ids = colData(x)[[g]], 
                                           assay.type = assay, statistics = "mean")
      e_ <- SummarizedExperiment::assay(e_, i = "mean")
      e[rowSums(e)==0,]=e_[rowSums(e)==0,]
      e <- scale(t(e))
    }
    
    e[e > 5] <- 5
    e[e < -5] <- -5
    if(auto){
      size.col=8
      cluster.size=as.numeric(unlist(strsplit(g,'_clust'))[2])
      if(cluster.size <= 20){
        size.row=8
        e.tmp=e
        if (min(e.tmp) < 0) {
          col_fun <- circlize::colorRamp2(
            breaks = c(min(e.tmp), 0, max(e.tmp)),
            colors = gradient_colors
          )
          hm <- ComplexHeatmap::Heatmap(e.tmp, name = "Expression", 
                                        row_title = g, 
                                        col = col_fun,
                                        cluster_rows = TRUE, 
                                        cluster_columns = TRUE,
                                        column_names_gp = grid::gpar(fontsize=size.col),
                                        row_names_gp = grid::gpar(fontsize=size.row))
        } else {
          hm <- ComplexHeatmap::Heatmap(e.tmp, name = "Expression", 
                                        row_title = g, 
                                        col = gradient_colors,
                                        cluster_rows = TRUE, 
                                        cluster_columns = TRUE,
                                        column_names_gp = grid::gpar(fontsize=size.col),
                                        row_names_gp = grid::gpar(fontsize=size.row))
        }
        if (is.null(out_dir)) {
          return(hm)
        } else {
          if(cluster.size > 10){
            out_base <- glue("{out_dir}/{g}-heatmap")
            png(filename = glue("{out_base}.png"), width = 10, 
                height = 5, units = "in", res = 300)
            ComplexHeatmap::draw(hm)
            dev.off()
          }
        }
        
      }else{
        if(segment){
          clust.increments=c(seq(1,cluster.size,20),cluster.size)
          
          if(cluster.size <= 30){
            size.row=3
          }else{
            size.row=2 
          }
          
          for(i in seq(2,length(clust.increments))){
            e.tmp=e[seq(clust.increments[i-1],clust.increments[i]),]
            if (min(e.tmp) < 0) {
              col_fun <- circlize::colorRamp2(
                breaks = c(min(e.tmp), 0, max(e.tmp)),
                colors = gradient_colors
              )
              hm <- ComplexHeatmap::Heatmap(e.tmp, name = "Expression",
                                            row_title = g,
                                            col = col_fun,
                                            cluster_rows = TRUE,
                                            cluster_columns = TRUE,
                                            column_names_gp = grid::gpar(fontsize=size.col),
                                            row_names_gp = grid::gpar(fontsize=size.row))
            } else {
              hm <- ComplexHeatmap::Heatmap(e.tmp, name = "Expression",
                                            row_title = g,
                                            col = gradient_colors,
                                            cluster_rows = TRUE,
                                            cluster_columns = TRUE,
                                            column_names_gp = grid::gpar(fontsize=size.col),
                                            row_names_gp = grid::gpar(fontsize=size.row))
            }
            if (is.null(out_dir)) {
              return(hm)
            } else {
              out_base <- glue("{out_dir}/{g}-{i-1}-heatmap")
              png(filename = glue("{out_base}.png"), width = 10, 
                  height = 5, units = "in", res = 300)
              ComplexHeatmap::draw(hm)
              dev.off()
            }
          }
        }else{
          if(cluster.size <= 30){
            size.row=3
          }else{
            size.row=2 
          }
          e.tmp=e
          if (min(e.tmp) < 0) {
            col_fun <- circlize::colorRamp2(
              breaks = c(min(e.tmp), 0, max(e.tmp)),
              colors = gradient_colors
            )
            hm <- ComplexHeatmap::Heatmap(e.tmp, name = "Expression",
                                          row_title = g,
                                          col = col_fun,
                                          cluster_rows = TRUE,
                                          cluster_columns = TRUE,
                                          column_names_gp = grid::gpar(fontsize=size.col),
                                          row_names_gp = grid::gpar(fontsize=size.row))
          } else {
            hm <- ComplexHeatmap::Heatmap(e.tmp, name = "Expression",
                                          row_title = g,
                                          col = gradient_colors,
                                          cluster_rows = TRUE,
                                          cluster_columns = TRUE,
                                          column_names_gp = grid::gpar(fontsize=size.col),
                                          row_names_gp = grid::gpar(fontsize=size.row))
          }
          if (is.null(out_dir)) {
            return(hm)
          } else {
            out_base <- glue("{out_dir}/{g}-heatmap")
            png(filename = glue("{out_base}.png"), width = 10,
                height = 5, units = "in", res = 300)
            ComplexHeatmap::draw(hm)
            dev.off()
          }
        }
      }
      
    }else{
      if (min(e) < 0) {
        col_fun <- circlize::colorRamp2(
          breaks = c(min(e), 0, max(e)),
          colors = gradient_colors
        )
        hm <- ComplexHeatmap::Heatmap(e, name = "Expression",
                                      row_title = g,
                                      col = col_fun,
                                      cluster_rows = TRUE,
                                      cluster_columns = TRUE,
                                      column_names_gp = grid::gpar(fontsize=size.col),
                                      row_names_gp = grid::gpar(fontsize=size.row))
      } else {
        hm <- ComplexHeatmap::Heatmap(e, name = "Expression",
                                      row_title = g,
                                      col = gradient_colors,
                                      cluster_rows = TRUE,
                                      cluster_columns = TRUE,
                                      column_names_gp = grid::gpar(fontsize=size.col),
                                      row_names_gp = grid::gpar(fontsize=size.row))
      }
      if (is.null(out_dir)) {
        return(hm)
      } else {
        out_base <- glue("{out_dir}/{g}-heatmap")
        png(filename = glue("{out_base}.png"), width = 10,
            height = 5, units = "in", res = 300)
        ComplexHeatmap::draw(hm)
        dev.off()
      }
    }
  }
}

plot_heatmap.mod_v1=function(x, group_by, assay = "logcounts", out_dir = NULL,size.col=4,size.row=4,auto=F,segment=F){
  if (!is(x, "SpatialExperiment")) {
    stop("input is not a SpatialExperiment object")
  }
  if (!is.character(group_by)) {
    stop("`group_by` is not a character string")
  }
  if (!all(group_by %in% names(colData(x)))) {
    stop("not all `group_by` values are present in the object")
  }
  if (!is.character(assay)) {
    stop("`assay` is not a character string")
  }
  if (!assay %in% assayNames(x)) {
    stop("input SpatialExperiment object does not have a `", 
         assay, "` assay")
  }
  if (!is.null(out_dir)) {
    if (!dir.exists(out_dir)) {
      stop("output directory `", out_dir, "` does not exist")
    }
  }
  gradient_colors <- rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
  for (g in group_by) {
    e <- scuttle::summarizeAssayByGroup(x, ids = colData(x)[[g]], 
                                        assay.type = assay, statistics = "median")
    e <- SummarizedExperiment::assay(e, i = "median")
    # with guard restored
    if (sum(rowSums(e) == 0) == 0) {
      e <- scale(t(e))
    } else {
      e_ <- scuttle::summarizeAssayByGroup(x, ids = colData(x)[[g]],
                                           assay.type = assay, statistics = "mean")
      e_ <- SummarizedExperiment::assay(e_, i = "mean")
      e[rowSums(e) == 0, ] <- e_[rowSums(e) == 0, ]
      e <- scale(t(e))
    }
    e[e > 5] <- 5
    e[e < -5] <- -5
    if(auto){
      size.col=4
      cluster.size=as.numeric(unlist(strsplit(g,'_clust'))[2])
      if(cluster.size <= 20){
        size.row=4
        hm <- ComplexHeatmap::Heatmap(e, name = "Expression", 
                                      row_title = g, col = gradient_colors, cluster_rows = TRUE, 
                                      cluster_columns = TRUE,column_names_gp = grid::gpar(fontsize=size.col),row_names_gp = grid::gpar(fontsize=size.row))
        if (is.null(out_dir)) {
          return(hm)
        } else {
          out_base <- glue("{out_dir}/{g}-heatmap")
          png(filename = glue("{out_base}.png"), width = 10, 
              height = 5, units = "in", res = 300)
          ComplexHeatmap::draw(hm)
          dev.off()
        }
        
      }else{
        if(segment){
          clust.increments=c(seq(1,cluster.size,20),cluster.size)
          
          if(cluster.size <= 30){
            size.row=3
          }else{
            size.row=2 
          }
          
          for(i in seq(2,length(clust.increments))){
            e.tmp=e[seq(clust.increments[i-1],clust.increments[i]),]
            hm <- ComplexHeatmap::Heatmap(e.tmp, name = "Expression", 
                                          row_title = g, col = gradient_colors, cluster_rows = TRUE, 
                                          cluster_columns = TRUE,column_names_gp = grid::gpar(fontsize=size.col),row_names_gp = grid::gpar(fontsize=size.row))
            if (is.null(out_dir)) {
              return(hm)
            } else {
              out_base <- glue("{out_dir}/{g}-{i-1}-heatmap")
              png(filename = glue("{out_base}.png"), width = 10, 
                  height = 5, units = "in", res = 300)
              ComplexHeatmap::draw(hm)
              dev.off()
            }
          }
        }else{
          if(cluster.size <= 30){
            size.row=3
          }else{
            size.row=2 
          }
          hm <- ComplexHeatmap::Heatmap(e, name = "Expression",
                                        row_title = g, col = gradient_colors, cluster_rows = TRUE,
                                        cluster_columns = TRUE,column_names_gp = grid::gpar(fontsize=size.col),row_names_gp = grid::gpar(fontsize=size.row))
          if (is.null(out_dir)) {
            return(hm)
          } else {
            out_base <- glue("{out_dir}/{g}-heatmap")
            png(filename = glue("{out_base}.png"), width = 10,
                height = 5, units = "in", res = 300)
            ComplexHeatmap::draw(hm)
            dev.off()
          }
        }
      }
      
    }else{
      hm <- ComplexHeatmap::Heatmap(e, name = "Expression",
                                    row_title = g, col = gradient_colors, cluster_rows = TRUE,
                                    cluster_columns = TRUE,column_names_gp = grid::gpar(fontsize=size.col),row_names_gp = grid::gpar(fontsize=size.row))
      if (is.null(out_dir)) {
        return(hm)
      } else {
        out_base <- glue("{out_dir}/{g}-heatmap")
        png(filename = glue("{out_base}.png"), width = 10,
            height = 5, units = "in", res = 300)
        ComplexHeatmap::draw(hm)
        dev.off()
      }
    }
  }
}

create_object.mod=function (x, expression_cols = NULL, metadata_cols = NULL, skip_cols = NULL, 
                            clean_names = TRUE, transformation = NULL, out_dir = NULL,down_sample_idx=NULL,classifier.label=NULL,plot_diagnostic=F,min.cells=10,max.cells=1e5){
  if (is.character(x)) {
    # -------------------------------------------------------------------------
    # FIX: two-pass read to avoid the 2^31-1 byte limit on very large files.
    #
    # Pass 1 — header only (nrows=0): determine which columns survive
    #           skip_cols and whether expression_cols are resolvable, so that
    #           Pass 2 can load only the columns we actually need.
    # -------------------------------------------------------------------------
    file_path <- x

    hdr <- data.table::fread(file_path, nrows=0, data.table=FALSE)
    all_cols <- names(hdr)

    # Resolve skip_cols against the header
    if (!is.null(skip_cols)) {
      if (length(skip_cols) == 1) {
        skip_cols_resolved <- str_subset(all_cols, pattern = skip_cols)
      } else {
        skip_cols_resolved <- skip_cols
      }
      keep_cols <- setdiff(all_cols, skip_cols_resolved)
    } else {
      skip_cols_resolved <- character(0)
      keep_cols <- all_cols
    }

    # Pass 2 — shell-pipe approach.
    # fread(select=) still decompresses the full .gz into an internal buffer
    # before filtering, which hits the 2^31-1 byte limit on large files.
    # Using cmd= pipes zcat|cut through the shell so column selection and
    # decompression happen outside R entirely — only the trimmed stream enters.
    col_indices  <- which(all_cols %in% keep_cols)
    col_idx_str  <- paste(col_indices, collapse = ",")
    is_gz        <- grepl("\\.gz$", file_path, ignore.case = TRUE)
    decomp_cmd   <- if (is_gz) paste0("zcat '", file_path, "'") else paste0("cat '", file_path, "'")
    shell_cmd    <- paste0(decomp_cmd, " | cut -f", col_idx_str)

    message("Reading file via shell pipe (", length(keep_cols), " of ",
            length(all_cols), " columns): ", file_path)
    x <- data.table::fread(cmd = shell_cmd, stringsAsFactors = FALSE, data.table = FALSE)

    # skip_cols already applied via select= above; clear it so the later
    # skip_cols block inside the function body is a no-op
    skip_cols <- NULL

    print(head(x))
    if(!is.null(classifier.label)){
      hits=rowSums(do.call('cbind',lapply(classifier.label,function(y) str_detect(x$`Classifier Label`,y))))
      x=x[hits==1,]
    }
    if(!is.null(down_sample_idx) || dim(x)[1] > max.cells){
      if(is.null(down_sample_idx)){
        down_sample_idx=sample(dim(x)[1],max.cells)
      }
      x=x[down_sample_idx,]
    }
  }
  if (!is.data.frame(x)) {
    stop("input is not a file or a data frame")
  }
  if (!is.null(out_dir)) {
    if (dir.exists(out_dir)) {
      stop("output directory already exists")
    }
  }
  if (nrow(x) < min.cells) {
    print(x)
    print(nrow(x))
    stop("data frame has too few rows/cells")
  }
  if (ncol(x) < 10) {
    stop("data frame has too few columns")
  }
  message("number of input table rows: ", nrow(x))
  message("number of input table columns: ", ncol(x))
  message("")
  if (!is.null(skip_cols)) {
    if (length(skip_cols) == 1) {
      skip_cols <- str_subset(names(x), pattern = skip_cols)
    }
    x <- x[, setdiff(names(x), skip_cols)]
    message("skipped columns: ", toString(skip_cols), "\n")
  }

  if (is.null(expression_cols)) {
    print('detecting columns')
    expression_cols <- detect_exprs_cols(x)
  }
  exp0=expression_cols
  if (length(expression_cols) == 1) {
    expression_cols <- str_subset(names(x), pattern = expression_cols)
  }else{
    expression_cols=paste0(expression_cols,collapse = "|")
    expression_cols <- str_subset(names(x), pattern = expression_cols)
    # the above eliminates qupath expression columns
    if(length(expression_cols) == 0){
      expression_cols=exp0
    }
  }
  if (length(setdiff(expression_cols, names(x))) > 0) {
    stop("missing markers: ", toString(setdiff(expression_cols, 
                                               names(x))))
  }
  exprs <- try(x[, expression_cols])
  if('try-error' %in% class(exprs)){
    exprs <- x[, ..expression_cols]
  }
  if (is.null(metadata_cols)) {
    metadata_cols <- setdiff(names(x), expression_cols)
  }
  if (length(metadata_cols) == 1) {
    metadata_cols <- str_subset(names(x), pattern = metadata_cols)
  }
  if (length(setdiff(metadata_cols, names(x))) > 0) {
    stop("missing metadata columns: ", toString(setdiff(metadata_cols, 
                                                        names(x))))
  }
  
  x <- try(x[, metadata_cols])
  if('try-error' %in% class(x)){
    x <- x[, ..metadata_cols]
  }
  if (clean_names) {
    
    clean_col_names <- function(x) {
      
      # check if the input is valid
      if (!is.data.frame(x)) {
        stop("input is not a data frame")
      }
      
      # Platform-specific column name adjustments
      # Only apply destructive colon-stripping for MAV data (identified by "Cyc_" pattern)
      if(any(grepl("Cyc_", names(x)))) {
        names(x) <- str_remove(names(x), ":Cyc_.*")
        names(x) <- str_remove(names(x), ".*:")
      }
      
      # adjust Visiopharm column names (only if MP pattern found)
      if(any(grepl("^MP", names(x)))) {
        names(x) <- str_remove(names(x), "MP.*\\(")
      }
      
      # adjust HALO column names
      names(x) <- str_remove(names(x), " Cell Intensity")      
      # clean column names
      x <- janitor::clean_names(x, case = "none")
      
      # adjust common column names
      #
      # QuPath exports centroid columns as space-delimited names with an optional
      # unit suffix: "Centroid X µm" | "Centroid X px" | "Centroid X". After
      # janitor::clean_names(case="none") the space delimiters become underscores
      # AND the micro symbol µ is mangled to 'm' (not 'u'), producing any of:
      #   Centroid_X, Centroid_X_px, Centroid_X_pixels, Centroid_X_mm (← µm),
      #   Centroid_X_um, Centroid_X_µm.
      # Match all variants with a single regex.
      if (!"x" %in% names(x)) {
        names(x)[names(x) == "X"] <- "x"
        names(x)[names(x) == "X_X"] <- "x"
        names(x)[names(x) == "CtrX"] <- "x"
        qupath_x <- grep("^Centroid_X(_[^_]+)?$", names(x), value = TRUE)
        if (length(qupath_x) > 0) {
          names(x)[names(x) == qupath_x[1]] <- "x"
        }

        if (all(c("XMin", "XMax") %in% names(x))) {
          x$x <- (x$XMin + x$XMax) / 2
        }
        if (all(c("x_min", "x_max") %in% names(x))) {
          x$x <- (x$x_min + x$x_max) / 2
        }
      }
      if (!"y" %in% names(x)) {
        names(x)[names(x) == "Y"] <- "y"
        names(x)[names(x) == "Y_Y"] <- "y"
        names(x)[names(x) == "CtrY"] <- "y"
        qupath_y <- grep("^Centroid_Y(_[^_]+)?$", names(x), value = TRUE)
        if (length(qupath_y) > 0) {
          names(x)[names(x) == qupath_y[1]] <- "y"
        }

        if (all(c("YMin", "YMax") %in% names(x))) {
          x$y <- (x$YMin + x$YMax) / 2
        }
        if (all(c("y_min", "y_max") %in% names(x))) {
          x$y <- (x$y_min + x$y_max) / 2
        }
      }
      
      if (!"z" %in% names(x)) {
        names(x)[names(x) == "Z"] <- "z"
        names(x)[names(x) == "Z_Z"] <- "z"
      }
      if (!"reg" %in% names(x)) {
        names(x)[names(x) == "Region"] <- "reg"
        names(x)[names(x) == "region"] <- "reg"
      }
      if (!"tile_num" %in% names(x)) {
        names(x)[names(x) == "tile_nr"] <- "tile_num"
        names(x)[names(x) == "tile_number"] <- "tile_num"
      }
      if (!"size" %in% names(x)) {
        names(x)[names(x) == "Cell_Area_mm2"] <- "size"
        names(x)[names(x) == "cell_area_um2"] <- "size"
      }
      
      return(x)
    }
    
    exprs <- clean_col_names(exprs)
    x <- clean_col_names(x)
    print(names(x))
        
    if (!"cell_id" %in% names(x)) {
      names(x)[names(x) == "CellID"] <- "cell_id"
      names(x)[names(x) == "label"] <- "cell_id"
      names(x)[names(x) == "Object_Id"] <- "cell_id"
    }
    if (!"cell_id" %in% names(x)) {
      x$cell_id <- rownames(x)
    }
    if (is.numeric(x$cell_id)) {
      x$cell_id <- str_pad(as.character(x$cell_id), width = 7, 
                           pad = "0")
      x$cell_id <- str_c("C", x$cell_id)
    }
    x$cell_id <- make.names(x$cell_id, unique = TRUE)
  }
  if (!"cell_id" %in% names(x)) {
    stop("data frame must contain `cell_id` column")
  }
  if (!"x" %in% names(x)) {
    stop("data frame must contain `x` column")
  }
  if (!"y" %in% names(x)) {
    stop("data frame must contain `y` column")
  }
  rownames(x) <- x$cell_id
  exprs <- as.matrix(exprs)
  exprs <- exprs[, sort(colnames(exprs))]
  rownames(exprs) <- rownames(x)
  message("number of expression columns: ", ncol(exprs))
  message("number of metadata columns: ", ncol(x))
  message("")
  message("expression columns: ", toString(colnames(exprs)), 
          "\n")
  message("metadata columns: ", toString(colnames(x)), "\n")
  if (!is.null(out_dir)) {
    dir.create(out_dir)
  }
  if(plot_diagnostic){
    expression_dir=glue('{out_dir}/expression-diagnostic-plots')
    for(i in seq(dim(exprs)[2])){
      pdf(glue('{out_dir}/{colnames(exprs)[i]}.pdf'))
      plot(density(exprs[,i],main=colnames(exprs)[i]))
      abline(v=quantile(exprs[,i]))
      dev.off()
    }
  }
  
  na.idx=seq(dim(t(exprs))[1])[rowVars(t(exprs))==0]
  genes0var=paste0(colnames(exprs)[na.idx],collapse=',')
  print('na.idx: ')
  print(na.idx)
  
  vars <- rowVars(t(exprs), na.rm = TRUE)
  na.idx <- which(vars == 0 | is.na(vars))
  print('na.idx: ')
  print(na.idx)

  if(length(na.omit(na.idx)) == length(na.idx)){  
    if(length(na.idx) > 0){
      exprs=exprs[,-na.idx]
      message(glue('removed {length(na.idx)} genes due to 0 variance'))
      message(glue("removed: {genes0var}"))
    }
  }
  # Remove rows (cells) with any remaining NAs
  na_rows <- !complete.cases(exprs)
  if(any(na_rows)){
    message(glue('removed {sum(na_rows)} cells with NA values'))
    exprs <- exprs[!na_rows, ]
    x <- x[!na_rows, ]
  }    
  s <- SpatialExperiment::SpatialExperiment(assay = list(counts = t(exprs)), 
                                            colData = x, spatialCoordsNames = c("x", "y"))
  if (!is.null(transformation)) {
    s <- transform(s, method = transformation, out_dir = out_dir)
    s <- run_umap(s, n_threads = 4)
  }
  if (!is.null(out_dir)) {
    message("saving object")
    saveRDS(s, paste0(out_dir, "/spe.rds"))
  }
  return(s)
}

cluster.mod=function (x, method = c("leiden"), resolution = 1, n_neighbors = 50, 
                      out_dir = NULL,max_clust=100) 
{
  method <- match.arg(method)
  if (!is(x, "SpatialExperiment")) {
    stop("input is not a SpatialExperiment object")
  }
  if (!is.numeric(resolution)) {
    stop("`resolution` is not a number")
  }
  if (!is.numeric(n_neighbors)) {
    stop("`n_neighbors` is not a number")
  }
  if (!"exprs" %in% assayNames(x)) {
    stop("input SpatialExperiment object does not have a `exprs` assay")
  }
  if (!is.null(out_dir)) {
    if (!dir.exists(out_dir)) {
      stop("output directory `", out_dir, "` does not exist")
    }
    clusters_dir <- glue("{out_dir}/clusters")
    dir.create(clusters_dir)
  }
  exprs_mat <- assay(x, "exprs")
  if (method == "leiden") {
    g <- scran::buildSNNGraph(exprs_mat, transposed = FALSE, 
                              k = n_neighbors)
    n_clust_prev <- 0
    for (res_num in resolution) {
      message(glue("clustering using resolution of {res_num}"))
      set.seed(99)
      clusters <- igraph::cluster_leiden(g, objective_function = "modularity", 
                                         resolution_parameter = res_num, n_iterations = 10)
      clusters <- clusters$membership
      res_str <- format(as.numeric(res_num), nsmall = 1)
      res_str <- stringr::str_pad(res_str, width = 3, side = "left", 
                                  pad = "0")
      res_str <- stringr::str_c("res", res_str)
      n_clust <- length(unique(clusters))
      if(n_clust <= max_clust){
        clusters_label <- glue("cluster_{method}_{res_str}_clust{n_clust}")
        clusters <- as.character(clusters)
        clusters <- stringr::str_pad(clusters, width = 2, 
                                     side = "left", pad = "0")
        clusters <- stringr::str_c("C", clusters)
        
        if (n_clust > n_clust_prev) {
          x[[clusters_label]] <- factor(clusters)
          n_clust_prev <- n_clust
          if (!is.null(out_dir)) {
            cluster_summary <- janitor::tabyl(x[[clusters_label]])
            colnames(cluster_summary) <- c("cluster", "cells_num", 
                                           "cells_freq")
            readr::write_csv(cluster_summary, glue("{clusters_dir}/{clusters_label}-summary.csv"))
            plot_heatmap.mod_v1(x, group_by = clusters_label, 
                         out_dir = clusters_dir)
            plot_spatial(x, color_by = clusters_label, 
                         out_dir = clusters_dir)
            plot_dr(x, dr = "UMAP", color_by = clusters_label, 
                    out_dir = clusters_dir)
          }
        }
      }else{
        message(glue('{n_clust} clusters found, ending clustering'))
        saveRDS(glue('clusters:{res_num}'),glue('{out_dir}/lock.rds'))
        break
      }
    }
  }
  if (!is.null(out_dir)) {
    message("saving object")
    saveRDS(x, paste0(out_dir, "/spe.rds"))
  }
  return(x)
}

# ---------------------------------------------------------------------------
# export_anndata.mod(spe, out_dir, filename, X_name)
#
# Writes a SingleCellExperiment / SpatialExperiment (as produced by
# create_object.mod / cluster.mod) out as an AnnData .h5ad file, alongside
# spe.rds, via the Bioconductor 'zellkonverter' package.
#
# SpatialExperiment keeps cell coordinates in a dedicated `spatialCoords`
# slot that zellkonverter does not know about, so they are copied into a
# `reducedDim` named "spatial" first — this lands them at `adata.obsm
# ['spatial']`, the convention scanpy/squidpy expect.
# ---------------------------------------------------------------------------
export_anndata.mod <- function(spe, out_dir, filename = "spe.h5ad", X_name = "counts") {
  if (!is(spe, "SingleCellExperiment")) {
    stop("input is not a SingleCellExperiment/SpatialExperiment object")
  }
  if (!requireNamespace("zellkonverter", quietly = TRUE)) {
    stop("the 'zellkonverter' package is required for AnnData export. Install with:\n",
         "  BiocManager::install(\"zellkonverter\")")
  }
  if (!dir.exists(out_dir)) {
    stop("output directory `", out_dir, "` does not exist")
  }
  if (is(spe, "SpatialExperiment") && !"spatial" %in% SingleCellExperiment::reducedDimNames(spe)) {
    SingleCellExperiment::reducedDim(spe, "spatial") <- SpatialExperiment::spatialCoords(spe)
  }
  x_name <- if (X_name %in% SummarizedExperiment::assayNames(spe)) {
    X_name
  } else {
    SummarizedExperiment::assayNames(spe)[1]
  }
  out_path <- file.path(out_dir, filename)
  message("exporting AnnData (X = '", x_name, "' assay): ", out_path)
  zellkonverter::writeH5AD(spe, out_path, X_name = x_name)
  invisible(out_path)
}

# ---------------------------------------------------------------------------
# .detect_segmentation_format(x)
#
# Inspects the column headers of a segmentation CSV (path or already-loaded
# data.frame) and returns one of:
#   "qupath" — QuPath-style table with space-delimited "Centroid X [µm|px]"
#              columns and/or colon-separated per-compartment mean columns
#              (e.g. "CD3: Cell: Mean").
#   "halo"   — HALO-style table with "<marker> Cell/Nucleus/Cytoplasm Intensity"
#              columns and/or HALO metadata columns (Classification, Completeness).
#   "mesmer" — Mesmer-style table with a (label, y, x[, size, *_min/*_max])
#              header block followed by bare marker columns.
#   "codex"  — Per-cell quantitative export (CODEX / PhenoCycler / Visiopharm /
#              CellProfiler-style) with explicit stat-suffix column names:
#              <marker>_Mean_intensity, <marker>_P25_intensity, etc. alongside
#              morphometric shape descriptors (Compactness, Elongation, …).
#              Only <marker>_Mean_intensity columns are used as expression; all
#              other stat suffixes and shape columns are pushed to metadata.
#
# QuPath is checked first because QuPath exports may contain substrings that
# would otherwise trip HALO heuristics (e.g. "Cell" in "CD3: Cell: Mean").
# CODEX is checked before Mesmer because CODEX tables may contain x/y columns
# that would otherwise trigger the Mesmer heuristic.
#
# If no fingerprint matches, emits a warning and falls back to "halo" (the
# historical default). Users can bypass detection by passing an explicit
# skip_cols regex to RunPhenomenalist().
# ---------------------------------------------------------------------------
.detect_segmentation_format <- function(x) {
  if (is.character(x) && length(x) == 1) {
    hdr <- data.table::fread(x, nrows = 0)
    cols <- names(hdr)
  } else if (is.data.frame(x)) {
    cols <- names(x)
  } else {
    cols <- as.character(x)
  }
  # QuPath: space-delimited "Centroid X [µm|px|<unit>]" or colon-separated means.
  qupath_centroid <- any(grepl("^Centroid\\s+[XY](\\s|$)", cols, perl = TRUE))
  qupath_mean     <- any(grepl(":\\s*(Cell|Nucleus|Cytoplasm|Membrane)\\s*:\\s*Mean",
                               cols, perl = TRUE))
  if (qupath_centroid || qupath_mean) return("qupath")
  # HALO: "<marker> Cell Intensity", etc., and/or HALO metadata.
  halo_hits <- grepl("(?i)(Cell|Nucleus|Cytoplasm|Membrane)[ .]Intensity|Classification|Completeness",
                     cols, perl = TRUE)
  if (any(halo_hits)) return("halo")
  # CODEX / quantitative stat-suffix exports: tables with explicit _Mean_intensity
  # columns alongside other per-pixel stat suffixes (_P25_intensity, _Max_intensity,
  # etc.) and/or morphometric shape descriptors (Compactness, Elongation, …).
  # Both conditions must be present to avoid false-positives on custom tables that
  # happen to have one but not the other.
  codex_intensity <- any(grepl("_Mean_intensity", cols, fixed = TRUE))
  codex_stats     <- any(grepl("_(Max|Min|Median|P\\d+|Std_dev)_intensity", cols, perl = TRUE))
  if (codex_intensity && codex_stats) return("codex")
  # Mesmer: bare label/y/x header block.
  mesmer_base <- c("label", "y", "x")
  if (all(mesmer_base %in% cols)) return("mesmer")
  warning("Could not auto-detect segmentation format from column names; ",
          "assuming HALO. Pass skip_cols=<regex> to override.")
  "halo"
}

phenomenalist.preprocess=function(x,failed.markers=NULL,nuclear.markers=NULL,else.cytoplasm=F){
  if(is.null(failed.markers) & is.null(nuclear.markers)){
    return(NULL)
  }else{
    HALO <- identical(.detect_segmentation_format(x), "halo")
    if(HALO){
      if(is.null(nuclear.markers)){
        # FIX: nrows=0 — only the header is needed to resolve column names.
        # The full file can exceed R's 2^31-1 byte string limit on large datasets.
        x=data.table::fread(x, nrows=0)
        
        names_original=names(x)
        names(x)=tolower(names(x))
        failed.markers=tolower(failed.markers)
        
        cell.intensity=names(x)[grep(tolower('Cell.Intensity'),names(x))]
        cell.intensity=cell.intensity[-unlist(lapply(failed.markers,function(marker) grep(marker,cell.intensity)))]
        expression.columns=cell.intensity
        expression.columns = names_original[match(expression.columns,names(x))]
        print(expression.columns)
      }else{
        
        nuclear_markers=nuclear.markers
        
        # FIX: nrows=0 — only the header is needed to resolve column names.
        x=data.table::fread(x, nrows=0)
        names_original=names(x)
        names(x)=tolower(names(x))
        nuclear_markers=tolower(nuclear_markers)
        failed.markers=tolower(failed.markers)
        qupath=F

        expression.columns=names(x)[unlist(lapply(nuclear_markers,function(m) grep(m,names(x))))]
        
        nuclear.intensity.cols=expression.columns[grep(tolower('Nucleus Intensity'),expression.columns)]
        if(length(nuclear.intensity.cols) == 0){
          nuclear.intensity.cols=expression.columns[grep(tolower('Nucleus.Intensity'),expression.columns)]
        }
        if(length(grep('mean',expression.columns)) > 0){
          qupath=T
          expression.columns=expression.columns[grep('mean',expression.columns)]
          if(length(nuclear.intensity.cols) == 0){
            nuclear.intensity.cols=expression.columns[grep(tolower('Nucleus'),expression.columns)]
          }
        }
        cell.intensity=names(x)[grep(tolower('Cell Intensity'),names(x))]
        if(length(cell.intensity) == 0){
          cell.intensity=names(x)[grep(tolower('Cell.Intensity'),names(x))]
        }
        if(qupath){
          cell.intensity=names(x)[grep(tolower('Cell'),names(x))]
          cell.intensity=cell.intensity[grep('mean',cell.intensity)]
        }
        cell.intensity=cell.intensity[-unlist(lapply(c(nuclear_markers,failed.markers),function(marker) grep(marker,cell.intensity)))]
        job2.selected_expression_cols=c(nuclear.intensity.cols,cell.intensity)
        if(else.cytoplasm){
          cytoplasm.intensity=names(x)[grep(tolower('Cytoplasm Intensity'),names(x))]
          cytoplasm.intensity.nuc_failed=lapply(c(nuclear_markers,failed.markers),function(marker) grep(marker,cytoplasm.intensity))
          if(length(unlist(cytoplasm.intensity.nuc_failed)) > 0){
            cytoplasm.intensity=cytoplasm.intensity[-unlist(cytoplasm.intensity.nuc_failed)]
          }
          
          job3.selected_expression_cols=c(nuclear.intensity.cols,cytoplasm.intensity)
          expression.columns=job3.selected_expression_cols
        }else{
          expression.columns=job2.selected_expression_cols
          print(expression.columns)
        }
        expression.columns = names_original[match(expression.columns,names(x))]
      }
      
    }else{
      # FIX: nrows=0 — only column names are needed here; ncol() replaces dim(x)[2].
      x=data.table::fread(x, nrows=0)
      if('size' %in% colnames(x)){
        init.idx=4
      }else{
        init.idx=3
      }
      init.idx=max(na.omit(match(c('label','y','x','size','y_min','x_min','y_max','x_max'),colnames(x))))+1

      cell.intensity=colnames(x)[seq(init.idx,ncol(x))]
      if(!is.null(failed.markers)){
        cell.intensity=cell.intensity[-unlist(lapply(failed.markers,function(marker) grep(marker,cell.intensity)))]
      }
      
      expression.columns=cell.intensity
    }
    
    return(expression.columns)
  }
}

prepare_mask_inputs=function(spe,mask.only,out_dir,res,failed.markers,label=NULL){
  library(glue)
  message(glue('generating mask-generation inputs for resolution {res} clusters'))
  
  spatial_obj=data.frame(spatialCoords(spe))
  res.hits=str_detect(names(colData(spe)),glue('res{res}'))
  hits.found=ifelse(sum(res.hits)>0,T,F)
  while(hits.found==F){
    res=res-1
    res.hits=str_detect(names(colData(spe)),glue('res{res}'))
    hits.found=ifelse(sum(res.hits)>0,T,F)
  }
  clust.tmp=colData(spe)[[names(colData(spe))[res.hits]]]
  spatial_obj$cluster=clust.tmp
  
  mask_inputs_dir = paste0(out_dir, "/mask-inputs/")
  dir.create(mask_inputs_dir, showWarnings = FALSE)
  if(!is.null(label)){
    out_name=glue('{mask_inputs_dir}/{label}-spatial-anno-res{res}.csv')
  }else{
    out_name=glue('{mask_inputs_dir}/spatial-anno-res{res}.csv')
  }
  write.csv(spatial_obj,out_name)
  
  message('generating marker metadata')
  all.markers=row.names(spe)
  all.markers=do.call('rbind',strsplit(all.markers,'*_Cytoplasm|_Nucleus'))[,1]
  if(!is.null(mask.only)){
    working.markers=c(all.markers,mask.only)
  }else{
    working.markers=all.markers
  }
  if(!is.null(label)){
    out_name=glue('{mask_inputs_dir}/{label}-marker-metadata.csv')
  }else{
    out_name=glue('{mask_inputs_dir}/marker-metadata.csv')
  }
  write.csv(working.markers,out_name)
}

#source('/Users/ee699/TRIC/phenomenalist/R/plot-scatter.R')
plot_dr.mod=function (x, dr, color_by, assay = "logcounts", smooth = FALSE, 
                      range = c(0.01, 0.99), out_dir = NULL,h=NULL,w=NULL,pdf=F) 
{
  if (!is(x, "SpatialExperiment")) {
    stop("input is not a SpatialExperiment object")
  }
  if (!is.character(dr)) {
    stop("`dr` is not a character string")
  }
  if (!dr %in% reducedDimNames(x)) {
    stop("SpatialExperiment object does not have a dimensionality reduction `", 
         dr, "`")
  }
  if (!is.null(out_dir)) {
    if (!dir.exists(out_dir)) {
      stop("output directory `", out_dir, "` does not exist")
    }
  }
  if (all(color_by %in% names(colData(x)))) {
    vals <- as.data.frame(colData(x), stringsAsFactors = FALSE)[color_by]
  }
  else if (all(color_by %in% rownames(x))) {
    if (!is.character(assay)) {
      stop("`assay` is not a character string")
    }
    if (!assay %in% assayNames(x)) {
      stop("input SpatialExperiment object does not have a `", 
           assay, "` assay")
    }
    vals <- assay(x, i = assay)
    vals <- as.data.frame(t(vals), stringsAsFactors = FALSE)[color_by]
  }
  else {
    stop("not all `color_by` values are present in the object")
  }
  coords <- reducedDim(x, type = dr)
  coords <- coords[, 1:2]
  colnames(coords) <- paste0(dr, 1:2)
  coords <- as.data.frame(coords)
  coords <- cbind(coords, vals[rownames(coords), , drop = FALSE])
  for (val in colnames(vals)) {
    p <- plot_scatter(data = coords, x = names(coords)[1], 
                      y = names(coords)[2], color_by = val, smooth = smooth, 
                      range = range, title = val)
    if (is.null(out_dir)) {
      return(p)
    } else {
      out_base <- glue("{out_dir}/{val}-{dr}")
      if (smooth) {
        out_base <- glue("{out_base}-smooth")
      }
      message(glue("generating {dr} plot for {val}"))
      ggsave(filename = ifelse(pdf,glue("{out_base}.pdf"),glue("{out_base}.png")), plot = p, 
             width = ifelse(is.null(w),8,w), height = ifelse(is.null(h),5,h))
    }
  }
  return(p)
}

plot_spatial.mod=function (x, color_by, assay = "logcounts", smooth = FALSE, range = c(0.01, 
                                                                                       0.99), out_dir = NULL,h=NULL,w=NULL,pdf=F,colors=NULL) 
{
  if (!is(x, "SpatialExperiment")) {
    stop("input is not a SpatialExperiment object")
  }
  if (!is.null(out_dir)) {
    if (!dir.exists(out_dir)) {
      stop("output directory `", out_dir, "` does not exist")
    }
  }
  if (all(color_by %in% names(colData(x)))) {
    vals <- as.data.frame(colData(x)[color_by], stringsAsFactors = FALSE)
  }
  else if (all(color_by %in% rownames(x))) {
    if (!is.character(assay)) {
      stop("`assay` is not a character string")
    }
    if (!assay %in% assayNames(x)) {
      stop("input SpatialExperiment object does not have a `", 
           assay, "` assay")
    }
    vals <- assay(x, i = assay)
    vals <- t(vals)[, color_by, drop = FALSE]
    vals <- as.data.frame(vals, stringsAsFactors = FALSE)
  }
  else {
    stop("not all `color_by` values are present in the object")
  }
  if (!is.null(out_dir)) {
    dir.create(out_dir, showWarnings = FALSE)
  }
  coords <- spatialCoords(x)
  coords <- as.data.frame(coords)
  coords$y <- coords$y * -1
  coords <- cbind(coords, vals[rownames(coords), , drop = FALSE])
  ratio <- max(coords$y)/max(coords$x)
  for (val in colnames(vals)) {
    if(is.null(colors)){
      p <- plot_scatter(data = coords, x = "x", y = "y", color_by = val, 
                        smooth = smooth, range = range, title = val, aspect_ratio = ratio)
    }else{
      p <- plot_scatter(data = coords, x = "x", y = "y", color_by = val, 
                        smooth = smooth, range = range, title = val, aspect_ratio = ratio)+scale_color_manual(values=colors)
    }
    if (is.null(out_dir)) {
      return(p)
    } else {
      out_base <- glue("{out_dir}/{val}-spatial")
      if (smooth) {
        out_base <- glue("{out_base}-smooth")
      }
      message(glue("generating spatial plot for {val}"))
      ggsave(filename = ifelse(pdf,glue("{out_base}.pdf"),glue("{out_base}.png")), plot = p, 
             width = ifelse(is.null(w),8,w), height = ifelse(is.null(h),5,h))
    }
  }
}

assign_celltype_with_template=function(obj,phenotyping_template,cluster,mclust=T){
  
  g=cluster
  multiply_vectors <- function(boolean_gates) {
    Reduce(`*`, boolean_gates)
  }
  assay='logcounts'
  e <- scuttle::summarizeAssayByGroup(obj, ids = colData(obj)[[g]], assay.type = assay, statistics = "median")
  e <- SummarizedExperiment::assay(e, i = "median")
  e <- scale(t(e))
  
  if(mclust){
    require(mclust)
    scaled_ = assay(obj,'exprs')
    mclust.cols=apply(e,2,function(x){Mclust(x,G=2)})
    # find higher mean group
    group_means=lapply(mclust.cols,function(x) x$parameters$mean)
  }
  
  markers=lapply(seq(nrow(phenotyping_template)),function(x) unlist(strsplit(phenotyping_template$MARKERS[x],'[,]')))
  
  message('generating marker gates')
  marker_gates=do.call('cbind',lapply(markers,function(x){
    tmp=x
    boolean_gates=lapply(tmp,function(y){
      
      # pull last character (sign):
      sign.tmp=substr(y,nchar(y),nchar(y))
      # strip sign, leave only marker:
      marker.tmp=unlist(strsplit(y,'[+]|[-]'))
      # ensure marker has no numbers or letters after: (e.g. CD3 --> CD3 not CD31)
      pattern <- glue("{marker.tmp}(?![a-zA-Z0-9])")
      if(sign.tmp == '+'){
        if(!mclust){
          return(e[,grep(pattern,colnames(e),perl=T)] >0)
        }else{
          marker.idx=grep(pattern,colnames(e),perl=T)
          return(mclust.cols[[marker.idx]]$classification == names(group_means[[marker.idx]])[group_means[[marker.idx]]==max(group_means[[marker.idx]])])
        }
      }else{
        if(!mclust){
          return(e[,grep(pattern,colnames(e),perl=T)] <0)
        }else{
          marker.idx=grep(pattern,colnames(e),perl=T)
          return(mclust.cols[[marker.idx]]$classification == names(group_means[[marker.idx]])[group_means[[marker.idx]]==min(group_means[[marker.idx]])])
        }
      }
    })
    
    # entry-wise multiplication of boolean decision vectors:
    decision_vector=multiply_vectors(boolean_gates = boolean_gates)
    return(decision_vector)
  }))
  colnames(marker_gates)=phenotyping_template$CELLTYPE
  message('assigning celltypes')
  print(head(marker_gates))
  assign_celltype=sapply(seq(nrow(marker_gates)),function(x){
    tmp=colnames(marker_gates)[marker_gates[x,]==1]
    if(length(tmp) == 0){
      return('Unannotated')
    }else{
      return(paste0(tmp,collapse = ','))
    }
  })
  mapping=cbind(row.names(marker_gates),assign_celltype)
  print(mapping)
  assigned_celltype=as.character(obj[[g]])
  for(i in seq_along(mapping[,1])){
    assigned_celltype[assigned_celltype==mapping[i,1]]=mapping[i,2]
  }
  obj[[glue('{g}_annotations_template')]]=assigned_celltype
  
  return(obj)
}
