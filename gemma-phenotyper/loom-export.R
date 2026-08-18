# Minimal loom writer (no LoomExperiment/rhdf5 dependency — just hdf5r,
# since loom is a plain HDF5 layout: a /matrix dataset plus /row_attrs and
# /col_attrs groups). See https://linnarssonlab.org/loompy/format/index.html
#
# Row/col attrs are named to match anndata.read_loom()'s defaults (var_names
# = "Gene", obs_names = "CellID") so files opened in scanpy/squidpy work with
# no extra arguments, even though "Gene" is a slight misnomer for mIF markers.

# matrix: numeric, rows = markers/genes, cols = cells (loom's required orientation).
# row_attrs / col_attrs: named lists of atomic vectors matching matrix dims.
#   row_attrs must include "Gene"; col_attrs must include "CellID".
write_loom <- function(matrix, row_attrs, col_attrs, file_path) {
  require(hdf5r)

  stopifnot(is.matrix(matrix))
  stopifnot("Gene" %in% names(row_attrs))
  stopifnot("CellID" %in% names(col_attrs))

  clean_attrs <- function(attrs, expected_len) {
    Filter(Negate(is.null), lapply(attrs, function(x) {
      if (is.factor(x)) x <- as.character(x)
      if (!is.atomic(x)) {
        message("write_loom: dropping non-atomic attribute")
        return(NULL)
      }
      if (length(x) != expected_len) {
        message("write_loom: dropping attribute with mismatched length")
        return(NULL)
      }
      x
    }))
  }
  row_attrs <- clean_attrs(row_attrs, nrow(matrix))
  col_attrs <- clean_attrs(col_attrs, ncol(matrix))

  if (file.exists(file_path)) file.remove(file_path)
  h5 <- H5File$new(file_path, mode = "w")
  on.exit(h5$close_all(), add = TRUE)

  # v3.0.0 spec stores global attributes as datasets under /attrs (not root
  # HDF5 attributes) — loompy's validator requires this group to exist.
  attrs_grp <- h5$create_group("attrs")
  attrs_grp[["LOOM_SPEC_VERSION"]] <- "3.0.0"

  # R matrices are column-major; loom/HDF5 readers (h5py, loompy) expect
  # row-major (genes x cells) shape/order, so the on-disk dataset must be
  # the transpose of what R's dim(matrix) reports, or shape AND values both
  # come out transposed when read back in Python.
  h5$create_dataset("matrix", t(matrix), chunk_dims = "auto", gzip_level = 4)

  row_grp <- h5$create_group("row_attrs")
  for (nm in names(row_attrs)) row_grp[[nm]] <- row_attrs[[nm]]

  col_grp <- h5$create_group("col_attrs")
  for (nm in names(col_attrs)) col_grp[[nm]] <- col_attrs[[nm]]

  invisible(file_path)
}
