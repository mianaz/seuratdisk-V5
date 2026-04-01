#' @include zzz.R
#' @include FormatRegistry.R
#' @include Connect.R
#' @include TestObject.R
#' @include Transpose.R
#' @include PadMatrix.R
#' @importFrom utils setTxtProgressBar
#' @importFrom hdf5r H5File h5attr H5S
#' @importFrom tools file_path_sans_ext
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Convert single-cell datasets between formats
#'
#' Universal converter between single-cell file formats and object types.
#' Supports arbitrary source/destination pairs by routing through Seurat as
#' a hub format. Direct HDF5-level paths (h5ad <-> h5seurat) are used when
#' available for memory efficiency.
#'
#' @param source Source dataset: a Seurat object, SingleCellExperiment, loom
#'   connection, filename path, or H5File connection
#' @param dest Name/path of destination file or format. Supported formats:
#'   h5seurat, h5ad, loom, rds. Also accepts \code{"sce"} to return
#'   a SingleCellExperiment object (in-memory, no file created).
#' @param assay For h5Seurat -> other formats: name of assay to convert.
#'   For other formats -> h5Seurat: name to assign to the assay.
#'   Default is "RNA".
#' @param overwrite Logical; if \code{TRUE}, overwrite an existing destination file.
#'   Default is \code{FALSE}.
#' @param verbose Logical; if \code{TRUE} (default), show progress updates
#' @param standardize Logical; if \code{TRUE}, convert Seurat-style metadata column names
#'   to scanpy/AnnData conventions when converting to h5ad format. For example,
#'   \code{nCount_RNA} becomes \code{n_counts}, \code{nFeature_RNA} becomes \code{n_genes}.
#'   Only applicable for conversions to h5ad format. Default is \code{FALSE}.
#' @param ... Arguments passed to specific conversion methods
#'
#' @return For file destinations, invisibly returns the destination filename.
#'   For \code{dest = "sce"}, returns a SingleCellExperiment object.
#'
#' @inheritSection H5ADToH5Seurat AnnData/H5AD to h5Seurat
#' @inheritSection H5SeuratToH5AD h5Seurat to AnnData/H5AD
#'
#' @details
#' \strong{Supported Formats:}
#' \itemize{
#'   \item \strong{R objects}: Seurat, SingleCellExperiment (requires
#'     \pkg{SingleCellExperiment}), loom (R6 connection)
#'   \item \strong{File formats}: h5seurat, h5ad, loom, rds
#' }
#'
#' Any source format can be converted to any destination format. Conversions
#' without a direct path go through Seurat as a universal hub:
#' \code{Source -> Seurat -> Destination}.
#'
#' \strong{Direct Paths} (memory-efficient, no full dataset loading):
#' \itemize{
#'   \item \code{h5ad <-> h5seurat}: Direct HDF5-level copy
#' }
#'
#' \strong{Key Features:}
#' \itemize{
#'   \item Preserves expression matrices, metadata, and dimensional reductions
#'   \item For Visium/spatial data: reconstructs images with scale factors
#'   \item Handles multiple data layers (V5 compatibility)
#' }
#'
#' @seealso
#' \code{\link{SeuratToH5AD}} for direct Seurat to h5ad convenience function
#' \code{\link{SaveH5Seurat}} to save Seurat objects
#' \code{\link{LoadH5Seurat}} to load h5Seurat files
#' \code{\link{LoadH5AD}} to directly load h5ad files
#' \code{\link{Connect}} to establish file connections
#'
#' @examples
#' \dontrun{
#' library(srtdisk)
#' library(Seurat)
#'
#' # --- Any format to any format ---
#' Convert("data.h5ad", dest = "data.h5seurat")     # h5ad -> h5seurat
#' Convert("data.h5ad", dest = "data.rds")           # h5ad -> RDS
#' Convert("data.h5ad", dest = "data.loom")           # h5ad -> loom
#' Convert("data.loom", dest = "data.h5seurat")       # loom -> h5seurat
#' Convert("data.rds",  dest = "data.h5ad")           # RDS -> h5ad
#'
#' # --- From R objects ---
#' Convert(seurat_obj, dest = "output.h5ad")          # Seurat -> h5ad
#' Convert(seurat_obj, dest = "output.loom")           # Seurat -> loom
#' Convert(seurat_obj, dest = "output.rds")            # Seurat -> RDS
#' Convert(sce_obj, dest = "output.h5ad")              # SCE -> h5ad
#' sce <- Convert(seurat_obj, dest = "sce")            # Seurat -> SCE
#' }
#'
#' @name Convert
#' @rdname Convert
#'
#' @export
#'
Convert <- function(source, dest, assay, overwrite = FALSE, verbose = TRUE, standardize = FALSE, ...) {
  if (!missing(x = dest) && !is.character(x = dest)) {
    stop("'dest' must be a filename or type", call. = FALSE)
  }
  UseMethod(generic = 'Convert', object = source)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom hdf5r H5File
#'
#' @rdname Convert
#' @method Convert character
#' @export
#'
Convert.character <- function(
  source,
  dest,
  assay,
  overwrite = FALSE,
  verbose = TRUE,
  standardize = FALSE,
  ...
) {
  if (!file.exists(source)) {
    stop("Source file not found: ", source, call. = FALSE)
  }
  stype <- FileType(file = source)

  # RDS: read the object and dispatch on its class
  if (stype == 'rds') {
    obj <- readRDS(file = source)
    if (!inherits(x = obj, what = 'Seurat')) {
      obj <- tryCatch(
        expr = as.Seurat(x = obj, verbose = verbose),
        error = function(e) {
          stop("RDS file does not contain a Seurat-coercible object: ",
               conditionMessage(e), call. = FALSE)
        }
      )
    }
    if (missing(x = assay)) assay <- DefaultAssay(object = obj)
    return(Convert(
      source = obj, dest = dest, assay = assay,
      overwrite = overwrite, verbose = verbose, standardize = standardize, ...
    ))
  }

  # Loom: open as loom R6 object and dispatch to Convert.loom
  if (stype == 'loom') {
    lfile <- Connect(filename = source, type = 'loom', mode = 'r')
    on.exit(expr = lfile$close_all(), add = TRUE)
    if (missing(x = assay)) {
      assay <- tryCatch(
        expr = DefaultAssay(object = lfile),
        error = function(...) 'RNA'
      )
    }
    return(Convert(
      source = lfile, dest = dest, assay = assay,
      overwrite = overwrite, verbose = verbose, standardize = standardize, ...
    ))
  }

  # HDF5-based formats: check for direct path first, otherwise use hub
  dtype <- FileType(file = dest)
  if (missing(x = assay)) {
    # Try to read default assay from the source file
    hfile_tmp <- tryCatch(
      expr = Connect(filename = source, force = TRUE),
      error = function(...) NULL
    )
    if (!is.null(x = hfile_tmp)) {
      assay <- tryCatch(
        expr = DefaultAssay(object = hfile_tmp),
        error = function(...) {
          warning("'assay' not set, setting to 'RNA'", call. = FALSE, immediate. = TRUE)
          "RNA"
        }
      )
      hfile_tmp$close_all()
    } else {
      assay <- 'RNA'
    }
  }

  # Check for a registered direct HDF5 path (memory-efficient)
  direct_fn <- GetDirectPath(stype = stype, dtype = dtype)
  if (!is.null(x = direct_fn)) {
    if (tolower(x = dest) == dtype) {
      dest <- paste(file_path_sans_ext(x = source), dtype, sep = '.')
    }
    hfile <- Connect(filename = source, force = TRUE)
    on.exit(expr = hfile$close_all())
    dfile <- direct_fn(
      source = hfile, dest = dest, assay = assay,
      overwrite = overwrite, verbose = verbose,
      standardize = standardize, ...
    )
    if (is.character(x = dfile)) {
      return(invisible(x = dfile))
    }
    dfile$close_all()
    return(invisible(x = dfile$filename))
  }

  # No direct path: use hub conversion (loader handles its own file I/O)
  if (tolower(x = dest) == dtype) {
    dest <- paste(file_path_sans_ext(x = source), dtype, sep = '.')
  }
  return(HubConvert(
    source_file = source,
    dest_file = dest,
    stype = stype,
    dtype = dtype,
    assay = assay,
    overwrite = overwrite,
    verbose = verbose,
    ...
  ))
}

#' @rdname Convert
#' @method Convert H5File
#' @export
#'
Convert.H5File <- function(
  source,
  dest = 'h5seurat',
  assay = 'RNA',
  overwrite = FALSE,
  verbose = TRUE,
  ...
) {
  stype <- FileType(file = source$filename)
  dtype <- FileType(file = dest)
  if (tolower(x = dest) == dtype) {
    dest <- paste(file_path_sans_ext(x = source$filename), dtype, sep = '.')
  }

  # 1. Try registered direct HDF5-level path (memory-efficient)
  direct_fn <- GetDirectPath(stype = stype, dtype = dtype)
  if (!is.null(x = direct_fn)) {
    if (verbose) {
      message("Using direct ", toupper(x = stype), " -> ", toupper(x = dtype), " conversion")
    }
    dfile <- direct_fn(
      source = source, dest = dest, assay = assay,
      overwrite = overwrite, verbose = verbose, ...
    )
    return(dfile)
  }

  # 2. Hub conversion: load source -> Seurat -> save dest
  return(HubConvert(
    source_file = source$filename,
    dest_file = dest,
    stype = stype,
    dtype = dtype,
    assay = assay,
    overwrite = overwrite,
    verbose = verbose,
    ...
  ))
}

#' @rdname Convert
#' @method Convert h5Seurat
#' @export
#'
Convert.h5Seurat <- function(
  source,
  dest = 'h5ad',
  assay = DefaultAssay(object = source),
  overwrite = FALSE,
  verbose = TRUE,
  standardize = FALSE,
  ...
) {
  dtype <- FileType(file = dest)
  if (tolower(x = dest) == dtype) {
    dest <- paste(file_path_sans_ext(x = source$filename), dtype, sep = '.')
  }

  # Try direct path first (h5seurat -> h5ad is registered)
  direct_fn <- GetDirectPath(stype = 'h5seurat', dtype = dtype)
  if (!is.null(x = direct_fn)) {
    return(direct_fn(
      source = source, dest = dest, assay = assay,
      overwrite = overwrite, verbose = verbose,
      standardize = standardize, ...
    ))
  }

  # Hub: load h5seurat -> Seurat -> save dest
  return(HubConvert(
    source_file = source$filename,
    dest_file = dest,
    stype = 'h5seurat',
    dtype = dtype,
    assay = assay,
    overwrite = overwrite,
    verbose = verbose,
    ...
  ))
}

#' @rdname Convert
#' @method Convert Seurat
#' @export
Convert.Seurat <- function(
  source,
  dest,
  assay = DefaultAssay(object = source),
  overwrite = FALSE,
  verbose = TRUE,
  standardize = FALSE,
  ...
) {
  if (missing(x = dest)) {
    stop("'dest' must be provided for Seurat object conversion", call. = FALSE)
  }
  type <- FileType(file = dest)

  # Handle in-memory SCE destination (no file created)
  if (type == 'sce') {
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
      stop("Package 'SingleCellExperiment' is required to convert to SCE", call. = FALSE)
    }
    if (verbose) message("Converting Seurat to SingleCellExperiment...")
    return(Seurat::as.SingleCellExperiment(x = source, assay = assay))
  }

  if (tolower(x = dest) == type) {
    dest <- paste(Project(object = source), type, sep = '.')
  }

  # Use registered saver if available
  saver <- GetSaver(ext = type)
  if (!is.null(x = saver)) {
    saver(object = source, filename = dest, overwrite = overwrite, verbose = verbose, ...)
    if (verbose) message("Successfully created: ", dest)
    return(invisible(x = dest))
  }

  stop("Cannot convert Seurat objects to '", type, "' format", call. = FALSE)
}

#' @rdname Convert
#' @method Convert loom
#' @export
#'
Convert.loom <- function(
  source,
  dest,
  assay = 'RNA',
  overwrite = FALSE,
  verbose = TRUE,
  standardize = FALSE,
  ...
) {
  if (missing(x = dest)) {
    stop("'dest' must be provided for loom conversion", call. = FALSE)
  }
  if (verbose) message("Loading loom file as Seurat object...")
  seurat_obj <- as.Seurat(x = source, verbose = verbose)
  return(Convert(
    source = seurat_obj, dest = dest, assay = assay,
    overwrite = overwrite, verbose = verbose, standardize = standardize, ...
  ))
}

#' @rdname Convert
#' @method Convert SingleCellExperiment
#' @export
#'
Convert.SingleCellExperiment <- function(
  source,
  dest,
  assay = NULL,
  overwrite = FALSE,
  verbose = TRUE,
  standardize = FALSE,
  ...
) {
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("Package 'SingleCellExperiment' is required", call. = FALSE)
  }
  if (missing(x = dest)) {
    stop("'dest' must be provided for SingleCellExperiment conversion", call. = FALSE)
  }
  if (verbose) message("Converting SingleCellExperiment to Seurat...")
  # Detect available SCE assays and set conversion params accordingly.
  # as.Seurat defaults to data="logcounts" which may not exist.
  sce_assays <- SummarizedExperiment::assayNames(source)
  sce_args <- list(x = source)
  if ("counts" %in% sce_assays) sce_args$counts <- "counts"
  if ("logcounts" %in% sce_assays) {
    sce_args$data <- "logcounts"
  } else {
    # Use alist to preserve NULL (list$key <- NULL removes the key in R)
    sce_args <- c(sce_args, alist(data = NULL))
  }
  extra_args <- list(...)
  sce_args <- c(sce_args, extra_args[!names(extra_args) %in% names(sce_args)])
  seurat_obj <- do.call(Seurat::as.Seurat, sce_args)
  effective_assay <- assay %||% DefaultAssay(object = seurat_obj)
  return(Convert(
    source = seurat_obj, dest = dest, assay = effective_assay,
    overwrite = overwrite, verbose = verbose, standardize = standardize
  ))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Hub Conversion Helper
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Internal hub converter: load any format -> Seurat -> save to any format.
# Used when no direct HDF5-level path exists between (stype, dtype).
#
# @param source_file Path to source file
# @param dest_file Path to destination file
# @param stype Source format extension (lowercase)
# @param dtype Destination format extension (lowercase)
# @param assay Assay name
# @param overwrite Whether to overwrite existing dest
# @param verbose Show progress messages
# @param ... Additional arguments passed to saver
#
# @return Invisible destination filepath
#
# @keywords internal
#
HubConvert <- function(source_file, dest_file, stype, dtype,
                       assay = 'RNA', overwrite = FALSE, verbose = TRUE, ...) {
  loader <- GetLoader(ext = stype)
  saver <- GetSaver(ext = dtype)

  if (is.null(x = loader)) {
    stop("No loader registered for source format '", stype,
         "'. Supported: ", paste(ListFormats()$loaders, collapse = ", "),
         call. = FALSE)
  }
  if (is.null(x = saver)) {
    stop("No saver registered for destination format '", dtype,
         "'. Supported: ", paste(ListFormats()$savers, collapse = ", "),
         call. = FALSE)
  }

  if (verbose) {
    message("Converting ", toupper(x = stype), " -> Seurat -> ", toupper(x = dtype))
  }

  seurat_obj <- loader(file = source_file, assay = assay, verbose = verbose)
  saver(object = seurat_obj, filename = dest_file,
        overwrite = overwrite, verbose = verbose, ...)

  return(invisible(x = dest_file))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Implementations
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# AnnData Preprocessing Helpers
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Sanitize a single column name to valid R/Seurat format
# Replaces /, spaces, commas, semicolons, colons, backslashes with underscore
# Removes consecutive and trailing underscores
# Ensures valid R name (can't start with number)
#
# @param x Character string to sanitize
#
# @return Sanitized character string
#
# @keywords internal
#
SanitizeColumnName <- function(x) {
  # Replace problematic characters with underscore
  x <- gsub(pattern = '[/\\s,;:\\\\]', replacement = '_', x = x, perl = TRUE)
  # Remove consecutive underscores
  x <- gsub(pattern = '_+', replacement = '_', x = x)
  # Remove leading/trailing underscores
  x <- gsub(pattern = '^_|_$', replacement = '', x = x)
  # Ensure valid R name (can't start with number)
  if (grepl('^[0-9]', x)) {
    x <- paste0('X', x)
  }
  return(x)
}

# Sanitize all column names, handling duplicates
# Returns a named vector mapping original -> sanitized names
#
# @param names Character vector of column names to sanitize
#
# @return Named character vector where names are original and values are sanitized
#
# @keywords internal
#
SanitizeColumnNames <- function(names) {
  sanitized <- vapply(names, SanitizeColumnName, character(1), USE.NAMES = FALSE)
  # Handle duplicates by appending __dupN
  if (anyDuplicated(sanitized)) {
    counts <- table(sanitized)
    dups <- names(counts[counts > 1])
    for (dup in dups) {
      idx <- which(sanitized == dup)
      sanitized[idx[-1]] <- paste0(dup, '__dup', seq_along(idx[-1]))
    }
  }
  names(sanitized) <- names  # Map original -> sanitized
  return(sanitized)
}

# Flatten nullable dtype (mask+values group) to simple vector
# h5ad stores nullable dtypes as groups with 'mask' (boolean) and 'values' arrays
# In h5ad nullable dtypes, mask=TRUE means the value is MISSING (NA)
#
# @param col_group H5Group object that may contain mask+values structure
#
# @return Flattened vector with NA values where mask is TRUE, or NULL if not a mask+values structure
#
# @keywords internal
#
FlattenNullable <- function(col_group) {
  if (!inherits(col_group, 'H5Group')) {
    return(NULL)  # Not a group, skip
  }

  # Check if this is a mask+values structure (not categorical which has categories/codes)
  if (col_group$exists('mask') && col_group$exists('values')) {
    # Skip if this looks like a categorical (has categories attribute or encoding-type)
    if (col_group$exists('categories') ||
        isTRUE(x = AttrExists(x = col_group, name = 'encoding-type'))) {
      return(NULL)
    }

    values <- col_group[['values']][]
    mask <- col_group[['mask']][]

    # In h5ad nullable dtypes, mask=TRUE means the value is MISSING (NA)
    values[mask] <- NA
    return(values)
  }

  return(NULL)  # Not a mask+values structure
}

#' Convert AnnData/H5AD files to h5Seurat files
#'
#' @inheritParams Convert
#'
#' @return Returns a handle to \code{dest} as an \code{\link{h5Seurat}} object
#'
#' @importFrom Seurat Project<- DefaultAssay<-
#'
#' @section AnnData/H5AD to h5Seurat:
#' The AnnData/H5AD to h5Seurat conversion will try to automatically fill in
#' datasets based on data presence. It works in the following manner:
#' \subsection{Expression data}{
#'  The expression matrices \code{counts}, \code{data}, and \code{scale.data}
#'  are filled by \code{/X} and \code{/raw/X} in the following manner:
#'  \itemize{
#'   \item \code{counts} will be filled with \code{/raw/X} if present;
#'   otherwise, it will be filled with \code{/X}
#'   \item \code{data} will be filled with \code{/raw/X} if \code{/raw/X} is
#'   present and \code{/X} is dense; otherwise, it will be filled with \code{/X}
#'   \item \code{scale.data} will be filled with \code{/X} if it dense;
#'   otherwise, it will be empty
#'  }
#'  Feature names are taken from the feature-level metadata
#' }
#' \subsection{Feature-level metadata}{
#' Feature-level metadata is added to the \code{meta.features} datasets in each
#' assay. Feature names are taken from the dataset specified by the
#' \dQuote{_index} attribute, the \dQuote{_index} dataset, or the \dQuote{index}
#' dataset, in that order. Metadata is populated with \code{/raw/var} if
#' present, otherwise with \code{/var}; if both \code{/raw/var} and \code{/var}
#' are present, then \code{meta.features} will be populated with \code{/raw/var}
#' first, then \code{/var} will be added to it. For columns present in both
#' \code{/raw/var} and \code{/var}, the values in \code{/var} will be used
#' instead. \strong{Note}: it is possible for \code{/var} to have fewer features
#' than \code{/raw/var}; if this is the case, then only the features present in
#' \code{/var} will be overwritten, with the metadata for features \emph{not}
#' present in \code{/var} remaining as they were in \code{/raw/var} or empty
#' }
#' \subsection{Cell-level metadata}{
#' Cell-level metadata is added to \code{meta.data}; the row names of the
#' metadata (as determined by the value of the \dQuote{_index} attribute, the
#' \dQuote{_index} dataset, or the \dQuote{index} dataset, in that order) are
#' added to the \dQuote{cell.names} dataset instead. If the
#' \dQuote{__categories} dataset is present, each dataset within
#' \dQuote{__categories} will be stored as a factor group. Cell-level metadata
#' will be added as an HDF5 group unless factors are \strong{not} present and
#' the \code{SeuratDisk.dtypes.dataframe_as_group} option is \code{FALSE}
#' }
#' \subsection{Dimensional reduction information:}{
#'  Cell embeddings are taken from \code{/obsm}; dimensional reductions are
#'  named based on their names from \code{obsm} by removing the preceding
#'  \dQuote{X_}.For example, if a dimensional reduction is named \dQuote{X_pca}
#'  in \code{/obsm}, the resulting dimensional reduction information will be
#'  named \dQuote{pca}. The key will be set to one of the following:
#'  \itemize{
#'   \item \dQuote{PC_} if \dQuote{pca} is present in the dimensional reduction
#'   name (\code{grepl("pca", reduction.name, ignore.case = TRUE)})
#'   \item \dQuote{tSNE_} if \dQuote{tsne} is present in the dimensional
#'   reduction name (\code{grepl("tsne", reduction.name, ignore.case = TRUE)})
#'   \item \code{reduction.name_} for all other reductions
#'  }
#'  Remember that the preceding \dQuote{X_} will be removed from the reduction
#'  name before converting to a key. Feature loadings are taken from
#'  \code{/varm} and placed in the associated dimensional reduction. The
#'  dimensional reduction is determine from the loadings name in \code{/varm}:
#'  \itemize{
#'   \item \dQuote{PCs} will be added to a dimensional reduction named
#'   \dQuote{pca}
#'   \item All other loadings in \code{/varm} will be added to a dimensional
#'   reduction named \code{tolower(loading)} (eg. a loading named \dQuote{ICA}
#'   will be added to a dimensional reduction named \dQuote{ica})
#'  }
#'  If a dimensional reduction cannot be found according to the rules above, the
#'  loading will not be taken from the AnnData/H5AD file. Miscellaneous
#'  information will be taken from \code{/uns/reduction} where \code{reduction}
#'  is the name of the reduction in \code{/obsm} without the preceding
#'  \dQuote{X_}; if no dimensional reduction information present, then
#'  miscellaneous information will not be taken from the AnnData/H5AD file.
#'  Standard deviations are taken from a dataset \code{/uns/reduction/variance};
#'  the variances will be converted to standard deviations and added to the
#'  \code{stdev} dataset of a dimensional reduction
#' }
#' \subsection{Nearest-neighbor graph}{
#'  If a nearest neighbor graph is present in \code{/uns/neighbors/distances},
#'  it will be added as a graph dataset in the h5Seurat file and associated with
#'  \code{assay}; if a value is present in \code{/uns/neighbors/params/method},
#'  the name of the graph will be \code{assay_method}, otherwise, it will be
#'  \code{assay_anndata}
#' }
#' \subsection{Layers}{
#'  TODO: add this
#' }
#' \subsection{Miscellaneous information}{
#'  All groups and datasets from \code{/uns} will be copied to \code{misc} in
#'  the h5Seurat file except for the following:
#'  \itemize{
#'   \item Any group or dataset named the same as a dimensional reduction (eg.
#'   \code{/uns/pca})
#'   \item \code{/uns/neighbors}
#'  }
#' }
#'
#' @keywords internal
#'
H5ADToH5Seurat <- function(
  source,
  dest,
  assay = 'RNA',
  overwrite = FALSE,
  verbose = TRUE
) {
  if (file.exists(dest)) {
    if (overwrite) {
      file.remove(dest)
    } else {
      stop("Destination h5Seurat file exists", call. = FALSE)
    }
  }
  dfile <- h5Seurat$new(filename = dest, mode = WriteMode(overwrite = FALSE))
  # Get rownames from an H5AD data frame
  #
  # @param dset Name of data frame
  #
  # @return Returns the name of the dataset that contains the rownames
  #
  GetRownames <- function(dset) {
    if (inherits(x = source[[dset]], what = 'H5Group')) {
      # rownames <- if (source[[dset]]$attr_exists(attr_name = '_index')) {
      rownames <- if (isTRUE(x = AttrExists(x = source[[dset]], name = '_index'))) {
        h5attr(x = source[[dset]], which = '_index')
      } else if (source[[dset]]$exists(name = '_index')) {
        '_index'
      } else if (source[[dset]]$exists(name = 'index')) {
        'index'
      } else {
        stop("Cannot find rownames in ", dset, call. = FALSE)
      }
    } else {
      # TODO: fix this
      stop("Don't know how to handle datasets", call. = FALSE)
      # rownames(x = source[[dset]])
    }
    return(rownames)
  }
  # Read categorical feature names from H5AD
  #
  # @param dset Name of var dataframe group (e.g., 'raw/var' or 'var')
  #
  # @return Returns feature names (gene symbols) if categorical encoding exists, NULL otherwise
  #
  # Cache for ReadCategoricalFeatures results to avoid redundant HDF5 reads
  .categorical_cache <- new.env(parent = emptyenv())

  ReadCategoricalFeatures <- function(dset) {
    # Return cached result if available
    cache_key <- dset
    if (exists(cache_key, envir = .categorical_cache)) {
      return(get(cache_key, envir = .categorical_cache))
    }

    # Check if feature_name exists and is a categorical group
    feature_name_path <- paste0(dset, '/feature_name')
    if (!source$exists(name = feature_name_path)) {
      assign(cache_key, NULL, envir = .categorical_cache)
      return(NULL)
    }

    feature_name_group <- source[[feature_name_path]]
    if (!inherits(x = feature_name_group, what = 'H5Group')) {
      assign(cache_key, NULL, envir = .categorical_cache)
      return(NULL)
    }

    # Check for categorical encoding using simple attr_exists (skip space check)
    if (!feature_name_group$attr_exists(attr_name = 'encoding-type')) {
      assign(cache_key, NULL, envir = .categorical_cache)
      return(NULL)
    }

    encoding_type <- h5attr(x = feature_name_group, which = 'encoding-type')
    if (encoding_type != 'categorical') {
      assign(cache_key, NULL, envir = .categorical_cache)
      return(NULL)
    }

    # Read categories and codes
    if (!feature_name_group$exists(name = 'categories') ||
        !feature_name_group$exists(name = 'codes')) {
      assign(cache_key, NULL, envir = .categorical_cache)
      return(NULL)
    }

    categories <- feature_name_group[['categories']]$read()
    codes <- feature_name_group[['codes']]$read()

    # Decode: categories[codes] (R uses 1-based indexing, codes are 0-based)
    feature_names <- as.character(categories[codes + 1])

    # Cache and return
    assign(cache_key, feature_names, envir = .categorical_cache)
    return(feature_names)
  }
  # Read compound HDF5 dataset using Python h5py
  #
  # Some h5ad files store obs/var as compound HDF5 datasets that hdf5r cannot read.
  # This function uses Python's h5py library as a fallback to read these datasets.
  #
  # @param h5file Path to the h5ad file
  # @param dataset_path Path to the dataset within the file (e.g., 'obs', 'var')
  #
  # @return A data.frame with the dataset contents, or NULL if Python is not available
  #
  ReadCompoundDataset <- function(h5file, dataset_path) {
    # Check if Python with h5py and pandas is available
    python_cmd <- tryCatch({
      system2("python3", c("-c", shQuote("import h5py, pandas")),
              stdout = FALSE, stderr = FALSE)
      "python3"
    }, error = function(e) {
      tryCatch({
        system2("python", c("-c", shQuote("import h5py, pandas")),
                stdout = FALSE, stderr = FALSE)
        "python"
      }, error = function(e) {
        if (verbose) {
          message("  Warning: Python with h5py/pandas not available, cannot read compound dataset")
        }
        return(NULL)
      })
    })

    if (is.null(python_cmd)) {
      return(NULL)
    }

    # Create a temporary file for the CSV output
    temp_csv <- tempfile(fileext = ".csv")
    on.exit(unlink(temp_csv), add = TRUE)

    # Python script to read compound dataset and save as CSV
    # This script handles categorical columns by decoding integer codes to category labels
    # sprintf arguments: h5file, dataset_path (for reading data), dataset_path (for __categories path), temp_csv
    python_script <- sprintf('
import h5py
import pandas as pd
import sys
import numpy as np
import json

try:
    with h5py.File("%s", "r") as f:
        data = f["%s"][:]
        df = pd.DataFrame(data)

        # Decode bytes to strings if needed
        for col in df.columns:
            if df[col].dtype == object:
                try:
                    df[col] = df[col].str.decode("utf-8")
                except:
                    pass

        # Handle categorical columns stored in __categories
        categories_path = "%s/__categories"
        if categories_path in f:
            categories_group = f[categories_path]
            for col in categories_group.keys():
                if col in df.columns:
                    # Read the categories (labels)
                    categories = categories_group[col][:]
                    # Decode bytes to strings if needed (handle bytes, ASCII, and Unicode strings)
                    if categories.dtype == object or categories.dtype.kind in ("S", "U"):
                        categories = np.array([c.decode("utf-8") if isinstance(c, bytes) else str(c) for c in categories])
                    else:
                        categories = categories.astype(str)

                    # Get the codes (integer indices) from the dataframe
                    codes = df[col].values

                    # Map codes to categories using vectorized operations (codes are 0-based indices)
                    # Handle -1 or out-of-range codes as missing/NA
                    valid_mask = (codes >= 0) & (codes < len(categories))
                    decoded = np.empty(len(codes), dtype=object)
                    decoded[valid_mask] = categories[codes[valid_mask]]
                    decoded[~valid_mask] = None  # Will be NA in R

                    df[col] = decoded

        # Handle list/dict columns by converting to JSON or semicolon-delimited strings
        for col in df.columns:
            if df[col].dtype == object:
                # Check if any non-null value is a list or dict
                non_null = df[col].dropna()
                if len(non_null) > 0:
                    sample = non_null.iloc[0]
                    if isinstance(sample, (list, tuple)):
                        # Convert lists to semicolon-delimited or JSON strings
                        def convert_list(x):
                            if x is None or (isinstance(x, float) and np.isnan(x)):
                                return ""
                            if isinstance(x, (list, tuple)):
                                # Simple list of primitives: semicolon-delimited
                                if all(isinstance(i, (str, int, float, type(None))) for i in x):
                                    return ";".join(str(i) for i in x if i is not None)
                                # Complex nested structure: JSON
                                return json.dumps(x)
                            return str(x)
                        df[col] = df[col].apply(convert_list)
                    elif isinstance(sample, dict):
                        # Convert dicts to JSON strings
                        def convert_dict(x):
                            if x is None or (isinstance(x, float) and np.isnan(x)):
                                return ""
                            if isinstance(x, dict):
                                return json.dumps(x)
                            return str(x)
                        df[col] = df[col].apply(convert_dict)

        df.to_csv("%s", index=False)
except Exception as e:
    sys.stderr.write(str(e))
    sys.exit(1)
', h5file, dataset_path, dataset_path, temp_csv)

    # Execute Python script
    result <- tryCatch({
      system2(python_cmd, c("-c", shQuote(python_script)),
              stdout = TRUE, stderr = TRUE)
    }, error = function(e) {
      if (verbose) {
        message("  Warning: Failed to read compound dataset with Python: ", conditionMessage(e))
      }
      return(NULL)
    })

    # Check if CSV was created
    if (!file.exists(temp_csv) || file.size(temp_csv) == 0) {
      if (verbose) {
        message("  Warning: Python failed to read compound dataset")
        if (length(result) > 0) {
          message("  Python error: ", paste(result, collapse = "\n"))
        }
      }
      return(NULL)
    }

    # Read the CSV
    tryCatch({
      df <- read.csv(temp_csv, stringsAsFactors = FALSE, check.names = FALSE)
      return(df)
    }, error = function(e) {
      if (verbose) {
        message("  Warning: Failed to read CSV from Python: ", conditionMessage(e))
      }
      return(NULL)
    })
  }
  # Copy LZF-compressed sparse matrix using Python to temporary file
  #
  # Some h5ad files use LZF compression which hdf5r cannot read.
  # This function uses Python to read LZF-compressed sparse matrices and
  # write them decompressed to a temporary file, which can then be copied.
  #
  # @param src_file Path to source h5ad file
  # @param src_path Path to sparse matrix group in source (e.g., 'X')
  #
  # @return Path to temporary H5 file with decompressed data, or NULL if failed
  #
  CopySparseMatrixDecompressed <- function(src_file, src_path) {
    # Check if Python with h5py is available
    python_cmd <- tryCatch({
      system2("python3", c("-c", shQuote("import h5py, numpy")),
              stdout = FALSE, stderr = FALSE)
      "python3"
    }, error = function(e) {
      tryCatch({
        system2("python", c("-c", shQuote("import h5py, numpy")),
                stdout = FALSE, stderr = FALSE)
        "python"
      }, error = function(e) {
        if (verbose) {
          message("  Warning: Python with h5py not available for LZF decompression")
        }
        return(NULL)
      })
    })

    if (is.null(python_cmd)) {
      return(NULL)
    }

    # Create temporary file
    temp_file <- tempfile(fileext = ".h5")

    # Python script to copy sparse matrix to temp file without LZF compression
    python_script <- sprintf('
import h5py
import numpy as np
import sys

try:
    # Open source file and create temp file
    with h5py.File("%s", "r") as src, h5py.File("%s", "w") as dst:
        # Read the sparse matrix components
        src_grp = src["%s"]

        # Create group at root
        dst_grp = dst.create_group("matrix")

        # Copy each component without LZF compression
        for comp in ["data", "indices", "indptr"]:
            if comp in src_grp:
                data = src_grp[comp][:]
                # Write without compression
                dst_grp.create_dataset(comp, data=data, compression=None)

        # Copy attributes
        for attr_name, attr_value in src_grp.attrs.items():
            dst_grp.attrs[attr_name] = attr_value

        sys.exit(0)
except Exception as e:
    sys.stderr.write(str(e))
    sys.exit(1)
', src_file, temp_file, src_path)

    # Execute Python script
    result <- system2(python_cmd, c("-c", shQuote(python_script)),
                     stdout = TRUE, stderr = TRUE)

    if (!is.null(attr(result, "status")) && attr(result, "status") != 0) {
      if (verbose) {
        message("  Warning: Python LZF decompression failed: ", paste(result, collapse = "\n"))
      }
      unlink(temp_file)
      return(NULL)
    }

    if (!file.exists(temp_file)) {
      return(NULL)
    }

    return(temp_file)
  }
  ColToFactor <- function(dfgroup) {
    if (dfgroup$exists(name = '__categories')) {
      for (i in names(x = dfgroup[['__categories']])) {
        tname <- basename(path = tempfile(tmpdir = ''))
        dfgroup$obj_copy_to(dst_loc = dfgroup, dst_name = tname, src_name = i)
        dfgroup$link_delete(name = i)
        # Because AnnData stores logicals as factors, but have too many levels
        # for factors
        bool.check <- dfgroup[['__categories']][[i]]$dims == 2
        if (isTRUE(x = bool.check)) {
          bool.check <- all(sort(x = dfgroup[['__categories']][[i]][]) == c('False', 'True'))
        }
        if (isTRUE(x = bool.check)) {
          dfgroup$create_dataset(
            name = i,
            robj = dfgroup[[tname]][] + 1L,
            dtype = dfgroup[[tname]]$get_type()
          )
        } else {
          dfgroup$create_group(name = i)
          dfgroup[[i]]$create_dataset(
            name = 'values',
            robj = dfgroup[[tname]][] + 1L,
            dtype = dfgroup[[tname]]$get_type()
          )
          if (IsDType(x = dfgroup[['__categories']][[i]], dtype = 'H5T_STRING')) {
            dfgroup$obj_copy_to(
              dst_loc = dfgroup,
              dst_name = paste0(i, '/levels'),
              src_name = paste0('__categories/', i)
            )
          } else {
            dfgroup[[i]]$create_dataset(
              name = 'levels',
              robj = as.character(x = dfgroup[[H5Path('__categories', i)]][]),
              dtype = StringType()
            )
          }
        }
        dfgroup$link_delete(name = tname)
      }
      dfgroup$link_delete(name = '__categories')
    }
    return(invisible(x = NULL))
  }
  # Helper to create fake cell names when metadata is unavailable
  CreateFakeCellNames <- function() {
    ncells <- if (inherits(x = assay.group[['data']], what = 'H5Group')) {
      assay.group[['data/indptr']]$dims - 1
    } else {
      assay.group[['data']]$dims[2]
    }
    dfile$create_group(name = 'meta.data')
    dfile$create_dataset(
      name = 'cell.names',
      robj = paste0('Cell', seq.default(from = 1, to = ncells)),
      dtype = GuessDType(x = 'Cell1')
    )
  }
  # Sanitize column names and handle nullable dtypes in HDF5 dataframes
  SanitizeH5DFColumns <- function(dfgroup) {
    col_names <- names(dfgroup)
    # Skip internal names
    col_names <- col_names[!col_names %in% c('__categories', '_index')]
    col_name_map <- SanitizeColumnNames(col_names)
    # Only process columns that need renaming or are nullable groups
    needs_rename <- vapply(names(col_name_map), function(n) n != col_name_map[[n]], logical(1))
    for (old_name in names(col_name_map)) {
      new_name <- col_name_map[[old_name]]
      col_obj <- dfgroup[[old_name]]

      # Handle nullable dtype (mask+values) groups
      if (inherits(col_obj, 'H5Group')) {
        flattened <- FlattenNullable(col_obj)
        if (!is.null(flattened)) {
          # Delete old group and create flattened dataset
          dfgroup$link_delete(name = old_name)
          dfgroup$create_dataset(
            name = new_name,
            robj = flattened,
            dtype = GuessDType(x = if (length(flattened) > 0) flattened[1] else flattened)
          )
          next
        }
        # Skip groups that don't need renaming (categorical groups handled elsewhere)
        if (old_name == new_name) next
      }

      # Rename if needed
      if (old_name != new_name) {
        dfgroup$obj_copy_from(
          src_loc = dfgroup,
          src_name = old_name,
          dst_name = new_name
        )
        dfgroup$link_delete(name = old_name)
      }
    }
  }
  ds.map <- c(
    scale.data = if (inherits(x = source[['X']], what = 'H5D')) {
      'X'
    } else {
      NULL
    },
    # Always use X for data slot - X contains normalized data when raw/X has counts
    data = 'X',
    counts = if (source$exists(name = 'raw')) {
      'raw/X'
    } else {
      'X'
    }
  )
  # Add assay data
  assay.group <- dfile[['assays']]$create_group(name = assay)
  for (i in seq_along(along.with = ds.map)) {
    if (verbose) {
      message("Adding ", ds.map[[i]], " as ", names(x = ds.map)[i])
    }
    dst <- names(x = ds.map)[i]

    # Check if this is a sparse matrix with LZF compression
    has_lzf <- FALSE
    if (inherits(x = source[[ds.map[[i]]]], what = 'H5Group')) {
      # Check if this sparse matrix group has LZF-compressed components
      if (source[[ds.map[[i]]]]$exists(name = 'data')) {
        tryCatch({
          # Try to read first element - will fail if LZF compressed
          test_read <- source[[ds.map[[i]]]][['data']][1]
        }, error = function(e) {
          if (grepl("can't synchronously read data|Read failed", conditionMessage(e))) {
            has_lzf <<- TRUE
          }
        })
      }
    }

    # Handle LZF compression with Python
    if (has_lzf) {
      if (verbose) {
        message("  Detected LZF compression; using Python to decompress")
      }

      # Get source file path
      src_file <- if (inherits(x = source, what = 'H5File')) {
        source$filename
      } else {
        as.character(source)
      }

      # Use Python to decompress to temporary file
      temp_file <- CopySparseMatrixDecompressed(
        src_file = src_file,
        src_path = ds.map[[i]]
      )

      if (!is.null(temp_file)) {
        # Copy from temporary file
        temp_h5 <- H5File$new(filename = temp_file, mode = 'r')
        assay.group$obj_copy_from(
          src_loc = temp_h5,
          src_name = 'matrix',
          dst_name = dst
        )
        temp_h5$close_all()
        unlink(temp_file)
      } else {
        if (verbose) {
          message("  Warning: LZF decompression failed, attempting normal copy")
        }
        # Fallback to normal copy (will likely fail on load, but conversion completes)
        assay.group$obj_copy_from(
          src_loc = source,
          src_name = ds.map[[i]],
          dst_name = dst
        )
      }
    } else {
      # Normal copy for non-LZF data
      assay.group$obj_copy_from(
        src_loc = source,
        src_name = ds.map[[i]],
        dst_name = dst
      )
    }

    dst_obj <- assay.group[[dst]]
    if (dst_obj$attr_exists(attr_name = 'shape')) {
      dims <- rev(x = h5attr(x = dst_obj, which = 'shape'))
      dst_obj$create_attr(
        attr_name = 'dims',
        robj = dims,
        dtype = GuessDType(x = dims)
      )
      dst_obj$attr_delete(attr_name = 'shape')
    }
  }
  features.source <- ifelse(
    test = source$exists(name = 'raw') && source$exists(name = 'raw/var'),
    yes = 'raw/var',
    no = 'var'
  )
  # Track if we've already handled compound var dataset
  var_compound_handled <- FALSE

  if (inherits(x = source[[features.source]], what = 'H5Group')) {
    # Try to read categorical feature_name (gene symbols) first
    categorical_features <- ReadCategoricalFeatures(dset = features.source)

    if (!is.null(categorical_features)) {
      # Use gene symbols from categorical encoding
      if (verbose) {
        message("Using gene symbols from categorical feature_name encoding")
      }
      assay.group$create_dataset(
        name = 'features',
        robj = categorical_features,
        dtype = GuessDType(x = categorical_features[1])
      )
    } else {
      # Fallback to index dataset (typically Ensembl IDs)
      features.dset <- GetRownames(dset = features.source)
      assay.group$obj_copy_from(
        src_loc = source,
        src_name = paste(features.source, features.dset, sep = '/'),
        dst_name = 'features'
      )
    }
  } else {
    # When var is an H5D dataset (e.g., compound datatype), try to read it with Python
    if (verbose) {
      message("Detected compound var dataset; attempting to read with Python")
    }

    # Get the source file path (need to handle both filename and H5File)
    h5file_path <- if (inherits(x = source, what = 'H5File')) {
      source$filename
    } else {
      as.character(source)
    }

    var_data <- ReadCompoundDataset(h5file_path, features.source)

    if (!is.null(var_data) && nrow(var_data) > 0) {
      # Successfully read compound dataset with Python
      if (verbose) {
        message("  Successfully read compound var dataset with ", nrow(var_data), " features")
      }

      # Extract feature names: prioritize feature_name (gene symbols) over index (Ensembl IDs)
      feature_names <- if ('feature_name' %in% names(var_data)) {
        as.character(var_data$feature_name)
      } else if ('index' %in% names(var_data)) {
        as.character(var_data$index)
      } else {
        as.character(var_data[[1]])
      }

      # Store the var data for later use in meta.features
      # Create meta.features group and add columns
      if (!assay.group$exists(name = 'meta.features')) {
        assay.group$create_group(name = 'meta.features')
      }
      # Sanitize column names for var metadata
      col_name_map <- SanitizeColumnNames(names(var_data))
      for (col_name in names(var_data)) {
        # Skip the column used for feature names (rownames)
        # If using feature_name for rownames, preserve index (Ensembl IDs) as metadata
        if ('feature_name' %in% names(var_data)) {
          if (col_name == 'feature_name') {
            next  # Skip feature_name as it's used for rownames
          }
        } else if (col_name == 'index' || col_name == names(var_data)[1]) {
          next  # Skip the index column as it's used for rownames
        }
        col_data <- var_data[[col_name]]
        sanitized_name <- col_name_map[[col_name]]
        assay.group[['meta.features']]$create_dataset(
          name = sanitized_name,
          robj = col_data,
          dtype = GuessDType(x = if (length(col_data) > 0) col_data[1] else col_data)
        )
      }
      # Add rownames to meta.features
      assay.group[['meta.features']]$create_dataset(
        name = '_index',
        robj = feature_names,
        dtype = GuessDType(x = feature_names[1])
      )
      assay.group[['meta.features']]$create_attr(
        attr_name = '_index',
        robj = '_index',
        dtype = GuessDType(x = '_index')
      )
      assay.group[['meta.features']]$create_attr(
        attr_name = 's4class',
        robj = 'DFrame',
        dtype = StringType()
      )
      # Set a flag to indicate we've handled meta.features from compound var
      var_compound_handled <- TRUE
    } else {
      # Fallback: create generic feature names
      feature_names <- tryCatch(
        expr = {
          # Try to get rownames directly
          rn <- rownames(x = source[[features.source]])
          # Check if rownames actually returned valid data
          if (is.null(rn) || length(rn) == 0) {
            stop("rownames returned NULL or empty")
          }
          rn
        },
        error = function(e1) {
          # If that fails, try to get dimensions and create generic names
          tryCatch(
            expr = {
              # Get number of features from var dataset dimensions
              nfeatures <- source[[features.source]]$dims
              if (is.null(nfeatures) || length(nfeatures) == 0) {
                stop("dims returned NULL or empty")
              }
              if (verbose) {
                message(
                  "  Warning: Could not read var with Python; using generic feature names ",
                  "(Feature1, Feature2, ..., Feature", nfeatures, ")"
                )
              }
              paste0("Feature", seq_len(nfeatures))
            },
            error = function(e2) {
              # Last resort: get from X matrix dimensions
              if (inherits(x = source[['X']], what = 'H5Group')) {
                # Sparse matrix - check shape attribute
                if (source[['X']]$attr_exists(attr_name = 'shape')) {
                  shape <- h5attr(x = source[['X']], which = 'shape')
                  nfeatures <- shape[1]
                } else {
                  stop("Cannot determine number of features from h5ad file", call. = FALSE)
                }
              } else {
                # Dense matrix
                nfeatures <- source[['X']]$dims[1]
              }
              if (verbose) {
                message(
                  "  Warning: Cannot extract feature names from var; using generic names ",
                  "(Feature1, Feature2, ..., Feature", nfeatures, ")"
                )
              }
              paste0("Feature", seq_len(nfeatures))
            }
          )
        }
      )
    }

    assay.group$create_dataset(
      name = 'features',
      robj = feature_names,
      dtype = GuessDType(x = feature_names[1])
    )
  }
  scaled <- !is.null(x = ds.map['scale.data']) && !is.na(x = ds.map['scale.data'])
  if (scaled) {
    if (inherits(x = source[['var']], what = 'H5Group')) {
      # Try to read categorical feature_name for scaled features
      scaled_categorical_features <- ReadCategoricalFeatures(dset = 'var')

      if (!is.null(scaled_categorical_features)) {
        # Use gene symbols from categorical encoding
        if (verbose) {
          message("Using gene symbols from categorical feature_name encoding for scaled features")
        }
        assay.group$create_dataset(
          name = 'scaled.features',
          robj = scaled_categorical_features,
          dtype = GuessDType(x = scaled_categorical_features[1])
        )
      } else {
        # Fallback to index dataset
        scaled.dset <- GetRownames(dset = 'var')
        assay.group$obj_copy_from(
          src_loc = source,
          src_name = paste0('var/', scaled.dset),
          dst_name = 'scaled.features'
        )
      }
    } else {
      # When var is an H5D dataset, use the same feature names as before
      scaled_feature_names <- tryCatch(
        expr = {
          rn <- rownames(x = source[['var']])
          # Check if rownames actually returned valid data
          if (is.null(rn) || length(rn) == 0) {
            stop("rownames returned NULL or empty")
          }
          rn
        },
        error = function(e) {
          # Get number of scaled features from var or X
          nfeatures <- if (!is.null(source[['var']]$dims) && length(source[['var']]$dims) > 0) {
            source[['var']]$dims
          } else if (inherits(x = source[['X']], what = 'H5Group')) {
            if (source[['X']]$attr_exists(attr_name = 'shape')) {
              h5attr(x = source[['X']], which = 'shape')[1]
            } else {
              stop("Cannot determine number of scaled features", call. = FALSE)
            }
          } else {
            source[['X']]$dims[1]
          }
          paste0("Feature", seq_len(nfeatures))
        }
      )

      assay.group$create_dataset(
        name = 'scaled.features',
        robj = scaled_feature_names,
        dtype = GuessDType(x = scaled_feature_names[1])
      )
    }
  }
  assay.group$create_attr(
    attr_name = 'key',
    robj = paste0(tolower(x = assay), '_'),
    dtype = GuessDType(x = assay)
  )
  # Set default assay
  DefaultAssay(object = dfile) <- assay
  # Add feature-level metadata
  if (!getOption(x = "SeuratDisk.dtypes.dataframe_as_group", default = FALSE)) {
    warning(
      "Adding feature-level metadata as a compound is not yet supported",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  # TODO: Support compound metafeatures
  if (Exists(x = source, name = 'raw/var')) {
    if (inherits(x = source[['raw/var']], what = 'H5Group')) {
      if (verbose) {
        message("Adding meta.features from raw/var")
      }
      assay.group$obj_copy_from(
        src_loc = source,
        src_name = 'raw/var',
        dst_name = 'meta.features'
      )
      if (scaled) {
        features.use <- SafeH5DRead(assay.group[['features']]) %in% SafeH5DRead(assay.group[['scaled.features']])
        features.use <- which(x = features.use)
        meta.scaled <- names(x = source[['var']])
        meta.scaled <- meta.scaled[!meta.scaled %in% c('__categories', scaled.dset)]
        for (mf in meta.scaled) {
          # Skip if not a dataset (e.g., skip groups like feature_types, genome)
          if (!inherits(x = source[['var']][[mf]], what = 'H5D')) {
            if (verbose) {
              message("Skipping ", mf, " (not a dataset)")
            }
            next
          }
          if (!mf %in% names(x = assay.group[['meta.features']])) {
            if (verbose) {
              message("Adding ", mf, " from scaled feature-level metadata")
            }
            assay.group[['meta.features']]$create_dataset(
              name = mf,
              dtype = source[['var']][[mf]]$get_type(),
              space = H5S$new(dims = assay.group[['features']]$dims)
            )
          } else if (verbose) {
            message("Merging ", mf, " from scaled feature-level metadata")
          }
          assay.group[['meta.features']][[mf]][features.use] <- source[['var']][[mf]]$read()
        }
      }
    } else if (!var_compound_handled) {
      warning(
        "Cannot yet add feature-level metadata from compound datasets",
        call. = FALSE,
        immediate. = TRUE
      )
      if (!assay.group$exists(name = 'meta.features')) {
        assay.group$create_group(name = 'meta.features')
      }
    }
  } else {
    if (inherits(x = source[['var']], what = 'H5Group')) {
      if (verbose) {
        message("Adding meta.features from var")
      }
      assay.group$obj_copy_from(
        src_loc = source,
        src_name = 'var',
        dst_name = 'meta.features'
      )
    } else if (!var_compound_handled) {
      # Only warn if we didn't already handle compound var with Python
      warning(
        "Cannot yet add feature-level metadata from compound datasets",
        call. = FALSE,
        immediate. = TRUE
      )
      if (!assay.group$exists(name = 'meta.features')) {
        assay.group$create_group(name = 'meta.features')
      }
    }
  }
  ColToFactor(dfgroup = assay.group[['meta.features']])
  # Sanitize column names and handle nullable dtypes in meta.features
  if (assay.group$exists(name = 'meta.features')) {
    SanitizeH5DFColumns(dfgroup = assay.group[['meta.features']])
  }
  mf_group <- assay.group[['meta.features']]
  if (mf_group$attr_exists(attr_name = 'column-order')) {
    colnames <- h5attr(x = mf_group, which = 'column-order')
    # Sanitize column names in the column-order attribute
    colnames <- vapply(colnames, SanitizeColumnName, character(1), USE.NAMES = FALSE)
    mf_group$create_attr(
      attr_name = 'colnames',
      robj = colnames,
      dtype = GuessDType(x = if (length(colnames) > 0) colnames[1] else colnames)
    )
  }
  if (inherits(x = source[['var']], what = 'H5Group')) {
    rownames_to_delete <- GetRownames(dset = 'var')
    # Only delete if the link exists in meta.features
    if (assay.group[['meta.features']]$exists(name = rownames_to_delete)) {
      assay.group[['meta.features']]$link_delete(name = rownames_to_delete)
    } else if (verbose) {
      message("Note: Rownames '", rownames_to_delete, "' not found in meta.features, skipping deletion")
    }
  }

  # Extract variable features from highly_variable column if available
  # Only check if var is an H5Group (compound datasets don't support sub-paths)
  if (inherits(x = source[['var']], what = 'H5Group') &&
      source$exists(name = 'var/highly_variable')) {
    if (verbose) {
      message("Converting highly_variable to variable.features")
    }

    tryCatch({
      # Check if highly_variable is a categorical (H5Group with categories/codes)
      # CellxGene stores as categorical with string "True"/"False" values
      hv_obj <- source[['var/highly_variable']]

      # Check object type
      is_h5group <- inherits(hv_obj, 'H5Group')
      has_categories <- is_h5group && hv_obj$exists('categories')
      has_codes <- is_h5group && hv_obj$exists('codes')

      if (is_h5group && has_categories && has_codes) {
        # Decode categorical: categories[codes + 1]
        categories <- hv_obj[['categories']][]
        codes <- hv_obj[['codes']][]
        decoded_values <- categories[codes + 1L]
        # Convert string "True"/"False" to logical
        highly_variable <- decoded_values == "True"
        if (verbose) {
          message("Decoded categorical highly_variable: ", sum(highly_variable), " of ", length(highly_variable), " are True")
        }
      } else if (inherits(hv_obj, 'H5D')) {
        # Direct read for non-categorical H5D dataset
        highly_variable <- hv_obj[]

        # Convert to logical if needed
        if (is.numeric(highly_variable)) {
          highly_variable <- as.logical(highly_variable)
        }
      } else {
        # Unsupported type
        if (verbose) {
          message("highly_variable is not a recognized type (H5Group categorical or H5D): ", class(hv_obj)[1])
        }
        highly_variable <- logical(0)
      }

      # Get feature names from /var (same dimension as highly_variable)
      # Use the same feature names as the Seurat object: prioritize feature_name (gene symbols)
      categorical_features <- ReadCategoricalFeatures(dset = 'var')
      if (!is.null(categorical_features)) {
        all_features <- categorical_features
      } else {
        # Fallback to index dataset (Ensembl IDs) or compound labels
        rownames_dset <- GetRownames(dset = 'var')
        if (source[['var']]$exists(rownames_dset)) {
          all_features <- source[[paste('var', rownames_dset, sep = '/')]][]
        } else {
          all_features <- source[['var']]$get_cpd_labels()
        }
      }

      # Get variable features
      if (is.logical(highly_variable) && length(highly_variable) == length(all_features)) {
        variable_features <- all_features[highly_variable]

        # Write to assay group
        if (length(variable_features) > 0) {
          assay.group$create_dataset(
            name = 'variable.features',
            robj = variable_features,
            dtype = GuessDType(x = variable_features)
          )
          if (verbose) {
            message("Added ", length(variable_features), " variable features")
          }
        }
      } else if (verbose) {
        message("Warning: highly_variable format incompatible with features (logical: ",
                is.logical(highly_variable), ", length: ", length(highly_variable),
                ", all_features length: ", length(all_features), ")")
      }
    }, error = function(e) {
      if (verbose) {
        message("Could not convert highly_variable to variable.features: ", conditionMessage(e))
      }
    })
  }

  # Add cell-level metadata
  if (source$exists(name = 'obs') && inherits(x = source[['obs']], what = 'H5Group')) {
    if (!source[['obs']]$exists(name = '__categories') && !getOption(x = "SeuratDisk.dtypes.dataframe_as_group", default = TRUE)) {
      warning(
        "Conversion from H5AD to h5Seurat allowing compound datasets is not yet implemented",
        call. = FALSE,
        immediate. = TRUE
      )
    }
    dfile$obj_copy_from(
      src_loc = source,
      src_name = 'obs',
      dst_name = 'meta.data'
    )
    # Normalize h5ad categorical format (categories/codes) to SeuratDisk format (levels/values)
    # h5ad uses 0-based indexing for codes, while R factors use 1-based indexing
    NormalizeH5ADCategorical <- function(dfgroup) {
      col_names <- names(dfgroup)
      # Skip internal/non-column entries
      col_names <- col_names[!col_names %in% c('__categories', '_index', 'index')]
      for (col_name in col_names) {
        col_obj <- dfgroup[[col_name]]
        # Check if this is an h5ad categorical group (has 'categories' and 'codes')
        if (!inherits(col_obj, 'H5Group')) next
        sub_names <- names(col_obj)
        if (!all(c('categories', 'codes') %in% sub_names)) next
        # Rename 'categories' to 'levels'
        if (!'levels' %in% sub_names) {
          col_obj$obj_copy_from(
            src_loc = col_obj,
            src_name = 'categories',
            dst_name = 'levels'
          )
          col_obj$link_delete(name = 'categories')
        }
        # Convert 'codes' to 'values' with 0-based to 1-based indexing conversion
        if (!'values' %in% sub_names) {
          codes_0based <- col_obj[['codes']]$read()
          codes_dtype <- col_obj[['codes']]$get_type()
          col_obj$create_dataset(
            name = 'values',
            robj = codes_0based + 1L,
            dtype = codes_dtype
          )
          col_obj$link_delete(name = 'codes')
        }
      }
    }
    NormalizeH5ADCategorical(dfgroup = dfile[['meta.data']])
    ColToFactor(dfgroup = dfile[['meta.data']])
    SanitizeH5DFColumns(dfgroup = dfile[['meta.data']])
    md_group <- dfile[['meta.data']]
    if (md_group$attr_exists(attr_name = 'column-order')) {
      colnames <- h5attr(x = md_group, which = 'column-order')
      # Sanitize column names in the column-order attribute
      colnames <- vapply(colnames, SanitizeColumnName, character(1), USE.NAMES = FALSE)
      md_group$create_attr(
        attr_name = 'colnames',
        robj = colnames,
        dtype = GuessDType(x = if (length(colnames) > 0) colnames[1] else colnames)
      )
    }
    rownames <- GetRownames(dset = 'obs')
    dfile$obj_copy_from(
      src_loc = dfile,
      src_name = paste0('meta.data/', rownames),
      dst_name = 'cell.names'
    )
    dfile[['meta.data']]$link_delete(name = rownames)
  } else if (source$exists(name = 'obs') && inherits(x = source[['obs']], what = 'H5D')) {
    # Handle compound obs dataset using Python
    if (verbose) {
      message("Detected compound obs dataset; attempting to read with Python")
    }

    # Get the source file path
    h5file_path <- if (inherits(x = source, what = 'H5File')) {
      source$filename
    } else {
      as.character(source)
    }

    obs_data <- ReadCompoundDataset(h5file_path, 'obs')

    if (!is.null(obs_data) && nrow(obs_data) > 0) {
      # Successfully read compound dataset with Python
      if (verbose) {
        message("  Successfully read compound obs dataset with ", nrow(obs_data), " cells")
      }

      # Extract cell names from the 'index' column or first column
      cell_names <- if ('index' %in% names(obs_data)) {
        as.character(obs_data$index)
      } else {
        as.character(obs_data[[1]])
      }

      # Create meta.data group and add columns
      dfile$create_group(name = 'meta.data')
      # Sanitize column names for obs metadata
      obs_col_name_map <- SanitizeColumnNames(names(obs_data))
      for (col_name in names(obs_data)) {
        if (col_name == 'index' || col_name == names(obs_data)[1]) {
          next  # Skip the index column as it's used for cell names
        }
        col_data <- obs_data[[col_name]]
        sanitized_name <- obs_col_name_map[[col_name]]
        dfile[['meta.data']]$create_dataset(
          name = sanitized_name,
          robj = col_data,
          dtype = GuessDType(x = if (length(col_data) > 0) col_data[1] else col_data)
        )
      }
      # Add rownames to meta.data
      dfile[['meta.data']]$create_dataset(
        name = '_index',
        robj = cell_names,
        dtype = GuessDType(x = cell_names[1])
      )
      dfile[['meta.data']]$create_attr(
        attr_name = '_index',
        robj = '_index',
        dtype = GuessDType(x = '_index')
      )
      dfile[['meta.data']]$create_attr(
        attr_name = 's4class',
        robj = 'DFrame',
        dtype = StringType()
      )

      # Create cell.names dataset
      dfile$create_dataset(
        name = 'cell.names',
        robj = cell_names,
        dtype = GuessDType(x = cell_names[1])
      )
    } else {
      # Fallback if Python reading fails
      warning(
        "Could not read compound obs dataset with Python; creating fake cell names",
        call. = FALSE,
        immediate. = TRUE
      )
      CreateFakeCellNames()
    }
  } else {
    warning(
      "No cell-level metadata present, creating fake cell names",
      call. = FALSE,
      immediate. = TRUE
    )
    CreateFakeCellNames()
  }

  # Restore Seurat-specific metadata from uns['seurat'] if available
  if (source$exists(name = 'uns/seurat')) {
    if (verbose) {
      message("Restoring Seurat-specific metadata from uns['seurat']")
    }

    seurat_group <- source[['uns/seurat']]

    # Check if this file was originally from Seurat
    if (seurat_group$exists(name = 'version')) {
      version_info <- seurat_group[['version']][]
      if (verbose) {
        message("  Found Seurat metadata from ", version_info)
      }
    }

    # Restore original assay name if different from current
    if (seurat_group$exists(name = 'assay_name')) {
      original_assay <- seurat_group[['assay_name']][]
      if (original_assay != assay && verbose) {
        message("  Original assay name was '", original_assay, "', now using '", assay, "'")
      }
    }

    # Note: Graph names are stored for reference, but graphs themselves
    # will be restored from uns/neighbors/distances in the existing code
    if (seurat_group$exists(name = 'graph_names') && verbose) {
      graph_names <- seurat_group[['graph_names']][]
      message("  Original graph names: ", paste(graph_names, collapse = ", "))
    }
  }

  # Add dimensional reduction information
  if (source$exists(name = 'obsm')) {
    # Add cell embeddings
    if (inherits(x = source[['obsm']], what = 'H5Group')) {
      for (reduc in names(x = source[['obsm']])) {
        sreduc <- gsub(pattern = '^X_', replacement = '', x = reduc)

        # Skip spatial - it will be stored in boundaries/centroids for VisiumV2
        if (sreduc == 'spatial') {
          message("Skipping ", reduc, " - spatial data stored in boundaries/centroids")
          next
        }

        reduc.group <- dfile[['reductions']]$create_group(name = sreduc)
        message("Adding ", reduc, " as cell embeddings for ", sreduc)
        Transpose(
          x = source[['obsm']][[reduc]],
          dest = reduc.group,
          dname = 'cell.embeddings',
          verbose = FALSE
        )
        reduc.group$create_group(name = 'misc')
        reduc.group$create_attr(
          attr_name = 'active.assay',
          robj = assay,
          dtype = GuessDType(x = assay)
        )
        key <- paste0(
          if (grepl(pattern = 'pca', x = sreduc, ignore.case = TRUE)) {
            'PC'
          } else if (grepl(pattern = 'tsne', x = sreduc, ignore.case = TRUE)) {
            'tSNE'
          } else {
            sreduc
          },
          '_'
        )
        reduc.group$create_attr(
          attr_name = 'key',
          robj = key,
          dtype = GuessDType(x = reduc)
        )
        global <- BoolToInt(x = grepl(
          pattern = 'tsne|umap',
          x = sreduc,
          ignore.case = TRUE
        ))
        reduc.group$create_attr(
          attr_name = 'global',
          robj = global,
          dtype = GuessDType(x = global)
        )
      }
    } else {
      warning(
        "Reading compound dimensional reductions not yet supported, please update your H5AD file",
        call. = FALSE,
        immediate. = TRUE
      )
    }
    # Add feature loadings
    if (source$exists(name = 'varm')) {
      if (inherits(x = source[['varm']], what = 'H5Group')) {
        for (reduc in names(x = source[['varm']])) {
          sreduc <- switch(EXPR = reduc, 'PCs' = 'pca', tolower(x = reduc))
          if (!isTRUE(x = sreduc %in% names(x = dfile[['reductions']]))) {
            warning(
              "Cannot find a reduction named ",
              sreduc,
              " (",
              reduc,
              " in varm)",
              call. = FALSE,
              immediate. = TRUE
            )
            next
          }
          if (isTRUE(x = verbose)) {
            message("Adding ", reduc, " as feature loadings fpr ", sreduc)
          }
          Transpose(
            x = source[['varm']][[reduc]],
            dest = dfile[['reductions']][[sreduc]],
            dname = 'feature.loadings',
            verbose = FALSE
          )
          reduc.features <- dfile[['reductions']][[sreduc]][['feature.loadings']]$dims[1]
          assay.features <- if (prod(assay.group[['features']]$dims) == reduc.features) {
            'features'
          } else if (assay.group$exists(name = 'scaled.features') && prod(assay.group[['scaled.features']]$dims) == reduc.features) {
            'scaled.features'
          } else {
            NULL
          }
          if (is.null(x = assay.features)) {
            warning(
              "Cannot find features for feature loadings, will not be able to load",
              call. = FALSE,
              immediate. = TRUE
            )
          } else {
            dfile[['reductions']][[sreduc]]$obj_copy_from(
              src_loc = assay.group,
              src_name = assay.features,
              dst_name = 'features'
            )
          }
        }
      } else {
        warning(
          "Reading compound dimensional reductions not yet supported",
          call. = FALSE,
          immediate. = TRUE
        )
      }
    }
    # Add miscellaneous information
    if (source$exists(name = 'uns')) {
      for (reduc in names(x = source[['uns']])) {
        if (!isTRUE(x = reduc %in% names(x = dfile[['reductions']]))) {
          next
        }
        if (verbose) {
          message("Adding miscellaneous information for ", reduc)
        }
        dfile[['reductions']][[reduc]]$link_delete(name = 'misc')
        dfile[['reductions']][[reduc]]$obj_copy_from(
          src_loc = source[['uns']],
          src_name = reduc,
          dst_name = 'misc'
        )
        if ('variance' %in% names(x = dfile[['reductions']][[reduc]][['misc']])) {
          if (verbose) {
            message("Adding standard deviations for ", reduc)
          }
          dfile[['reductions']][[reduc]]$create_dataset(
            name = 'stdev',
            robj = sqrt(x = dfile[['reductions']][[reduc]][['misc']][['variance']][]),
            dtype = GuessDType(x = 1.0)
          )
        }
      }
    }
  }
  # Add project and cell identities
  Project(object = dfile) <- 'AnnData'
  idents <- dfile$create_group(name = 'active.ident')
  idents$create_dataset(
    name = 'values',
    dtype = GuessDType(x = 1L),
    space = H5S$new(dims = dfile[['cell.names']]$dims)
  )
  idents$create_dataset(
    name = 'levels',
    robj = 'AnnData',
    dtype = GuessDType(x = 'AnnData')
  )
  idents[['values']]$write(
    args = list(seq.default(from = 1, to = idents[['values']]$dims)),
    value = 1L
  )
  # Add nearest-neighbor graphs
  # In AnnData: uns/neighbors/distances = k-NN distance matrix
  #             uns/neighbors/connectivities = weighted adjacency/connectivity matrix (like SNN)
  # In Seurat:  assay_nn = nearest neighbor graph (distances)
  #             assay_snn = shared nearest neighbor graph (connectivities)

  if (Exists(x = source, name = 'uns/neighbors')) {
    # Save distances as NN graph
    if (source[['uns/neighbors']]$exists('distances')) {
      nn.graph.name <- paste(assay, 'nn', sep = '_')
      if (verbose) {
        message("Saving nearest-neighbor distances as ", nn.graph.name)
      }
      dfile[['graphs']]$obj_copy_from(
        src_loc = source,
        src_name = 'uns/neighbors/distances',
        dst_name = nn.graph.name
      )
      # Convert 'shape' attribute to 'dims' for compatibility
      if (isTRUE(x = AttrExists(x = dfile[['graphs']][[nn.graph.name]], name = 'shape'))) {
        # Only create 'dims' if it doesn't already exist
        if (!isTRUE(x = AttrExists(x = dfile[['graphs']][[nn.graph.name]], name = 'dims'))) {
          dfile[['graphs']][[nn.graph.name]]$create_attr(
            attr_name = 'dims',
            robj = h5attr(x = dfile[['graphs']][[nn.graph.name]], which = 'shape'),
            dtype = GuessDType(x = h5attr(
              x = dfile[['graphs']][[nn.graph.name]],
              which = 'shape'
            ))
          )
        }
        dfile[['graphs']][[nn.graph.name]]$attr_delete(attr_name = 'shape')
      }
      # Only create assay.used attribute if it doesn't already exist
      if (!isTRUE(x = AttrExists(x = dfile[['graphs']][[nn.graph.name]], name = 'assay.used'))) {
        dfile[['graphs']][[nn.graph.name]]$create_attr(
          attr_name = 'assay.used',
          robj = assay,
          dtype = GuessDType(x = assay)
        )
      }
    }

    # Save connectivities as SNN graph
    if (source[['uns/neighbors']]$exists('connectivities')) {
      snn.graph.name <- paste(assay, 'snn', sep = '_')
      if (verbose) {
        message("Saving nearest-neighbor connectivities as ", snn.graph.name)
      }
      dfile[['graphs']]$obj_copy_from(
        src_loc = source,
        src_name = 'uns/neighbors/connectivities',
        dst_name = snn.graph.name
      )
      # Convert 'shape' attribute to 'dims' for compatibility
      if (isTRUE(x = AttrExists(x = dfile[['graphs']][[snn.graph.name]], name = 'shape'))) {
        # Only create 'dims' if it doesn't already exist
        if (!isTRUE(x = AttrExists(x = dfile[['graphs']][[snn.graph.name]], name = 'dims'))) {
          dfile[['graphs']][[snn.graph.name]]$create_attr(
            attr_name = 'dims',
            robj = h5attr(x = dfile[['graphs']][[snn.graph.name]], which = 'shape'),
            dtype = GuessDType(x = h5attr(
              x = dfile[['graphs']][[snn.graph.name]],
              which = 'shape'
            ))
          )
        }
        dfile[['graphs']][[snn.graph.name]]$attr_delete(attr_name = 'shape')
      }
      # Only create assay.used attribute if it doesn't already exist
      if (!isTRUE(x = AttrExists(x = dfile[['graphs']][[snn.graph.name]], name = 'assay.used'))) {
        dfile[['graphs']][[snn.graph.name]]$create_attr(
          attr_name = 'assay.used',
          robj = assay,
          dtype = GuessDType(x = assay)
        )
      }
    }
  }

  # Also check obsp for neighbor graphs (scanpy default location)
  # This handles h5ad files where neighbors are stored directly in obsp
  # rather than in uns/neighbors
  if (source$exists(name = 'obsp')) {
    obsp_names <- names(x = source[['obsp']])
    nn.graph.name <- paste(assay, 'nn', sep = '_')
    snn.graph.name <- paste(assay, 'snn', sep = '_')

    # Transfer distances as NN graph if not already done from uns/neighbors
    if ('distances' %in% obsp_names && !dfile[['graphs']]$exists(nn.graph.name)) {
      if (verbose) {
        message("Saving obsp/distances as ", nn.graph.name)
      }
      dfile[['graphs']]$obj_copy_from(
        src_loc = source,
        src_name = 'obsp/distances',
        dst_name = nn.graph.name
      )
      # Convert 'shape' attribute to 'dims' for compatibility
      if (isTRUE(x = AttrExists(x = dfile[['graphs']][[nn.graph.name]], name = 'shape'))) {
        if (!isTRUE(x = AttrExists(x = dfile[['graphs']][[nn.graph.name]], name = 'dims'))) {
          dfile[['graphs']][[nn.graph.name]]$create_attr(
            attr_name = 'dims',
            robj = h5attr(x = dfile[['graphs']][[nn.graph.name]], which = 'shape'),
            dtype = GuessDType(x = h5attr(
              x = dfile[['graphs']][[nn.graph.name]],
              which = 'shape'
            ))
          )
        }
        dfile[['graphs']][[nn.graph.name]]$attr_delete(attr_name = 'shape')
      }
      if (!isTRUE(x = AttrExists(x = dfile[['graphs']][[nn.graph.name]], name = 'assay.used'))) {
        dfile[['graphs']][[nn.graph.name]]$create_attr(
          attr_name = 'assay.used',
          robj = assay,
          dtype = GuessDType(x = assay)
        )
      }
    }

    # Transfer connectivities as SNN graph if not already done from uns/neighbors
    if ('connectivities' %in% obsp_names && !dfile[['graphs']]$exists(snn.graph.name)) {
      if (verbose) {
        message("Saving obsp/connectivities as ", snn.graph.name)
      }
      dfile[['graphs']]$obj_copy_from(
        src_loc = source,
        src_name = 'obsp/connectivities',
        dst_name = snn.graph.name
      )
      # Convert 'shape' attribute to 'dims' for compatibility
      if (isTRUE(x = AttrExists(x = dfile[['graphs']][[snn.graph.name]], name = 'shape'))) {
        if (!isTRUE(x = AttrExists(x = dfile[['graphs']][[snn.graph.name]], name = 'dims'))) {
          dfile[['graphs']][[snn.graph.name]]$create_attr(
            attr_name = 'dims',
            robj = h5attr(x = dfile[['graphs']][[snn.graph.name]], which = 'shape'),
            dtype = GuessDType(x = h5attr(
              x = dfile[['graphs']][[snn.graph.name]],
              which = 'shape'
            ))
          )
        }
        dfile[['graphs']][[snn.graph.name]]$attr_delete(attr_name = 'shape')
      }
      if (!isTRUE(x = AttrExists(x = dfile[['graphs']][[snn.graph.name]], name = 'assay.used'))) {
        dfile[['graphs']][[snn.graph.name]]$create_attr(
          attr_name = 'assay.used',
          robj = assay,
          dtype = GuessDType(x = assay)
        )
      }
    }
  }

  # Add miscellaneous information
  if (source$exists(name = 'uns')) {
    misc <- setdiff(
      x = names(x = source[['uns']]),
      y = c('neighbors', names(x = dfile[['reductions']]))
    )
    for (i in misc) {
      if (verbose) {
        message("Adding ", i, " to miscellaneous data")
      }
      tryCatch({
        dfile[['misc']]$obj_copy_from(
          src_loc = source[['uns']],
          src_name = i,
          dst_name = i
        )
      }, error = function(e) {
        if (verbose) {
          message("  Warning: Could not copy ", i, " to miscellaneous data: ", conditionMessage(e))
        }
      })
    }
  }
  # Add spatial images from uns/spatial
  if (source$exists(name = 'uns/spatial')) {
    if (verbose) {
      message("Adding spatial images from uns/spatial")
    }
    spatial_uns <- source[['uns/spatial']]
    library_ids <- names(x = spatial_uns)

    # Detect coordinate format early to determine if image needs transformation
    # Native h5ad (scanpy/CellXGene): shape is (n_cells, 2)
    # h5ad from Seurat conversion: shape is (2, n_cells)
    is_native_h5ad <- FALSE
    spatial_key <- if (source$exists(name = 'obsm/X_spatial')) {
      'obsm/X_spatial'
    } else if (source$exists(name = 'obsm/spatial')) {
      'obsm/spatial'
    } else {
      NULL
    }
    if (!is.null(spatial_key)) {
      coords_dims <- dim(source[[spatial_key]]$read())
      is_native_h5ad <- coords_dims[2] == 2  # (n_cells, 2) format
      if (verbose) {
        message("  Detected ", if(is_native_h5ad) "native h5ad" else "Seurat-converted", " coordinate format")
      }
    }

    for (lib_id in library_ids) {
      lib_group <- spatial_uns[[lib_id]]

      # Skip non-group items (e.g., 'is_single' dataset)
      if (!inherits(lib_group, 'H5Group')) {
        next
      }

      # Check if images exist for this library
      if (lib_group$exists(name = 'images')) {
        images_group <- lib_group[['images']]
        image_types <- names(x = images_group)

        # Create image group in h5seurat (using library_id as image name)
        # Remove illegal characters from library_id for use as image name
        img_name <- gsub(pattern = '[^A-Za-z0-9_]', replacement = '_', x = lib_id)

        if (verbose) {
          message("  Processing library '", lib_id, "' as image '", img_name, "'")
        }

        # Create image group if it doesn't exist
        if (!dfile[['images']]$exists(name = img_name)) {
          img_h5 <- dfile[['images']]$create_group(name = img_name)

          # Add required attributes for image indexing
          # assay attribute
          img_h5$create_attr(
            attr_name = 'assay',
            robj = assay,  # Use the assay parameter
            dtype = GuessDType(x = assay)
          )

          # s4class attribute (VisiumV2 for spatial data)
          img_h5$create_attr(
            attr_name = 's4class',
            robj = 'VisiumV2',
            dtype = GuessDType(x = 'VisiumV2')
          )

          # global attribute (make image globally available)
          img_h5$create_attr(
            attr_name = 'global',
            robj = 1L,  # Use integer 1 for TRUE
            dtype = GuessDType(x = 1L)
          )
        } else {
          img_h5 <- dfile[['images']][[img_name]]
        }

        # Add key dataset (required by Seurat)
        if (!img_h5$exists(name = 'key')) {
          img_h5$create_dataset(
            name = 'key',
            robj = paste0(img_name, '_'),
            dtype = GuessDType(x = paste0(img_name, '_'))
          )
        }

        # Store lowres as the primary 'image' dataset (used for visualization)
        # Note: hires can be added later if needed
        if ('lowres' %in% image_types) {
          if (verbose) {
            message("    Adding lowres as primary image")
          }

          # Read lowres image data from h5ad
          img_data <- images_group[['lowres']]$read()

          # Check if this is a 3D array
          if (length(dim(img_data)) == 3L) {
            # h5ad/anndata stores images as (channels, width, height) in HDF5
            # Seurat expects (height, width, channels)
            # aperm(c(3,2,1)) converts (C,W,H) -> (H,W,C)
            img_arr <- aperm(img_data, c(3L, 2L, 1L))

            # Write to h5seurat as 'image' dataset
            if (img_h5$exists(name = 'image')) {
              img_h5$link_delete(name = 'image')
            }
            img_h5$create_dataset(
              name = 'image',
              robj = img_arr,
              dtype = GuessDType(x = img_arr)
            )
          } else {
            if (verbose) {
              message("      Warning: Lowres image is not a 3D array, skipping")
            }
          }
        } else if ('hires' %in% image_types) {
          # Fallback to hires if lowres not available
          if (verbose) {
            message("    Adding hires as primary image (lowres not available)")
          }

          img_data <- images_group[['hires']]$read()
          if (length(dim(img_data)) == 3L) {
            # h5ad/anndata stores images as (channels, width, height) in HDF5
            # Seurat expects (height, width, channels)
            # aperm(c(3,2,1)) converts (C,W,H) -> (H,W,C)
            img_arr <- aperm(img_data, c(3L, 2L, 1L))
            if (img_h5$exists(name = 'image')) {
              img_h5$link_delete(name = 'image')
            }
            img_h5$create_dataset(
              name = 'image',
              robj = img_arr,
              dtype = GuessDType(x = img_arr)
            )
          }
        }

        # Add scale factors if available
        if (lib_group$exists(name = 'scalefactors')) {
          if (verbose) {
            message("    Adding scale factors")
          }

          # Create scale.factors group
          if (img_h5$exists(name = 'scale.factors')) {
            img_h5$link_delete(name = 'scale.factors')
          }
          sf_group <- img_h5$create_group(name = 'scale.factors')

          # Map from h5ad names to h5seurat names
          sf_map <- c(
            tissue_hires_scalef = 'hires',
            tissue_lowres_scalef = 'lowres',
            spot_diameter_fullres = 'spot',
            fiducial_diameter_fullres = 'fiducial'
          )

          sf_src <- lib_group[['scalefactors']]
          for (sf_h5ad_name in names(x = sf_src)) {
            # Find corresponding h5seurat name from the map
            # sf_map has h5ad names as keys and h5seurat names as values
            if (sf_h5ad_name %in% names(sf_map)) {
              sf_h5seurat_name <- unname(sf_map[sf_h5ad_name])
            } else {
              # If no mapping found, use the h5ad name directly
              sf_h5seurat_name <- sf_h5ad_name
            }

            sf_value <- sf_src[[sf_h5ad_name]][]

            # Convert diameter to radius for spot and fiducial
            if (sf_h5ad_name %in% c('spot_diameter_fullres', 'fiducial_diameter_fullres')) {
              sf_value <- sf_value / 2.0
            }

            sf_group$create_dataset(
              name = sf_h5seurat_name,
              robj = sf_value,
              dtype = GuessDType(x = sf_value)
            )
          }

          # If lowres scale factor is missing but hires exists, set lowres = hires
          # This is needed when h5ad only has hires image (no lowres)
          # Seurat's SpatialFeaturePlot defaults to using lowres scale
          if (!sf_group$exists('lowres') && sf_group$exists('hires')) {
            hires_val <- sf_group[['hires']][]
            sf_group$create_dataset(
              name = 'lowres',
              robj = hires_val,
              dtype = GuessDType(x = hires_val)
            )
            if (verbose) {
              message("      Setting lowres scale factor = hires (", hires_val, ")")
            }
          }
        }

        # Create boundaries/centroids structure for VisiumV2
        # This requires spatial coordinates from obsm/X_spatial or obsm/spatial
        spatial_key <- if (source$exists(name = 'obsm/X_spatial')) {
          'obsm/X_spatial'
        } else if (source$exists(name = 'obsm/spatial')) {
          'obsm/spatial'
        } else {
          NULL
        }

        if (!is.null(spatial_key)) {
          if (verbose) {
            message("    Creating boundaries/centroids for VisiumV2 (using ", spatial_key, ")")
          }

          # Read spatial coordinates from h5ad
          # Use $read() instead of [] to avoid HDF5 indexing issues
          coords_h5ad <- source[[spatial_key]]$read()

          # Detect coordinate orientation:
          # - Native h5ad (scanpy/CellXGene): shape is (n_cells, 2) where rows are cells
          # - h5ad from Seurat conversion: shape is (2, n_cells) where columns are cells
          # We detect by checking which dimension is 2
          coords_dims <- dim(coords_h5ad)
          is_native_h5ad <- coords_dims[2] == 2  # (n_cells, 2) format

          if (verbose) {
            message("    Coordinate array shape: ", paste(coords_dims, collapse=" x "),
                    " (", if(is_native_h5ad) "native h5ad format" else "Seurat-converted format", ")")
          }

          # Get cell names from the h5Seurat file (already loaded earlier)
          cell_names <- as.character(dfile[['cell.names']][])

          # Get library IDs to filter cells for this library
          # Try common library ID column names
          lib_id_col <- NULL
          lib_ids_vec <- NULL
          for (col_name in c('sangerID', 'library_id', 'sample', 'batch')) {
            if (source$exists(paste0('obs/', col_name))) {
              lib_id_obj <- source[[paste0('obs/', col_name)]]

              # Handle categorical/factor structure
              if (inherits(lib_id_obj, 'H5Group')) {
                if (lib_id_obj$exists('codes') && lib_id_obj$exists('categories')) {
                  codes <- lib_id_obj[['codes']][]
                  categories <- lib_id_obj[['categories']][]
                  lib_ids_vec <- categories[codes + 1]  # Convert to 1-indexed
                  lib_id_col <- col_name
                  break
                }
              } else {
                lib_ids_vec <- as.character(lib_id_obj[])
                lib_id_col <- col_name
                break
              }
            }
          }

          # Filter cells for this library
          if (!is.null(lib_ids_vec)) {
            cell_mask <- lib_ids_vec == lib_id
            if (sum(cell_mask) == 0) {
              if (verbose) {
                message("    Warning: No cells found for library ", lib_id, ", skipping centroids")
              }
            } else {
              # Filter coordinates based on format
              if (is_native_h5ad) {
                # Native format (n_cells, 2): filter rows
                filtered_coords <- coords_h5ad[cell_mask, , drop = FALSE]
              } else {
                # Seurat-converted format (2, n_cells): filter columns
                filtered_coords <- coords_h5ad[, cell_mask, drop = FALSE]
              }
              filtered_cells <- cell_names[cell_mask]

              if (verbose) {
                message("    Filtered to ", length(filtered_cells), " cells for library ", lib_id)
              }

              # Convert to (n, 2) format for centroids
              if (is_native_h5ad) {
                # Native h5ad: already in (n, 2) format
                coords_matrix <- filtered_coords
              } else {
                # Seurat-converted: transpose from (2, n) to (n, 2)
                coords_matrix <- t(filtered_coords)
              }

              # Note: Coordinates are kept in full resolution
              # Seurat's GetTissueCoordinates() will handle scaling via scale.factors
              # when called with scale = "lowres" or scale = "hires"

              # Calculate radius from scale factors (in full resolution)
              radius <- if (img_h5$exists('scale.factors/spot')) {
                img_h5[['scale.factors/spot']][]
              } else {
                # Default Visium spot diameter is ~143 pixels in full res
                # (based on spot_diameter_fullres from the data)
                71.5  # radius = diameter / 2
              }

              # Create boundaries group
              if (!img_h5$exists('boundaries')) {
                boundaries_group <- img_h5$create_group('boundaries')
              } else {
                boundaries_group <- img_h5[['boundaries']]
              }

              # Create centroids group
              if (boundaries_group$exists('centroids')) {
                boundaries_group$link_delete('centroids')
              }
              centroids_group <- boundaries_group$create_group('centroids')

              # Store coords (n x 2 matrix)
              centroids_group$create_dataset(
                name = 'coords',
                robj = coords_matrix,
                dtype = GuessDType(x = coords_matrix)
              )

              # Store cell names
              centroids_group$create_dataset(
                name = 'cells',
                robj = filtered_cells,
                dtype = GuessDType(x = filtered_cells)
              )

              # Store radius (spot radius for visualization)
              centroids_group$create_dataset(
                name = 'radius',
                robj = radius,
                dtype = GuessDType(x = radius)
              )

              # Store nsides (0 = infinite sides/circle, returns centroids not vertices)
              centroids_group$create_dataset(
                name = 'nsides',
                robj = 0L,
                dtype = GuessDType(x = 0L)
              )

              # Store theta (rotation angle, 0 for Visium)
              centroids_group$create_dataset(
                name = 'theta',
                robj = 0.0,
                dtype = GuessDType(x = 0.0)
              )

              if (verbose) {
                message("    Created centroids with ", nrow(coords_matrix), " spots")
              }
            }
          } else {
            if (verbose) {
              message("    Warning: No library ID column found, using all cells for centroids")
            }

            # Use all cells if no library ID found
            # Convert to (n, 2) format for centroids based on detected format
            if (is_native_h5ad) {
              # Native h5ad: already in (n, 2) format
              coords_matrix <- coords_h5ad
            } else {
              # Seurat-converted: transpose from (2, n) to (n, 2)
              coords_matrix <- t(coords_h5ad)
            }

            radius <- if (img_h5$exists('scale.factors/spot')) {
              img_h5[['scale.factors/spot']][]
            } else {
              71.5  # Default radius in full resolution
            }

            if (!img_h5$exists('boundaries')) {
              boundaries_group <- img_h5$create_group('boundaries')
            } else {
              boundaries_group <- img_h5[['boundaries']]
            }

            if (boundaries_group$exists('centroids')) {
              boundaries_group$link_delete('centroids')
            }
            centroids_group <- boundaries_group$create_group('centroids')

            centroids_group$create_dataset(
              name = 'coords',
              robj = coords_matrix,
              dtype = GuessDType(x = coords_matrix)
            )

            centroids_group$create_dataset(
              name = 'cells',
              robj = cell_names,
              dtype = GuessDType(x = cell_names)
            )

            centroids_group$create_dataset(
              name = 'radius',
              robj = radius,
              dtype = GuessDType(x = radius)
            )

            centroids_group$create_dataset(
              name = 'nsides',
              robj = 0L,
              dtype = GuessDType(x = 0L)
            )

            centroids_group$create_dataset(
              name = 'theta',
              robj = 0.0,
              dtype = GuessDType(x = 0.0)
            )
          }
        }
      }
    }
  }
  # Add layers to the RNA assay (V5 style)
  if (Exists(x = source, name = 'layers')) {
    # Get the RNA assay
    rna.assay <- dfile[['assays']][[assay]]

    # For V5, we add layers as named slots within the assay
    # Common layer names: counts, data, scale.data
    for (layer in names(x = source[['layers']])) {
      # Map layer names to appropriate slots
      # Check if this layer should be added (avoid duplicates)
      layer.name <- layer

      # Skip if this slot already exists (was filled from raw/X or X)
      if (rna.assay$exists(name = layer.name)) {
        if (verbose) {
          message("Skipping layer ", layer, " (already exists in RNA assay)")
        }
        next
      }

      if (verbose) {
        message("Adding layer ", layer, " to RNA assay")
      }

      # Copy the layer data to the appropriate slot in RNA assay
      rna.assay$obj_copy_from(
        src_loc = source[['layers']],
        src_name = layer,
        dst_name = layer.name
      )

      # Fix dimensions attributes if needed
      if (isTRUE(x = AttrExists(x = rna.assay[[layer.name]], name = 'shape'))) {
        dims <- rev(x = h5attr(x = rna.assay[[layer.name]], which = 'shape'))
        rna.assay[[layer.name]]$create_attr(
          attr_name = 'dims',
          robj = dims,
          dtype = GuessDType(x = dims)
        )
        rna.assay[[layer.name]]$attr_delete(attr_name = 'shape')
      }
    }
  }
  return(dfile)
}

#' Convert h5Seurat files to H5AD files
#'
#' @inheritParams Convert
#'
#' @return Returns a handle to \code{dest} as an \code{\link[hdf5r]{H5File}}
#' object
#'
#' @section h5Seurat to AnnData/H5AD:
#' The h5Seurat to AnnData/H5AD conversion will try to automatically fill in
#' datasets based on data presence. Data presense is determined by the h5Seurat
#' index (\code{source$index()}). It works in the following manner:
#' \subsection{Assay data}{
#'  \itemize{
#'   \item \code{X} will be filled with \code{scale.data} if \code{scale.data}
#'   is present; otherwise, it will be filled with \code{data}
#'   \item \code{var} will be filled with \code{meta.features} \strong{only} for
#'   the features present in \code{X}; for example, if \code{X} is filled with
#'   \code{scale.data}, then \code{var} will contain only features that have
#'   been scaled
#'   \item \code{raw.X} will be filled with \code{data} if \code{X} is filled
#'   with \code{scale.data}; otherwise, it will be filled with \code{counts}. If
#'   \code{counts} is not present, then \code{raw} will not be filled
#'   \item \code{raw.var} will be filled with \code{meta.features} with the
#'   features present in \code{raw.X}; if \code{raw.X} is not filled, then
#'   \code{raw.var} will not be filled
#'  }
#' }
#' \subsection{Cell-level metadata}{
#'  Cell-level metadata is added to \code{obs}
#' }
#' \subsection{Dimensional reduction information}{
#'  Only dimensional reductions associated with \code{assay} or marked as
#'  \link[SeuratObject:IsGlobal]{global} will be transfered to the H5AD file. For
#'  every reduction \code{reduc}:
#'  \itemize{
#'   \item cell embeddings are placed in \code{obsm} and renamed to
#'   \code{X_reduc}
#'   \item feature loadings, if present, are placed in \code{varm} and renamed
#'   to either \dQuote{PCs} if \code{reduc} is \dQuote{pca} otherwise
#'   \code{reduc} in all caps
#'  }
#'  For example, if \code{reduc} is \dQuote{ica}, then cell embeddings will be
#'  \dQuote{X_ica} in \code{obsm} and feature loaodings, if present, will be
#'  \dQuote{ICA} in \code{varm}
#' }
#' \subsection{Nearest-neighbor graphs}{
#'  If a nearest-neighbor graph is associated with \code{assay}, it will be
#'  added to \code{uns/neighbors/distances}; if more than one graph is present,
#'  then \strong{only} the last graph according to the index will be added.
#' }
#' \subsection{Layers}{
#'  Data from other assays can be added to \code{layers} if they have the same
#'  shape as \code{X} (same number of cells and features). To determine this,
#'  the shape of each alternate assays's \code{scale.data} and \code{data} slots
#'  are determined. If they are the same shape as \code{X}, then that slot
#'  (\code{scale.data} is given priority over \code{data}) will be added as a
#'  layer named the name of the assay (eg. \dQuote{SCT}). In addition, the
#'  features names will be added to \code{var} as \code{assay_features}
#'  (eg. \dQuote{SCT_features}).
#' }
#'
#' @keywords internal
#'
H5SeuratToH5AD <- function(
  source,
  dest,
  assay = DefaultAssay(object = source),
  overwrite = FALSE,
  verbose = TRUE,
  standardize = FALSE
) {
  if (file.exists(dest)) {
    if (overwrite) {
      file.remove(dest)
    } else {
      stop("Destination H5AD file exists", call. = FALSE)
    }
  }
  rownames <- '_index'
  dfile <- H5File$new(filename = dest, mode = WriteMode(overwrite = FALSE))

  # Define mapping tables for Seurat to scanpy naming conventions
  if (standardize) {
    if (verbose) {
      message("Standardization enabled: will convert Seurat names to scanpy conventions")
    }
    # Cell-level metadata (obs) mappings
    obs_name_map <- c(
      "nCount_RNA" = "n_counts",
      "nFeature_RNA" = "n_genes",
      "percent.mt" = "percent_mito",
      "seurat_clusters" = "clusters"
    )

    # Feature-level metadata (var) mappings
    var_name_map <- c(
      "vst.variable" = "highly_variable",
      "vst.mean" = "means",
      "vst.variance" = "dispersions",
      "vst.variance.standardized" = "dispersions_norm",
      "vst.variance.expected" = "dispersions_expected"
    )
  } else {
    obs_name_map <- character(0)
    var_name_map <- character(0)
  }

  # Transfer data frames from h5Seurat files to H5AD files
  #
  # @param src Source dataset
  # @param dname Name of destination
  # @param index Integer values of rows to take
  #
  # @return Invisibly returns \code{NULL}
  #
  TransferDF <- function(src, dname, index) {
    # Helper function to apply name mapping based on destination type
    apply_name_map <- function(col_name, dest_type) {
      if (!standardize) return(col_name)

      name_map <- if (dest_type == "obs") obs_name_map else var_name_map

      # If column name has a mapping, use it; otherwise keep original
      if (col_name %in% names(name_map)) {
        return(unname(name_map[col_name]))
      }
      return(col_name)
    }

    # Handle V5 environment-wrapped metadata
    if (dname == "obs" && is.environment(x = src)) {
      # Extract from environment if possible
      if (exists("hgroup", envir = src, inherits = FALSE)) {
        src <- get("hgroup", envir = src, inherits = FALSE)
      } else if (source$exists(name = 'meta.data')) {
        # Direct metadata transfer for V5
        if (verbose) {
          message("Handling V5 metadata transfer for obs")
        }

        # Create obs group
        if (!dfile$exists(name = 'obs')) {
          dfile$create_group(name = 'obs')
        }

        # Transfer metadata columns
        # Access the underlying HDF5 group directly to avoid V5 environment wrappers
        meta_group <- source[['meta.data']]
        if (is.environment(x = meta_group) && exists("hgroup", envir = meta_group, inherits = FALSE)) {
          meta_group <- get("hgroup", envir = meta_group, inherits = FALSE)
        }
        if (inherits(x = meta_group, what = 'H5Group')) {
          for (col in names(x = meta_group)) {
            if (col == '__categories' || col == '_index') next

            # Apply name mapping if standardize is enabled
            mapped_col <- apply_name_map(col, "obs")

            if (IsFactor(x = meta_group[[col]])) {
              # Use newer anndata format: each categorical is a group with categories/codes
              dfile[['obs']]$create_group(name = mapped_col)

              # Convert codes to 0-based indices, handling NA values
              # R factor codes are 1-based; h5ad uses 0-based with -1 for missing
              raw_codes <- meta_group[[col]][['values']][index]
              codes_values <- ifelse(is.na(raw_codes), -1L, raw_codes - 1L)
              n_categories <- length(meta_group[[col]][['levels']][])
              has_na <- any(is.na(raw_codes))

              # Choose appropriate integer dtype based on number of categories
              # Must use signed integers when NA values present (need -1)
              codes_dtype <- if (has_na || n_categories <= 127) {
                hdf5r::h5types$H5T_NATIVE_INT8
              } else if (n_categories <= 255) {
                hdf5r::h5types$H5T_NATIVE_UINT8
              } else if (n_categories <= 32767) {
                hdf5r::h5types$H5T_NATIVE_INT16
              } else if (n_categories <= 65535) {
                hdf5r::h5types$H5T_NATIVE_UINT16
              } else {
                hdf5r::h5types$H5T_NATIVE_INT32
              }

              # Add codes dataset with proper dtype for scanpy compatibility
              dfile[['obs']][[mapped_col]]$create_dataset(
                name = 'codes',
                robj = codes_values,
                dtype = codes_dtype
              )

              # Add encoding attributes to codes
              dfile[['obs']][[mapped_col]][['codes']]$create_attr(
                attr_name = 'encoding-type',
                robj = 'array',
                dtype = GuessDType(x = 'array'),
                space = H5S$new(type = "scalar")
              )
              dfile[['obs']][[mapped_col]][['codes']]$create_attr(
                attr_name = 'encoding-version',
                robj = '0.2.0',
                dtype = GuessDType(x = '0.2.0'),
                space = H5S$new(type = "scalar")
              )

              # Add categories dataset
              dfile[['obs']][[mapped_col]]$create_dataset(
                name = 'categories',
                robj = as.character(meta_group[[col]][['levels']][]),
                dtype = StringType('utf8')
              )

              # Add encoding attributes to categories
              AddAnndataEncoding(dfile[['obs']][[mapped_col]][['categories']], encoding_type = 'string-array')

              # Add attributes to the categorical group
              dfile[['obs']][[mapped_col]]$create_attr(
                attr_name = 'encoding-type',
                robj = 'categorical',
                dtype = GuessDType(x = 'categorical'),
                space = H5S$new(type = "scalar")
              )
              dfile[['obs']][[mapped_col]]$create_attr(
                attr_name = 'encoding-version',
                robj = '0.2.0',
                dtype = GuessDType(x = '0.2.0'),
                space = H5S$new(type = "scalar")
              )
              dfile[['obs']][[mapped_col]]$create_attr(
                attr_name = 'ordered',
                robj = FALSE,
                dtype = GuessDType(x = FALSE),
                space = H5S$new(type = "scalar")
              )
            } else {
              # Regular column
              dfile[['obs']]$create_dataset(
                name = mapped_col,
                robj = if (is.null(index)) meta_group[[col]][] else meta_group[[col]][index],
                dtype = meta_group[[col]]$get_type()
              )
            }
          }

          # Add column order with mapped names
          if (meta_group$attr_exists(attr_name = 'colnames')) {
            original_colnames <- h5attr(x = meta_group, which = 'colnames')
            if (standardize) {
              # Apply mapping to column names
              mapped_colnames <- sapply(original_colnames, function(cn) apply_name_map(cn, "obs"))
            } else {
              mapped_colnames <- original_colnames
            }
            dfile[['obs']]$create_attr(
              attr_name = 'column-order',
              robj = mapped_colnames,
              dtype = GuessDType(x = mapped_colnames)
            )
          }

          # Add encoding attributes
          encoding.info <- c('type' = 'dataframe', 'version' = '0.2.0')
          names(x = encoding.info) <- paste0('encoding-', names(x = encoding.info))
          for (i in seq_along(along.with = encoding.info)) {
            attr.name <- names(x = encoding.info)[i]
            if (!dfile[['obs']]$attr_exists(attr_name = attr.name)) {
              dfile[['obs']]$create_attr(
                attr_name = attr.name,
                robj = encoding.info[i],
                dtype = GuessDType(x = encoding.info[i]),
                space = H5S$new(type = "scalar")
              )
            }
          }
        }
        return(invisible(x = NULL))
      }
    }

    # Validate source for standard transfer
    if (is.null(x = src) || !inherits(x = src, what = c('H5D', 'H5Group'))) {
      warning(paste("TransferDF: invalid source for", dname), immediate. = TRUE)
      return(invisible(x = NULL))
    }

    if (verbose) {
      message("Transfering ", basename(path = src$get_obj_name()), " to ", dname)
    }

    if (inherits(x = src, what = 'H5D')) {
      # Create a name mapping function for this dname
      map_fn <- if (standardize) {
        function(col_name) apply_name_map(col_name, dname)
      } else {
        NULL
      }
      CompoundToGroup(
        src = src,
        dst = dfile,
        dname = dname,
        order = 'column-order',
        index = index,
        name_map_fn = map_fn
      )
      # Add dataframe encoding attributes (required by anndata)
      encoding.info <- c('type' = 'dataframe', 'version' = '0.2.0')
      names(x = encoding.info) <- paste0('encoding-', names(x = encoding.info))
      for (i in seq_along(along.with = encoding.info)) {
        attr.name <- names(x = encoding.info)[i]
        attr.value <- encoding.info[i]
        if (dfile[[dname]]$attr_exists(attr_name = attr.name)) {
          dfile[[dname]]$attr_delete(attr_name = attr.name)
        }
        dfile[[dname]]$create_attr(
          attr_name = attr.name,
          robj = attr.value,
          dtype = GuessDType(x = attr.value),
          space = H5S$new(type = "scalar")
        )
      }
    } else if (inherits(x = src, what = 'H5Group')) {
      dfile$create_group(name = dname)
      for (i in src$names) {
        # Skip internal HDF5 entries
        if (i %in% c('__categories', '_index')) next
        # Apply name mapping if standardize is enabled
        mapped_i <- apply_name_map(i, dname)

        if (IsFactor(x = src[[i]])) {
          # Use newer anndata format: each categorical is a group with categories/codes
          dfile[[dname]]$create_group(name = mapped_i)

          # Convert codes to 0-based indices, handling NA values
          # R factor codes are 1-based; h5ad uses 0-based with -1 for missing
          raw_codes <- src[[i]][['values']][index]
          codes_values <- ifelse(is.na(raw_codes), -1L, raw_codes - 1L)
          n_categories <- length(src[[i]][['levels']][])
          has_na <- any(is.na(raw_codes))

          # Choose appropriate integer dtype based on number of categories
          # Must use signed integers when NA values present (need -1)
          codes_dtype <- if (has_na || n_categories <= 127) {
            hdf5r::h5types$H5T_NATIVE_INT8
          } else if (n_categories <= 255) {
            hdf5r::h5types$H5T_NATIVE_UINT8
          } else if (n_categories <= 32767) {
            hdf5r::h5types$H5T_NATIVE_INT16
          } else if (n_categories <= 65535) {
            hdf5r::h5types$H5T_NATIVE_UINT16
          } else {
            hdf5r::h5types$H5T_NATIVE_INT32
          }

          # Add codes dataset with proper dtype for scanpy compatibility
          dfile[[dname]][[mapped_i]]$create_dataset(
            name = 'codes',
            robj = codes_values,
            dtype = codes_dtype
          )

          # Add encoding attributes to codes
          dfile[[dname]][[mapped_i]][['codes']]$create_attr(
            attr_name = 'encoding-type',
            robj = 'array',
            dtype = GuessDType(x = 'array'),
            space = H5S$new(type = "scalar")
          )
          dfile[[dname]][[mapped_i]][['codes']]$create_attr(
            attr_name = 'encoding-version',
            robj = '0.2.0',
            dtype = GuessDType(x = '0.2.0'),
            space = H5S$new(type = "scalar")
          )

          # Add categories dataset
          dfile[[dname]][[mapped_i]]$create_dataset(
            name = 'categories',
            robj = as.character(src[[i]][['levels']][]),
            dtype = StringType('utf8')
          )

          # Add encoding attributes to categories
          AddAnndataEncoding(dfile[[dname]][[mapped_i]][['categories']], encoding_type = 'string-array')

          # Add attributes to the categorical group
          AddAnndataEncoding(dfile[[dname]][[mapped_i]], encoding_type = 'categorical')
          dfile[[dname]][[mapped_i]]$create_attr(
            attr_name = 'ordered',
            robj = FALSE,
            dtype = GuessDType(x = FALSE),
            space = H5S$new(type = "scalar")
          )
        } else {
          # Regular dataset
          dfile[[dname]]$create_dataset(
            name = mapped_i,
            robj = src[[i]][index],
            dtype = src[[i]]$get_type()
          )

          # Add encoding attributes
          dfile[[dname]][[mapped_i]]$create_attr(
            attr_name = 'encoding-type',
            robj = 'array',
            dtype = GuessDType(x = 'array'),
            space = H5S$new(type = "scalar")
          )
          dfile[[dname]][[mapped_i]]$create_attr(
            attr_name = 'encoding-version',
            robj = '0.2.0',
            dtype = GuessDType(x = '0.2.0'),
            space = H5S$new(type = "scalar")
          )
        }
      }
      # Add column order with mapped names
      if (src$attr_exists(attr_name = 'colnames')) {
        original_colnames <- h5attr(x = src, which = 'colnames')
        if (standardize) {
          # Apply mapping to column names
          mapped_colnames <- sapply(original_colnames, function(cn) apply_name_map(cn, dname))
        } else {
          mapped_colnames <- original_colnames
        }
        dfile[[dname]]$create_attr(
          attr_name = 'column-order',
          robj = mapped_colnames,
          dtype = GuessDType(x = mapped_colnames)
        )
      }
      encoding.info <- c('type' = 'dataframe', 'version' = '0.2.0')
      names(x = encoding.info) <- paste0('encoding-', names(x = encoding.info))
      for (i in seq_along(along.with = encoding.info)) {
        attr.name <- names(x = encoding.info)[i]
        attr.value <- encoding.info[i]
        if (dfile[[dname]]$attr_exists(attr_name = attr.name)) {
          dfile[[dname]]$attr_delete(attr_name = attr.name)
        }
        dfile[[dname]]$create_attr(
          attr_name = attr.name,
          robj = attr.value,
          dtype = GuessDType(x = attr.value),
          space = H5S$new(type = "scalar")
        )
      }
    }
    return(invisible(x = NULL))
  }
  # Because AnnData can't figure out that sparse matrices are stored as groups
  AddEncoding <- function(dname) {
    encoding.info <- c('type' = 'csr_matrix', 'version' = '0.1.0')
    names(x = encoding.info) <- paste0('encoding-', names(x = encoding.info))
    if (inherits(x = dfile[[dname]], what = 'H5Group')) {
      # Add encoding type and version
      for (i in seq_along(along.with = encoding.info)) {
        attr.name <- names(x = encoding.info)[i]
        attr.value <- encoding.info[i]
        if (dfile[[dname]]$attr_exists(attr_name = attr.name)) {
          dfile[[dname]]$attr_delete(attr_name = attr.name)
        }
        # Flush to clear hdf5r's attribute cache after delete
        dfile$flush()
        dfile[[dname]]$create_attr(
          attr_name = attr.name,
          robj = attr.value,
          dtype = GuessDType(x = attr.value),
          space = H5S$new(type = "scalar")
        )
      }

      # Recompute shape for AnnData convention (n_obs, n_vars).
      # Seurat stores CSC matrices (indptr = column pointers for cells),
      # which we reinterpret as CSR for h5ad (indptr = row pointers for obs).
      if (dfile[[dname]]$exists(name = 'indptr') &&
          dfile[[dname]]$exists(name = 'indices')) {
        indptr_len <- prod(dfile[[dname]][['indptr']]$dims)
        nrows <- as.integer(indptr_len - 1L)

        all_indices <- dfile[[dname]][['indices']][]
        ncols <- if (length(all_indices) > 0) {
          as.integer(max(all_indices) + 1L)
        } else {
          0L
        }

        shape <- c(nrows, ncols)
        if (dfile[[dname]]$attr_exists(attr_name = 'shape')) {
          dfile[[dname]]$attr_delete(attr_name = 'shape')
        }
        # Flush to clear hdf5r's attribute cache after delete
        dfile$flush()
        dfile[[dname]]$create_attr(
          attr_name = 'shape',
          robj = shape,
          dtype = GuessDType(x = shape)
        )
      }
    }
    return(invisible(x = NULL))
  }
  # Add assay data
  assay.group <- source[['assays']][[assay]]

  # Check if this is a V5 h5Seurat file with layers structure
  has_layers <- assay.group$exists(name = 'layers') && inherits(assay.group[['layers']], 'H5Group')

  # Determine the appropriate data paths based on structure
  if (has_layers) {
    # V5 structure: data is under layers/
    # Prioritize layers/data (all genes) over scale.data (variable features only)
    if (assay.group[['layers']]$exists(name = 'data')) {
      x.data <- 'layers/data'
      raw.data <- if (assay.group[['layers']]$exists(name = 'counts')) {
        'layers/counts'
      } else {
        NULL
      }
    } else if (assay.group[['layers']]$exists(name = 'counts')) {
      x.data <- 'layers/counts'
      raw.data <- NULL
    } else if (assay.group$exists(name = 'scale.data')) {
      # Fallback: only scale.data available (unusual)
      x.data <- 'scale.data'
      raw.data <- NULL
    } else {
      stop("Cannot find data or counts in V5 h5Seurat file", call. = FALSE)
    }
  } else {
    # Legacy structure
    # Prioritize data (all genes) over scale.data (variable features only)
    if (source$index()[[assay]]$slots[['data']]) {
      x.data <- 'data'
      raw.data <- if (source$index()[[assay]]$slots[['counts']]) {
        'counts'
      } else {
        NULL
      }
    } else if (source$index()[[assay]]$slots[['counts']]) {
      x.data <- 'counts'
      raw.data <- NULL
    } else if (source$index()[[assay]]$slots[['scale.data']]) {
      # Fallback: only scale.data available (unusual)
      x.data <- 'scale.data'
      raw.data <- NULL
    } else {
      stop("Cannot find data or counts in h5Seurat file", call. = FALSE)
    }
  }

  if (verbose) {
    message("Adding ", x.data, " from ", assay, " as X")
  }
  # Get shape from source BEFORE copying (most reliable)
  source_dims <- NULL
  if (grepl("/", x.data)) {
    # Handle nested paths like "layers/data"
    parts <- strsplit(x.data, "/")[[1]]
    src_obj <- assay.group
    for (part in parts) {
      src_obj <- src_obj[[part]]
    }
    if (src_obj$attr_exists(attr_name = 'dims')) {
      source_dims <- h5attr(x = src_obj, which = 'dims')
    }
  } else {
    # Simple path
    if (assay.group$exists(name = x.data) &&
        assay.group[[x.data]]$attr_exists(attr_name = 'dims')) {
      source_dims <- h5attr(x = assay.group[[x.data]], which = 'dims')
    }
  }

  # Copy the data
  CopyH5MatrixData(
    src_group = assay.group,
    src_path = x.data,
    dst_loc = dfile,
    dst_name = 'X',
    verbose = verbose
  )

  # Ensure shape attribute exists (required by anndata)
  if (!dfile[['X']]$attr_exists(attr_name = 'shape')) {
    # Try to get from copied dims attribute first
    if (dfile[['X']]$attr_exists(attr_name = 'dims')) {
      dims <- h5attr(x = dfile[['X']], which = 'dims')
      dfile[['X']]$create_attr(
        attr_name = 'shape',
        robj = rev(x = dims),
        dtype = GuessDType(x = dims)
      )
      dfile[['X']]$attr_delete(attr_name = 'dims')
    } else if (!is.null(source_dims)) {
      # Use dims from source
      dfile[['X']]$create_attr(
        attr_name = 'shape',
        robj = rev(x = source_dims),
        dtype = GuessDType(x = source_dims)
      )
    }
  }

  AddEncoding(dname = 'X')
  # Determine n_features from the SOURCE metadata (reliable).
  # Note: reading shape from dfile[['X']] after AddEncoding is unreliable
  # due to hdf5r attribute caching after attr_delete + create_attr.
  if (has_layers && assay.group$exists(name = 'meta.data/_index')) {
    n_features <- length(SafeH5DRead(assay.group[['meta.data/_index']]))
  } else if (assay.group$exists(name = 'features')) {
    n_features <- prod(assay.group[['features']]$dims)
  } else {
    # Last resort: close and reopen dfile to clear cache, then read shape
    dfile$close_all()
    dfile <- H5File$new(filename = dest, mode = 'r+')
    x_shape <- h5attr(x = dfile[['X']], which = 'shape')
    n_features <- x_shape[2]
  }

  x.features <- switch(
    EXPR = x.data,
    'scale.data' = {
      if (has_layers && assay.group$exists(name = 'meta.data/_index')) {
        # V5: match meta.data/_index against scaled.features
        which(x = SafeH5DRead(assay.group[['meta.data/_index']]) %in% SafeH5DRead(assay.group[['scaled.features']]))
      } else {
        # Legacy: match features against scaled.features
        which(x = SafeH5DRead(assay.group[['features']]) %in% SafeH5DRead(assay.group[['scaled.features']]))
      }
    },
    seq.default(from = 1, to = n_features)
  )
  # Add meta.features with validation
  if (assay.group$exists(name = 'meta.features')) {
    meta.features.src <- assay.group[['meta.features']]
    if (inherits(x = meta.features.src, what = c('H5D', 'H5Group'))) {
      TransferDF(
        src = meta.features.src,
        dname = 'var',
        index = x.features
      )
    } else {
      warning("meta.features is not a valid H5D or H5Group, creating empty var", immediate. = TRUE)
      dfile$create_group(name = 'var')
    }
  } else {
    dfile$create_group(name = 'var')
  }
  # Add feature names
  if (Exists(x = dfile[['var']], name = rownames)) {
    dfile[['var']]$link_delete(name = rownames)
  }
  # Get feature names from correct location based on structure
  if (has_layers && assay.group$exists(name = 'meta.data/_index')) {
    # V5 structure: feature names are in meta.data/_index
    features_data <- SafeH5DRead(assay.group[['meta.data/_index']])
  } else {
    # Legacy structure: feature names are in features dataset
    features_data <- SafeH5DRead(assay.group[['features']])
  }
  # Ensure features are character strings
  features_subset <- as.character(features_data[x.features])
  dfile[['var']]$create_dataset(
    name = rownames,
    robj = features_subset,
    dtype = StringType('utf8')
  )

  # Add encoding attributes to _index
  dfile[['var']][[rownames]]$create_attr(
    attr_name = 'encoding-type',
    robj = 'string-array',
    dtype = GuessDType(x = 'string-array'),
    space = H5S$new(type = "scalar")
  )
  dfile[['var']][[rownames]]$create_attr(
    attr_name = 'encoding-version',
    robj = '0.2.0',
    dtype = GuessDType(x = '0.2.0'),
    space = H5S$new(type = "scalar")
  )

  # Add _index attribute pointing to itself
  dfile[['var']]$create_attr(
    attr_name = rownames,
    robj = rownames,
    dtype = GuessDType(x = rownames),
    space = H5S$new(type = "scalar")
  )

  # Add highly variable genes information if available
  if (assay.group$exists(name = 'variable.features')) {
    if (verbose) {
      message("Adding highly variable gene information")
    }
    variable.features <- assay.group[['variable.features']][]
    all.features <- dfile[['var']][[rownames]][]

    # Create boolean array for highly variable status
    is.variable <- all.features %in% variable.features

    dfile[['var']]$create_dataset(
      name = 'highly_variable',
      robj = is.variable,
      dtype = GuessDType(x = is.variable)
    )

    # Add number of variable features to uns
    if (!dfile$exists(name = 'uns')) {
      dfile$create_group(name = 'uns')
    }
    dfile[['uns']]$create_dataset(
      name = 'n_variable_features',
      robj = sum(is.variable),
      dtype = GuessDType(x = sum(is.variable))
    )
  }
  # Because AnnData requries meta.features and can't build an empty data frame
  if (!dfile[['var']]$attr_exists(attr_name = 'column-order')) {
    var.cols <- setdiff(
      x = names(x = dfile[['var']]),
      y = c(rownames, '__categories')
    )
    if (!length(x = var.cols)) {
      var.cols <- 'features'
      dfile[['var']]$obj_copy_to(
        dst_loc = dfile[['var']],
        dst_name = var.cols,
        src_name = rownames
      )
    }
    dfile[['var']]$create_attr(
      attr_name = 'column-order',
      robj = var.cols,
      dtype = GuessDType(x = var.cols)
    )
  }
  
  # Add encoding, to ensure compatibility with python's anndata > 0.8.0:
  encoding.info <- c('type' = 'dataframe', 'version' = '0.2.0')
  names(x = encoding.info) <- paste0('encoding-', names(x = encoding.info))
  for (i in seq_along(along.with = encoding.info)) {
    attr.name <- names(x = encoding.info)[i]
    attr.value <- encoding.info[i]
    if (dfile[['var']]$attr_exists(attr_name = attr.name)) {
      dfile[['var']]$attr_delete(attr_name = attr.name)
    }
    dfile[['var']]$create_attr(
      attr_name = attr.name,
      robj = attr.value,
      dtype = GuessDType(x = attr.value),
      space = H5S$new(type = "scalar")
    )
  }
  
  # Add raw
  if (!is.null(x = raw.data)) {
    if (verbose) {
      message("Adding ", raw.data, " from ", assay, " as raw")
    }

    # Get shape from source BEFORE copying (most reliable)
    raw_source_dims <- NULL
    if (grepl("/", raw.data)) {
      # Handle nested paths like "layers/counts"
      parts <- strsplit(raw.data, "/")[[1]]
      src_obj <- assay.group
      for (part in parts) {
        src_obj <- src_obj[[part]]
      }
      if (src_obj$attr_exists(attr_name = 'dims')) {
        raw_source_dims <- h5attr(x = src_obj, which = 'dims')
      }
    } else {
      # Simple path
      if (assay.group$exists(name = raw.data) &&
          assay.group[[raw.data]]$attr_exists(attr_name = 'dims')) {
        raw_source_dims <- h5attr(x = assay.group[[raw.data]], which = 'dims')
      }
    }

    dfile$create_group(name = 'raw')
    CopyH5MatrixData(
      src_group = assay.group,
      src_path = raw.data,
      dst_loc = dfile[['raw']],
      dst_name = 'X',
      verbose = verbose
    )

    # Ensure shape attribute exists for raw/X (required by anndata)
    if (!dfile[['raw/X']]$attr_exists(attr_name = 'shape')) {
      # Try to get from copied dims attribute first
      if (dfile[['raw/X']]$attr_exists(attr_name = 'dims')) {
        dims <- h5attr(x = dfile[['raw/X']], which = 'dims')
        dfile[['raw/X']]$create_attr(
          attr_name = 'shape',
          robj = rev(x = dims),
          dtype = GuessDType(x = dims)
        )
        dfile[['raw/X']]$attr_delete(attr_name = 'dims')
      } else if (!is.null(raw_source_dims)) {
        # Use dims from source
        dfile[['raw/X']]$create_attr(
          attr_name = 'shape',
          robj = rev(x = raw_source_dims),
          dtype = GuessDType(x = raw_source_dims)
        )
      }
    }

    AddEncoding(dname = 'raw/X')
    # Get number of features from the ACTUAL raw/X matrix dimensions
    raw_x_shape <- NULL
    if (dfile[['raw/X']]$attr_exists(attr_name = 'shape')) {
      raw_x_shape <- h5attr(x = dfile[['raw/X']], which = 'shape')
    } else if (dfile[['raw/X']]$attr_exists(attr_name = 'dims')) {
      raw_x_shape <- rev(h5attr(x = dfile[['raw/X']], which = 'dims'))
    }
    if (!is.null(raw_x_shape) && length(raw_x_shape) >= 2) {
      n_raw_features <- raw_x_shape[2]
    } else {
      n_raw_features <- prod(assay.group[['features']]$dims)
    }
    # Add meta.features with validation
    if (assay.group$exists(name = 'meta.features')) {
      meta.features.src <- assay.group[['meta.features']]
      if (inherits(x = meta.features.src, what = c('H5D', 'H5Group'))) {
        TransferDF(
          src = meta.features.src,
          dname = 'raw/var',
          index = seq.default(from = 1, to = n_raw_features)
        )
      } else {
        warning("meta.features is not a valid H5D or H5Group in raw, creating empty var", immediate. = TRUE)
        dfile[['raw']]$create_group(name = 'var')
      }
    } else {
      dfile[['raw']]$create_group(name = 'var')
    }
    # Add feature names
    if (Exists(x = dfile[['raw/var']], name = rownames)) {
      dfile[['raw/var']]$link_delete(name = rownames)
    }
    # Get feature names from correct location based on structure
    if (has_layers && assay.group$exists(name = 'meta.data/_index')) {
      # V5 structure: feature names are in meta.data/_index
      features_data_raw <- SafeH5DRead(assay.group[['meta.data/_index']])
    } else {
      # Legacy structure: feature names are in features dataset
      features_data_raw <- SafeH5DRead(assay.group[['features']])
    }
    # Subset to match actual raw/X dimensions
    features_data_raw <- as.character(features_data_raw[seq_len(n_raw_features)])
    dfile[['raw/var']]$create_dataset(
      name = rownames,
      robj = features_data_raw,
      dtype = GuessDType(x = features_data_raw[1])
    )
    # Add encoding attributes to raw/var/_index (required for AnnData to read gene names correctly)
    dfile[['raw/var']][[rownames]]$create_attr(
      attr_name = 'encoding-type',
      robj = 'string-array',
      dtype = GuessDType(x = 'string-array'),
      space = H5S$new(type = "scalar")
    )
    dfile[['raw/var']][[rownames]]$create_attr(
      attr_name = 'encoding-version',
      robj = '0.2.0',
      dtype = GuessDType(x = '0.2.0'),
      space = H5S$new(type = "scalar")
    )
    dfile[['raw/var']]$create_attr(
      attr_name = rownames,
      robj = rownames,
      dtype = GuessDType(x = rownames),
      space = H5S$new(type = "scalar")
    )
    # Add dataframe encoding attributes to raw/var group (required for AnnData)
    # Only add column-order if it doesn't already exist (TransferDF may have added it)
    if (!dfile[['raw/var']]$attr_exists(attr_name = 'column-order')) {
      # Get list of columns (everything except _index)
      raw_var_cols <- setdiff(dfile[['raw/var']]$names, rownames)
      if (length(raw_var_cols) > 0) {
        dfile[['raw/var']]$create_attr(
          attr_name = 'column-order',
          robj = raw_var_cols,
          dtype = GuessDType(x = raw_var_cols)
        )
      } else {
        # Empty column-order for dataframes with only _index
        dfile[['raw/var']]$create_attr(
          attr_name = 'column-order',
          robj = character(0),
          dtype = StringType('utf8')
        )
      }
    }
    # Add encoding-type and encoding-version
    encoding.info <- c('type' = 'dataframe', 'version' = '0.2.0')
    names(x = encoding.info) <- paste0('encoding-', names(x = encoding.info))
    for (i in seq_along(along.with = encoding.info)) {
      attr.name <- names(x = encoding.info)[i]
      attr.value <- encoding.info[i]
      if (dfile[['raw/var']]$attr_exists(attr_name = attr.name)) {
        dfile[['raw/var']]$attr_delete(attr_name = attr.name)
      }
      dfile[['raw/var']]$create_attr(
        attr_name = attr.name,
        robj = attr.value,
        dtype = GuessDType(x = attr.value),
        space = H5S$new(type = "scalar")
      )
    }
  }
  # Add cell-level metadata with validation
  if (source$exists(name = 'meta.data')) {
    meta.src <- source[['meta.data']]

    # First attempt normal transfer
    TransferDF(
      src = meta.src,
      dname = 'obs',
      index = seq.default(from = 1, to = length(x = Cells(x = source)))
    )

    # If obs was not created successfully, try direct extraction
    if (!dfile$exists(name = 'obs')) {
      if (verbose) {
        message("Attempting direct metadata extraction for V5 compatibility")
      }

      # Create obs group
      dfile$create_group(name = 'obs')

      # Try to read metadata directly from h5Seurat structure
      # Unwrap V5 environment wrappers if present
      meta.src.unwrapped <- meta.src
      if (is.environment(x = meta.src.unwrapped) && exists("hgroup", envir = meta.src.unwrapped, inherits = FALSE)) {
        meta.src.unwrapped <- get("hgroup", envir = meta.src.unwrapped, inherits = FALSE)
      }
      if (inherits(x = meta.src.unwrapped, what = 'H5Group')) {
        # Get all column names except internal ones
        meta_cols <- setdiff(names(x = meta.src.unwrapped), c('__categories', '_index'))

        for (col in meta_cols) {
          tryCatch({
            if (IsFactor(x = meta.src.unwrapped[[col]])) {
              # Use newer anndata format: each categorical is a group with categories/codes
              dfile[['obs']]$create_group(name = col)

              # Convert codes to 0-based indices, handling NA values
              # R factor codes are 1-based; h5ad uses 0-based with -1 for missing
              raw_codes <- meta.src.unwrapped[[col]][['values']][]
              codes_values <- ifelse(is.na(raw_codes), -1L, raw_codes - 1L)
              n_categories <- length(meta.src.unwrapped[[col]][['levels']][])
              has_na <- any(is.na(raw_codes))

              # Choose appropriate integer dtype based on number of categories
              # Must use signed integers when NA values present (need -1)
              codes_dtype <- if (has_na || n_categories <= 127) {
                hdf5r::h5types$H5T_NATIVE_INT8
              } else if (n_categories <= 255) {
                hdf5r::h5types$H5T_NATIVE_UINT8
              } else if (n_categories <= 32767) {
                hdf5r::h5types$H5T_NATIVE_INT16
              } else if (n_categories <= 65535) {
                hdf5r::h5types$H5T_NATIVE_UINT16
              } else {
                hdf5r::h5types$H5T_NATIVE_INT32
              }

              # Add codes dataset with proper dtype for scanpy compatibility
              dfile[['obs']][[col]]$create_dataset(
                name = 'codes',
                robj = codes_values,
                dtype = codes_dtype
              )

              # Add encoding attributes to codes
              dfile[['obs']][[col]][['codes']]$create_attr(
                attr_name = 'encoding-type',
                robj = 'array',
                dtype = GuessDType(x = 'array'),
                space = H5S$new(type = "scalar")
              )
              dfile[['obs']][[col]][['codes']]$create_attr(
                attr_name = 'encoding-version',
                robj = '0.2.0',
                dtype = GuessDType(x = '0.2.0'),
                space = H5S$new(type = "scalar")
              )

              # Add categories dataset
              dfile[['obs']][[col]]$create_dataset(
                name = 'categories',
                robj = as.character(meta.src.unwrapped[[col]][['levels']][]),
                dtype = StringType('utf8')
              )

              # Add encoding attributes to categories
              dfile[['obs']][[col]][['categories']]$create_attr(
                attr_name = 'encoding-type',
                robj = 'string-array',
                dtype = GuessDType(x = 'string-array'),
                space = H5S$new(type = "scalar")
              )
              dfile[['obs']][[col]][['categories']]$create_attr(
                attr_name = 'encoding-version',
                robj = '0.2.0',
                dtype = GuessDType(x = '0.2.0'),
                space = H5S$new(type = "scalar")
              )

              # Add attributes to the categorical group
              dfile[['obs']][[col]]$create_attr(
                attr_name = 'encoding-type',
                robj = 'categorical',
                dtype = GuessDType(x = 'categorical'),
                space = H5S$new(type = "scalar")
              )
              dfile[['obs']][[col]]$create_attr(
                attr_name = 'encoding-version',
                robj = '0.2.0',
                dtype = GuessDType(x = '0.2.0'),
                space = H5S$new(type = "scalar")
              )
              dfile[['obs']][[col]]$create_attr(
                attr_name = 'ordered',
                robj = FALSE,
                dtype = GuessDType(x = FALSE),
                space = H5S$new(type = "scalar")
              )
            } else {
              # Handle regular columns
              dfile[['obs']]$create_dataset(
                name = col,
                robj = meta.src.unwrapped[[col]][],
                dtype = meta.src.unwrapped[[col]]$get_type()
              )

              # Add encoding attributes
              dfile[['obs']][[col]]$create_attr(
                attr_name = 'encoding-type',
                robj = 'array',
                dtype = GuessDType(x = 'array'),
                space = H5S$new(type = "scalar")
              )
              dfile[['obs']][[col]]$create_attr(
                attr_name = 'encoding-version',
                robj = '0.2.0',
                dtype = GuessDType(x = '0.2.0'),
                space = H5S$new(type = "scalar")
              )
            }
          }, error = function(e) {
            if (verbose) {
              message("Could not transfer metadata column ", col, ": ", e$message)
            }
          })
        }

        # Add column-order attribute
        if (length(meta_cols) > 0) {
          if (!dfile[['obs']]$attr_exists(attr_name = 'column-order')) {
            dfile[['obs']]$create_attr(
              attr_name = 'column-order',
              robj = meta_cols,
              dtype = GuessDType(x = meta_cols)
            )
          }
        }

        # Add encoding attributes
        encoding.info <- c('type' = 'dataframe', 'version' = '0.1.0')
        names(x = encoding.info) <- paste0('encoding-', names(x = encoding.info))
        for (i in seq_along(along.with = encoding.info)) {
          attr.name <- names(x = encoding.info)[i]
          if (!dfile[['obs']]$attr_exists(attr_name = attr.name)) {
            dfile[['obs']]$create_attr(
              attr_name = attr.name,
              robj = encoding.info[i],
              dtype = GuessDType(x = encoding.info[i]),
              space = H5S$new(type = "scalar")
            )
          }
        }
      } else {
        warning("meta.data structure not recognized, creating minimal obs", immediate. = TRUE)
      }
    }
  } else {
    # Create empty obs group if meta.data doesn't exist
    dfile$create_group(name = 'obs')
  }

  # Add cell names if obs was created
  if (dfile$exists(name = 'obs')) {
    if (Exists(x = dfile[['obs']], name = rownames)) {
      dfile[['obs']]$link_delete(name = rownames)
    }
    cell_names <- as.character(Cells(x = source))
    dfile[['obs']]$create_dataset(
      name = rownames,
      robj = cell_names,
      dtype = StringType('utf8')
    )

    # Add encoding attributes to _index
    dfile[['obs']][[rownames]]$create_attr(
      attr_name = 'encoding-type',
      robj = 'string-array',
      dtype = GuessDType(x = 'string-array'),
      space = H5S$new(type = "scalar")
    )
    dfile[['obs']][[rownames]]$create_attr(
      attr_name = 'encoding-version',
      robj = '0.2.0',
      dtype = GuessDType(x = '0.2.0'),
      space = H5S$new(type = "scalar")
    )

    # Add _index attribute pointing to itself
    dfile[['obs']]$create_attr(
      attr_name = rownames,
      robj = rownames,
      dtype = GuessDType(x = rownames),
      space = H5S$new(type = "scalar")
    )

    # Ensure encoding attributes exist for anndata compatibility
    if (!dfile[['obs']]$attr_exists(attr_name = 'encoding-type')) {
      encoding.info <- c('type' = 'dataframe', 'version' = '0.2.0')
      names(x = encoding.info) <- paste0('encoding-', names(x = encoding.info))
      for (i in seq_along(along.with = encoding.info)) {
        attr.name <- names(x = encoding.info)[i]
        if (!dfile[['obs']]$attr_exists(attr_name = attr.name)) {
          dfile[['obs']]$create_attr(
            attr_name = attr.name,
            robj = encoding.info[i],
            dtype = GuessDType(x = encoding.info[i]),
            space = H5S$new(type = "scalar")
          )
        }
      }
    }
  }
  # Add dimensional reduction information
  if (!dfile$exists(name = 'obsm')) {
    obsm <- dfile$create_group(name = 'obsm')
  } else {
    obsm <- dfile[['obsm']]
  }
  if (!dfile$exists(name = 'varm')) {
    varm <- dfile$create_group(name = 'varm')
  } else {
    varm <- dfile[['varm']]
  }
  reductions <- source$index()[[assay]]$reductions
  for (reduc in names(x = reductions)) {
    if (verbose) {
      message("Adding dimensional reduction information for ", reduc)
    }
    Transpose(
      x = source[[H5Path('reductions', reduc, 'cell.embeddings')]],
      dest = obsm,
      dname = paste0('X_', reduc),
      verbose = FALSE
    )
    if (reductions[[reduc]]['feature.loadings']) {
      if (verbose) {
        message("Adding feature loadings for ", reduc)
      }
      loadings <- source[['reductions']][[reduc]][['feature.loadings']]
      reduc.features <- loadings$dims[1]
      x.features <- dfile[['var']][[rownames]]$dims
      varm.name <- switch(EXPR = reduc, 'pca' = 'PCs', toupper(x = reduc))
      # Because apparently AnnData requires nPCs == nrow(X)
      if (reduc.features < x.features) {
        pad <- paste0('pad_', varm.name)
        # Get feature indices and filter out NAs
        feature_indices <- match(
          x = SafeH5DRead(source[['reductions']][[reduc]][['features']]),
          table = SafeH5DRead(dfile[['var']][[rownames]])
        )
        # Replace NA values with new rows (append to end)
        na_indices <- which(is.na(feature_indices))
        if (length(na_indices) > 0) {
          # For NA indices, assign new row positions at the end
          feature_indices[na_indices] <- (reduc.features + 1L):(reduc.features + length(na_indices))
        }

        PadMatrix(
          src = loadings,
          dest = dfile[['varm']],
          dname = pad,
          dims = c(x.features, loadings$dims[2]),
          index = list(
            feature_indices,
            seq.default(from = 1, to = loadings$dims[2])
          )
        )
        loadings <- dfile[['varm']][[pad]]
      }
      Transpose(x = loadings, dest = varm, dname = varm.name, verbose = FALSE)
      if (reduc.features < x.features) {
        dfile$link_delete(name = loadings$get_obj_name())
      }
    }
  }
  # Add global dimensional reduction information
  global.reduc <- source$index()[['global']][['reductions']]
  for (reduc in global.reduc) {
    if (reduc %in% names(x = reductions)) {
      next
    } else if (verbose) {
      message("Adding dimensional reduction information for ", reduc, " (global)")
    }
    Transpose(
      x = source[[H5Path('reductions', reduc, 'cell.embeddings')]],
      dest = obsm,
      dname = paste0('X_', reduc),
      verbose = FALSE
    )
  }
  # Add spatial coordinates if available
  # Check for spatial data in images slot
  if (source$exists(name = 'images')) {
    cell_names <- SafeH5DRead(source[['cell.names']])
    images <- names(x = source[['images']])
    if (length(images) > 0 && verbose) {
      message("Processing spatial data from images")
    }

    # Process first image for now (extend for multiple images later)
    if (length(images) > 0) {
      img_name <- images[1]
      img_group <- source[['images']][[img_name]]

      # Check for coordinates
      spatial_matrix <- NULL
      if (img_group$exists(name = 'coordinates')) {
        if (verbose) {
          message("Adding spatial coordinates to obsm['spatial']")
        }

        # Read coordinates
        coord_group <- img_group[['coordinates']]

        # Get x and y coordinates
        if (coord_group$exists('x') && coord_group$exists('y')) {
          x_coords <- coord_group[['x']][]
          y_coords <- coord_group[['y']][]

          # Create spatial matrix (cells x 2) as [X, Y] for scanpy/squidpy
          spatial_matrix <- cbind(x_coords, y_coords)
        } else if (coord_group$exists('imagerow') && coord_group$exists('imagecol')) {
          # Visium-style coordinates: imagerow=Y, imagecol=X
          row_coords <- coord_group[['imagerow']][]
          col_coords <- coord_group[['imagecol']][]

          # Create spatial matrix (cells x 2) as [X, Y] = [col, row]
          spatial_matrix <- cbind(col_coords, row_coords)
        }
      } else if (img_group$exists(name = 'boundaries') &&
                 img_group[['boundaries']]$exists(name = 'centroids')) {
        centroid_group <- img_group[['boundaries/centroids']]
        if (centroid_group$exists(name = 'coords')) {
          coords_mat <- centroid_group[['coords']]$read()
          row_index <- seq_len(nrow(coords_mat))
          if (centroid_group$exists(name = 'cells')) {
            centroid_cells <- centroid_group[['cells']][]
            if (!is.null(cell_names)) {
              matched <- match(cell_names, centroid_cells)
              if (any(is.na(matched))) {
                warning("Unable to align all centroid coordinates with cell names")
              } else {
                row_index <- matched
              }
            }
          }
          coords_mat <- coords_mat[row_index, , drop = FALSE]
          # coords_mat is stored as [X, Y] in Seurat h5Seurat, keep as-is for scanpy
          spatial_matrix <- coords_mat
        }
      }
      if (!is.null(spatial_matrix)) {
        spatial_matrix <- as.matrix(spatial_matrix)
        storage.mode(spatial_matrix) <- 'double'
        obsm$create_dataset(
          name = 'spatial',
          robj = spatial_matrix,
          dtype = h5types$H5T_NATIVE_DOUBLE
        )
        Transpose(
          x = obsm[['spatial']],
          dest = obsm,
          dname = 'spatial',
          overwrite = TRUE,
          verbose = FALSE
        )
        obsm[['spatial']]$create_attr(
          attr_name = 'encoding-type',
          robj = 'array',
          dtype = GuessDType(x = 'array'),
          space = H5S$new(type = "scalar")
        )
        obsm[['spatial']]$create_attr(
          attr_name = 'encoding-version',
          robj = '0.2.0',
          dtype = GuessDType(x = '0.2.0'),
          space = H5S$new(type = "scalar")
        )
      }
    }
  }

  # Create uns if it doesn't exist
  if (!dfile$exists(name = 'uns')) {
    dfile$create_group(name = 'uns')
  }

  # Preserve Seurat-specific metadata in uns['seurat']
  if (verbose) {
    message("Preserving Seurat-specific metadata in uns['seurat']")
  }
  seurat_group <- dfile[['uns']]$create_group(name = 'seurat')

  # Store assay name
  seurat_group$create_dataset(
    name = 'assay_name',
    robj = assay,
    dtype = StringType('utf8')
  )

  # Store SeuratDisk version
  seurat_group$create_dataset(
    name = 'version',
    robj = 'SeuratDisk-V5',
    dtype = StringType('utf8')
  )

  # Preserve graph names if they exist
  if (source$exists(name = 'graphs')) {
    graph_names <- names(x = source[['graphs']])
    if (length(graph_names) > 0) {
      seurat_group$create_dataset(
        name = 'graph_names',
        robj = graph_names,
        dtype = StringType('utf8')
      )
    }
  }

  # Preserve obs column names to track which metadata came from Seurat
  if (dfile$exists(name = 'obs')) {
    if (dfile[['obs']]$attr_exists(attr_name = 'column-order')) {
      obs_cols <- h5attr(x = dfile[['obs']], which = 'column-order')
      seurat_group$create_dataset(
        name = 'obs_columns',
        robj = obs_cols,
        dtype = StringType('utf8')
      )
    }
  }

  # Add spatial metadata to uns if available
  if (source$exists(name = 'images')) {
    images <- names(x = source[['images']])

    if (length(images) > 0) {
      if (verbose) {
        message("Adding spatial metadata to uns['spatial']")
      }

      # Create spatial group in uns
      spatial_uns <- dfile[['uns']]$create_group(name = 'spatial')

      # Each Seurat image becomes a separate library in h5ad format
      # This matches the multi-library structure expected by scanpy/squidpy
      for (img_name in images) {
        img_group <- source[['images']][[img_name]]

        # Use the Seurat image name as the library ID (e.g., 'anterior1')
        lib_id <- img_name
        lib_group <- spatial_uns$create_group(name = lib_id)

        # Add scale factors for this library
        if (img_group$exists(name = 'scale.factors')) {
          sf_group <- lib_group$create_group(name = 'scalefactors')
          sf_src <- img_group[['scale.factors']]
          sf_map <- c(
            hires = 'tissue_hires_scalef',
            lowres = 'tissue_lowres_scalef',
            spot = 'spot_diameter_fullres',
            fiducial = 'fiducial_diameter_fullres'
          )
          for (sf_name in sf_src$names) {
            dst_name <- sf_map[[sf_name]]
            if (is.null(dst_name)) {
              dst_name <- sf_name
            }
            sf_value <- sf_src[[sf_name]][]
            sf_group$create_dataset(
              name = dst_name,
              robj = sf_value,
              dtype = h5types$H5T_NATIVE_DOUBLE,
              space = H5S$new(type = "scalar"),
              chunk_dims = NULL
            )
          }
        }

        # Add metadata
        meta_group <- lib_group$create_group(name = 'metadata')

        # Create images group with standard h5ad keys ('hires', 'lowres')
        images_group <- lib_group$create_group(name = 'images')

        # Add primary image as 'lowres' and 'hires'
        if (img_group$exists(name = 'image')) {
          img_data <- img_group[['image']]$read()
          if (length(dim(img_data)) == 3L) {
            img_arr <- aperm(img_data, c(3L, 2L, 1L))

            # Check for high-resolution image in attributes
            has_hires <- FALSE
            if (!is.null(img_group[['image']]$attr_names) &&
                'hires.image' %in% img_group[['image']]$attr_names) {
              hires_data <- h5attr(x = img_group[['image']], which = 'hires.image')
              if (!is.null(hires_data) && length(dim(hires_data)) == 3L) {
                # hires.image is already in correct dimensions (channels x width x height)
                images_group$create_dataset(
                  name = 'hires',
                  robj = hires_data,
                  dtype = h5types$H5T_NATIVE_DOUBLE
                )
                has_hires <- TRUE
              }
            }

            # If no separate hires image, use the primary image for both hires and lowres
            # This ensures compatibility with squidpy which defaults to 'hires'
            if (!has_hires) {
              images_group$create_dataset(
                name = 'hires',
                robj = img_arr,
                dtype = h5types$H5T_NATIVE_DOUBLE
              )
            }

            # Always create lowres
            images_group$create_dataset(
              name = 'lowres',
              robj = img_arr,
              dtype = h5types$H5T_NATIVE_DOUBLE
            )
          }
        }
      }
    }
  }

  # Add graphs to obsp (pairwise observations)
  graphs.available <- source$index()[[assay]]$graphs
  if (length(graphs.available) > 0 && source$exists('graphs')) {
    if (verbose) {
      message("Adding graph information to obsp")
    }
    if (!dfile$exists(name = 'obsp')) {
      dfile$create_group(name = 'obsp')
    }

    # Transfer all available graphs
    for (graph.name in graphs.available) {
      if (source[['graphs']]$exists(name = graph.name)) {
        # Map Seurat graph names to anndata conventions
        # e.g., RNA_nn -> connectivities, RNA_snn -> distances
        obsp.name <- if (grepl("_nn$", graph.name)) {
          "connectivities"
        } else if (grepl("_snn$", graph.name)) {
          "distances"
        } else {
          gsub(paste0("^", assay, "_"), "", graph.name)
        }

        if (verbose) {
          message("  - Adding ", graph.name, " as obsp/", obsp.name)
        }

        source[['graphs']]$obj_copy_to(
          dst_loc = dfile[['obsp']],
          dst_name = obsp.name,
          src_name = graph.name
        )

        # Add shape attribute for anndata compatibility
        if (source[['graphs']][[graph.name]]$attr_exists(attr_name = 'dims')) {
          dims <- h5attr(x = source[['graphs']][[graph.name]], which = 'dims')
          dfile[['obsp']][[obsp.name]]$create_attr(
            attr_name = 'shape',
            robj = rev(x = dims),
            dtype = GuessDType(x = dims)
          )
        }

        # Add encoding type
        AddEncoding(dname = paste0('obsp/', obsp.name))
      }
    }
  }

  # Add graph (legacy - keep for backward compatibility)
  graph <- source$index()[[assay]]$graphs
  graph <- graph[length(x = graph)]
  if (!is.null(x = graph)) {
    if (verbose) {
      message("Adding ", graph, " as neighbors")
    }
    dgraph <- dfile[['uns']]$create_group(name = 'neighbors')
    source[['graphs']]$obj_copy_to(
      dst_loc = dgraph,
      dst_name = 'distances',
      src_name = graph
    )
    if (source[['graphs']][[graph]]$attr_exists(attr_name = 'dims')) {
      dims <- h5attr(x = source[['graphs']][[graph]], which = 'dims')
      dgraph[['distances']]$create_attr(
        attr_name = 'shape',
        robj = rev(x = dims),
        dtype = GuessDType(x = dims)
      )
    }
    AddEncoding(dname = 'uns/neighbors/distances')
    # Add parameters
    dgraph$create_group(name = 'params')
    dgraph[['params']]$create_dataset(
      name = 'method',
      robj = gsub(pattern = paste0('^', assay, '_'), replacement = '', x = graph),
      dtype = GuessDType(x = graph)
    )
    cmdlog <- paste(
      paste0('FindNeighbors.', assay),
      unique(x = c(names(x = reductions), source$index()$global$reductions)),
      sep = '.',
      collapse = '|'
    )
    cmdlog <- grep(
      pattern = cmdlog,
      x = names(x = source[['commands']]),
      value = TRUE
    )
    if (length(x = cmdlog) > 1) {
      timestamps <- sapply(
        X = cmdlog,
        FUN = function(cmd) {
          ts <- if (source[['commands']][[cmd]]$attr_exists(attr_name = 'time.stamp')) {
            h5attr(x = source[['commands']][[cmd]], which = 'time.stamp')
          } else {
            NULL
          }
          return(ts)
        },
        simplify = TRUE,
        USE.NAMES = FALSE
      )
      timestamps <- Filter(f = Negate(f = is.null), x = timestamps)
      cmdlog <- cmdlog[order(timestamps, decreasing = TRUE)][1]
    }
    if (length(x = cmdlog) && !is.na(x = cmdlog)) {
      cmdlog <- source[['commands']][[cmdlog]]
      if ('k.param' %in% names(x = cmdlog)) {
        dgraph[['params']]$obj_copy_from(
          src_loc = cmdlog,
          src_name = 'k.param',
          dst_name = 'n_neighbors'
        )
      }
    }
  }
  # Add layers
  other.assays <- setdiff(
    x = names(x = source$index()),
    y = c(assay, 'global', 'no.assay')
  )
  if (length(x = other.assays)) {
    x.dims <- Dims(x = dfile[['X']])
    layers <- dfile$create_group(name = 'layers')
    for (other in other.assays) {
      layer.slot <- NULL
      other.group <- source[['assays']][[other]]

      # Check if this assay has V5 structure with layers
      other.has_layers <- other.group$exists(name = 'layers') && inherits(other.group[['layers']], 'H5Group')

      for (slot in c('scale.data', 'data')) {
        # Determine the actual path for the slot
        actual_path <- if (other.has_layers && slot %in% c('data', 'counts')) {
          paste0('layers/', slot)
        } else {
          slot
        }

        # Check if the slot exists at the determined path
        slot.exists <- if (other.has_layers && slot %in% c('data', 'counts')) {
          other.group[['layers']]$exists(name = slot)
        } else {
          other.group$exists(name = slot)
        }

        if (slot.exists) {
          slot.dims <- Dims(x = other.group[[actual_path]])
          if (isTRUE(all.equal(slot.dims, x.dims))) {
            layer.slot <- actual_path
            break
          }
        }
      }

      if (!is.null(x = layer.slot)) {
        if (verbose) {
          message("Adding ", layer.slot, " from ", other, " as a layer")
        }
        layers$obj_copy_from(
          src_loc = source[['assays']][[other]],
          src_name = layer.slot,
          dst_name = other
        )
        if (layers[[other]]$attr_exists(attr_name = 'dims')) {
          dims <- h5attr(x = layers[[other]], which = 'dims')
          layers[[other]]$create_attr(
            attr_name = 'shape',
            robj = rev(x = dims),
            dtype = GuessDType(x = dims)
          )
          layers[[other]]$attr_delete(attr_name = 'dims')
        }
        AddEncoding(dname = paste('layers', other, sep = '/'))
        layer.features <- switch(
          EXPR = layer.slot,
          'scale.data' = 'scaled.features',
          'features'
        )
        var.name <- paste0(other, '_features')
        dfile[['var']]$obj_copy_from(
          src_loc = source[['assays']][[other]],
          src_name = layer.features,
          dst_name = var.name
        )
        col.order <- h5attr(x = dfile[['var']], which = 'column-order')
        col.order <- c(col.order, var.name)
        dfile[['var']]$attr_rename(
          old_attr_name = 'column-order',
          new_attr_name = 'old-column-order'
        )
        dfile[['var']]$create_attr(
          attr_name = 'column-order',
          robj = col.order,
          dtype = GuessDType(x = col.order)
        )
        dfile[['var']]$attr_delete(attr_name = 'old-column-order')
      }
    }
  }

  # Ensure all required anndata groups exist (even if empty)
  required_groups <- c('obsm', 'obsp', 'varm', 'varp', 'layers', 'uns')
  for (group_name in required_groups) {
    if (!dfile$exists(name = group_name)) {
      dfile$create_group(name = group_name)
    }
  }

  # Add encoding attributes to all top-level groups
  groups_to_encode <- c('obsm', 'obsp', 'varm', 'varp', 'layers', 'uns')
  for (group_name in groups_to_encode) {
    if (dfile$exists(name = group_name) && inherits(dfile[[group_name]], 'H5Group')) {
      # Add encoding-type
      if (!dfile[[group_name]]$attr_exists(attr_name = 'encoding-type')) {
        dfile[[group_name]]$create_attr(
          attr_name = 'encoding-type',
          robj = 'dict',
          dtype = GuessDType(x = 'dict'),
          space = H5S$new(type = "scalar")
        )
      }
      # Add encoding-version
      if (!dfile[[group_name]]$attr_exists(attr_name = 'encoding-version')) {
        dfile[[group_name]]$create_attr(
          attr_name = 'encoding-version',
          robj = '0.1.0',
          dtype = GuessDType(x = '0.1.0'),
          space = H5S$new(type = "scalar")
        )
      }
    }
  }

  # Add encoding attributes to dimensional reductions in obsm
  if (dfile$exists(name = 'obsm')) {
    for (reduc_name in names(dfile[['obsm']])) {
      if (!dfile[['obsm']][[reduc_name]]$attr_exists(attr_name = 'encoding-type')) {
        dfile[['obsm']][[reduc_name]]$create_attr(
          attr_name = 'encoding-type',
          robj = 'array',
          dtype = GuessDType(x = 'array'),
          space = H5S$new(type = "scalar")
        )
      }
      if (!dfile[['obsm']][[reduc_name]]$attr_exists(attr_name = 'encoding-version')) {
        dfile[['obsm']][[reduc_name]]$create_attr(
          attr_name = 'encoding-version',
          robj = '0.2.0',
          dtype = GuessDType(x = '0.2.0'),
          space = H5S$new(type = "scalar")
        )
      }
    }
  }

  dfile$flush()
  return(dfile)
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Helper Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Read observation metadata from an H5AD file
#'
#' Read the cell-level metadata (obs) from an H5AD file and return it as a
#' data.frame. Handles both categorical (factor) and standard data types.
#'
#' @param file Path to an H5AD file
#'
#' @return A \code{data.frame} with cell metadata, where row names are cell
#'   barcodes from the \code{_index} dataset and columns correspond to
#'   observation annotations
#'
#' @export
#'
readH5AD_obs <- function(file) {
  suppressWarnings(expr = hfile <- Connect(filename = file, force = TRUE))
  on.exit(hfile$close_all())
  hfile_obs <- hfile[['obs']]
  obs_groups <- setdiff(names(hfile_obs), c('__categories', '_index'))
  matrix <- as.data.frame(
    x = matrix(data = NA,
               nrow = hfile_obs[['_index']]$dims[1],
               ncol = length(obs_groups))
    )
  colnames(matrix) <- obs_groups
  rownames(matrix) <- hfile_obs[['_index']][]
  if ('__categories' %in% names(x = hfile_obs)) {
    # Legacy h5ad format: __categories group with per-column category lists
    hfile_cate <- hfile_obs[['__categories']]
    for (i in seq_along(obs_groups)) {
      obs.i <- obs_groups[i]
      obs_value_i <- hfile_obs[[obs.i]][]
      if (obs.i %in% names(x = hfile_cate)) {
        categories <- as.character(hfile_cate[[obs.i]][])
        codes <- obs_value_i
        codes[codes == -1L] <- NA_integer_
        obs_value_i <- factor(categories[codes + 1L], levels = categories)
      }
      matrix[, i] <- obs_value_i
    }
  } else {
    # Modern h5ad format: categorical groups with categories/codes sub-datasets
    for (i in seq_along(obs_groups)) {
      obs.i <- obs_groups[i]
      col_obj <- hfile_obs[[obs.i]]
      if (inherits(col_obj, 'H5Group') &&
          all(c("categories", "codes") %in% names(col_obj))) {
        categories <- as.character(col_obj[['categories']]$read())
        codes <- col_obj[['codes']]$read()
        codes[codes == -1L] <- NA_integer_
        obs_value_i <- factor(categories[codes + 1L], levels = categories)
      } else if (inherits(col_obj, 'H5D')) {
        obs_value_i <- col_obj$read()
      } else {
        obs_value_i <- NA
      }
      matrix[, i] <- obs_value_i
    }
  }
  return(matrix)
}

#' Read observation embeddings from an H5AD file
#'
#' Read the cell embeddings (obsm) from an H5AD file and return them as a
#' named list of matrices. Each entry corresponds to a dimensional reduction
#' (e.g., PCA, UMAP, t-SNE).
#'
#' @param file Path to an H5AD file
#'
#' @return A named list of matrices, where each matrix has cells as rows and
#'   embedding dimensions as columns. Names are derived from the obsm keys
#'   with the \dQuote{X_} prefix removed.
#'
#' @export
#'
readH5AD_obsm <- function(file) {
  suppressWarnings(hfile <- Connect(filename = file, force = TRUE))
  on.exit(hfile$close_all())
  hfile_obsm <- hfile[['obsm']]
  if (length(names(hfile_obsm)) == 0) {
    message('No obsm found in this object')
    return(list())
  }
  obsm_set <- names(hfile_obsm)
  cells.name <- hfile[['obs']][['_index']][]
  n_cells <- length(cells.name)
  obsm.list <- lapply(obsm_set, function(x) {
    emb <- hfile_obsm[[x]][,]
    # hdf5r may return transposed (n_dims x n_cells) due to row/column major difference
    if (nrow(emb) != n_cells && ncol(emb) == n_cells) {
      emb <- t(emb)
    }
    rownames(emb) <- cells.name
    key.name <- gsub('X_', '', x)
    colnames(emb) <- paste0(key.name, "_", seq_len(ncol(emb)))
    return(emb)
  })
  names(obsm.list) <- gsub('X_', '', obsm_set)
  return(obsm.list)
}

#' Direct Seurat to H5AD Conversion
#'
#' Converts a Seurat object directly to an H5AD file, handling the intermediate
#' h5Seurat file automatically. The intermediate file is created in a temporary
#' location and removed after conversion.
#'
#' @param object A Seurat object to convert
#' @param filename Output H5AD filename. If not specified, uses the project name
#'   with .h5ad extension.
#' @param assay Name of assay to convert. Default is the default assay.
#' @param overwrite Logical; overwrite existing file. Default FALSE.
#' @param verbose Logical; show progress messages. Default TRUE.
#' @param standardize Logical; convert Seurat metadata names to scanpy conventions.
#'   Default FALSE.
#' @param ... Additional arguments passed to SaveH5Seurat and Convert.
#'
#' @return Invisibly returns the path to the created H5AD file.
#'
#' @details
#' This function provides a convenient one-step conversion from Seurat objects
#' to H5AD format (used by Python's scanpy/anndata). Internally, it:
#' \enumerate{
#'   \item Saves the Seurat object to a temporary h5Seurat file
#'   \item Converts the h5Seurat file to H5AD format
#'   \item Removes the intermediate h5Seurat file
#' }
#'
#' This is useful when you want to export data for Python analysis without
#' keeping the intermediate h5Seurat file.
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(srtdisk)
#'
#' # Load or create a Seurat object
#' pbmc <- pbmc_small
#'
#' # Convert directly to H5AD
#' SeuratToH5AD(pbmc, filename = "pbmc.h5ad")
#'
#' # The file can now be loaded in Python:
#' # import scanpy as sc
#' # adata = sc.read_h5ad("pbmc.h5ad")
#' }
#'
#' @seealso \code{\link{SaveH5Seurat}}, \code{\link{Convert}}, \code{\link{LoadH5AD}}
#'
#' @export
#'
SeuratToH5AD <- function(
  object,
  filename = NULL,
  assay = DefaultAssay(object = object),
  overwrite = FALSE,
  verbose = TRUE,
  standardize = FALSE,
  ...
) {
  if (!inherits(x = object, what = 'Seurat')) {
    stop("'object' must be a Seurat object", call. = FALSE)
  }
  if (is.null(x = filename)) {
    filename <- paste0(Project(object = object), '.h5ad')
  }
  if (!grepl(pattern = '\\.h5ad$', x = filename, ignore.case = TRUE)) {
    filename <- paste0(filename, '.h5ad')
  }
  Convert(source = object, dest = filename, assay = assay,
          overwrite = overwrite, verbose = verbose, standardize = standardize, ...)
}

