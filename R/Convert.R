#' @include zzz.R
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

#' Convert an on-disk single-cell dataset to another format
#'
#' HDF5-based single-cell datasets can be converted from one format to another
#' using minimal memory. Details about conversion formats implemented are
#' provided below
#'
#' @inheritSection H5ADToH5Seurat AnnData/H5AD to h5Seurat
#' @inheritSection H5SeuratToH5AD h5Seurat to AnnData/H5AD
#'
#' @param source Source dataset
#' @param dest Name of destination dataset
#' @param assay Converting from \code{\link{h5Seurat}}: name of assay to write
#' out; converting to \code{\link{h5Seurat}}: name to store assay data as
#' @param overwrite Overwrite existing \code{dest}
#' @param verbose Show progress updates
#' @param ... Arguments passed to other methods
#'
#' @return If \code{source} is a \code{character}, invisibly returns
#' \code{dest}; otherwise, returns an \code{\link[hdf5r]{H5File}}, or
#' filetype-specific subclass of \code{H5File} (eg. \code{\link{h5Seurat}}),
#' connection to \code{dest}
#'
#' @name Convert
#' @rdname Convert
#'
#' @export
#'
Convert <- function(source, dest, assay, overwrite = FALSE, verbose = TRUE, ...) {
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
  ...
) {
  hfile <- Connect(filename = source, force = TRUE)
  if (missing(x = assay)) {
    assay <- tryCatch(
      expr = DefaultAssay(object = hfile),
      error = function(...) {
        warning(
          "'assay' not set, setting to 'RNA'",
          call. = FALSE,
          immediate. = TRUE
        )
        "RNA"
      }
    )
  }
  on.exit(expr = hfile$close_all())
  dfile <- Convert(
    source = hfile,
    dest = dest,
    assay = assay,
    overwrite = overwrite,
    verbose = verbose,
    ...
  )
  dfile$close_all()
  return(invisible(x = dfile$filename))
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
  dfile <- switch(
    EXPR = stype,
    'h5ad' = switch(
      EXPR = dtype,
      'h5seurat' = H5ADToH5Seurat(
        source = source,
        dest = dest,
        assay = assay,
        overwrite = overwrite,
        verbose = verbose
      ),
      stop("Unable to convert H5AD files to ", dtype, " files", call. = FALSE)
    ),
    'h5mu' = switch(
      EXPR = dtype,
      'h5seurat' = H5MUToH5Seurat(
        source = source,
        dest = dest,
        assay = assay,
        overwrite = overwrite,
        verbose = verbose
      ),
      'h5ad' = H5MUToH5AD(
        source = source,
        dest = dest,
        modality = assay,
        overwrite = overwrite,
        verbose = verbose
      ),
      stop("Unable to convert H5MU files to ", dtype, " files", call. = FALSE)
    ),
    'h5seurat' = switch(
      EXPR = dtype,
      'h5mu' = H5SeuratToH5MU(
        source = source,
        dest = dest,
        assay = assay,
        overwrite = overwrite,
        verbose = verbose
      ),
      stop("Unable to convert h5Seurat files to ", dtype, " files", call. = FALSE)
    ),
    stop("Unknown file type: ", stype, call. = FALSE)
  )
  return(dfile)
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
  ...
) {
  type <- FileType(file = dest)
  if (tolower(x = dest) == type) {
    dest <- paste(file_path_sans_ext(x = source$filename), type, sep = '.')
  }
  dfile <- switch(
    EXPR = type,
    'h5ad' = H5SeuratToH5AD(
      source = source,
      dest = dest,
      assay = assay,
      overwrite = overwrite,
      verbose = verbose
    ),
    stop("Unable to convert h5Seurat files to ", type, " files", call. = FALSE)
  )
  return(dfile)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Implementations
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
#' \code{\link[SeuratDisk]{SeuratDisk.dtype.dataframe_as_group}} is \code{FALSE}
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
        # col.order <- h5attr(x = dfile[['var']], which = 'column-order')
        # col.order <- c(col.order, var.name)
        # dfile[['var']]$attr_rename(
        #   old_attr_name = 'column-order',
        #   new_attr_name = 'old-column-order'
        # )
        # dfile[['var']]$create_attr(
        #   attr_name = 'column-order',
        #   robj = col.order,
        #   dtype = GuessDType(x = col.order)
        # )
        # dfile[['var']]$attr_delete(attr_name = 'old-column-order')
      }
      dfgroup$link_delete(name = '__categories')
    }
    return(invisible(x = NULL))
  }
  ds.map <- c(
    scale.data = if (inherits(x = source[['X']], what = 'H5D')) {
      'X'
    } else {
      NULL
    },
    data = if (inherits(x = source[['X']], what = 'H5D') && source$exists(name = 'raw')) {
      'raw/X'
    } else {
      'X'
    },
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
    assay.group$obj_copy_from(
      src_loc = source,
      src_name = ds.map[[i]],
      dst_name = dst
    )
    # if (assay.group[[dst]]$attr_exists(attr_name = 'shape')) {
    if (isTRUE(x = AttrExists(x = assay.group[[dst]], name = 'shape'))) {
      dims <- rev(x = h5attr(x = assay.group[[dst]], which = 'shape'))
      assay.group[[dst]]$create_attr(
        attr_name = 'dims',
        robj = dims,
        dtype = GuessDType(x = dims)
      )
      assay.group[[dst]]$attr_delete(attr_name = 'shape')
    }
  }
  features.source <- ifelse(
    test = source$exists(name = 'raw') && source$exists(name = 'raw/var'),
    yes = 'raw/var',
    no = 'var'
  )
  if (inherits(x = source[[features.source]], what = 'H5Group')) {
    features.dset <- GetRownames(dset = features.source)
    assay.group$obj_copy_from(
      src_loc = source,
      src_name = paste(features.source, features.dset, sep = '/'),
      dst_name = 'features'
    )
  } else {
    tryCatch(
      expr = assay.group$create_dataset(
        name = 'features',
        robj = rownames(x = source[[features.source]]),
        dtype = GuessDType(x = "")
      ),
      error = function(...) {
        stop("Cannot find feature names in this H5AD file", call. = FALSE)
      }
    )
  }
  scaled <- !is.null(x = ds.map['scale.data']) && !is.na(x = ds.map['scale.data'])
  if (scaled) {
    if (inherits(x = source[['var']], what = 'H5Group')) {
      scaled.dset <- GetRownames(dset = 'var')
      assay.group$obj_copy_from(
        src_loc = source,
        src_name = paste0('var/', scaled.dset),
        dst_name = 'scaled.features'
      )
    } else {
      tryCatch(
        expr = assay.group$create_dataset(
          name = 'scaled.features',
          robj = rownames(x = source[['var']]),
          dtype = GuessDType(x = "")
        ),
        error = function(...) {
          stop("Cannot find scaled features in this H5AD file", call. = FALSE)
        }
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
    } else {
      warning(
        "Cannot yet add feature-level metadata from compound datasets",
        call. = FALSE,
        immediate. = TRUE
      )
      assay.group$create_group(name = 'meta.features')
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
    } else {
      warning(
        "Cannot yet add feature-level metadata from compound datasets",
        call. = FALSE,
        immediate. = TRUE
      )
      assay.group$create_group(name = 'meta.features')
    }
  }
  ColToFactor(dfgroup = assay.group[['meta.features']])
  # if (assay.group[['meta.features']]$attr_exists(attr_name = 'column-order')) {
  if (isTRUE(x = AttrExists(x = assay.group[['meta.features']], name = 'column-order'))) {
    colnames <- h5attr(
      x = assay.group[['meta.features']],
      which = 'column-order'
    )
    assay.group[['meta.features']]$create_attr(
      attr_name = 'colnames',
      robj = colnames,
      dtype = GuessDType(x = colnames)
    )
  }
  if (inherits(x = source[['var']], what = 'H5Group')) {
    assay.group[['meta.features']]$link_delete(name = GetRownames(dset = 'var'))
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
    NormalizeH5ADCategorical <- function(dfgroup) {
      for (col_name in names(dfgroup)) {
        col_obj <- dfgroup[[col_name]]
        # Check if this is an h5ad categorical group (has 'categories' and 'codes')
        if (inherits(col_obj, 'H5Group') &&
            all(c('categories', 'codes') %in% names(col_obj))) {
          # Rename 'categories' to 'levels'
          if (!col_obj$exists('levels')) {
            col_obj$obj_copy_from(
              src_loc = col_obj,
              src_name = 'categories',
              dst_name = 'levels'
            )
            col_obj$link_delete(name = 'categories')
          }
          # Rename 'codes' to 'values'
          if (!col_obj$exists('values')) {
            col_obj$obj_copy_from(
              src_loc = col_obj,
              src_name = 'codes',
              dst_name = 'values'
            )
            col_obj$link_delete(name = 'codes')
          }
        }
      }
    }
    NormalizeH5ADCategorical(dfgroup = dfile[['meta.data']])
    ColToFactor(dfgroup = dfile[['meta.data']])
    # if (dfile[['meta.data']]$attr_exists(attr_name = 'column-order')) {
    if (isTRUE(x = AttrExists(x = dfile[['meta.data']], name = 'column-order'))) {
      colnames <- h5attr(x = dfile[['meta.data']], which = 'column-order')
      dfile[['meta.data']]$create_attr(
        attr_name = 'colnames',
        robj = colnames,
        dtype = GuessDType(x = colnames)
      )
    }
    rownames <- GetRownames(dset = 'obs')
    dfile$obj_copy_from(
      src_loc = dfile,
      src_name = paste0('meta.data/', rownames),
      dst_name = 'cell.names'
    )
    dfile[['meta.data']]$link_delete(name = rownames)
  } else {
    warning(
      "No cell-level metadata present, creating fake cell names",
      call. = FALSE,
      immediate. = TRUE
    )
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
  # Add dimensional reduction information
  if (source$exists(name = 'obsm')) {
    # Add cell embeddings
    if (inherits(x = source[['obsm']], what = 'H5Group')) {
      for (reduc in names(x = source[['obsm']])) {
        sreduc <- gsub(pattern = '^X_', replacement = '', x = reduc)
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
  # Add nearest-neighbor graph
  if (Exists(x = source, name = 'uns/neighbors/distances')) {
    graph.name <- paste(
      assay,
      ifelse(
        test = source$exists(name = 'uns/neighbors/params/method'),
        yes = source[['uns/neighbors/params/method']][1],
        no = 'anndata'
      ),
      sep = '_'
    )
    if (verbose) {
      message("Saving nearest-neighbor graph as ", graph.name)
    }
    dfile[['graphs']]$obj_copy_from(
      src_loc = source,
      src_name = 'uns/neighbors/distances',
      dst_name = graph.name
    )
    # if (dfile[['graphs']][[graph.name]]$attr_exists(attr_name = 'shape')) {
    if (isTRUE(x = AttrExists(x = dfile[['graphs']], name = 'shape'))) {
      dfile[['graphs']][[graph.name]]$create_attr(
        attr_name = 'dims',
        robj = h5attr(x = dfile[['graphs']][[graph.name]], which = 'shape'),
        dtype = GuessDType(x = h5attr(
          x = dfile[['graphs']][[graph.name]],
          which = 'shape'
        ))
      )
      dfile[['graphs']][[graph.name]]$attr_delete(attr_name = 'shape')
    }
    dfile[['graphs']][[graph.name]]$create_attr(
      attr_name = 'assay.used',
      robj = assay,
      dtype = GuessDType(x = assay)
    )
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
      dfile[['misc']]$obj_copy_from(
        src_loc = source[['uns']],
        src_name = i,
        dst_name = i
      )
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
#'  \link[Seurat:IsGlobal]{global} will be transfered to the H5AD file. For
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
  verbose = TRUE
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
  # Transfer data frames from h5Seurat files to H5AD files
  #
  # @param src Source dataset
  # @param dname Name of destination
  # @param index Integer values of rows to take
  #
  # @return Invisibly returns \code{NULL}
  #
  TransferDF <- function(src, dname, index) {
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
        meta_group <- source[['meta.data']]
        if (inherits(x = meta_group, what = 'H5Group')) {
          for (col in names(x = meta_group)) {
            if (col == '__categories' || col == '_index') next

            if (IsFactor(x = meta_group[[col]])) {
              # Handle factors
              if (!dfile[['obs']]$exists(name = '__categories')) {
                dfile[['obs']]$create_group(name = '__categories')
              }
              dfile[['obs']]$create_dataset(
                name = col,
                robj = meta_group[[col]][['values']][index] - 1L,
                dtype = meta_group[[col]][['values']]$get_type()
              )
              dfile[['obs']][['__categories']]$create_dataset(
                name = col,
                robj = meta_group[[col]][['levels']][],
                dtype = meta_group[[col]][['levels']]$get_type()
              )
            } else {
              # Regular column
              dfile[['obs']]$create_dataset(
                name = col,
                robj = if (is.null(index)) meta_group[[col]][] else meta_group[[col]][index],
                dtype = meta_group[[col]]$get_type()
              )
            }
          }

          # Add column order
          if (meta_group$attr_exists(attr_name = 'colnames')) {
            dfile[['obs']]$create_attr(
              attr_name = 'column-order',
              robj = h5attr(x = meta_group, which = 'colnames'),
              dtype = GuessDType(x = h5attr(x = meta_group, which = 'colnames'))
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
                space = Scalar()
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
      CompoundToGroup(
        src = src,
        dst = dfile,
        dname = dname,
        order = 'column-order',
        index = index
      )
    } else if (inherits(x = src, what = 'H5Group')) {
      dfile$create_group(name = dname)
      for (i in src$names) {
        if (IsFactor(x = src[[i]])) {
          # Use newer anndata format: each categorical is a group with categories/codes
          dfile[[dname]]$create_group(name = i)

          # Add codes dataset (0-based indices)
          dfile[[dname]][[i]]$create_dataset(
            name = 'codes',
            robj = src[[i]][['values']][index] - 1L,
            dtype = src[[H5Path(i, 'values')]]$get_type()
          )

          # Add encoding attributes to codes
          dfile[[dname]][[i]][['codes']]$create_attr(
            attr_name = 'encoding-type',
            robj = 'array',
            dtype = GuessDType(x = 'array'),
            space = Scalar()
          )
          dfile[[dname]][[i]][['codes']]$create_attr(
            attr_name = 'encoding-version',
            robj = '0.2.0',
            dtype = GuessDType(x = '0.2.0'),
            space = Scalar()
          )

          # Add categories dataset
          dfile[[dname]][[i]]$create_dataset(
            name = 'categories',
            robj = src[[i]][['levels']][],
            dtype = src[[H5Path(i, 'levels')]]$get_type()
          )

          # Add encoding attributes to categories
          dfile[[dname]][[i]][['categories']]$create_attr(
            attr_name = 'encoding-type',
            robj = 'array',
            dtype = GuessDType(x = 'array'),
            space = Scalar()
          )
          dfile[[dname]][[i]][['categories']]$create_attr(
            attr_name = 'encoding-version',
            robj = '0.2.0',
            dtype = GuessDType(x = '0.2.0'),
            space = Scalar()
          )

          # Add attributes to the categorical group
          dfile[[dname]][[i]]$create_attr(
            attr_name = 'encoding-type',
            robj = 'categorical',
            dtype = GuessDType(x = 'categorical'),
            space = Scalar()
          )
          dfile[[dname]][[i]]$create_attr(
            attr_name = 'encoding-version',
            robj = '0.2.0',
            dtype = GuessDType(x = '0.2.0'),
            space = Scalar()
          )
          dfile[[dname]][[i]]$create_attr(
            attr_name = 'ordered',
            robj = FALSE,
            dtype = GuessDType(x = FALSE),
            space = Scalar()
          )
        } else {
          # Regular dataset
          dfile[[dname]]$create_dataset(
            name = i,
            robj = src[[i]][index],
            dtype = src[[i]]$get_type()
          )

          # Add encoding attributes
          dfile[[dname]][[i]]$create_attr(
            attr_name = 'encoding-type',
            robj = 'array',
            dtype = GuessDType(x = 'array'),
            space = Scalar()
          )
          dfile[[dname]][[i]]$create_attr(
            attr_name = 'encoding-version',
            robj = '0.2.0',
            dtype = GuessDType(x = '0.2.0'),
            space = Scalar()
          )
        }
      }
      if (src$attr_exists(attr_name = 'colnames')) {
        dfile[[dname]]$create_attr(
          attr_name = 'column-order',
          robj = h5attr(x = src, which = 'colnames'),
          dtype = GuessDType(x = h5attr(x = src, which = 'colnames'))
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
          space = Scalar()
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
        dfile[[dname]]$create_attr(
          attr_name = attr.name,
          robj = attr.value,
          dtype = GuessDType(x = attr.value),
          space = Scalar()
        )
      }

      # Add shape attribute for sparse matrices (required by anndata)
      if (dfile[[dname]]$exists(name = 'indptr') &&
          dfile[[dname]]$exists(name = 'indices')) {
        # For CSR matrix: nrows = length(indptr) - 1
        indptr_len <- prod(dfile[[dname]][['indptr']]$dims)
        nrows <- as.integer(indptr_len - 1L)

        # Get ncols from max index + 1 (most reliable for sparse matrices)
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
    if (assay.group$exists(name = 'scale.data')) {
      x.data <- 'scale.data'
      raw.data <- if (assay.group[['layers']]$exists(name = 'data')) {
        'layers/data'
      } else if (assay.group[['layers']]$exists(name = 'counts')) {
        'layers/counts'
      } else {
        NULL
      }
    } else if (assay.group[['layers']]$exists(name = 'data')) {
      x.data <- 'layers/data'
      raw.data <- if (assay.group[['layers']]$exists(name = 'counts')) {
        'layers/counts'
      } else {
        NULL
      }
    } else if (assay.group[['layers']]$exists(name = 'counts')) {
      # Only counts available - use it for X and no raw
      x.data <- 'layers/counts'
      raw.data <- NULL
    } else {
      stop("Cannot find data or counts in V5 h5Seurat file", call. = FALSE)
    }
  } else {
    # Legacy structure
    if (source$index()[[assay]]$slots[['scale.data']]) {
      x.data <- 'scale.data'
      raw.data <- 'data'
    } else {
      x.data <- 'data'
      raw.data <- if (source$index()[[assay]]$slots[['counts']]) {
        'counts'
      } else {
        NULL
      }
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
  assay.group$obj_copy_to(dst_loc = dfile, dst_name = 'X', src_name = x.data)

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
  # Get number of features based on structure
  if (has_layers && assay.group$exists(name = 'meta.data/_index')) {
    # V5 structure: get count from meta.data/_index
    n_features <- length(SafeH5DRead(assay.group[['meta.data/_index']]))
  } else {
    # Legacy structure: use features dims
    n_features <- prod(assay.group[['features']]$dims)
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
    dtype = h5types$H5T_STRING$new(size = Inf)$set_cset(cset = h5const$H5T_CSET_UTF8)$set_strpad(strpad = h5const$H5T_STR_NULLTERM)
  )

  # Add encoding attributes to _index
  dfile[['var']][[rownames]]$create_attr(
    attr_name = 'encoding-type',
    robj = 'string-array',
    dtype = GuessDType(x = 'string-array'),
    space = Scalar()
  )
  dfile[['var']][[rownames]]$create_attr(
    attr_name = 'encoding-version',
    robj = '0.2.0',
    dtype = GuessDType(x = '0.2.0'),
    space = Scalar()
  )

  # Add _index attribute pointing to itself
  dfile[['var']]$create_attr(
    attr_name = rownames,
    robj = rownames,
    dtype = GuessDType(x = rownames),
    space = Scalar()
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
      space = Scalar()
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
    assay.group$obj_copy_to(
      dst_loc = dfile[['raw']],
      dst_name = 'X',
      src_name = raw.data
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
    # Add meta.features with validation
    if (assay.group$exists(name = 'meta.features')) {
      meta.features.src <- assay.group[['meta.features']]
      if (inherits(x = meta.features.src, what = c('H5D', 'H5Group'))) {
        TransferDF(
          src = meta.features.src,
          dname = 'raw/var',
          index = seq.default(from = 1, to = prod(assay.group[['features']]$dims))
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
    dfile[['raw/var']]$create_dataset(
      name = rownames,
      robj = features_data_raw,
      dtype = GuessDType(x = features_data_raw[1])
    )
    dfile[['raw/var']]$create_attr(
      attr_name = rownames,
      robj = rownames,
      dtype = GuessDType(x = rownames),
      space = Scalar()
    )
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
      if (inherits(x = meta.src, what = 'H5Group')) {
        # Get all column names except internal ones
        meta_cols <- setdiff(names(x = meta.src), c('__categories', '_index'))

        for (col in meta_cols) {
          tryCatch({
            if (IsFactor(x = meta.src[[col]])) {
              # Handle factors
              if (!dfile[['obs']]$exists(name = '__categories')) {
                dfile[['obs']]$create_group(name = '__categories')
              }
              dfile[['obs']]$create_dataset(
                name = col,
                robj = meta.src[[col]][['values']][] - 1L,
                dtype = meta.src[[col]][['values']]$get_type()
              )
              dfile[['obs']][['__categories']]$create_dataset(
                name = col,
                robj = meta.src[[col]][['levels']][],
                dtype = meta.src[[col]][['levels']]$get_type()
              )
            } else {
              # Handle regular columns
              dfile[['obs']]$create_dataset(
                name = col,
                robj = meta.src[[col]][],
                dtype = meta.src[[col]]$get_type()
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
              space = Scalar()
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
      dtype = h5types$H5T_STRING$new(size = Inf)$set_cset(cset = h5const$H5T_CSET_UTF8)$set_strpad(strpad = h5const$H5T_STR_NULLTERM)
    )

    # Add encoding attributes to _index
    dfile[['obs']][[rownames]]$create_attr(
      attr_name = 'encoding-type',
      robj = 'string-array',
      dtype = GuessDType(x = 'string-array'),
      space = Scalar()
    )
    dfile[['obs']][[rownames]]$create_attr(
      attr_name = 'encoding-version',
      robj = '0.2.0',
      dtype = GuessDType(x = '0.2.0'),
      space = Scalar()
    )

    # Add _index attribute pointing to itself
    dfile[['obs']]$create_attr(
      attr_name = rownames,
      robj = rownames,
      dtype = GuessDType(x = rownames),
      space = Scalar()
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
            space = Scalar()
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

          # Create spatial matrix (cells x 2)
          spatial_matrix <- cbind(y_coords, x_coords)
        } else if (coord_group$exists('imagerow') && coord_group$exists('imagecol')) {
          # Visium-style coordinates
          row_coords <- coord_group[['imagerow']][]
          col_coords <- coord_group[['imagecol']][]

          # Create spatial matrix (cells x 2)
          spatial_matrix <- cbind(row_coords, col_coords)
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
          spatial_matrix <- cbind(coords_mat[, 2], coords_mat[, 1])
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
          space = Scalar()
        )
        obsm[['spatial']]$create_attr(
          attr_name = 'encoding-version',
          robj = '0.2.0',
          dtype = GuessDType(x = '0.2.0'),
          space = Scalar()
        )
      }
    }
  }

  # Create uns
  dfile$create_group(name = 'uns')

  # Add spatial metadata to uns if available
  if (source$exists(name = 'images')) {
    images <- names(x = source[['images']])

    if (length(images) > 0) {
      if (verbose) {
        message("Adding spatial metadata to uns['spatial']")
      }

      # Create spatial group in uns
      spatial_uns <- dfile[['uns']]$create_group(name = 'spatial')

      # Process each image (using first as primary library)
      for (i in seq_along(images)) {
        img_name <- images[i]
        img_group <- source[['images']][[img_name]]

        # Use library_1, library_2, etc. as keys
        lib_id <- paste0('library_', i)
        lib_group <- spatial_uns$create_group(name = lib_id)

        # Add scale factors if available
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
              space = Scalar(),
              chunk_dims = NULL
            )
          }
        }

        # Add metadata
        meta_group <- lib_group$create_group(name = 'metadata')
        meta_group$create_dataset(
          name = 'image_name',
          robj = img_name
        )

        # Add images if available
        if (img_group$exists(name = 'image')) {
          images_group <- lib_group$create_group(name = 'images')
          img_data <- img_group[['image']]$read()
          if (length(dim(img_data)) == 3L) {
            img_arr <- aperm(img_data, c(3L, 2L, 1L))
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
          space = Scalar()
        )
      }
      # Add encoding-version
      if (!dfile[[group_name]]$attr_exists(attr_name = 'encoding-version')) {
        dfile[[group_name]]$create_attr(
          attr_name = 'encoding-version',
          robj = '0.1.0',
          dtype = GuessDType(x = '0.1.0'),
          space = Scalar()
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
          space = Scalar()
        )
      }
      if (!dfile[['obsm']][[reduc_name]]$attr_exists(attr_name = 'encoding-version')) {
        dfile[['obsm']][[reduc_name]]$create_attr(
          attr_name = 'encoding-version',
          robj = '0.2.0',
          dtype = GuessDType(x = '0.2.0'),
          space = Scalar()
        )
      }
    }
  }

  dfile$flush()
  return(dfile)
}


#' Convert H5MU files to h5Seurat files
#'
#' @inheritParams Convert
#'
#' @return Returns a handle to \code{dest} as an \code{\link{h5Seurat}} object
#'
#' @keywords internal
#'
H5MUToH5Seurat <- function(
  source,
  dest,
  assay = 'RNA',
  overwrite = FALSE,
  verbose = TRUE
) {
  if (verbose) {
    message("Converting H5MU to h5Seurat via Seurat object...")
  }

  # Load h5mu file as Seurat object
  seurat_obj <- LoadH5MU(
    file = source$filename,
    verbose = verbose
  )

  # Save as h5Seurat
  h5seurat_file <- SaveH5Seurat(
    object = seurat_obj,
    filename = dest,
    overwrite = overwrite,
    verbose = verbose
  )

  # Return h5Seurat connection
  dfile <- h5Seurat$new(filename = dest, mode = 'r')
  return(dfile)
}


#' Convert h5Seurat files to H5MU files
#'
#' @inheritParams Convert
#'
#' @return Returns a handle to \code{dest} as an \code{\link[hdf5r]{H5File}} object
#'
#' @keywords internal
#'
H5SeuratToH5MU <- function(
  source,
  dest,
  assay = DefaultAssay(object = source),
  overwrite = FALSE,
  verbose = TRUE
) {
  if (verbose) {
    message("Converting h5Seurat to H5MU via Seurat object...")
  }

  # Load h5Seurat as Seurat object
  seurat_obj <- LoadH5Seurat(
    file = source$filename,
    verbose = verbose
  )

  # Save as h5mu
  h5mu_file <- SaveH5MU(
    object = seurat_obj,
    filename = dest,
    overwrite = overwrite,
    verbose = verbose
  )

  # Return H5File connection to h5mu
  dfile <- H5File$new(filename = dest, mode = 'r')
  return(dfile)
}


#' Convert H5MU files to H5AD files (extract single modality)
#'
#' @inheritParams Convert
#' @param modality Name of modality to extract from h5mu file
#'
#' @return Returns a handle to \code{dest} as an \code{\link[hdf5r]{H5File}} object
#'
#' @keywords internal
#'
H5MUToH5AD <- function(
  source,
  dest,
  modality = 'rna',
  overwrite = FALSE,
  verbose = TRUE
) {
  if (verbose) {
    message("Extracting modality '", modality, "' from H5MU to H5AD...")
  }

  # Load h5mu file
  seurat_obj <- LoadH5MU(
    file = source$filename,
    modalities = modality,
    verbose = verbose
  )

  # Get the corresponding assay name
  assay_names <- Assays(seurat_obj)
  if (length(assay_names) == 0) {
    stop("No assays found in converted object", call. = FALSE)
  }

  # Use first assay (should be the only one if modality was specified)
  target_assay <- assay_names[1]

  # Save as h5Seurat first
  temp_h5seurat <- tempfile(fileext = ".h5seurat")
  on.exit(file.remove(temp_h5seurat), add = TRUE)

  SaveH5Seurat(
    object = seurat_obj,
    filename = temp_h5seurat,
    overwrite = TRUE,
    verbose = FALSE
  )

  # Convert h5Seurat to h5ad
  temp_h5seurat_conn <- Connect(filename = temp_h5seurat, force = TRUE)
  dfile <- H5SeuratToH5AD(
    source = temp_h5seurat_conn,
    dest = dest,
    assay = target_assay,
    overwrite = overwrite,
    verbose = verbose
  )

  temp_h5seurat_conn$close_all()

  return(dfile)
}


# Read obs meta data from h5ad file and return a data.frame
#' @export
#'
readH5AD_obs <- function(file) {
  suppressWarnings(expr = hfile <- SeuratDisk:: Connect(filename = file, force = TRUE))
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
    hfile_cate <- hfile_obs[['__categories']]
    for (i in seq_along(obs_groups)) {
      obs.i <- obs_groups[i]
      obs_value_i <- hfile_obs[[obs.i]][]
      if (obs.i %in% names(x = hfile_cate)){
        obs_value_i <- factor(x = obs_value_i, labels =  hfile_cate[[obs.i]][])
      }
      matrix[,i] <- obs_value_i
    }
  } else {
    for (i in seq_along(obs_groups)) {
      obs.i <- obs_groups[i]
      if (all(names(hfile_obs[[obs.i]]) == c("categories", "codes"))) {
        if (
          length(unique(hfile_obs[[obs.i]][['codes']][])) == length(hfile_obs[[obs.i]][['categories']][])
          ) {
          obs_value_i <- factor(
            x = hfile_obs[[obs.i]][['codes']][],
            labels =  hfile_obs[[obs.i]][['categories']][]
            )
        } else {
          obs_value_i <- hfile_obs[[obs.i]][['codes']][]
        }
      
      } else {
        obs_value_i <- hfile_obs[[obs.i]][]
      }
      matrix[,i] <- obs_value_i
    }
  }
  hfile$close_all()
  return(matrix)
}

# Read obsm from h5ad file and return a list of embeddings
#' @export
#' 
readH5AD_obsm <-  function(file) {
  hfile <- SeuratDisk:: Connect(filename = file, force = TRUE)
  hfile_obsm <- hfile[['obsm']]
  if (length(names(hfile_obsm)) == 0) {
    message('No obsm if found in this object')
    return (list())
  }
  obsm_set <- names(hfile_obsm)
  cells.name <- hfile[['obs']][['_index']][]
  obsm.list <- lapply(obsm_set, function(x) {
    emb <- t(hfile_obsm[[x]][,])
    rownames(emb) <- cells.name
    key.name <- gsub('X_', '', x)
    colnames(emb) <- paste0(key.name, "_", 1:ncol(emb))
    return(emb)
  })
  names(obsm.list) <- gsub('X_', '',obsm_set)
  hfile$close_all()
  return(obsm.list)
}
