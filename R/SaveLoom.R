#' @importFrom Seurat Project
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Save a \code{\link[Seurat]{Seurat}} object to a loom file
#'
#' @inheritParams SaveH5Seurat
#'
#' @return \code{SaveLoom}: Invisibly returns \code{filename}
#'
#' @name SaveLoom
#' @rdname SaveLoom
#'
#' @export
#'
SaveLoom <- function(object, filename, overwrite = FALSE, verbose = TRUE, ...) {
  UseMethod(generic = 'SaveLoom', object = object)
}

#' @return \code{as.loom}: A \code{\link{loom}} object
#'
#' @name SaveLoom
#' @rdname SaveLoom
#'
#' @export
#'
as.loom <- function(x, ...) {
  UseMethod(generic = 'as.loom', object = x)
}

#' @keywords internal
#' @noRd
setGeneric(
  name = 'WriteMatrix',
  def = function(x, name, lfile, transpose = TRUE, verbose = TRUE) {
    if (!inherits(x = lfile, what = 'loom')) {
      stop("'lfile' must be a loom object", call. = FALSE)
    }
    standardGeneric(f = 'WriteMatrix')
  },
  signature = c('x')
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname SaveLoom
#' @method SaveLoom default
#' @export
#'
SaveLoom.default <- function(
  object,
  filename,
  overwrite = FALSE,
  verbose = TRUE,
  ...
) {
  tryCatch(
    expr = object <- as.Seurat(object = object, verbose = verbose, ...),
    error = function(...) {
      stop(
        "Unable to coerce an object of class ",
        paste(class(x = object), collapse = ', '),
        " to a Seurat object",
        call. = FALSE
      )
    }
  )
  if (missing(x = filename)) {
    filename <- paste0(Project(object = object), '.loom')
  }
  return(invisible(x = SaveLoom(
    object = object,
    filename = filename,
    overwrite = overwrite,
    verbose = verbose,
    ...
  )))
}

#' @rdname SaveLoom
#' @method SaveLoom Seurat
#' @export
#'
SaveLoom.Seurat <- function(
  object,
  filename = paste0(Project(object = object), '.loom'),
  overwrite = FALSE,
  verbose = TRUE,
  ...
) {
  loom <- as.loom(
    x = object,
    filename = filename,
    overwrite = overwrite,
    verbose = verbose,
    ...
  )
  loom$close_all()
  return(invisible(x = loom$filename))
}

#' @rdname SaveLoom
#' @method as.loom default
#' @export
#'
as.loom.default <- function(x, filename, overwrite = FALSE, verbose = TRUE) {
  tryCatch(
    expr = x <- as.Seurat(object = x, verbose = verbose, ...),
    error = function(...) {
      stop(
        "Unable to coerce an object of class ",
        paste(class(x = x), collapse = ', '),
        " to a Seurat object",
        call. = FALSE
      )
    }
  )
  if (missing(x = filename)) {
    filename <- paste0(Project(object = x), '.loom')
  }
  return(as.loom(
    object = x,
    filename = filename,
    overwrite = overwrite,
    verbose = verbose,
    ...
  ))
}

#' @rdname SaveLoom
#' @method as.loom H5File
#' @export
#'
as.loom.H5File <- function(x, ...) {
  # Validate that this is a loom file
  if (!inherits(x = x, what = 'H5File')) {
    stop("'x' must be an H5File object", call. = FALSE)
  }
  # Check for required loom structure
  if (!x$exists(name = 'matrix')) {
    stop("H5File does not contain required '/matrix' dataset", call. = FALSE)
  }
  required.groups <- c('layers', 'row_attrs', 'col_attrs')
  for (group in required.groups) {
    if (!x$exists(name = group)) {
      stop("H5File does not contain required '", group, "' group", call. = FALSE)
    }
  }
  # Wrap the H5File in a loom object
  loom.obj <- loom$new(
    filename = x$filename,
    mode = ifelse(test = x$get_file_id()$is_writeable(), yes = 'r+', no = 'r')
  )
  return(loom.obj)
}

#' @importFrom tools file_ext
#' @importFrom Seurat Assays DefaultAssay GetAssayData Reductions Embeddings Loadings Stdev
#'
#' @rdname SaveLoom
#' @method as.loom Seurat
#' @export
#'
as.loom.Seurat <- function(
  x,
  filename = paste0(Project(object = x), '.loom'),
  overwrite = FALSE,
  verbose = TRUE,
  ...
) {
  if (!grepl(pattern = '^loom$', x = file_ext(x = filename))) {
    filename <- paste0(filename, '.loom')
  }
  if (file.exists(filename)) {
    if (isTRUE(x = overwrite)) {
      warning(
        "Overwriting previous file ",
        filename,
        call. = FALSE,
        immediate. = TRUE
      )
      file.remove(filename)
    } else {
      stop("Loom file at ", filename, " already exists", call. = FALSE)
    }
  }
  lfile <- loom$new(filename = filename, mode = 'w')
  AddSlots <- function(assay) {
    for (slot in c('counts', 'scale.data')) {
      mat <- GetAssayData(object = x, slot = slot, assay = assay)
      if (identical(x = dim(x = mat), y = dim(x = x))) {
        if (isTRUE(x = verbose)) {
          message("Adding slot ", slot, " for assay ", assay)
        }
        name <- ifelse(
          test = assay == DefaultAssay(object = x),
          yes = slot,
          no = paste(assay, slot, sep = '_')
        )
        lfile$add_layer(
          x = mat,
          name = ifelse(
            test = assay == DefaultAssay(object = x),
            yes = slot,
            no = paste(assay, slot, sep = '_')
          ),
          verbose = verbose
        )
      }
    }
    return(invisible(x = NULL))
  }
  # Add assays
  assays.write <- lapply(
    X = Assays(object = x),
    FUN = function(i) {
      return(if (identical(x = dim(x = x[[i]]), y = dim(x = x))) {
        i
      } else {
        NULL
      })
    }
  )
  assays.write <- unlist(x = Filter(f = Negate(f = is.null), x = assays.write))
  assays.write <- setdiff(x = DefaultAssay(object = x), assays.write)
  if (verbose) {
    message("Saving data from ", DefaultAssay(object = x), " as /matrix")
  }
  WriteMatrix(
    x = GetAssayData(object = x, slot = 'data'),
    name = 'matrix',
    lfile = lfile,
    verbose = verbose
  )
  AddSlots(assay = DefaultAssay(object = x))
  for (assay in assays.write) {
    if (isTRUE(x = verbose)) {
      message("Adding data for ", assay)
    }
    lfile$add_layer(
      x = GetAssayData(object = x, slot = data, assay = assay),
      name = assay,
      verbose = verbose
    )
    AddSlots(assay = assay)
  }
  # Add dimensional reduction information
  if (length(Reductions(object = x)) > 0) {
    # Create parent reductions group first
    if (!lfile$exists('reductions')) {
      lfile$create_group(name = 'reductions')
    }
  }
  for (reduc in Reductions(object = x)) {
    tryCatch({
      if (verbose) {
        message("Adding dimensional reduction ", reduc)
      }
      reduc.obj <- x[[reduc]]
      # Create a group for this reduction
      reduc.group <- lfile$create_group(name = H5Path('reductions', reduc))
      # Save embeddings (transpose to match loom convention: components x cells)
      embeddings <- Embeddings(object = reduc.obj)
      lfile$create_dataset(
        name = H5Path('reductions', reduc, 'embeddings'),
        robj = t(x = embeddings),
        dtype = GuessDType(x = embeddings[1, 1])
      )
      # Save loadings if available
      loadings <- tryCatch(Loadings(object = reduc.obj), error = function(e) NULL)
      if (!is.null(x = loadings) && !IsMatrixEmpty(x = loadings)) {
        lfile$create_dataset(
          name = H5Path('reductions', reduc, 'loadings'),
          robj = loadings,
          dtype = GuessDType(x = loadings[1, 1])
        )
      }
      # Save standard deviation if available
      stdev <- tryCatch(Stdev(object = reduc.obj), error = function(e) NULL)
      if (!is.null(x = stdev) && length(x = stdev) > 0) {
        lfile$create_dataset(
          name = H5Path('reductions', reduc, 'stdev'),
          robj = stdev,
          dtype = GuessDType(x = stdev)
        )
      }
    }, error = function(e) {
      warning(
        "Failed to save dimensional reduction ", reduc, ": ", e$message,
        call. = FALSE,
        immediate. = TRUE
      )
    })
  }
  # Add graphs
  for (graph in SafeGraphs(object = x)) {
    lfile$add_graph(x = x[[graph]], name = graph, type = 'col')
  }
  # Add metadata
  lfile$add_attribute(x = colnames(x = x), name = 'CellID', type = 'col')
  for (md in names(x = x[[]])) {
    lfile$add_attribute(x = x[[md, drop = TRUE]], name = md, type = 'col')
  }
  # Add feature-level metadata
  lfile$add_attribute(x = rownames(x = x), name = 'Gene', type = 'row')
  # Add feature metadata from the default assay
  default.assay <- DefaultAssay(object = x)
  if (!is.null(x = default.assay) && default.assay %in% Assays(object = x)) {
    feature.meta <- x[[default.assay]][[]]
    if (ncol(x = feature.meta) > 0) {
      for (feat.md in colnames(x = feature.meta)) {
        tryCatch({
          lfile$add_attribute(
            x = feature.meta[[feat.md]],
            name = feat.md,
            type = 'row'
          )
        }, error = function(e) {
          if (verbose) {
            warning(
              "Failed to save feature metadata '", feat.md, "': ", e$message,
              call. = FALSE,
              immediate. = TRUE
            )
          }
        })
      }
    }
  }
  return(lfile)
}

#' @importFrom Matrix t
#' @importFrom hdf5r H5S
#' @importFrom methods slot
#' @importFrom utils setTxtProgressBar
#' @keywords internal
#' @noRd
setMethod(
  f = 'WriteMatrix',
  signature = c('x' = 'dgCMatrix'),
  definition = function(x, name, lfile, transpose = TRUE, verbose = TRUE) {
    if (isTRUE(x = transpose)) {
      x <- Matrix::t(x = x)
    }
    lfile$create_dataset(
      name = name,
      dtype = GuessDType(x = x[1, 1]),
      space = H5S$new(dims = dim(x = x))
    )
    MARGIN <- GetMargin(dims = lfile[[name]]$dims)
    chunk.points <- ChunkPoints(
      dsize = lfile[[name]]$dims[MARGIN],
      csize = lfile[[name]]$chunk_dims[MARGIN]
    )
    dims <- vector(mode = 'list', length = 2L)
    dims[[-MARGIN]] <- seq_len(length.out = lfile[[name]]$dims[-MARGIN])
    if (isTRUE(x = verbose)) {
      pb <- PB()
    }
    for (i in seq_len(length.out = nrow(x = chunk.points))) {
      dims[[MARGIN]] <- seq.default(
        from = chunk.points[i, 'start'],
        to = chunk.points[i, 'end']
      )
      lfile[[name]]$write(
        args = dims,
        value = as.matrix(x = x[dims[[1]], dims[[2]]])
      )
      if (isTRUE(x = verbose)) {
        setTxtProgressBar(pb = pb, value = i / nrow(x = chunk.points))
      }
    }
    if (isTRUE(x = verbose)) {
      close(con = pb)
    }
    return(invisible(x = NULL))
  }
)

#' @keywords internal
#' @noRd
setMethod(
  f = 'WriteMatrix',
  signature = c('x' = 'matrix'),
  definition = function(x, name, lfile, transpose = TRUE, verbose = TRUE) {
    if (isTRUE(x = transpose)) {
      x <- t(x = x)
    }
    lfile$create_dataset(
      name = name,
      robj = x,
      dtype = GuessDType(x = x[1, 1])
    )
    return(invisible(x = NULL))
  }
)
