#' Read HDF5 Files
#'
#' Read data from HDF5 files
#'
#' @param x An HDF5 dataset (H5D) or group (H5Group) object
#' @param ... Arguments passed to other methods
#'
#' @return Varies depending on the method being called
#'
#' @name ReadH5
#' @rdname ReadH5
#'
#' @seealso \code{\link{as.sparse}}
#'
NULL

#' @param x An HDF5 dataset (H5D) object
#'
#' @return \code{as.array}: an array with the data from the HDF5 dataset
#'
#' @importFrom hdf5r h5attr
#' @importMethodsFrom Matrix t
#'
#' @aliases as.array
#'
#' @rdname ReadH5
#' @method as.array H5D
#' @export
#'
as.array.H5D <- function(x, ...) {
  # Simplified V5-compatible reading
  if (length(x$dims) >= 2) {
    # For multidimensional arrays, try comma indexing first (most reliable)
    tryCatch(
      expr = as.array(x = x[,], ...),
      error = function(e) {
        # Fallback to direct read
        as.array(x = x$read(), ...)
      }
    )
  } else {
    # For 1D arrays, use direct read
    as.array(x = x$read(), ...)
  }
}

#' @inheritParams base::as.data.frame
#'
#' @return \code{as.data.frame}: returns a \code{\link[base]{data.frame}} with
#' the data from the HDF5 dataset
#'
#' @rdname ReadH5
#' @method as.data.frame H5D
#' @export
#'
as.data.frame.H5D <- function(x, row.names = NULL, optional = FALSE, ...) {
  # Simplified V5-compatible reading
  df <- if (length(x$dims) >= 2) {
    # For multidimensional datasets, try comma indexing first
    tryCatch(
      expr = x[,],
      error = function(e) x$read()
    )
  } else {
    # For 1D datasets, use bracket notation or fallback to read
    tryCatch(
      expr = x[],
      error = function(e) x$read()
    )
  }

  # Handle row names if provided
  if (!is.null(x = row.names)) {
    row.names(x = df) <- row.names
  }

  # Handle column names unless optional
  if (!optional) {
    tryCatch({
      colnames(x = df) <- make.names(names = x$get_type()$get_cpd_labels())
    }, error = function(e) {
      # Column names may not be available for all datasets
    })
  }

  # Handle logical columns properly
  if (x$attr_exists(attr_name = 'logicals')) {
    bool.cols <- intersect(
      x = h5attr(x = x, which = 'logicals'),
      y = colnames(x = df)
    )
    for (i in bool.cols) {
      df[[i]] <- as.logical(x = df[[i]])
    }
  }

  return(df)
}

#' @rdname ReadH5
#' @method as.data.frame H5Group
#' @export
#'
as.data.frame.H5Group <- function(x, row.names = NULL, optional = FALSE, ...) {
  df <- NULL
  idx <- NULL
  colnames <- NULL
  if ('index' %in% names(x = x)) {
    idx <- x[['index']][]
  }
  if (x$attr_exists(attr_name = 'colnames')) {
    colnames <- h5attr(x = x, which = 'colnames')
  }

  for (i in names(x = x)) {
    if (i == 'index') {
      next
    }

    # Check the type of the object first
    obj_class <- class(x[[i]])

    # Handle H5Groups (which might be factors)
    if (inherits(x[[i]], "H5Group")) {
      # Check if this is a factor (has levels and values)
      if (IsFactor(x = x[[i]])) {
        # Read as factor
        tryCatch({
          values <- as.integer(x = x[[i]][['values']][])
          levels <- x[[i]][['levels']][]

          # Convert 0-based indices (h5ad format) to 1-based indices (R format)
          values <- values + 1L

          if (!is.null(x = idx)) {
            values <- values[order(idx)]
          }

          if (is.null(x = df)) {
            df <- data.frame(row.names = row.names %||% seq_along(along.with = values))
          }

          # Check if ordered
          if (x[[i]]$attr_exists(attr_name = 'ordered') &&
              h5attr(x = x[[i]], which = 'ordered')) {
            df[[i]] <- ordered(x = levels[values], levels = levels)
          } else {
            df[[i]] <- factor(x = levels[values], levels = levels)
          }
        }, error = function(e) {
          # Skip if can't read as factor
          NULL
        })
      }
      next
    }

    # Skip environment objects
    if (inherits(x[[i]], "environment")) {
      next
    }

    # Safe reading with error handling
    tryCatch({
      dset <- as.vector(x = x[[i]][])

      if (!is.null(x = idx)) {
        dset <- dset[order(idx)]
      }

      if (is.null(x = df)) {
        df <- data.frame(row.names = row.names %||% seq_along(along.with = dset))
      }

      # Handle factors
      if (IsFactor(x = x[[i]])) {
        if (x[[i]]$attr_exists(attr_name = 'ordered') &&
            h5attr(x = x[[i]], which = 'ordered')) {
          df[[i]] <- ordered(
            x = x[[i]][['levels']][dset],
            levels = x[[i]][['levels']][]
          )
        } else {
          df[[i]] <- factor(
            x = x[[i]][['levels']][dset],
            levels = x[[i]][['levels']][]
          )
        }
      } else {
        df[[i]] <- dset
      }
    }, error = function(e) {
      # Skip columns that can't be read
      NULL  # Silently skip problematic columns
    })
  }

  # Ensure we have a valid data frame
  if (is.null(df)) {
    df <- data.frame(row.names = row.names %||% character(0))
  }

  # Validate columns before subsetting
  if (!is.null(x = colnames) && ncol(df) > 0) {
    valid_cols <- intersect(colnames, names(df))
    if (length(valid_cols) > 0) {
      df <- df[, valid_cols, drop = FALSE]
    }
  }

  return(df)
}

#' @return \code{as.list}: a list with the data from the HDF5 dataset
#'
#' @rdname ReadH5
#' @method as.list H5D
#' @export
#'
as.list.H5D <- function(x, ...) {
  return(lapply(X = x[], FUN = identity))
}

#' @rdname ReadH5
#' @method as.list H5Group
#' @export
#'
as.list.H5Group <- function(x, recursive = TRUE, ...) {
  lst <- list()
  for (i in names(x = x)) {
    # Check if this is a subgroup
    if (inherits(x[[i]], "H5Group")) {
      if (recursive) {
        # Recursively convert subgroups to lists
        lst[[i]] <- as.list(x[[i]], recursive = recursive, ...)
      } else {
        # Keep as H5Group object if not recursive
        lst[[i]] <- x[[i]]
      }
    } else {
      # Handle datasets
      stype <- GetStringType(stype = x[[i]])
      # Ensure stype is a single value
      if (length(stype) != 1 || is.null(stype)) {
        # Default handling if stype is not valid
        lst[[i]] <- x[[i]][]
      } else {
        lst[[i]] <- switch(
          EXPR = stype,
          'H5T_STRING' = as.character(x = x[[i]][]),
          'numeric' = as.vector(x = x[[i]][]),
          'factor' = as.factor.H5D(x = x[[i]]),
          x[[i]][]  # Default case
        )
      }
    }
  }
  return(lst)
}

#' @param transpose Transpose the matrix before returning
#'
#' @return \code{as.matrix}: a matrix with the data from the HDF5 dataset
#'
#' @aliases as.matrix
#'
#' @rdname ReadH5
#' @method as.matrix H5D
#' @export
#'
as.matrix.H5D <- function(x, transpose = FALSE, ...) {
  # Simplified V5-compatible reading
  obj <- if (length(x$dims) >= 2) {
    # For multidimensional arrays, try comma indexing first
    tryCatch(
      expr = x[,],
      error = function(e) x$read()
    )
  } else {
    # For 1D arrays, use direct read
    x$read()
  }

  # Apply transpose if needed
  if (transpose) {
    obj <- t(x = obj)
  }

  # Convert to matrix
  return(as.matrix(x = obj))
}

#' @rdname ReadH5
#' @method as.matrix H5Group
#' @export
#'
as.matrix.H5Group <- function(x, ...) {
  return(as.sparse(x = x, ...))
}

#' Convert an HDF5 dataset to a sparse matrix
#'
#' @param x An HDF5 dataset or group
#'
#' @return A sparse matrix
#'
#' @importFrom Matrix sparseMatrix
#'
#' @rdname as.sparse
#' @method as.sparse H5D
#' @export
#'
as.sparse.H5D <- function(x) {
  if (!x$attr_exists(attr_name = 'dims')) {
    stop("Cannot create a sparse matrix without dimensions", call. = FALSE)
  }

  return(sparseMatrix(
    i = x$read()[, 1],
    j = x$read()[, 2],
    x = x$read()[, 3],
    dims = rev(x = h5attr(x = x, which = 'dims'))
  ))
}

#' @rdname as.sparse
#' @method as.sparse H5Group
#' @export
#'
as.sparse.H5Group <- function(x) {
  if (!all(c('indices', 'indptr', 'data') %in% names(x = x))) {
    stop("Not a sparse matrix", call. = FALSE)
  }
  if (!x$attr_exists(attr_name = 'dims')) {
    stop("Cannot create a sparse matrix without dimensions", call. = FALSE)
  }

  return(sparseMatrix(
    i = as.integer(x = x[['indices']][] + 1),
    p = as.integer(x = x[['indptr']][]),
    x = as.vector(x = x[['data']][]),
    dims = rev(x = h5attr(x = x, which = 'dims'))
  ))
}

#' Convert an HDF5 dataset to a factor
#'
#' @param x An HDF5 dataset
#'
#' @return A factor
#'
#' @keywords internal
#'
as.factor.H5D <- function(x) {
  values <- as.integer(x = x[['values']][])
  if (x$attr_exists(attr_name = 'ordered') &&
      h5attr(x = x, which = 'ordered')) {
    return(ordered(
      x = x[['levels']][values],
      levels = x[['levels']][]
    ))
  } else {
    return(factor(
      x = x[['levels']][values],
      levels = x[['levels']][]
    ))
  }
}

#' Check if an HDF5 dataset is a factor
#'
#' @param x An HDF5 dataset
#'
#' @return TRUE if the dataset represents a factor
#'
#' @keywords internal
#'
IsFactor <- function(x) {
  return(is(object = x, class2 = 'H5Group') &&
         all(c('levels', 'values') %in% names(x = x)))
}

#' Get string type from HDF5 dataset
#'
#' @param stype An HDF5 dataset or type
#'
#' @return String representation of the type
#'
#' @keywords internal
#'
GetStringType <- function(stype) {
  if (inherits(x = stype, what = c('H5T', 'H5T_COMPOUND'))) {
    stype <- stype$get_class()
  } else if (inherits(x = stype, what = c('H5D', 'H5Group'))) {
    if (IsFactor(x = stype)) {
      return('factor')
    }
    stype <- stype$get_type()$get_class()
  }

  return(switch(
    EXPR = stype,
    'H5T_INTEGER' = 'numeric',
    'H5T_FLOAT' = 'numeric',
    as.character(x = stype)
  ))
}

#' Check for compound type in HDF5 dataset
#'
#' @param x An HDF5 dataset
#'
#' @return TRUE if the dataset has a compound type
#'
#' @keywords internal
#'
IsCompound <- function(x) {
  dtype <- tryCatch(
    expr = x$get_type(),
    error = function(...) return(NULL)
  )

  if (is.null(x = dtype)) {
    return(FALSE)
  }

  return(inherits(x = dtype, what = 'H5T_COMPOUND'))
}

#' Convert compound dataset to group
#'
#' @param src Source compound dataset
#' @param dst Destination file
#' @param dname Destination name
#' @param order Order attribute name
#' @param index Index for subsetting
#'
#' @keywords internal
#'
CompoundToGroup <- function(src, dst, dname, order, index = NULL) {
  dst$create_group(name = dname)
  dgroup <- dst[[dname]]

  for (i in src$get_type()$get_cpd_labels()) {
    values <- src$read(fields = i)
    if (!is.null(x = index)) {
      values <- values[index]
    }
    dgroup$create_dataset(name = i, robj = values)
  }

  if (src$attr_exists(attr_name = order)) {
    dgroup$create_attr(
      attr_name = order,
      robj = h5attr(x = src, which = order)
    )
  }

  return(invisible(x = NULL))
}