#' @title Safe H5D Dataset Access
#' @description Functions for safely accessing HDF5 datasets regardless of dimensionality
#' @keywords internal
#'
NULL

#' Safely read an H5D dataset
#'
#' @param dataset An H5D dataset object
#' @param subset Optional subset indices
#' @return The data from the dataset
#' @keywords internal
#'
SafeH5DRead <- function(dataset, subset = NULL) {
  if (!inherits(x = dataset, what = 'H5D')) {
    stop("Input must be an H5D dataset", call. = FALSE)
  }

  ndims <- length(dataset$dims)

  # If subset is provided, use it
  if (!is.null(subset)) {
    if (ndims == 1) {
      return(dataset[subset])
    } else if (ndims == 2) {
      if (is.list(subset)) {
        return(dataset[subset[[1]], subset[[2]]])
      } else {
        # Assume subset is for first dimension, take all of second
        return(dataset[subset, ])
      }
    } else {
      stop("Cannot subset datasets with >2 dimensions", call. = FALSE)
    }
  }

  # No subset - read all data
  if (ndims == 0) {
    # Scalar
    return(dataset$read())
  } else if (ndims == 1) {
    # 1D array
    return(dataset[])
  } else if (ndims == 2) {
    # 2D array - use appropriate access method
    # Try different methods to find one that works
    result <- tryCatch(
      expr = dataset[,],
      error = function(e) {
        tryCatch(
          expr = dataset$read(),
          error = function(e2) {
            # Last resort: explicit indexing
            dataset[1:dataset$dims[1], 1:dataset$dims[2]]
          }
        )
      }
    )

    # Handle special case where 2nd dimension is 1 (should be 1D vector)
    if (!is.null(dim(result)) && ncol(result) == 1) {
      result <- drop(result)
    }

    return(result)
  } else {
    # Higher dimensional - use read()
    return(dataset$read())
  }
}

#' Check if H5D dataset is effectively 1D
#'
#' @param dataset An H5D dataset object
#' @return Logical indicating if dataset is effectively 1D
#' @keywords internal
#'
IsEffectively1D <- function(dataset) {
  if (!inherits(x = dataset, what = 'H5D')) {
    return(FALSE)
  }

  ndims <- length(dataset$dims)

  if (ndims == 1) {
    return(TRUE)
  } else if (ndims == 2) {
    # Check if one dimension is 1
    return(dataset$dims[1] == 1 || dataset$dims[2] == 1)
  } else {
    return(FALSE)
  }
}

#' Safely convert H5Group to list, handling 3D+ arrays
#'
#' This function converts an HDF5 group to an R list, properly handling
#' datasets with 3 or more dimensions (which fail with standard as.list).
#' This is needed for Squidpy UMAP data that stores 3D arrays.
#'
#' @param h5obj An H5Group or H5D object
#' @param recursive Whether to recursively convert subgroups
#' @return A list representation of the HDF5 structure
#' @keywords internal
#'
SafeH5GroupToList <- function(h5obj, recursive = TRUE) {
  # If it's a dataset, read it directly
  if (inherits(h5obj, "H5D")) {
    return(SafeH5DRead(h5obj))
  }

  # If it's not a group, fall back to standard conversion
  if (!inherits(h5obj, "H5Group")) {
    return(tryCatch(
      as.list(h5obj, recursive = recursive),
      error = function(e) {
        warning("Failed to convert HDF5 object: ", conditionMessage(e))
        return(NULL)
      }
    ))
  }

  # It's a group - convert manually to handle 3D+ datasets
  result <- list()
  obj_names <- names(h5obj)

  if (length(obj_names) == 0) {
    return(result)
  }

  for (name in obj_names) {
    item <- tryCatch(h5obj[[name]], error = function(e) NULL)

    # Skip if we couldn't access the item
    if (is.null(item)) {
      next
    }

    if (inherits(item, "H5D")) {
      # It's a dataset - check dimensionality and dataspace
      # Wrap in tryCatch to handle datasets with special attributes
      result[[name]] <- tryCatch({
        ndims <- length(item$dims)

        if (ndims <= 2) {
          # 1D or 2D - use standard SafeH5DRead
          SafeH5DRead(item)
        } else {
          # 3D or higher - read directly with dataset$read()
          # This preserves the array structure
          item$read()
        }
      }, error = function(e) {
        # Skip datasets that can't be read (e.g., NULL dataspace, special attributes)
        # Silently skip rather than warning for common edge cases
        msg <- conditionMessage(e)
        if (!grepl("Dataspace has to be simple", msg)) {
          warning("Skipping dataset '", name, "': ", msg)
        }
        NULL
      })

      # Remove NULL entries
      if (is.null(result[[name]])) {
        result[[name]] <- NULL
      }
    } else if (inherits(item, "H5Group") && recursive) {
      # It's a subgroup - recurse
      result[[name]] <- SafeH5GroupToList(item, recursive = TRUE)
    } else {
      # Other type - try standard conversion
      result[[name]] <- tryCatch(
        as.list(item, recursive = recursive),
        error = function(e) {
          warning("Failed to convert HDF5 item '", name, "': ", conditionMessage(e))
          NULL
        }
      )

      # Remove NULL entries
      if (is.null(result[[name]])) {
        result[[name]] <- NULL
      }
    }
  }

  return(result)
}

