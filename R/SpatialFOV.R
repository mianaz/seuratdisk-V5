#' @include zzz.R
#'
NULL

#' Functions for handling FOV (Field of View) spatial objects
#'
#' FOV objects are used for imaging-based spatial transcriptomics
#' (Xenium, MERFISH, CosMx, etc.)
#'

#' Extract centroids from FOV object
#'
#' @param fov_obj FOV object from Seurat
#'
#' @return Matrix of centroids (cells x 2) with x,y coordinates
#' @keywords internal
ExtractFOVCentroids <- function(fov_obj) {

  centroids <- NULL

  tryCatch({
    if (inherits(fov_obj, "FOV")) {
      # FOV objects have a centroids slot
      if (!is.null(fov_obj@centroids)) {
        cents <- fov_obj@centroids

        # Centroids may be stored as a Centroids object
        if (inherits(cents, "Centroids")) {
          # Extract coordinates
          coords <- Seurat::GetTissueCoordinates(cents)
          if (all(c("x", "y") %in% colnames(coords))) {
            centroids <- as.matrix(coords[, c("x", "y")])
          }
        } else if (is.matrix(cents) || is.data.frame(cents)) {
          # Already in matrix/df format
          if (all(c("x", "y") %in% colnames(cents))) {
            centroids <- as.matrix(cents[, c("x", "y")])
          }
        }
      }
    }

    # Validate
    if (!is.null(centroids)) {
      if (!is.matrix(centroids)) {
        centroids <- as.matrix(centroids)
      }
      if (ncol(centroids) != 2) {
        warning("Centroids should have exactly 2 columns (x, y)", call. = FALSE)
        return(NULL)
      }
    }

  }, error = function(e) {
    warning("Failed to extract FOV centroids: ", e$message, call. = FALSE)
    return(NULL)
  })

  return(centroids)
}

#' Extract segmentation boundaries from FOV object
#'
#' @param fov_obj FOV object from Seurat
#'
#' @return List with segmentation data or NULL
#' @keywords internal
ExtractFOVSegmentation <- function(fov_obj) {

  segmentation <- NULL

  tryCatch({
    if (inherits(fov_obj, "FOV")) {
      # Check for boundaries/segmentation slot
      if (!is.null(fov_obj@boundaries)) {
        bounds <- fov_obj@boundaries

        # Boundaries may be a Segmentation object
        if (inherits(bounds, "Segmentation")) {
          # Extract polygon data
          # This is complex - for now, just note that it exists
          segmentation <- list(
            type = "Segmentation",
            available = TRUE,
            note = "Full polygon extraction not yet implemented"
          )
        }
      }
    }
  }, error = function(e) {
    warning("Failed to extract FOV segmentation: ", e$message, call. = FALSE)
    return(NULL)
  })

  return(segmentation)
}

#' Extract molecule positions from FOV object
#'
#' @param fov_obj FOV object from Seurat
#'
#' @return Data frame with molecule positions or NULL
#' @keywords internal
ExtractFOVMolecules <- function(fov_obj) {

  molecules <- NULL

  tryCatch({
    if (inherits(fov_obj, "FOV")) {
      # Check for molecules slot
      if (!is.null(fov_obj@molecules)) {
        mols <- fov_obj@molecules

        # Molecules may contain x, y, and gene information
        if (is.data.frame(mols)) {
          molecules <- mols
        }
      }
    }
  }, error = function(e) {
    warning("Failed to extract FOV molecules: ", e$message, call. = FALSE)
    return(NULL)
  })

  return(molecules)
}

#' Convert FOV object to h5ad-compatible structure
#'
#' @param fov_obj FOV object from Seurat
#' @param verbose Show progress messages
#'
#' @return List with spatial data in h5ad-compatible format
#' @keywords internal
ConvertFOVToH5AD <- function(fov_obj, verbose = TRUE) {

  if (!inherits(fov_obj, "FOV")) {
    stop("Object is not an FOV object", call. = FALSE)
  }

  result <- list()

  # Extract centroids
  if (verbose) message("  Extracting FOV centroids...")
  centroids <- ExtractFOVCentroids(fov_obj)
  if (!is.null(centroids)) {
    result$centroids <- centroids
    if (verbose) message("    Found ", nrow(centroids), " centroids")
  }

  # Extract segmentation
  if (verbose) message("  Checking for segmentation data...")
  segmentation <- ExtractFOVSegmentation(fov_obj)
  if (!is.null(segmentation)) {
    result$segmentation <- segmentation
    if (verbose) message("    Segmentation data available")
  }

  # Extract molecules
  if (verbose) message("  Checking for molecule data...")
  molecules <- ExtractFOVMolecules(fov_obj)
  if (!is.null(molecules)) {
    result$molecules <- molecules
    if (verbose) message("    Found ", nrow(molecules), " molecules")
  }

  # Add metadata
  result$technology <- "imaging"  # vs "array-based" for Visium
  result$fov_name <- fov_obj@key %||% "fov"

  return(result)
}

#' Convert h5ad spatial data to FOV object
#'
#' @param spatial_data List with spatial coordinates and optionally segmentation
#' @param key FOV key/name
#'
#' @return FOV object or NULL
#' @keywords internal
#' @note This is a placeholder - full FOV reconstruction requires SeuratObject functions
H5ADToFOV <- function(spatial_data, key = "fov") {

  # This would require CreateFOV(), CreateCentroids(), CreateSegmentation()
  # from SeuratObject, which may not be available in all versions

  # For now, return NULL and store data in misc
  warning(
    "Full FOV object reconstruction not yet implemented. ",
    "Spatial data will be stored in @misc slot instead.",
    call. = FALSE
  )

  return(NULL)
}
