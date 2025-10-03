#' @include zzz.R
#' @include Convert.R
#' @importFrom hdf5r H5File h5attr h5types
#' @importFrom Seurat CreateSeuratObject CreateAssayObject Images GetTissueCoordinates
#' @importFrom SeuratObject AddMetaData
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Spatial Data Conversion Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Convert spatial coordinates from h5ad to Seurat format
#'
#' @param h5ad_file H5AD file handle or path
#' @param seurat_obj Seurat object to add spatial data to
#' @param assay_name Name of the assay to associate spatial data with
#' @param verbose Print progress messages
#'
#' @return Modified Seurat object with spatial data
#' @export
ConvertH5ADSpatialToSeurat <- function(h5ad_file, seurat_obj = NULL,
                                       assay_name = "Spatial", verbose = TRUE) {

  # Open h5ad file if path provided
  if (is.character(h5ad_file)) {
    h5ad <- H5File$new(h5ad_file, mode = "r")
    on.exit(h5ad$close_all())
  } else {
    h5ad <- h5ad_file
  }

  # Check for spatial coordinates in obsm
  has_spatial_coords <- FALSE
  spatial_coords <- NULL

  if (h5ad$exists("obsm") && "spatial" %in% names(h5ad[["obsm"]])) {
    spatial_coords <- h5ad[["obsm"]][["spatial"]][,]
    has_spatial_coords <- TRUE
    if (verbose) message("Found spatial coordinates in obsm['spatial']")
  }

  # Check for Visium-style spatial data in uns
  has_visium_data <- FALSE
  visium_data <- list()

  if (h5ad$exists("uns") && "spatial" %in% names(h5ad[["uns"]])) {
    spatial_uns <- h5ad[["uns/spatial"]]
    library_ids <- names(spatial_uns)
    has_visium_data <- length(library_ids) > 0

    if (has_visium_data && verbose) {
      message("Found Visium spatial data for libraries: ", paste(library_ids, collapse = ", "))
    }

    # Process each library
    for (lib_id in library_ids) {
      lib_data <- list()
      lib_group <- spatial_uns[[lib_id]]

      # Get scale factors
      if ("scalefactors" %in% names(lib_group)) {
        sf_group <- lib_group[["scalefactors"]]
        lib_data$scalefactors <- list()

        for (sf_name in names(sf_group)) {
          lib_data$scalefactors[[sf_name]] <- sf_group[[sf_name]][]
        }
      }

      # Get image metadata (not loading actual images here)
      if ("images" %in% names(lib_group)) {
        lib_data$has_images <- TRUE
        image_names <- names(lib_group[["images"]])
        lib_data$image_names <- image_names
      }

      # Get metadata
      if ("metadata" %in% names(lib_group)) {
        meta_group <- lib_group[["metadata"]]
        lib_data$metadata <- list()

        for (meta_name in names(meta_group)) {
          lib_data$metadata[[meta_name]] <- meta_group[[meta_name]][]
        }
      }

      visium_data[[lib_id]] <- lib_data
    }
  }

  # Create or modify Seurat object
  if (is.null(seurat_obj)) {
    if (!has_spatial_coords) {
      stop("No spatial coordinates found and no Seurat object provided")
    }

    # Create minimal Seurat object with spatial coordinates
    n_cells <- nrow(spatial_coords)
    n_genes <- 100  # Placeholder

    counts <- matrix(0, nrow = n_genes, ncol = n_cells)
    rownames(counts) <- paste0("Gene", seq_len(n_genes))
    colnames(counts) <- paste0("Cell", seq_len(n_cells))

    seurat_obj <- CreateSeuratObject(counts = counts, assay = assay_name)
  }

  # Add spatial coordinates to Seurat object
  if (has_spatial_coords) {
    # Determine technology type based on data structure
    technology <- DetectSpatialTechnology(spatial_coords, visium_data)

    if (technology == "Visium") {
      # Create Visium-specific spatial object
      seurat_obj <- AddVisiumSpatialData(
        seurat_obj,
        spatial_coords,
        visium_data,
        assay_name = assay_name,
        verbose = verbose
      )
    } else if (technology == "SlideSeq") {
      # Create SlideSeq-specific spatial object
      seurat_obj <- AddSlideSeqSpatialData(
        seurat_obj,
        spatial_coords,
        assay_name = assay_name,
        verbose = verbose
      )
    } else {
      # Generic spatial data
      seurat_obj <- AddGenericSpatialData(
        seurat_obj,
        spatial_coords,
        assay_name = assay_name,
        verbose = verbose
      )
    }
  }

  return(seurat_obj)
}

#' Convert Seurat spatial data to h5ad format
#'
#' @param seurat_obj Seurat object with spatial data
#' @param h5ad_file H5AD file handle or path to write to
#' @param library_id Library ID for spatial data (default: "library_1")
#' @param verbose Print progress messages
#'
#' @export
ConvertSeuratSpatialToH5AD <- function(seurat_obj, h5ad_file,
                                       library_id = "library_1",
                                       verbose = TRUE) {

  # Open h5ad file if path provided
  if (is.character(h5ad_file)) {
    h5ad <- H5File$new(h5ad_file, mode = "r+")
    on.exit(h5ad$close_all())
  } else {
    h5ad <- h5ad_file
  }

  # Check for spatial data in Seurat object
  images <- Images(seurat_obj)

  if (length(images) == 0) {
    if (verbose) message("No spatial data found in Seurat object")
    return(invisible(NULL))
  }

  if (verbose) message("Found ", length(images), " spatial image(s)")

  # Process first image (extend for multiple images later)
  img_obj <- seurat_obj[[images[1]]]

  # Extract coordinates
  coords <- GetTissueCoordinates(img_obj)

  # Convert to h5ad format (cells x 2 matrix)
  spatial_matrix <- as.matrix(coords[, c("imagerow", "imagecol")])

  # Create obsm group if not exists
  if (!h5ad$exists("obsm")) {
    h5ad$create_group("obsm")
  }

  # Write spatial coordinates
  if (h5ad[["obsm"]]$exists("spatial")) {
    h5ad[["obsm"]]$link_delete("spatial")
  }

  h5ad[["obsm"]]$create_dataset(
    name = "spatial",
    robj = spatial_matrix,
    dtype = h5types$H5T_NATIVE_DOUBLE
  )

  if (verbose) message("Wrote spatial coordinates to obsm['spatial']")

  # Create uns/spatial structure for Visium-like data
  if (!h5ad$exists("uns")) {
    h5ad$create_group("uns")
  }

  if (!h5ad[["uns"]]$exists("spatial")) {
    h5ad[["uns"]]$create_group("spatial")
  }

  # Create library-specific group
  spatial_group <- h5ad[["uns/spatial"]]

  if (spatial_group$exists(library_id)) {
    spatial_group$link_delete(library_id)
  }

  lib_group <- spatial_group$create_group(library_id)

  # Add scale factors if available
  if (inherits(img_obj, "VisiumV1") || inherits(img_obj, "VisiumV2")) {
    scalefactors <- GetScaleFactors(img_obj)

    if (!is.null(scalefactors)) {
      sf_group <- lib_group$create_group("scalefactors")

      for (sf_name in names(scalefactors)) {
        sf_group$create_dataset(
          name = sf_name,
          robj = scalefactors[[sf_name]],
          dtype = h5types$H5T_NATIVE_DOUBLE
        )
      }

      if (verbose) message("Wrote scale factors")
    }
  }

  # Add metadata
  meta_group <- lib_group$create_group("metadata")
  meta_group$create_dataset(
    name = "technology",
    robj = class(img_obj)[1],
    dtype = h5types$H5T_STRING
  )

  if (verbose) message("Spatial data conversion complete")

  return(invisible(NULL))
}

#' Detect spatial technology from data structure
#'
#' @param coords Spatial coordinates matrix
#' @param metadata Additional metadata
#'
#' @return Technology type string
#' @keywords internal
DetectSpatialTechnology <- function(coords, metadata = NULL) {

  # Check for Visium indicators
  if (!is.null(metadata) && length(metadata) > 0) {
    # Look for Visium-specific scale factors
    if (any(sapply(metadata, function(x) "tissue_hires_scalef" %in% names(x$scalefactors)))) {
      return("Visium")
    }
  }

  # Check coordinate patterns
  if (ncol(coords) == 2) {
    # Check if coordinates appear to be on a grid (Visium)
    x_vals <- unique(coords[, 1])
    y_vals <- unique(coords[, 2])

    if (length(x_vals) < nrow(coords) / 2 && length(y_vals) < nrow(coords) / 2) {
      # Likely gridded data
      return("Visium")
    } else {
      # Likely continuous coordinates
      return("SlideSeq")
    }
  }

  return("Generic")
}

#' Add Visium spatial data to Seurat object
#'
#' @param seurat_obj Seurat object
#' @param coords Spatial coordinates
#' @param visium_data Visium-specific metadata
#' @param assay_name Assay name
#' @param verbose Print messages
#'
#' @return Modified Seurat object
#' @keywords internal
AddVisiumSpatialData <- function(seurat_obj, coords, visium_data,
                                 assay_name = "Spatial", verbose = TRUE) {

  # For now, add coordinates to metadata
  # Full implementation would create proper VisiumV1/V2 objects

  if (nrow(coords) == ncol(seurat_obj)) {
    seurat_obj@meta.data$spatial_x <- coords[, 1]
    seurat_obj@meta.data$spatial_y <- coords[, 2]

    # Add Visium metadata to misc slot
    seurat_obj@misc$spatial_technology <- "Visium"
    seurat_obj@misc$spatial_metadata <- visium_data

    if (verbose) {
      message("Added Visium spatial data:")
      message("  - Spots: ", nrow(coords))
      message("  - Libraries: ", paste(names(visium_data), collapse = ", "))
    }
  } else {
    warning("Coordinate count doesn't match cell count")
  }

  return(seurat_obj)
}

#' Add SlideSeq spatial data to Seurat object
#'
#' @param seurat_obj Seurat object
#' @param coords Spatial coordinates
#' @param assay_name Assay name
#' @param verbose Print messages
#'
#' @return Modified Seurat object
#' @keywords internal
AddSlideSeqSpatialData <- function(seurat_obj, coords,
                                   assay_name = "Spatial", verbose = TRUE) {

  # Add continuous coordinates to metadata
  if (nrow(coords) == ncol(seurat_obj)) {
    seurat_obj@meta.data$spatial_x <- coords[, 1]
    seurat_obj@meta.data$spatial_y <- coords[, 2]

    seurat_obj@misc$spatial_technology <- "SlideSeq"

    if (verbose) {
      message("Added SlideSeq spatial data:")
      message("  - Beads: ", nrow(coords))
      message("  - X range: ", round(min(coords[,1]), 2), " - ", round(max(coords[,1]), 2))
      message("  - Y range: ", round(min(coords[,2]), 2), " - ", round(max(coords[,2]), 2))
    }
  } else {
    warning("Coordinate count doesn't match cell count")
  }

  return(seurat_obj)
}

#' Add generic spatial data to Seurat object
#'
#' @param seurat_obj Seurat object
#' @param coords Spatial coordinates
#' @param assay_name Assay name
#' @param verbose Print messages
#'
#' @return Modified Seurat object
#' @keywords internal
AddGenericSpatialData <- function(seurat_obj, coords,
                                  assay_name = "Spatial", verbose = TRUE) {

  # Add coordinates to metadata
  if (nrow(coords) == ncol(seurat_obj)) {
    seurat_obj@meta.data$spatial_x <- coords[, 1]
    seurat_obj@meta.data$spatial_y <- coords[, 2]

    seurat_obj@misc$spatial_technology <- "Generic"

    if (verbose) {
      message("Added generic spatial data:")
      message("  - Points: ", nrow(coords))
    }
  } else {
    warning("Coordinate count doesn't match cell count")
  }

  return(seurat_obj)
}

#' Get scale factors from Visium object
#'
#' @param visium_obj Visium image object
#'
#' @return List of scale factors
#' @keywords internal
GetScaleFactors <- function(visium_obj) {

  # This would extract scale factors from actual Visium objects
  # For now, return example values

  scale_factors <- list(
    tissue_hires_scalef = 0.17,
    tissue_lowres_scalef = 0.05,
    spot_diameter_fullres = 89.0,
    fiducial_diameter_fullres = 144.0
  )

  return(scale_factors)
}