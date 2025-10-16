#' @include zzz.R
#' @include Convert.R
#' @importFrom hdf5r H5File h5attr h5types
#' @importFrom Seurat CreateSeuratObject CreateAssayObject Images GetTissueCoordinates scalefactors
#' @importFrom SeuratObject AddMetaData Cells CreateFOV
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

  cell_names <- NULL

  if (!is.null(seurat_obj)) {
    cell_names <- Cells(seurat_obj)
  } else if (h5ad$exists("obs")) {
    obs_group <- h5ad[["obs"]]
    if (obs_group$exists("_index")) {
      cell_names <- as.character(obs_group[["_index"]][])
    } else if (obs_group$exists("index")) {
      cell_names <- as.character(obs_group[["index"]][])
    }
  }

  if (h5ad$exists("obsm") && "spatial" %in% names(h5ad[["obsm"]])) {
    spatial_coords <- h5ad[["obsm"]][["spatial"]][,]
    spatial_coords <- as.matrix(spatial_coords)
    if (!is.null(dim(spatial_coords))) {
      if (!is.null(cell_names)) {
        if (nrow(spatial_coords) == length(cell_names) && ncol(spatial_coords) == 2L) {
          rownames(spatial_coords) <- cell_names
        } else if (ncol(spatial_coords) == length(cell_names) && nrow(spatial_coords) == 2L) {
          spatial_coords <- t(spatial_coords)
          rownames(spatial_coords) <- cell_names
        }
      } else if (nrow(spatial_coords) == 2L) {
        spatial_coords <- t(spatial_coords)
      }
    }

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
    cell_names <- Cells(seurat_obj)
  }

  # Add spatial coordinates to Seurat object
  if (has_spatial_coords) {
    if (!is.null(cell_names) && nrow(spatial_coords) == length(cell_names)) {
      spatial_coords <- spatial_coords[cell_names, , drop = FALSE]
    }
    if (ncol(spatial_coords) == 2L) {
      colnames(spatial_coords) <- c('imagerow', 'imagecol')
    }

    # Determine technology type based on data structure
    technology <- DetectSpatialTechnology(spatial_coords, visium_data)

    if (technology == "Visium") {
      # Create Visium-specific spatial object
      seurat_obj <- AddVisiumSpatialData(
        seurat_obj,
        spatial_coords,
        visium_data,
        h5ad,
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
    robj = class(img_obj)[1]
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
AddVisiumSpatialData <- function(seurat_obj, coords, visium_data, h5ad,
                                 assay_name = "Spatial", verbose = TRUE) {

  coords <- as.matrix(coords)
  cell_names <- Cells(seurat_obj)

  if (nrow(coords) != length(cell_names) || ncol(coords) != 2L) {
    warning("Coordinate count doesn't match cell count")
    return(seurat_obj)
  }

  coords <- coords[cell_names, , drop = FALSE]
  colnames(coords) <- c('imagerow', 'imagecol')

  seurat_obj@meta.data$spatial_x <- coords[, 'imagecol']
  seurat_obj@meta.data$spatial_y <- coords[, 'imagerow']

  seurat_obj@misc$spatial_technology <- "Visium"
  seurat_obj@misc$spatial_metadata <- visium_data

  if (is.null(visium_data) || length(visium_data) == 0L || !h5ad$exists("uns") ||
      !h5ad[["uns"]]$exists("spatial")) {
    if (verbose) {
      message("Stored Visium coordinates; no image data available in h5ad")
    }
    return(seurat_obj)
  }

  read_image_dataset <- function(dataset) {
    arr <- dataset$read()
    dims <- dataset$dims
    if (length(dims) == 3L && dims[1] == 3L) {
      arr <- aperm(arr, c(3L, 2L, 1L))
    } else if (length(dims) == 3L && dims[3] == 3L) {
      arr <- aperm(arr, c(1L, 2L, 3L))
    }
    arr[arr < 0] <- 0
    arr[arr > 1] <- 1
    storage.mode(arr) <- 'double'
    arr
  }

  sanitize_key <- function(x) {
    key <- gsub(pattern = '\\.+', replacement = '_', x = make.names(x))
    key <- gsub(pattern = '_+', replacement = '_', x = key)
    key <- gsub(pattern = '^_', replacement = '', x = key)
    key <- gsub(pattern = '_$', replacement = '', x = key)
    if (nchar(key) == 0L) {
      key <- 'spatial'
    }
    paste0(key, '_')
  }

  spatial_uns <- tryCatch(h5ad[["uns/spatial"]], error = function(e) NULL)
  if (is.null(spatial_uns)) {
    return(seurat_obj)
  }

  for (lib_id in names(visium_data)) {
    lib_info <- visium_data[[lib_id]]
    lib_group <- tryCatch(spatial_uns[[lib_id]], error = function(e) NULL)
    if (is.null(lib_group) || !lib_group$exists("images")) {
      next
    }

    images_group <- lib_group[["images"]]
    lowres_img <- if (images_group$exists("lowres")) {
      read_image_dataset(images_group[["lowres"]])
    } else {
      NULL
    }
    if (is.null(lowres_img)) {
      next
    }

    hires_img <- if (images_group$exists("hires")) {
      read_image_dataset(images_group[["hires"]])
    } else {
      NULL
    }

    sf_values <- lib_info$scalefactors %||% list()
    spot <- sf_values[['spot_diameter_fullres']] %||% sf_values[['spot']]
    fiducial <- sf_values[['fiducial_diameter_fullres']] %||% sf_values[['fiducial']]
    hires_sf <- sf_values[['tissue_hires_scalef']] %||% sf_values[['hires']]
    lowres_sf <- sf_values[['tissue_lowres_scalef']] %||% sf_values[['lowres']]

    scales <- scalefactors(
      spot = as.numeric(spot %||% NA_real_),
      fiducial = as.numeric(fiducial %||% NA_real_),
      hires = as.numeric(hires_sf %||% NA_real_),
      lowres = as.numeric(lowres_sf %||% NA_real_)
    )

    radius <- as.numeric(scales[['spot']])
    if (!is.finite(radius)) {
      radius <- 1
    }

    coord_df <- data.frame(
      imagerow = coords[, 'imagerow'],
      imagecol = coords[, 'imagecol'],
      row.names = cell_names,
      stringsAsFactors = FALSE
    )

    image_key <- sanitize_key(lib_id)
    fov <- CreateFOV(
      coords = coord_df[, c('imagerow', 'imagecol'), drop = FALSE],
      type = 'centroids',
      radius = radius,
      assay = assay_name,
      key = image_key
    )

    visium_image <- new(
      Class = 'VisiumV2',
      image = lowres_img,
      scale.factors = scales,
      coordinates = coord_df,
      molecules = fov@molecules,
      boundaries = fov@boundaries,
      assay = assay_name,
      key = image_key
    )
    visium_image <- visium_image[cell_names]
    if (!is.null(hires_img)) {
      attr(visium_image@image, 'hires.image') <- hires_img
    }
    seurat_obj[[lib_id]] <- visium_image

    if (verbose) {
      message("  Added Visium image for library '", lib_id, "'")
    }
  }

  if (verbose) {
    message("Added Visium spatial data:")
    message("  - Spots: ", nrow(coords))
    message("  - Libraries: ", paste(names(visium_data), collapse = ", "))
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
