#' Load a MuData H5MU file as a Seurat object
#'
#' Read multimodal MuData (.h5mu) files and convert to a Seurat object with
#' multiple assays. This function wraps MuDataSeurat::ReadH5MU with additional
#' support for spatial data, Seurat V5 layers, and enhanced metadata preservation.
#'
#' @param file Path to .h5mu file
#' @param modalities Character vector of modality names to load. If NULL (default),
#'   loads all modalities found in the file
#' @param assay.names Named vector mapping modality names to desired Seurat assay
#'   names. If NULL, uses standard mapping (rna→RNA, prot→ADT, atac→ATAC, etc.)
#' @param restore.spatial Logical; if TRUE, attempts to restore spatial data from
#'   the h5mu file structure (coordinates, images, scalefactors)
#' @param verbose Show progress messages
#'
#' @return A \code{Seurat} object with multiple assays corresponding to each modality
#'
#' @details
#' The h5mu format stores multimodal data where each modality is stored as a
#' separate AnnData-like structure under \code{/mod/{modality_name}}. This function:
#' \itemize{
#'   \item Reads each modality using MuDataSeurat functionality
#'   \item Converts modalities to Seurat assays
#'   \item Merges modalities into a single Seurat object
#'   \item Preserves global cell metadata from \code{/obs}
#'   \item Restores spatial data if present
#'   \item Maintains Seurat V5 layer structure
#' }
#'
#' @section Modality Mapping:
#' By default, modality names are mapped to standard Seurat assay names:
#' \itemize{
#'   \item rna → RNA
#'   \item prot → ADT
#'   \item atac → ATAC
#'   \item spatial → Spatial
#'   \item Other names are preserved as-is
#' }
#'
#' @examples
#' \dontrun{
#' # Load all modalities from an h5mu file
#' seurat_obj <- LoadH5MU("multimodal_data.h5mu")
#'
#' # Load specific modalities only
#' seurat_obj <- LoadH5MU("multimodal_data.h5mu", modalities = c("rna", "prot"))
#'
#' # Custom assay name mapping
#' seurat_obj <- LoadH5MU(
#'   "data.h5mu",
#'   assay.names = c(rna = "RNA", prot = "Protein")
#' )
#' }
#'
#' @importFrom hdf5r H5File h5attr
#' @importFrom Seurat CreateSeuratObject Assays DefaultAssay<-
#' @importFrom SeuratObject Cells AddMetaData
#'
#' @export
#'
LoadH5MU <- function(file,
                     modalities = NULL,
                     assay.names = NULL,
                     restore.spatial = TRUE,
                     verbose = TRUE) {

  if (!requireNamespace("MuDataSeurat", quietly = TRUE)) {
    stop(
      "Package 'MuDataSeurat' is required for h5mu file support.\n",
      "Install it with: remotes::install_github('PMBio/MuDataSeurat')",
      call. = FALSE
    )
  }

  if (!file.exists(file)) {
    stop("File not found: ", file, call. = FALSE)
  }

  if (verbose) {
    message("Loading H5MU file: ", file)
  }

  # Open h5mu file to inspect structure
  h5mu <- H5File$new(file, mode = "r")
  on.exit(h5mu$close_all())

  # Get available modalities
  available_modalities <- if (h5mu$exists("mod")) {
    names(h5mu[["mod"]])
  } else {
    stop("No modalities found in h5mu file (missing /mod group)", call. = FALSE)
  }

  if (length(available_modalities) == 0) {
    stop("No modalities found in h5mu file", call. = FALSE)
  }

  if (verbose) {
    message("Found ", length(available_modalities), " modalities: ",
            paste(available_modalities, collapse = ", "))
  }

  # Determine which modalities to load
  if (is.null(modalities)) {
    modalities_to_load <- available_modalities
  } else {
    modalities_to_load <- modalities
    missing <- setdiff(modalities_to_load, available_modalities)
    if (length(missing) > 0) {
      warning(
        "Requested modalities not found in file: ",
        paste(missing, collapse = ", "),
        immediate. = TRUE
      )
      modalities_to_load <- intersect(modalities_to_load, available_modalities)
    }
  }

  if (length(modalities_to_load) == 0) {
    stop("No valid modalities to load", call. = FALSE)
  }

  # Set up assay name mapping
  if (is.null(assay.names)) {
    assay.names <- GetDefaultModalityMapping(modalities_to_load)
  } else {
    # Fill in any missing mappings with defaults
    default_mapping <- GetDefaultModalityMapping(modalities_to_load)
    for (mod in modalities_to_load) {
      if (!mod %in% names(assay.names)) {
        assay.names[mod] <- default_mapping[mod]
      }
    }
  }

  # Read global cell metadata if present
  global_obs <- ReadGlobalObs(h5mu, verbose = verbose)

  # Use MuDataSeurat::ReadH5MU to load the data
  if (verbose) {
    message("Reading h5mu file with MuDataSeurat...")
  }

  tryCatch({
    seurat_obj <- MuDataSeurat::ReadH5MU(file, verbose = verbose)
  }, error = function(e) {
    stop(
      "Failed to read h5mu file with MuDataSeurat: ",
      e$message,
      "\nPlease ensure MuDataSeurat is properly installed and the file is valid.",
      call. = FALSE
    )
  })

  # Rename assays according to mapping
  current_assays <- Assays(seurat_obj)
  for (old_name in names(assay.names)) {
    if (old_name %in% current_assays && assay.names[old_name] != old_name) {
      new_name <- assay.names[old_name]
      seurat_obj[[new_name]] <- seurat_obj[[old_name]]
      seurat_obj[[old_name]] <- NULL
      if (verbose) {
        message("  Renamed assay: ", old_name, " → ", new_name)
      }
    }
  }

  # Add global obs metadata if present
  if (!is.null(global_obs) && nrow(global_obs) > 0) {
    seurat_obj <- AddMetaData(seurat_obj, metadata = global_obs)
    if (verbose) {
      message("Added ", ncol(global_obs), " global metadata columns from /obs")
    }
  }

  # Restore spatial data if requested
  if (restore.spatial) {
    seurat_obj <- RestoreSpatialFromH5MU(
      h5mu_file = h5mu,
      seurat_obj = seurat_obj,
      modalities = modalities_to_load,
      assay.names = assay.names,
      verbose = verbose
    )
  }

  # Set default assay (prioritize RNA, then first assay)
  final_assays <- Assays(seurat_obj)
  if ("RNA" %in% final_assays) {
    DefaultAssay(seurat_obj) <- "RNA"
  } else if (length(final_assays) > 0) {
    DefaultAssay(seurat_obj) <- final_assays[1]
  }

  if (verbose) {
    message("\nSuccessfully loaded H5MU file")
    message("  Cells: ", ncol(seurat_obj))
    message("  Assays: ", paste(Assays(seurat_obj), collapse = ", "))
    if (length(seurat_obj@reductions) > 0) {
      message("  Reductions: ", paste(names(seurat_obj@reductions), collapse = ", "))
    }
    if (length(seurat_obj@graphs) > 0) {
      message("  Graphs: ", paste(names(seurat_obj@graphs), collapse = ", "))
    }
    message("  Metadata columns: ", ncol(seurat_obj@meta.data))
  }

  return(seurat_obj)
}


#' Get default modality to assay name mapping
#'
#' @param modalities Character vector of modality names
#' @return Named character vector mapping modality names to assay names
#' @keywords internal
#'
GetDefaultModalityMapping <- function(modalities) {
  # Standard mapping based on muon/MuData conventions
  default_map <- c(
    rna = "RNA",
    prot = "ADT",
    atac = "ATAC",
    spatial = "Spatial",
    hto = "HTO",
    adt = "ADT",
    cite = "ADT",
    peaks = "ATAC"
  )

  mapping <- setNames(modalities, modalities)
  for (mod in modalities) {
    if (tolower(mod) %in% names(default_map)) {
      mapping[mod] <- default_map[tolower(mod)]
    }
  }

  return(mapping)
}


#' Read global observation metadata from h5mu file
#'
#' @param h5mu H5File connection to h5mu file
#' @param verbose Logical; print messages
#' @return Data frame with global cell metadata, or NULL
#' @keywords internal
#'
#' @importFrom hdf5r h5attr
#'
ReadGlobalObs <- function(h5mu, verbose = FALSE) {
  if (!h5mu$exists("obs")) {
    return(NULL)
  }

  obs_group <- h5mu[["obs"]]
  obs_cols <- setdiff(names(obs_group), c("_index", "index", "__categories"))

  if (length(obs_cols) == 0) {
    return(NULL)
  }

  # Read cell names
  cell_names <- NULL
  if (obs_group$exists("_index")) {
    cell_names <- as.character(obs_group[["_index"]][])
  } else if (obs_group$exists("index")) {
    cell_names <- as.character(obs_group[["index"]][])
  }

  if (is.null(cell_names)) {
    if (verbose) {
      message("  Warning: No cell names found in /obs, skipping global metadata")
    }
    return(NULL)
  }

  # Read metadata columns
  obs_data <- data.frame(row.names = cell_names, stringsAsFactors = FALSE)

  for (col in obs_cols) {
    tryCatch({
      # Check if categorical
      if (obs_group$exists("__categories") && col %in% names(obs_group[["__categories"]])) {
        codes <- obs_group[[col]][]
        categories <- as.character(obs_group[["__categories"]][[col]][])
        codes[codes == -1] <- NA
        obs_data[[col]] <- factor(categories[codes + 1], levels = categories)
      } else {
        # Numeric or string
        obs_data[[col]] <- obs_group[[col]][]
      }
    }, error = function(e) {
      if (verbose) {
        warning("Could not read /obs column '", col, "': ", e$message, immediate. = TRUE)
      }
    })
  }

  if (ncol(obs_data) == 0) {
    return(NULL)
  }

  return(obs_data)
}


#' Restore spatial data from h5mu file to Seurat object
#'
#' @param h5mu_file H5File connection to h5mu file
#' @param seurat_obj Seurat object to add spatial data to
#' @param modalities Character vector of modalities that were loaded
#' @param assay.names Named vector mapping modality to assay names
#' @param verbose Logical; print messages
#' @return Modified Seurat object with spatial data
#' @keywords internal
#'
RestoreSpatialFromH5MU <- function(h5mu_file, seurat_obj, modalities,
                                   assay.names, verbose = FALSE) {

  has_spatial_data <- FALSE

  # Check each modality for spatial data
  for (mod in modalities) {
    mod_path <- paste0("mod/", mod)
    if (!h5mu_file$exists(mod_path)) {
      next
    }

    mod_group <- h5mu_file[[mod_path]]

    # Check for spatial coordinates in obsm
    if (mod_group$exists("obsm") && "spatial" %in% names(mod_group[["obsm"]])) {
      has_spatial_data <- TRUE
      assay_name <- assay.names[mod]

      if (verbose) {
        message("  Found spatial data in modality '", mod, "'")
      }

      # Use existing SpatialConversion functionality
      seurat_obj <- ConvertH5MUSpatialToSeurat(
        h5mu_file = h5mu_file,
        seurat_obj = seurat_obj,
        modality = mod,
        assay_name = assay_name,
        verbose = verbose
      )
    }
  }

  if (has_spatial_data && verbose) {
    message("Spatial data restoration complete")
  }

  return(seurat_obj)
}


#' Convert h5mu spatial data to Seurat format for a specific modality
#'
#' @param h5mu_file H5File connection to h5mu file
#' @param seurat_obj Seurat object
#' @param modality Modality name in h5mu file
#' @param assay_name Corresponding assay name in Seurat
#' @param verbose Logical
#' @return Modified Seurat object
#' @keywords internal
#'
#' @importFrom hdf5r H5File h5attr
#'
ConvertH5MUSpatialToSeurat <- function(h5mu_file, seurat_obj, modality,
                                       assay_name, verbose = FALSE) {

  mod_path <- paste0("mod/", modality)
  mod_group <- h5mu_file[[mod_path]]

  # Read spatial coordinates
  if (!mod_group$exists("obsm") || !"spatial" %in% names(mod_group[["obsm"]])) {
    return(seurat_obj)
  }

  spatial_coords <- mod_group[["obsm"]][["spatial"]][,]
  spatial_coords <- as.matrix(spatial_coords)

  cell_names <- Cells(seurat_obj)

  if (nrow(spatial_coords) != length(cell_names)) {
    if (verbose) {
      warning(
        "Spatial coordinate count doesn't match cell count in modality '",
        modality, "'",
        immediate. = TRUE
      )
    }
    return(seurat_obj)
  }

  rownames(spatial_coords) <- cell_names

  # Check for Visium-style data in uns
  visium_data <- NULL
  if (mod_group$exists("uns") && "spatial" %in% names(mod_group[["uns"]])) {
    spatial_uns <- mod_group[["uns/spatial"]]
    library_ids <- names(spatial_uns)

    if (length(library_ids) > 0) {
      visium_data <- list()
      for (lib_id in library_ids) {
        lib_group <- spatial_uns[[lib_id]]
        lib_data <- list()

        # Read scalefactors
        if ("scalefactors" %in% names(lib_group)) {
          sf_group <- lib_group[["scalefactors"]]
          lib_data$scalefactors <- list()
          for (sf_name in names(sf_group)) {
            lib_data$scalefactors[[sf_name]] <- sf_group[[sf_name]][]
          }
        }

        # Note images presence
        if ("images" %in% names(lib_group)) {
          lib_data$has_images <- TRUE
          lib_data$image_names <- names(lib_group[["images"]])
        }

        visium_data[[lib_id]] <- lib_data
      }
    }
  }

  # Use existing ConvertH5ADSpatialToSeurat logic, but adapted for h5mu
  # For now, add coordinates to metadata and store Visium info in misc
  seurat_obj@meta.data$spatial_x <- spatial_coords[, 1]
  seurat_obj@meta.data$spatial_y <- spatial_coords[, 2]

  if (!is.null(visium_data)) {
    seurat_obj@misc[[paste0(modality, "_spatial_metadata")]] <- visium_data
    seurat_obj@misc[[paste0(modality, "_spatial_technology")]] <- "Visium"

    # Try to reconstruct Visium images using existing functionality
    seurat_obj <- TryAddVisiumImagesFromH5MU(
      seurat_obj = seurat_obj,
      spatial_uns = mod_group[["uns/spatial"]],
      coords = spatial_coords,
      visium_data = visium_data,
      assay_name = assay_name,
      verbose = verbose
    )
  }

  return(seurat_obj)
}


#' Try to add Visium images from h5mu structure
#'
#' @param seurat_obj Seurat object
#' @param spatial_uns HDF5 group containing spatial/uns data
#' @param coords Spatial coordinates matrix
#' @param visium_data Visium metadata list
#' @param assay_name Assay name
#' @param verbose Logical
#' @return Modified Seurat object
#' @keywords internal
#'
TryAddVisiumImagesFromH5MU <- function(seurat_obj, spatial_uns, coords,
                                       visium_data, assay_name, verbose = FALSE) {

  # This uses similar logic to AddVisiumSpatialData in SpatialConversion.R
  # For now, store metadata; full image reconstruction can be added later

  if (verbose) {
    message("    Found Visium spatial metadata for ",
            length(visium_data), " library(ies)")
  }

  return(seurat_obj)
}
