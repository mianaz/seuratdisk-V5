#' Load Spatial Demo Datasets for Testing
#'
#' This script demonstrates loading various spatial datasets from SeuratData
#' for testing spatial data conversion functionality
#'
#' @examples
#' \dontrun{
#' # Load required libraries
#' library(Seurat)
#' library(SeuratData)
#' library(SeuratDisk)
#'
#' # Install and load spatial datasets
#' load_spatial_demo_data()
#' }
#'
#' @export

load_spatial_demo_data <- function(datasets = c("stxBrain", "stxKidney", "ssHippo"),
                                   install_missing = TRUE,
                                   verbose = TRUE) {

  # Check if SeuratData is available
  if (!requireNamespace("SeuratData", quietly = TRUE)) {
    stop("SeuratData is required. Install with: remotes::install_github('satijalab/seurat-data')")
  }

  library(SeuratData)

  # Get list of available datasets
  available_data <- AvailableData()
  spatial_datasets <- available_data[grep("spatial|visium|slide",
                                          available_data$Summary,
                                          ignore.case = TRUE), ]

  if (verbose) {
    message("Available spatial datasets in SeuratData:")
    print(spatial_datasets[, c("Dataset", "Species", "Technology")])
  }

  # Dictionary of known spatial datasets and their properties
  spatial_data_info <- list(
    stxBrain = list(
      name = "stxBrain",
      desc = "Mouse Brain Serial Section 1 (Sagittal-Anterior)",
      technology = "Visium",
      species = "Mouse",
      install_name = "stxBrain.SeuratData"
    ),
    stxKidney = list(
      name = "stxKidney",
      desc = "Mouse Kidney Section (Coronal)",
      technology = "Visium",
      species = "Mouse",
      install_name = "stxKidney.SeuratData"
    ),
    ssHippo = list(
      name = "ssHippo",
      desc = "Slide-seq v2 mouse hippocampus",
      technology = "Slide-seqV2",
      species = "Mouse",
      install_name = "ssHippo.SeuratData"
    )
  )

  loaded_data <- list()

  for (dataset in datasets) {
    if (dataset %in% names(spatial_data_info)) {
      info <- spatial_data_info[[dataset]]

      if (verbose) {
        message("\n", paste(rep("=", 50), collapse = ""))
        message("Processing: ", info$name)
        message("Description: ", info$desc)
        message("Technology: ", info$technology)
      }

      # Check if dataset is installed
      if (!info$name %in% InstalledData()$Dataset) {
        if (install_missing) {
          if (verbose) message("Installing ", info$name, "...")
          InstallData(info$name)
        } else {
          warning(info$name, " not installed. Skipping...")
          next
        }
      }

      # Load the dataset
      if (verbose) message("Loading ", info$name, "...")
      data(list = info$name)
      loaded_data[[dataset]] <- get(info$name)

      # Print basic information
      if (verbose) {
        obj <- loaded_data[[dataset]]
        message("  Cells: ", ncol(obj))
        message("  Features: ", nrow(obj))
        message("  Assays: ", paste(Assays(obj), collapse = ", "))

        # Check for spatial information
        if (length(Images(obj)) > 0) {
          message("  Images: ", paste(Images(obj), collapse = ", "))

          # Get coordinates for first image
          coords <- GetTissueCoordinates(obj[[Images(obj)[1]]])
          message("  Spatial spots: ", nrow(coords))
        }
      }

    } else {
      warning("Unknown dataset: ", dataset)
    }
  }

  return(loaded_data)
}

#' Test Spatial Data Conversion Round-trip
#'
#' Tests converting spatial data from Seurat to h5ad and back
#'
#' @param seurat_obj Seurat object with spatial data
#' @param temp_dir Directory for temporary files
#' @param verbose Print progress messages
#'
#' @return List with original and converted objects for comparison
#' @export

test_spatial_roundtrip <- function(seurat_obj,
                                   temp_dir = tempdir(),
                                   verbose = TRUE) {

  # Paths for conversion
  h5seurat_file <- file.path(temp_dir, "spatial_test.h5seurat")
  h5ad_file <- file.path(temp_dir, "spatial_test.h5ad")
  h5seurat_back <- file.path(temp_dir, "spatial_test_back.h5seurat")

  if (verbose) {
    message("Testing spatial data round-trip conversion...")
    message("Original object summary:")
    message("  Cells: ", ncol(seurat_obj))
    message("  Features: ", nrow(seurat_obj))
    message("  Images: ", paste(Images(seurat_obj), collapse = ", "))
  }

  # Step 1: Save as h5Seurat
  if (verbose) message("\n1. Saving as h5Seurat...")
  SaveH5Seurat(seurat_obj, filename = h5seurat_file, overwrite = TRUE)

  # Step 2: Convert to h5ad
  if (verbose) message("2. Converting to h5ad...")
  Convert(h5seurat_file, dest = h5ad_file, overwrite = TRUE)

  # Step 3: Convert back to h5Seurat
  if (verbose) message("3. Converting back to h5Seurat...")
  Convert(h5ad_file, dest = h5seurat_back, overwrite = TRUE)

  # Step 4: Load the round-tripped data
  if (verbose) message("4. Loading round-tripped data...")
  seurat_roundtrip <- LoadH5Seurat(h5seurat_back)

  # Step 5: Compare spatial data
  if (verbose) {
    message("\nComparison:")
    message("  Original cells: ", ncol(seurat_obj))
    message("  Roundtrip cells: ", ncol(seurat_roundtrip))
    message("  Original features: ", nrow(seurat_obj))
    message("  Roundtrip features: ", nrow(seurat_roundtrip))

    # Check spatial coordinates
    if (length(Images(seurat_obj)) > 0) {
      orig_coords <- GetTissueCoordinates(seurat_obj[[Images(seurat_obj)[1]]])

      if (length(Images(seurat_roundtrip)) > 0) {
        rt_coords <- GetTissueCoordinates(seurat_roundtrip[[Images(seurat_roundtrip)[1]]])

        message("  Original spatial spots: ", nrow(orig_coords))
        message("  Roundtrip spatial spots: ", nrow(rt_coords))

        # Check if coordinates match
        if (all(dim(orig_coords) == dim(rt_coords))) {
          coord_match <- all.equal(orig_coords, rt_coords, check.attributes = FALSE)
          if (isTRUE(coord_match)) {
            message("  ✓ Spatial coordinates preserved!")
          } else {
            warning("  ✗ Spatial coordinates differ: ", coord_match)
          }
        }
      } else {
        warning("  ✗ No images in roundtripped object!")
      }
    }
  }

  # Clean up temp files
  if (verbose) message("\nCleaning up temporary files...")
  file.remove(c(h5seurat_file, h5ad_file, h5seurat_back))

  return(list(
    original = seurat_obj,
    roundtrip = seurat_roundtrip
  ))
}

#' Create Synthetic Spatial Data for Testing
#'
#' Creates a minimal Seurat object with spatial information for testing
#'
#' @param n_spots Number of spatial spots
#' @param n_genes Number of genes
#' @param technology Type of spatial technology to simulate
#'
#' @return Seurat object with simulated spatial data
#' @export

create_synthetic_spatial <- function(n_spots = 100,
                                    n_genes = 200,
                                    technology = c("Visium", "SlideSeq")) {

  technology <- match.arg(technology)

  # Create expression matrix
  counts <- matrix(
    rpois(n_spots * n_genes, lambda = 5),
    nrow = n_genes,
    ncol = n_spots
  )

  rownames(counts) <- paste0("Gene", seq_len(n_genes))
  colnames(counts) <- paste0("Spot", seq_len(n_spots))

  # Create Seurat object
  obj <- CreateSeuratObject(counts = counts, assay = "Spatial")

  # Add spatial coordinates
  if (technology == "Visium") {
    # Visium has hexagonal grid pattern
    coords <- data.frame(
      imagerow = sample(1:50, n_spots, replace = TRUE),
      imagecol = sample(1:50, n_spots, replace = TRUE),
      row.names = colnames(counts)
    )

    # Add tissue detection
    coords$tissue <- sample(c(0, 1), n_spots, replace = TRUE, prob = c(0.2, 0.8))

  } else {
    # SlideSeq has more continuous coordinates
    coords <- data.frame(
      x = runif(n_spots, 0, 5000),
      y = runif(n_spots, 0, 5000),
      row.names = colnames(counts)
    )
  }

  # Add to metadata
  obj <- AddMetaData(obj, coords)

  message("Created synthetic ", technology, " data:")
  message("  Spots: ", n_spots)
  message("  Genes: ", n_genes)

  return(obj)
}