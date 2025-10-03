#!/usr/bin/env Rscript

#' Test Spatial Data Conversion
#'
#' This script tests the spatial data conversion functionality
#' between h5Seurat and h5ad formats

# Load required libraries
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(hdf5r)

# Set up test environment
test_dir <- tempdir()
message("Test directory: ", test_dir)

# Test 1: Create synthetic spatial data and test round-trip conversion
message("\n", paste(rep("=", 60), collapse = ""))
message("TEST 1: Synthetic Spatial Data Round-trip")
message(paste(rep("=", 60), collapse = ""))

# Create synthetic spatial data
create_test_spatial_data <- function(n_cells = 100, n_genes = 200) {
  # Create expression matrix
  counts <- matrix(
    rpois(n_cells * n_genes, lambda = 5),
    nrow = n_genes,
    ncol = n_cells,
    dimnames = list(
      paste0("Gene", seq_len(n_genes)),
      paste0("Cell", seq_len(n_cells))
    )
  )

  # Create Seurat object
  obj <- CreateSeuratObject(counts = counts, assay = "Spatial")

  # Add spatial coordinates to metadata
  obj@meta.data$spatial_x <- runif(n_cells, 0, 100)
  obj@meta.data$spatial_y <- runif(n_cells, 0, 100)

  # Add some additional metadata
  obj@meta.data$tissue <- sample(c(0, 1), n_cells, replace = TRUE, prob = c(0.2, 0.8))

  return(obj)
}

# Create test object
test_obj <- create_test_spatial_data(n_cells = 50, n_genes = 100)
message("Created test object with ", ncol(test_obj), " cells and ", nrow(test_obj), " genes")

# Define file paths
h5seurat_file <- file.path(test_dir, "test_spatial.h5seurat")
h5ad_file <- file.path(test_dir, "test_spatial.h5ad")
h5seurat_back <- file.path(test_dir, "test_spatial_back.h5seurat")

# Test conversion pipeline
tryCatch({
  # Step 1: Save as h5Seurat
  message("\n1. Saving as h5Seurat...")
  SaveH5Seurat(test_obj, filename = h5seurat_file, overwrite = TRUE)
  message("   ✓ Saved to ", h5seurat_file)

  # Step 2: Convert to h5ad
  message("\n2. Converting to h5ad...")
  Convert(h5seurat_file, dest = h5ad_file, overwrite = TRUE, verbose = FALSE)
  message("   ✓ Converted to ", h5ad_file)

  # Step 3: Validate h5ad structure
  message("\n3. Validating h5ad structure...")
  h5ad <- H5File$new(h5ad_file, mode = "r")

  # Check for obsm['spatial']
  has_obsm_spatial <- FALSE
  if (h5ad$exists("obsm") && "spatial" %in% names(h5ad[["obsm"]])) {
    has_obsm_spatial <- TRUE
    spatial_coords <- h5ad[["obsm"]][["spatial"]][,]
    message("   ✓ Found spatial coordinates in obsm['spatial']")
    message("     Shape: ", paste(dim(spatial_coords), collapse = " x "))
  } else {
    message("   ✗ No spatial coordinates in obsm['spatial']")
  }

  # Check for uns['spatial']
  has_uns_spatial <- FALSE
  if (h5ad$exists("uns") && "spatial" %in% names(h5ad[["uns"]])) {
    has_uns_spatial <- TRUE
    message("   ✓ Found spatial metadata in uns['spatial']")

    spatial_uns <- h5ad[["uns/spatial"]]
    lib_ids <- names(spatial_uns)
    message("     Libraries: ", paste(lib_ids, collapse = ", "))
  } else {
    message("   ✗ No spatial metadata in uns['spatial']")
  }

  h5ad$close_all()

  # Step 4: Load h5ad back
  message("\n4. Loading h5ad file...")
  loaded_obj <- LoadH5AD(h5ad_file, verbose = FALSE)
  message("   ✓ Loaded object with ", ncol(loaded_obj), " cells")

  # Check if spatial coordinates were preserved
  if ("spatial_x" %in% colnames(loaded_obj@meta.data) &&
      "spatial_y" %in% colnames(loaded_obj@meta.data)) {
    message("   ✓ Spatial coordinates preserved in metadata")

    # Compare coordinates
    orig_x <- test_obj@meta.data$spatial_x
    loaded_x <- loaded_obj@meta.data$spatial_x[rownames(test_obj@meta.data)]

    coord_match <- all.equal(orig_x, loaded_x, tolerance = 1e-6)
    if (isTRUE(coord_match)) {
      message("   ✓ Coordinates match original data")
    } else {
      message("   ✗ Coordinates differ from original")
    }
  } else {
    message("   ✗ Spatial coordinates not found in loaded object")
  }

  message("\n✓ Test 1 PASSED")

}, error = function(e) {
  message("\n✗ Test 1 FAILED: ", e$message)
})

# Test 2: Test with SeuratData spatial datasets (if available)
message("\n", paste(rep("=", 60), collapse = ""))
message("TEST 2: SeuratData Spatial Dataset (if available)")
message(paste(rep("=", 60), collapse = ""))

if (requireNamespace("SeuratData", quietly = TRUE)) {
  library(SeuratData)

  # Check for available spatial datasets
  available_data <- AvailableData()
  spatial_datasets <- available_data[grep("spatial|visium|slide",
                                          available_data$Summary,
                                          ignore.case = TRUE), ]

  if (nrow(spatial_datasets) > 0) {
    message("Found ", nrow(spatial_datasets), " spatial dataset(s)")

    # Try to load stxBrain if available
    if ("stxBrain" %in% InstalledData()$Dataset) {
      message("\nTesting with stxBrain dataset...")

      data("stxBrain")
      brain_obj <- stxBrain

      # Test conversion
      brain_h5seurat <- file.path(test_dir, "brain.h5seurat")
      brain_h5ad <- file.path(test_dir, "brain.h5ad")

      tryCatch({
        # Save and convert
        SaveH5Seurat(brain_obj, filename = brain_h5seurat, overwrite = TRUE)
        Convert(brain_h5seurat, dest = brain_h5ad, overwrite = TRUE, verbose = FALSE)

        # Validate
        h5ad <- H5File$new(brain_h5ad, mode = "r")

        if (h5ad$exists("obsm") && "spatial" %in% names(h5ad[["obsm"]])) {
          message("   ✓ Spatial data successfully converted for Visium dataset")
        } else {
          message("   ✗ Spatial data not found in converted file")
        }

        h5ad$close_all()
        message("\n✓ Test 2 PASSED")

      }, error = function(e) {
        message("\n✗ Test 2 FAILED: ", e$message)
      })

    } else {
      message("stxBrain dataset not installed. Install with: InstallData('stxBrain')")
      message("Skipping Test 2")
    }
  } else {
    message("No spatial datasets found in SeuratData")
    message("Skipping Test 2")
  }
} else {
  message("SeuratData not available")
  message("Install with: remotes::install_github('satijalab/seurat-data')")
  message("Skipping Test 2")
}

# Test 3: Test multimodal spatial data
message("\n", paste(rep("=", 60), collapse = ""))
message("TEST 3: Multimodal Spatial Data")
message(paste(rep("=", 60), collapse = ""))

# Create multimodal spatial object
create_multimodal_spatial <- function(n_cells = 50) {
  # RNA assay
  rna_counts <- matrix(
    rpois(n_cells * 100, lambda = 5),
    nrow = 100,
    ncol = n_cells,
    dimnames = list(paste0("Gene", 1:100), paste0("Cell", 1:n_cells))
  )

  # ADT assay (protein)
  adt_counts <- matrix(
    rpois(n_cells * 20, lambda = 10),
    nrow = 20,
    ncol = n_cells,
    dimnames = list(paste0("Protein", 1:20), paste0("Cell", 1:n_cells))
  )

  # Create object
  obj <- CreateSeuratObject(counts = rna_counts, assay = "RNA")
  obj[["ADT"]] <- CreateAssayObject(counts = adt_counts)

  # Add spatial coordinates
  obj@meta.data$spatial_x <- runif(n_cells, 0, 100)
  obj@meta.data$spatial_y <- runif(n_cells, 0, 100)

  return(obj)
}

mm_obj <- create_multimodal_spatial(n_cells = 30)
message("Created multimodal object with RNA and ADT assays")

# Test conversion
mm_h5seurat <- file.path(test_dir, "multimodal.h5seurat")
mm_h5ad <- file.path(test_dir, "multimodal.h5ad")

tryCatch({
  SaveH5Seurat(mm_obj, filename = mm_h5seurat, overwrite = TRUE)
  Convert(mm_h5seurat, dest = mm_h5ad, assay = "RNA", overwrite = TRUE, verbose = FALSE)

  # Check result
  h5ad <- H5File$new(mm_h5ad, mode = "r")

  has_spatial <- h5ad$exists("obsm") && "spatial" %in% names(h5ad[["obsm"]])
  has_x <- h5ad$exists("X")

  h5ad$close_all()

  if (has_spatial && has_x) {
    message("   ✓ Multimodal spatial data successfully converted")
    message("\n✓ Test 3 PASSED")
  } else {
    message("   ✗ Missing data in conversion")
    message("\n✗ Test 3 FAILED")
  }

}, error = function(e) {
  message("\n✗ Test 3 FAILED: ", e$message)
})

# Summary
message("\n", paste(rep("=", 60), collapse = ""))
message("SPATIAL CONVERSION TEST SUMMARY")
message(paste(rep("=", 60), collapse = ""))
message("All basic spatial conversion tests completed")
message("Check above for individual test results")

# Clean up
message("\nCleaning up test files...")
test_files <- list.files(test_dir, pattern = "\\.(h5seurat|h5ad)$", full.names = TRUE)
file.remove(test_files)
message("Removed ", length(test_files), " test files")