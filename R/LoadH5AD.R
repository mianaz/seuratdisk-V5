#' Load an AnnData H5AD file as a Seurat object
#'
#' Direct conversion from H5AD format to Seurat object without intermediate h5Seurat
#'
#' @param file Path to H5AD file
#' @param assay.name Name for the primary assay (default: "RNA")
#' @param verbose Show progress messages
#'
#' @return A \code{Seurat} object
#'
#' @importFrom hdf5r H5File h5attr h5attr_names
#' @importFrom Matrix sparseMatrix
#' @importFrom Seurat CreateSeuratObject SetAssayData CreateDimReducObject
#' @importFrom SeuratObject Cells
#'
#' @export
#'
LoadH5AD <- function(file, assay.name = "RNA", verbose = TRUE) {
  if (!file.exists(file)) {
    stop("File not found: ", file, call. = FALSE)
  }

  h5ad <- H5File$new(file, mode = "r")
  on.exit(h5ad$close_all())

  if (verbose) {
    message("Loading H5AD file: ", file)
  }

  # Helper function to read H5AD sparse or dense matrix
  ReadH5ADMatrix <- function(h5_obj, transpose = TRUE) {
    if (inherits(h5_obj, "H5Group")) {
      # Sparse matrix (CSR format in h5ad)
      if (h5_obj$exists("data") && h5_obj$exists("indices") && h5_obj$exists("indptr")) {
        data <- h5_obj[["data"]][]
        indices <- h5_obj[["indices"]][] + 1L  # Convert to 1-based indexing
        indptr <- h5_obj[["indptr"]][]

        # Get dimensions from shape attribute
        if (h5_obj$attr_exists("shape")) {
          shape <- h5attr(h5_obj, "shape")
          n_rows <- shape[1]
          n_cols <- shape[2]
        } else {
          n_rows <- length(indptr) - 1
          n_cols <- max(indices)
        }

        # Create sparse matrix from CSR format
        row_indices <- rep(seq_len(n_rows), diff(indptr))
        mat <- sparseMatrix(
          i = row_indices,
          j = indices,
          x = data,
          dims = c(n_rows, n_cols),
          index1 = TRUE
        )

        if (transpose) {
          # h5ad stores as cells x genes, Seurat needs genes x cells
          mat <- t(mat)
        }
        return(mat)
      }
    } else if (inherits(h5_obj, "H5D")) {
      # Dense matrix
      mat <- h5_obj[,]
      if (transpose) {
        mat <- t(mat)
      }
      return(mat)
    }
    stop("Unknown matrix format", call. = FALSE)
  }

  # 1. Read cell names
  if (verbose) message("Reading cell names...")
  cell.names <- NULL
  if (h5ad$exists("obs")) {
    obs_group <- h5ad[["obs"]]
    if (obs_group$exists("_index")) {
      cell.names <- as.character(obs_group[["_index"]][])
    } else if (obs_group$exists("index")) {
      cell.names <- as.character(obs_group[["index"]][])
    }
  }

  # Generate cell names if not found
  if (is.null(cell.names)) {
    n_cells <- if (h5ad$exists("X")) {
      if (inherits(h5ad[["X"]], "H5Group") && h5ad[["X"]]$attr_exists("shape")) {
        h5attr(h5ad[["X"]], "shape")[1]
      } else if (inherits(h5ad[["X"]], "H5D")) {
        h5ad[["X"]]$dims[1]
      }
    } else {
      stop("Cannot determine number of cells", call. = FALSE)
    }
    cell.names <- paste0("Cell", seq_len(n_cells))
  }

  # 2. Read feature names
  if (verbose) message("Reading feature names...")
  feature.names <- NULL
  if (h5ad$exists("var")) {
    var_group <- h5ad[["var"]]
    if (var_group$exists("_index")) {
      feature.names <- as.character(var_group[["_index"]][])
    } else if (var_group$exists("index")) {
      feature.names <- as.character(var_group[["index"]][])
    }
  }

  # Generate feature names if not found
  if (is.null(feature.names)) {
    n_features <- if (h5ad$exists("X")) {
      if (inherits(h5ad[["X"]], "H5Group") && h5ad[["X"]]$attr_exists("shape")) {
        h5attr(h5ad[["X"]], "shape")[2]
      } else if (inherits(h5ad[["X"]], "H5D")) {
        h5ad[["X"]]$dims[2]
      }
    } else {
      stop("Cannot determine number of features", call. = FALSE)
    }
    feature.names <- paste0("Gene", seq_len(n_features))
  }

  # 3. Read main expression matrix
  if (verbose) message("Reading expression matrix...")
  if (!h5ad$exists("X")) {
    stop("No expression matrix (X) found in h5ad file", call. = FALSE)
  }

  # Determine if we need to transpose by checking matrix orientation
  # Standard h5ad: X is (n_obs x n_var) = (cells x features)
  # Seurat needs: (n_features x n_cells) = (features x cells)
  X_obj <- h5ad[["X"]]

  # Get actual matrix dimensions
  if (inherits(X_obj, "H5Group") && X_obj$attr_exists("shape")) {
    # Sparse matrix
    X_dims <- h5attr(X_obj, "shape")
    X_nrows <- X_dims[1]
    X_ncols <- X_dims[2]
  } else if (inherits(X_obj, "H5D")) {
    # Dense matrix
    X_nrows <- X_obj$dims[1]
    X_ncols <- X_obj$dims[2]
  } else {
    stop("Cannot determine X matrix dimensions", call. = FALSE)
  }

  # Compare with obs/var counts to determine orientation
  n_obs <- length(cell.names)
  n_var <- length(feature.names)

  # Decide whether to transpose
  # If X is already (features x cells), don't transpose
  # If X is (cells x features), transpose
  needs_transpose <- TRUE
  if (X_nrows == n_var && X_ncols == n_obs) {
    # Matrix is already (features x cells)
    needs_transpose <- FALSE
    if (verbose) message("  Matrix is already in (features x cells) format, no transpose needed")
  } else if (X_nrows == n_obs && X_ncols == n_var) {
    # Matrix is (cells x features), needs transpose
    needs_transpose <- TRUE
    if (verbose) message("  Matrix is in (cells x features) format, transposing to (features x cells)")
  } else {
    warning("Matrix dimensions don't match obs/var counts. Using default transpose.", immediate. = TRUE)
  }

  expr_matrix <- ReadH5ADMatrix(X_obj, transpose = needs_transpose)

  # Handle dimension mismatches
  if (nrow(expr_matrix) != length(feature.names)) {
    warning("Adjusting feature names to match matrix dimensions", immediate. = TRUE)
    feature.names <- feature.names[seq_len(nrow(expr_matrix))]
  }
  if (ncol(expr_matrix) != length(cell.names)) {
    warning("Adjusting cell names to match matrix dimensions", immediate. = TRUE)
    cell.names <- cell.names[seq_len(ncol(expr_matrix))]
  }

  # Handle NA values in names
  if (any(is.na(feature.names))) {
    warning("NA values detected in feature names. Replacing with generic names.", immediate. = TRUE)
    na_idx <- which(is.na(feature.names))
    feature.names[na_idx] <- paste0("Feature", na_idx)
  }
  if (any(is.na(cell.names))) {
    warning("NA values detected in cell names. Replacing with generic names.", immediate. = TRUE)
    na_idx <- which(is.na(cell.names))
    cell.names[na_idx] <- paste0("Cell", na_idx)
  }

  # Make feature and cell names unique to avoid rowname/colname conflicts
  if (any(duplicated(feature.names))) {
    warning("Duplicate feature names detected. Making unique.", immediate. = TRUE)
    feature.names <- make.unique(feature.names)
  }
  if (any(duplicated(cell.names))) {
    warning("Duplicate cell names detected. Making unique.", immediate. = TRUE)
    cell.names <- make.unique(cell.names)
  }

  rownames(expr_matrix) <- feature.names
  colnames(expr_matrix) <- cell.names

  # 4. Create Seurat object
  if (verbose) message("Creating Seurat object...")
  seurat_obj <- CreateSeuratObject(
    counts = expr_matrix,
    project = "H5AD",
    assay = assay.name,
    min.cells = 0,
    min.features = 0
  )

  # 5. Add raw counts if present
  if (h5ad$exists("raw") && h5ad[["raw"]]$exists("X")) {
    if (verbose) message("Adding raw counts...")

    # Get raw cell names
    raw_cells <- NULL
    if (h5ad[["raw"]]$exists("obs")) {
      raw_obs <- h5ad[["raw/obs"]]
      if (raw_obs$exists("_index")) {
        raw_cells <- as.character(raw_obs[["_index"]][])
      } else if (raw_obs$exists("index")) {
        raw_cells <- as.character(raw_obs[["index"]][])
      }
    }

    raw_features <- NULL
    if (h5ad[["raw"]]$exists("var")) {
      raw_var <- h5ad[["raw/var"]]
      if (raw_var$exists("_index")) {
        raw_features <- as.character(raw_var[["_index"]][])
      } else if (raw_var$exists("index")) {
        raw_features <- as.character(raw_var[["index"]][])
      }
    }

    if (!is.null(raw_features)) {
      raw_matrix <- ReadH5ADMatrix(h5ad[["raw/X"]], transpose = TRUE)

      # Match dimensions
      raw_features <- raw_features[seq_len(min(length(raw_features), nrow(raw_matrix)))]
      rownames(raw_matrix) <- raw_features

      # Set column names based on raw_cells or match to cell.names
      if (!is.null(raw_cells) && ncol(raw_matrix) == length(raw_cells)) {
        colnames(raw_matrix) <- raw_cells
        # Find common cells between raw and filtered
        common_cells <- intersect(cell.names, raw_cells)
        if (length(common_cells) > 0) {
          # Subset raw matrix to match current cells
          raw_matrix <- raw_matrix[, common_cells, drop = FALSE]
        } else {
          # No overlap, skip raw counts
          if (verbose) message("  No overlapping cells between raw and main matrix. Skipping raw counts.")
          raw_matrix <- NULL
        }
      } else if (ncol(raw_matrix) == length(cell.names)) {
        # Same number of cells, assume same order
        colnames(raw_matrix) <- cell.names
      } else {
        # Dimension mismatch and no cell names - skip
        if (verbose) message("  Raw matrix dimensions don't match. Skipping raw counts.")
        raw_matrix <- NULL
      }

      # Add raw counts to assay if we have a valid matrix
      if (!is.null(raw_matrix)) {
        # Find common features
        common_features <- intersect(feature.names, rownames(raw_matrix))
        if (length(common_features) > 0) {
          # Ensure we're only using cells that exist in the Seurat object
          common_cells_final <- intersect(cell.names, colnames(raw_matrix))
          if (length(common_cells_final) > 0) {
            raw_subset <- raw_matrix[common_features, common_cells_final, drop = FALSE]
            seurat_obj[[assay.name]] <- SetAssayData(
              object = seurat_obj[[assay.name]],
              layer = "counts",
              new.data = raw_subset
            )
          }
        }
      }
    }
  }

  # 6. Add layers if present
  if (h5ad$exists("layers")) {
    if (verbose) message("Adding layers...")
    layer_names <- names(h5ad[["layers"]])

    for (layer_name in layer_names) {
      if (verbose) message("  Adding layer: ", layer_name)

      layer_matrix <- ReadH5ADMatrix(h5ad[["layers"]][[layer_name]], transpose = TRUE)

      # Ensure dimensions match
      if (nrow(layer_matrix) == nrow(expr_matrix) && ncol(layer_matrix) == ncol(expr_matrix)) {
        rownames(layer_matrix) <- feature.names
        colnames(layer_matrix) <- cell.names

        # Map layer names to Seurat slots
        seurat_slot <- switch(layer_name,
          "counts" = "counts",
          "data" = "data",
          "log_normalized" = "data",
          "scale.data" = "scale.data",
          "scaled" = "scale.data",
          layer_name
        )

        tryCatch({
          seurat_obj[[assay.name]] <- SetAssayData(
            object = seurat_obj[[assay.name]],
            layer = seurat_slot,
            new.data = layer_matrix
          )
        }, error = function(e) {
          if (verbose) warning("Could not add layer ", layer_name, ": ", e$message, immediate. = TRUE)
        })
      }
    }
  }

  # 7. Add cell metadata
  if (h5ad$exists("obs")) {
    if (verbose) message("Adding cell metadata...")
    obs_group <- h5ad[["obs"]]
    obs_cols <- setdiff(names(obs_group), c("_index", "index", "__categories"))

    for (col in obs_cols) {
      if (verbose) message("  Adding metadata: ", col)

      tryCatch({
        # Skip if it's a group (not a dataset)
        if (!inherits(obs_group[[col]], "H5D")) {
          if (verbose) message("    Skipping (not a dataset)")
          next
        }

        # Check if categorical
        is_categorical <- FALSE
        if (obs_group$exists("__categories") && col %in% names(obs_group[["__categories"]])) {
          is_categorical <- TRUE
          codes <- obs_group[[col]][]
          categories <- as.character(obs_group[["__categories"]][[col]][])
          # h5ad uses -1 for NA in categorical
          codes[codes == -1] <- NA
          meta_values <- factor(categories[codes + 1], levels = categories)
        } else {
          # Numeric or string
          meta_values <- obs_group[[col]][]
          if (is.character(meta_values)) {
            meta_values <- as.character(meta_values)
          }
        }

        # Add to Seurat object
        seurat_obj[[col]] <- meta_values
      }, error = function(e) {
        if (verbose) warning("  Could not add metadata ", col, ": ", e$message, immediate. = TRUE)
      })
    }
  }

  # 8. Add dimensional reductions
  if (h5ad$exists("obsm")) {
    if (verbose) message("Adding dimensional reductions...")

    for (reduc_name in names(h5ad[["obsm"]])) {
      clean_name <- gsub("^X_", "", reduc_name)
      if (verbose) message("  Adding reduction: ", clean_name)

      embeddings <- h5ad[["obsm"]][[reduc_name]][,]

      # Check dimensions
      if (nrow(embeddings) != length(cell.names)) {
        warning("Skipping reduction ", clean_name, " - dimension mismatch", immediate. = TRUE)
        next
      }

      rownames(embeddings) <- cell.names
      colnames(embeddings) <- paste0(toupper(clean_name), "_", seq_len(ncol(embeddings)))

      # Determine key
      key <- switch(clean_name,
        "pca" = "PC_",
        "tsne" = "tSNE_",
        "umap" = "UMAP_",
        paste0(toupper(clean_name), "_")
      )

      reduc_obj <- CreateDimReducObject(
        embeddings = embeddings,
        key = key,
        assay = assay.name
      )

      seurat_obj[[clean_name]] <- reduc_obj
    }
  }

  # 9. Add feature metadata
  if (h5ad$exists("var")) {
    if (verbose) message("Adding feature metadata...")
    var_group <- h5ad[["var"]]
    var_cols <- setdiff(names(var_group), c("_index", "index", "__categories"))

    for (col in var_cols) {
      if (verbose) message("  Adding feature metadata: ", col)

      tryCatch({
        # Check if categorical
        if (var_group$exists("__categories") && col %in% names(var_group[["__categories"]])) {
          codes <- var_group[[col]][]
          categories <- as.character(var_group[["__categories"]][[col]][])
          codes[codes == -1] <- NA
          meta_values <- factor(categories[codes + 1], levels = categories)
        } else {
          meta_values <- var_group[[col]][]
        }

        # Special handling for highly_variable
        if (col == "highly_variable" && is.numeric(meta_values)) {
          # Convert to logical if stored as numeric
          meta_values <- as.logical(meta_values)
        }

        # Ensure length matches
        if (length(meta_values) == nrow(seurat_obj)) {
          seurat_obj[[assay.name]][[col]] <- meta_values

          # Set variable features if highly_variable column exists
          if (col == "highly_variable" && is.logical(meta_values)) {
            VariableFeatures(seurat_obj) <- feature.names[meta_values]
          }
        }
      }, error = function(e) {
        if (verbose) warning("Could not add feature metadata ", col, ": ", e$message, immediate. = TRUE)
      })
    }
  }

  # 10. Add neighbor graphs from obsp
  if (h5ad$exists("obsp")) {
    if (verbose) message("Adding neighbor graphs...")

    for (graph_name in names(h5ad[["obsp"]])) {
      if (verbose) message("  Adding graph: ", graph_name)

      graph_matrix <- ReadH5ADMatrix(h5ad[["obsp"]][[graph_name]], transpose = FALSE)

      # Ensure square matrix with correct dimensions
      if (nrow(graph_matrix) == ncol(graph_matrix) && nrow(graph_matrix) == ncol(seurat_obj)) {
        rownames(graph_matrix) <- cell.names
        colnames(graph_matrix) <- cell.names

        # Map to Seurat graph names
        seurat_graph_name <- switch(graph_name,
          "connectivities" = paste0(assay.name, "_snn"),
          "distances" = paste0(assay.name, "_nn"),
          graph_name
        )

        seurat_obj@graphs[[seurat_graph_name]] <- as.Graph(graph_matrix)
      }
    }
  }

  # 11. Add uns (unstructured) data to misc
  if (h5ad$exists("uns")) {
    if (verbose) message("Adding unstructured data...")
    uns_group <- h5ad[["uns"]]

    for (item in names(uns_group)) {
      tryCatch({
        if (inherits(uns_group[[item]], "H5D")) {
          # Simple dataset
          seurat_obj@misc[[item]] <- uns_group[[item]][]
        } else if (inherits(uns_group[[item]], "H5Group")) {
          # Complex group - store as list
          if (verbose) message("  Storing complex uns item: ", item)
          # For now, just note it exists
          seurat_obj@misc[[paste0(item, "_present")]] <- TRUE
        }
      }, error = function(e) {
        if (verbose) warning("Could not add uns item ", item, ": ", e$message, immediate. = TRUE)
      })
    }
  }

  # 12. Add spatial data support
  if (h5ad$exists("obsm") && "spatial" %in% names(h5ad[["obsm"]])) {
    seurat_obj <- ConvertH5ADSpatialToSeurat(
      h5ad_file = h5ad,
      seurat_obj = seurat_obj,
      assay_name = assay.name,
      verbose = verbose
    )
  }

  if (verbose) {
    message("\nSuccessfully loaded H5AD file")
    message("  Cells: ", ncol(seurat_obj))
    message("  Features: ", nrow(seurat_obj))
    message("  Assays: ", paste(names(seurat_obj@assays), collapse = ", "))
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
