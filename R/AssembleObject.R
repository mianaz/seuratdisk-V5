#' @include zzz.R
#' @importFrom hdf5r h5attr
#' @importFrom methods slot<- new
#' @importFrom Seurat Cells Key<- Key Cells scalefactors
#' @importFrom SeuratObject CreateFOV CreateCentroids
#'
NULL

#' Read sparse matrix from HDF5
#'
#' Internal function to read sparse matrix data from h5Seurat files,
#' handling both V4 and V5 formats
#'
#' @param h5_group The HDF5 group containing sparse matrix data
#' @param verbose Show progress updates
#'
#' @return A matrix object, or NULL if reading fails
#'
#' @keywords internal
#'
ReadSparseMatrix <- function(h5_group, verbose = FALSE) {
  # Helper for safe existence checking
  safe_exists_local <- function(group, path) {
    tryCatch(
      expr = group$exists(name = path),
      error = function(e) path %in% names(group)
    )
  }

  # Check for CSR/CSC sparse matrix components
  if (safe_exists_local(h5_group, "data") && safe_exists_local(h5_group, "indices") && safe_exists_local(h5_group, "indptr")) {
    if (verbose) {
      message("Reading sparse matrix in CSR/CSC format")
    }

    # Read components
    data_vals <- h5_group[["data"]][]
    indices <- h5_group[["indices"]][]
    indptr <- h5_group[["indptr"]][]

    # Get dimensions - check for shape attribute first (V5 style)
    if (h5_group$attr_exists("shape")) {
      shape <- h5attr(x = h5_group, which = "shape")
      nrows <- shape[1]
      ncols <- shape[2]
    } else if (h5_group$attr_exists("dims")) {
      # Fall back to dims attribute (older format)
      dims <- h5attr(x = h5_group, which = "dims")
      nrows <- dims[1]
      ncols <- dims[2]
    } else {
      # Infer dimensions from the data
      ncols <- length(indptr) - 1L
      nrows <- max(indices) + 1L  # 0-based indexing
      if (verbose) {
        message("Warning: No shape/dims attribute found, inferring dimensions: ", nrows, " x ", ncols)
      }
    }

    # Detect format (CSR vs CSC) based on indptr length
    if (length(indptr) - 1 == ncols) {
      # CSC format: indptr indexes columns
      if (verbose) {
        message("Detected CSC format")
      }
      sparse_mat <- Matrix::sparseMatrix(
        i = indices + 1L,  # Convert from 0-based to 1-based
        p = indptr,
        x = data_vals,
        dims = c(nrows, ncols),
        index1 = TRUE
      )
    } else if (length(indptr) - 1 == nrows) {
      # CSR format: indptr indexes rows - need to transpose
      if (verbose) {
        message("Detected CSR format, will transpose")
      }
      sparse_mat <- Matrix::sparseMatrix(
        j = indices + 1L,  # Convert from 0-based to 1-based
        p = indptr,
        x = data_vals,
        dims = c(ncols, nrows),  # Note: swapped dimensions for transpose
        index1 = TRUE
      )
      sparse_mat <- Matrix::t(sparse_mat)
    } else {
      stop("Cannot determine sparse matrix format from indptr length")
    }

    return(sparse_mat)  # Keep as sparse!
  }

  # Check for dgCMatrix components (older SeuratDisk format)
  if (safe_exists_local(h5_group, "i") && safe_exists_local(h5_group, "p") && safe_exists_local(h5_group, "x")) {
    if (verbose) {
      message("Reading sparse matrix in dgCMatrix format")
    }

    i_vals <- h5_group[["i"]][]
    p_vals <- h5_group[["p"]][]
    x_vals <- h5_group[["x"]][]

    # Get dimensions
    if (h5_group$attr_exists("Dim")) {
      dims <- h5attr(x = h5_group, which = "Dim")
    } else {
      dims <- c(max(i_vals) + 1L, length(p_vals) - 1L)
    }

    sparse_mat <- Matrix::sparseMatrix(
      i = i_vals + 1L,  # Convert to 1-based if needed
      p = p_vals,
      x = x_vals,
      dims = dims,
      index1 = TRUE
    )

    return(sparse_mat)  # Keep as sparse!
  }

  # If it's not a sparse matrix group, try reading as dense matrix
  if (verbose) {
    message("No sparse matrix components found, attempting dense matrix read")
  }

  tryCatch({
    # Read the matrix data
    mat_data <- h5_group[,]

    # Check if we should convert dense to sparse
    if (is.matrix(mat_data)) {
      zero_proportion <- sum(mat_data == 0) / length(mat_data)
      if (zero_proportion > 0.5) {
        if (verbose) {
          message("Converting dense to sparse (", round(zero_proportion * 100, 1), "% zeros)")
        }
        return(as(mat_data, "dgCMatrix"))
      }
    }

    return(mat_data)
  }, error = function(e) {
    if (verbose) {
      message("Failed to read matrix: ", conditionMessage(e))
    }
    return(NULL)
  })
}

#' Assemble an object from an h5Seurat file
#'
#' @param assay,reduction,graph,image,neighbor,cmd Name of assay, reduction,
#' graph, image, neighbor, or command to load
#' @param file A connected h5Seurat file to pull the data from
#' @param verbose Show progress updates
#'
#' @name AssembleObject
#' @rdname AssembleObject
#'
#' @keywords internal
#'
NULL

#' @param slots Optional vector of assay slots to load, defaults to all slots
#' present in assay
#'
#' @importFrom Seurat CreateAssayObject GetAssayData SetAssayData AddMetaData
#' VariableFeatures<-
#'
#' @return \code{AssembleAssay}: An \code{Assay} object
#'
#' @rdname AssembleObject
#'
#' @aliases AssembleAssay
#'
AssembleAssay <- function(assay, file, slots = NULL, verbose = TRUE) {
  index <- file$index()
  if (!assay %in% names(x = index)) {
    stop("Cannot find assay ", assay, " in this h5Seurat file", call. = FALSE)
  }
  slots.assay <- names(x = Filter(f = isTRUE, x = index[[assay]]$slots))
  slots <- slots %||% slots.assay
  slots <- match.arg(arg = slots, choices = slots.assay, several.ok = TRUE)
  if (!any(c('counts', 'data') %in% slots)) {
    stop("At least one of 'counts' or 'data' must be loaded", call. = FALSE)
  }
  assay.group <- file[['assays']][[assay]]

  # Helper function for safe path checking
  safe_exists <- function(group, path) {
    tryCatch(
      expr = {
        group$exists(name = path)
      },
      error = function(e) {
        # If exists() fails, try checking names directly
        path %in% names(group)
      }
    )
  }

  # Helper function to read datasets that might be 2D (V5 compatibility)
  safe_read_dataset <- function(dataset, dataset_name = "") {
    if (length(dataset$dims) > 1) {
      # V5: 2D datasets are often dummy data, check for real data elsewhere
      if (dataset_name == "features" && safe_exists(assay.group, "meta.data/_index")) {
        # V5 feature names are in meta.data/_index
        return(as.character(assay.group[["meta.data/_index"]][]))
      } else {
        # Fallback: read first column and ensure it's character
        result <- SafeH5DRead(dataset)
        return(as.character(result))
      }
    } else {
      # 1D dataset - standard read
      return(as.character(dataset[]))
    }
  }

  # CRITICAL: Get features FIRST before they're used in read_matrix_data
  if (safe_exists(assay.group, "meta.data/_index")) {
    # V5: Use the full feature index for proper sparse matrix reconstruction
    features <- FixFeatures(features = as.character(assay.group[["meta.data/_index"]][]))
    if (verbose) {
      message("Using V5 full feature space: ", length(features), " features")
    }
  } else {
    # V4 or direct features
    features <- FixFeatures(features = safe_read_dataset(assay.group[['features']], "features"))
    if (verbose) {
      message("Using direct features: ", length(features), " features")
    }
  }

  # Helper function to read matrices, handling both direct and V5 sparse layers
  read_matrix_data <- function(slot_name) {
    # Use ReadV5Layer if available (handles both V4 and V5 structures)
    if (exists("ReadV5Layer", envir = asNamespace("SeuratDisk"))) {
      mat <- ReadV5Layer(
        h5_group = assay.group,
        layer_name = slot_name,
        features = features,
        cells = Cells(x = file),
        verbose = verbose
      )
      if (!is.null(mat)) {
        return(mat)
      }
    }

    # Fallback to direct reading if ReadV5Layer not available or returns NULL
    # First try direct matrix (V4 style)
    if (safe_exists(assay.group, slot_name)) {
      if (verbose) {
        message("Reading direct matrix for slot '", slot_name, "'")
      }
      # Use ReadSparseMatrix which handles both sparse and dense formats
      mat <- ReadSparseMatrix(h5_group = assay.group[[slot_name]], verbose = verbose)
      if (!is.null(mat)) {
        # Ensure proper dimensions for V5 compatibility
        if (ncol(mat) != length(Cells(x = file)) || nrow(mat) != length(features)) {
          # Check if transposed
          if (nrow(mat) == length(Cells(x = file)) && ncol(mat) == length(features)) {
            if (verbose) {
              message("Transposing matrix to match expected dimensions")
            }
            mat <- t(mat)
          }
        }
      }
      return(mat)
    }

    # Try V5 sparse layers structure
    layers_path <- paste0("layers/", slot_name)
    if (safe_exists(assay.group, layers_path)) {
      if (verbose) {
        message("Reading V5 sparse layers for slot '", slot_name, "'")
      }
      layer_group <- assay.group[[layers_path]]

      # Use ReadSparseMatrix for V5 layers
      mat <- ReadSparseMatrix(h5_group = layer_group, verbose = verbose)

      # V5 sparse matrices might need special dimension handling
      if (!is.null(mat) && length(features) > 0) {
        expected_nrows <- length(features)
        expected_ncols <- length(Cells(x = file))

        # Check dimensions and adjust if needed
        if (nrow(mat) != expected_nrows || ncol(mat) != expected_ncols) {
          # If transposed
          if (nrow(mat) == expected_ncols && ncol(mat) == expected_nrows) {
            if (verbose) {
              message("Transposing V5 matrix to match expected dimensions")
            }
            mat <- t(mat)
          } else if (nrow(mat) < expected_nrows) {
            # V5 might have filtered features, need to expand
            if (verbose) {
              message("Expanding sparse matrix to full feature space: ", nrow(mat), " -> ", expected_nrows)
            }
            # Create full matrix with zeros for missing features
            full_mat <- matrix(0, nrow = expected_nrows, ncol = ncol(mat))
            # Copy existing data (assuming first features match)
            full_mat[1:nrow(mat), ] <- mat
            mat <- full_mat
          }
        }
      }

      return(mat)
    }

    # If neither direct nor sparse layers exist, return NULL instead of error
    # This allows more graceful fallback behavior
    if (verbose) {
      message("Warning: Cannot find matrix data for slot '", slot_name, "' in assay '", assay, "'")
    }
    return(NULL)
  }
  # Load the appropriate slots - prioritize counts for initialization
  # We need to load counts first if available, then add other layers
  if ('counts' %in% slots) {
    if (verbose) {
      message("Initializing ", assay, " with counts")
    }
    counts <- read_matrix_data('counts')

    # Check if matrix was successfully read
    if (is.null(counts)) {
      stop("Failed to read counts matrix for assay '", assay, "'", call. = FALSE)
    }

    if (verbose) { message("Matrix loaded, setting row/col names") }
    rownames(x = counts) <- features
    if (verbose) { message("Rownames set") }
    colnames(x = counts) <- Cells(x = file)
    if (verbose) { message("Colnames set, creating assay object") }

    # Create V5-compatible assay by using CreateSeuratObject and extracting the assay
    if (verbose) { message("Creating V5-compatible assay object") }
    temp_seurat <- CreateSeuratObject(counts = counts, min.cells = -1, min.features = -1)
    obj <- temp_seurat[['RNA']]
    if (verbose) { message("V5 assay object extracted, class:", class(obj)) }

  } else if ('data' %in% slots) {
    # Only if counts is not available, initialize with data
    if (verbose) {
      message("Initializing ", assay, " with data (no counts available)")
    }
    data <- read_matrix_data('data')

    # Check if matrix was successfully read
    if (is.null(data)) {
      stop("Failed to read data matrix for assay '", assay, "'", call. = FALSE)
    }

    if (verbose) { message("Data matrix loaded, setting row/col names") }
    rownames(x = data) <- features
    if (verbose) { message("Rownames set") }
    colnames(x = data) <- Cells(x = file)
    if (verbose) { message("Colnames set, creating assay object") }

    # Create V5-compatible assay by using CreateSeuratObject and extracting the assay
    if (verbose) { message("Creating V5-compatible assay object") }
    temp_seurat <- CreateSeuratObject(counts = data, min.cells = -1, min.features = -1)
    obj <- temp_seurat[['RNA']]
    if (verbose) { message("V5 assay object extracted, class:", class(obj)) }
  }
  if (verbose) { message("Setting assay key") }
  tryCatch({
    Key(object = obj) <- Key(object = assay.group)
    if (verbose) { message("Key set successfully") }
  }, error = function(e) {
    if (verbose) { message("Key setting failed: ", conditionMessage(e)) }
    # Continue without setting key if it fails
  })
  # Add remaining slots/layers (V5 compatibility)
  for (slot in slots) {
    # Skip slots that were used for initial object creation
    # Only skip counts if it was used for initialization (which happens when 'counts' is in slots)
    if (slot == 'counts' && 'counts' %in% slots) {
      if (verbose) {
        message("Skipping slot '", slot, "' - already used for object initialization")
      }
      next
    }
    # Only skip data if it was used for initialization (which only happens when counts is NOT in slots)
    if (slot == 'data' && !('counts' %in% slots) && 'data' %in% slots) {
      if (verbose) {
        message("Skipping slot '", slot, "' - already used for object initialization")
      }
      next
    }
    
    # Check if this slot/layer is already populated - use layer parameter for V5
    slot_empty <- tryCatch({
      if (inherits(obj, "Assay5")) {
        # V5: Use layer parameter
        IsMatrixEmpty(x = GetAssayData(object = obj, layer = slot))
      } else {
        # V4: Use slot parameter  
        IsMatrixEmpty(x = GetAssayData(object = obj, slot = slot))
      }
    }, error = function(e) {
      # If we can't check, assume empty
      TRUE
    })
    
    if (slot_empty) {
      if (verbose) {
        message("Adding ", slot, " for ", assay)
      }
      # Try to read the matrix data, skip if it doesn't exist
      tryCatch({
        dat <- read_matrix_data(slot)

        # Skip if matrix couldn't be read
        if (is.null(dat)) {
          if (verbose) {
            message("Skipping slot '", slot, "' - not found in file")
          }
          next
        }

        colnames(x = dat) <- Cells(x = file)
        rownames(x = dat) <- if (slot == 'scale.data') {
          FixFeatures(features = safe_read_dataset(assay.group[['scaled.features']], "scaled.features"))
        } else {
          features
        }

        # Set data using V5-compatible approach
        if (inherits(obj, "Assay5")) {
          # V5: Use layer parameter
          obj <- SetAssayData(object = obj, layer = slot, new.data = dat)
        } else {
          # V4: Use slot parameter
          obj <- SetAssayData(object = obj, slot = slot, new.data = dat)
        }

      }, error = function(e) {
        if (verbose) {
          message("Skipping slot '", slot, "' - failed to set: ", conditionMessage(e))
        }
      })
    } else {
      if (verbose) {
        message("Slot '", slot, "' already populated, skipping")
      }
    }
  }
  # Add meta features
  if (safe_exists(assay.group, 'meta.features')) {
    if (verbose) {
      message("Adding feature-level metadata for ", assay)
    }
    meta.data <- as.data.frame(
      x = assay.group[['meta.features']],
      row.names = features
    )
    if (ncol(x = meta.data)) {
      obj <- AddMetaData(
        object = obj,
        metadata = meta.data
      )
    }
  }
  # Add variable feature information
  if (safe_exists(assay.group, 'variable.features')) {
    if (verbose) {
      message("Adding variable feature information for ", assay)
    }
    VariableFeatures(object = obj) <- safe_read_dataset(assay.group[['variable.features']], "variable.features")
  }
  # Add miscellaneous information
  if (safe_exists(assay.group, 'misc')) {
    if (verbose) {
      message("Adding miscellaneous information for ", assay)
    }
    tryCatch({
      # Try to set misc slot if it exists
      if ('misc' %in% slotNames(obj)) {
        # Use SafeH5GroupToList to handle 3D+ arrays
        slot(object = obj, name = 'misc') <- SafeH5GroupToList(h5obj = assay.group[['misc']], recursive = TRUE)
      } else {
        if (verbose) {
          message("Misc slot not available in this assay type (Assay5), skipping")
        }
      }
    }, error = function(e) {
      if (verbose) {
        message("Failed to set misc information: ", conditionMessage(e))
      }
    })
  }
  # Handle S4 class reconstruction - completely skip for Assay5 objects
  # CRITICAL: Check object type FIRST before any S4 operations
  if (verbose) {
    message("Object class check: ", class(obj))
    message("Is Assay5: ", inherits(obj, "Assay5"))
    message("Has s4class attr: ", assay.group$attr_exists(attr_name = 's4class'))
  }
  
  if (inherits(obj, "Assay5")) {
    if (verbose) {
      message("Detected Assay5 object - skipping all S4 reconstruction to preserve object integrity")
    }
    # For Assay5, we're completely done - do not touch the object at all
  } else if (assay.group$attr_exists(attr_name = 's4class')) {
    # Only do S4 reconstruction for older Assay objects (V4 and below)
    if (verbose) {
      message("Processing legacy Assay object with S4 reconstruction")
    }
    tryCatch({
      classdef <- unlist(x = strsplit(
        x = h5attr(x = assay.group, which = 's4class'),
        split = ':'
      ))
      pkg <- classdef[1]
      cls <- classdef[2]
      
      formal <- methods::getClassDef(Class = cls, package = pkg, inherits = FALSE)
      missing <- setdiff(
        x = slotNames(x = formal),
        y = slotNames(x = methods::getClass(Class = 'Assay'))
      )
      missing <- intersect(x = missing, y = names(x = assay.group))
      missing <- sapply(
        X = missing,
        FUN = function(x) {
          return(as.list(x = assay.group[[x]], recursive = TRUE))
        },
        simplify = FALSE
      )
      obj <- c(SeuratObject::S4ToList(object = obj), missing)
      attr(x = obj, which = 'classDef') <- paste(classdef, collapse = ':')
      obj <- SeuratObject::ListToS4(x = obj)
    }, error = function(e) {
      if (verbose) {
        message("S4 class reconstruction failed for legacy assay: ", conditionMessage(e))
      }
    })
  }
  return(obj)
}

#' @importClassesFrom Seurat JackStrawData
#' @importFrom Seurat CreateDimReducObject Stdev IsGlobal JS<-
#'
#' @rdname AssembleObject
#'
AssembleDimReduc <- function(reduction, file, verbose = TRUE) {
  # Helper function to read datasets that might be 2D (V5 compatibility)
  safe_read_dataset <- function(dataset) {
    if (length(dataset$dims) > 1) {
      # 2D dataset - read first column and ensure it's character
      result <- dataset[, 1]
      return(as.character(result))
    } else {
      # 1D dataset - standard read
      return(as.character(dataset[]))
    }
  }
  
  index <- file$index()
  index.check <- vapply(
    X = setdiff(x = names(x = index), y = c('global', 'no.assay')),
    FUN = function(x) {
      return(reduction %in% names(x = index[[x]]$reductions))
    },
    FUN.VALUE = logical(length = 1L)
  )
  if (!any(index.check)) {
    stop(
      "Cannot find reduction ",
      reduction,
      " in this h5Seurat file",
      call. = FALSE
    )
  } else if (sum(index.check) > 1) {
    stop("Multiple reductions named ", reduction, call. = FALSE)
  }
  assay <- names(x = which(x = index.check))
  reduc.group <- file[['reductions']][[reduction]]
  key <- Key(object = reduc.group)
  # Pull cell embeddings
  if (index[[assay]]$reductions[[reduction]][['cell.embeddings']]) {
    if (verbose) {
      message("Adding cell embeddings for ", reduction)
    }
    embeddings <- as.matrix(x = reduc.group[['cell.embeddings']])
    rownames(x = embeddings) <- Cells(x = file)
    colnames(x = embeddings) <- paste0(key, 1:ncol(x = embeddings))
  } else {
    if (verbose) {
      warning(
        "No cell embeddings for ",
        reduction,
        call. = FALSE,
        immediate. = TRUE
      )
    }
    embeddings <- new(Class = 'matrix')
  }
  # Pull feature loadings
  if (index[[assay]]$reductions[[reduction]][['feature.loadings']]) {
    if (verbose) {
      message("Adding feature loadings for ", reduction)
    }
    loadings <- as.matrix(x = reduc.group[['feature.loadings']])
    rownames(x = loadings) <- safe_read_dataset(reduc.group[['features']])
    colnames(x = loadings) <- paste0(key, 1:ncol(x = loadings))
  } else {
    loadings <- new(Class = 'matrix')
  }
  # Pull projected loadings
  if (index[[assay]]$reductions[[reduction]][['feature.loadings.projected']]) {
    if (verbose) {
      message("Adding projected loadings for ", reduction)
    }
    projected <- as.matrix(x = reduc.group[['feature.loadings.projected']])
    rownames(x = projected) <- safe_read_dataset(reduc.group[['projected.features']])
    colnames(x = projected) <- paste0(key, 1:ncol(x = projected))
  } else {
    projected <- new(Class = 'matrix')
  }
  # Build the object
  obj <- CreateDimReducObject(
    embeddings = embeddings,
    loadings = loadings,
    projected = projected,
    assay = assay,
    stdev = Stdev(object = file, reduction = reduction),
    key = key,
    global = IsGlobal(object = reduc.group)
  )
  # Add misc
  if (reduc.group$exists(name = 'misc')) {
    if (verbose) {
      message("Adding miscellaneous information for ", reduction)
    }
    # Use SafeH5GroupToList to handle 3D+ arrays (e.g., from Squidpy UMAP)
    slot(object = obj, name = 'misc') <- SafeH5GroupToList(h5obj = reduc.group[['misc']], recursive = TRUE)
  }
  # Add jackstraw
  if (index[[assay]]$reductions[[reduction]][['jackstraw']]) {
    if (verbose) {
      message("Loading JackStraw data for ", reduction)
    }
    js <- new(Class = 'JackStrawData')
    for (slot in names(x = reduc.group[['jackstraw']])) {
      JS(object = js, slot = slot) <- as.matrix(x = reduc.group[['jackstraw']][[slot]])
    }
    JS(object = obj) <- js
  }
  return(obj)
}

#' @importFrom Seurat as.sparse as.Graph
#'
#' @rdname AssembleObject
#'
AssembleGraph <- function(graph, file, verbose = TRUE) {
  index <- file$index()
  obj <- as.sparse(x = file[['graphs']][[graph]])
  rownames(x = obj) <- colnames(x = obj) <- Cells(x = file)
  obj <- as.Graph(x = obj)
  if (file[['graphs']][[graph]]$attr_exists(attr_name = 'assay.used')) {
    assay <- h5attr(x = file[['graphs']][[graph]], which = 'assay.used')
    if (graph %in% index[[assay]]$graphs) {
      DefaultAssay(object = obj) <- assay
    }
  }
  return(obj)
}

#' @rdname AssembleObject
#'
AssembleImage <- function(image, file, verbose = TRUE) {
  index <- file$index()
  img_group <- file[['images']][[image]]

  # Get assay information
  assay <- if (img_group$attr_exists(attr_name = 'assay')) {
    h5attr(x = img_group, which = 'assay')
  } else {
    'Spatial'  # Default assay name for spatial data
  }

  # Get s4class to determine object type
  s4class <- if (img_group$attr_exists(attr_name = 's4class')) {
    h5attr(x = img_group, which = 's4class')
  } else {
    NULL
  }

  # Try to construct a proper Seurat spatial image object (VisiumV1/VisiumV2)
  # if we have the necessary components
  if (!is.null(s4class) && s4class %in% c('VisiumV1', 'VisiumV2')) {
    tryCatch({
      # Get image data
      image_data <- NULL
      if (img_group$exists(name = 'image')) {
        image_data <- img_group[['image']]$read()
        # Ensure proper format (should already be channels x width x height)
        if (length(dim(image_data)) == 3L) {
          # Normalize to 0-1 range if needed
          if (max(image_data) > 1) {
            image_data <- image_data / 255
          }
          storage.mode(image_data) <- 'double'
        }
      }

      # Get scale factors
      scale_factors <- NULL
      if (img_group$exists(name = 'scale.factors')) {
        sf_group <- img_group[['scale.factors']]

        # Read scale factor values
        spot <- if (sf_group$exists('spot')) sf_group[['spot']][] else NA_real_
        fiducial <- if (sf_group$exists('fiducial')) sf_group[['fiducial']][] else NA_real_
        hires <- if (sf_group$exists('hires')) sf_group[['hires']][] else NA_real_
        lowres <- if (sf_group$exists('lowres')) sf_group[['lowres']][] else NA_real_

        # Create scalefactors object
        scale_factors <- tryCatch({
          scalefactors(
            spot = as.numeric(spot),
            fiducial = as.numeric(fiducial),
            hires = as.numeric(hires),
            lowres = as.numeric(lowres)
          )
        }, error = function(e) {
          if (verbose) {
            message("Could not create scalefactors object: ", conditionMessage(e))
          }
          NULL
        })
      }

      # Get spatial coordinates and boundaries - different for VisiumV1 vs VisiumV2
      coordinates <- NULL
      boundaries_list <- list()

      # Try boundaries/centroids (VisiumV2 style)
      if (img_group$exists(name = 'boundaries')) {
        boundaries_group <- img_group[['boundaries']]
        if (boundaries_group$exists(name = 'centroids')) {
          centroids_group <- boundaries_group[['centroids']]

          # Read all centroids data
          if (centroids_group$exists(name = 'coords')) {
            if (verbose) {
              message("Reading spatial coordinates from boundaries/centroids")
            }

            coords_mat <- centroids_group[['coords']][,]
            cell_names <- if (centroids_group$exists('cells')) {
              as.character(centroids_group[['cells']][])
            } else {
              Cells(x = file)
            }

            # Read centroid parameters
            radius_val <- if (centroids_group$exists('radius')) {
              centroids_group[['radius']][]
            } else {
              as.numeric(scale_factors[['spot']])
            }

            theta_val <- if (centroids_group$exists('theta')) {
              centroids_group[['theta']][]
            } else {
              0
            }

            nsides_val <- if (centroids_group$exists('nsides')) {
              centroids_group[['nsides']][]
            } else {
              8L
            }

            # Create Centroids object using SeuratObject function
            tryCatch({
              # Ensure coords has proper column names and rownames
              if (is.null(colnames(coords_mat))) {
                colnames(coords_mat) <- c('x', 'y')
              }

              # Set rownames to cell names - CreateCentroids uses rownames for cells
              if (length(cell_names) == nrow(coords_mat)) {
                rownames(coords_mat) <- cell_names
              }

              centroids_obj <- SeuratObject::CreateCentroids(
                coords = coords_mat,
                nsides = as.integer(nsides_val),
                radius = as.numeric(radius_val),
                theta = as.numeric(theta_val)
              )

              boundaries_list[['centroids']] <- centroids_obj

              if (verbose) {
                message("Created Centroids object with ", nrow(coords_mat), " cells")
              }
            }, error = function(e) {
              if (verbose) {
                message("Could not create Centroids object: ", conditionMessage(e))
              }
            })
          }
        }
      }

      # For VisiumV1, try reductions/spatial for coordinates dataframe
      if (s4class == 'VisiumV1' && file$exists(name = 'reductions/spatial')) {
        if (verbose) {
          message("Reading spatial coordinates from reduction for VisiumV1")
        }
        spatial_reduc <- file[['reductions/spatial']]
        if (spatial_reduc$exists(name = 'cell.embeddings')) {
          coords_mat <- as.matrix(spatial_reduc[['cell.embeddings']])

          # Create coordinate dataframe for VisiumV1
          coordinates <- data.frame(
            imagerow = coords_mat[, 2],  # Y coordinate
            imagecol = coords_mat[, 1],  # X coordinate
            row.names = Cells(x = file),
            stringsAsFactors = FALSE
          )
        }
      }

      # Get key
      image_key <- if (img_group$exists(name = 'key')) {
        img_group[['key']][]
      } else {
        paste0(image, '_')
      }

      # Only attempt to create a full Visium object if we have the required components
      can_create_v2 <- !is.null(image_data) && !is.null(scale_factors) && length(boundaries_list) > 0
      can_create_v1 <- !is.null(image_data) && !is.null(scale_factors) && !is.null(coordinates)

      if (s4class == 'VisiumV2' && can_create_v2) {
        # VisiumV2 object
        if (verbose) {
          n_cells <- if (!is.null(boundaries_list$centroids)) {
            length(slot(boundaries_list$centroids, 'cells'))
          } else {
            0
          }
          message("Creating VisiumV2 object with ", n_cells, " spots")
        }

        tryCatch({
          obj <- new(
            Class = 'VisiumV2',
            image = image_data,
            scale.factors = scale_factors,
            molecules = list(),  # Empty molecules
            boundaries = boundaries_list,
            assay = assay,
            key = image_key
          )
          return(obj)
        }, error = function(e) {
          if (verbose) {
            message("Error creating VisiumV2 object: ", conditionMessage(e))
          }
        })
      } else if (s4class == 'VisiumV1' && can_create_v1) {
        # VisiumV1 - simpler structure with coordinates dataframe
        if (verbose) {
          message("Creating VisiumV1 object with ", nrow(coordinates), " spots")
        }

        tryCatch({
          obj <- new(
            Class = 'VisiumV1',
            image = image_data,
            scale.factors = scale_factors,
            coordinates = coordinates,
            spot.radius = as.numeric(scale_factors[['spot']]) / 2,
            assay = assay,
            key = image_key
          )
          return(obj)
        }, error = function(e) {
          if (verbose) {
            message("Error creating VisiumV1 object: ", conditionMessage(e))
          }
        })
      }

      # If we get here, we couldn't create a full object
      if (verbose) {
        message("Could not create full ", s4class, " object, falling back to list")
      }
    }, error = function(e) {
      if (verbose) {
        message("Error creating Visium object: ", conditionMessage(e))
        message("Falling back to simple list structure")
      }
    })
  }

  # Fallback: Use simple list structure (original behavior)
  obj <- tryCatch({
    as.list(x = img_group, recursive = TRUE, row.names = Cells(x = file))
  }, error = function(e) {
    if (verbose) {
      message("Standard as.list failed for image, using safe conversion: ", conditionMessage(e))
    }
    SafeH5GroupToList(h5obj = img_group, recursive = TRUE)
  })

  # Try to set assay attribute
  tryCatch({
    # Check if image is in the index - wrap in tryCatch to handle structured objects
    if (!is.null(index[[assay]]$images) && image %in% index[[assay]]$images) {
      tryCatch({
        DefaultAssay(object = obj) <- assay
      }, error = function(e) {
        if (verbose) {
          message("Note: Setting assay as attribute instead of DefaultAssay")
        }
        attr(obj, 'assay') <- assay
      })
    } else {
      # If not in index or index check failed, just set assay attribute
      if (verbose) {
        message("Setting assay as attribute for image")
      }
      attr(obj, 'assay') <- assay
    }
  }, error = function(e) {
    # If the index check itself fails, just set the assay attribute
    if (verbose) {
      message("Could not check image index, setting assay as attribute: ", conditionMessage(e))
    }
    attr(obj, 'assay') <- assay
  })

  return(obj)
}

#' @importClassesFrom Seurat Neighbor
#'
#' @rdname AssembleObject
#'
AssembleNeighbor <- function(neighbor, file, verbose = TRUE) {
  neighbor.group <- file[['neighbors']][[neighbor]]
  obj <- new(
    Class = 'Neighbor',
    nn.idx =  as.matrix(x = neighbor.group[["nn.idx"]]),
    nn.dist = as.matrix(x = neighbor.group[["nn.dist"]]),
    cell.names =  as.matrix(x = neighbor.group[["cell.names"]])[,1]
  )
  return(obj)
}

#' @importClassesFrom Seurat SeuratCommand
#'
#' @rdname AssembleObject
#'
AssembleSeuratCommand <- function(cmd, file, verbose = TRUE) {
  index <- file$index()
  index.check <- vapply(
    X = setdiff(x = names(x = index), y = 'global'),
    FUN = function(x) {
      return(cmd %in% index[[x]]$commands)
    },
    FUN.VALUE = logical(length = 1L)
  )
  if (!any(index.check)) {
    stop("Cannot find command ", cmd, " in this h5Seurat file", call. = FALSE)
  } else if (sum(index.check) > 1) {
    stop("Multiple commands named ", cmd, call. = FALSE)
  }
  cmd.group <- file[['commands']][[cmd]]
  cmdlog <- new(
    Class = 'SeuratCommand',
    name = h5attr(x = cmd.group, which = 'name'),
    time.stamp = as.POSIXct(x = h5attr(x = cmd.group, which = 'time.stamp')),
    call.string = h5attr(x = cmd.group, 'call.string'),
    params = as.list(x = cmd.group, recursive = TRUE)
  )
  assay <- names(x = which(x = index.check))
  if (assay != 'no.assay') {
    slot(object = cmdlog, name = 'assay.used') <- assay
  }
  return(cmdlog)
}
