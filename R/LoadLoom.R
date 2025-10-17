#' @include zzz.R
#' @include loom.R
#' @include loom_bindings.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Loom-file Loading
#'
#' Load data from a loom file into a \code{\link[Seurat]{Seurat}} object
#'
#' @param file,x Name of loom file or a \code{loom} object to load data from
#' @param assay Name of assay to store expression data as; if \code{NULL}, will
#' search for an HDF5 attribute named \code{SEURAT_ASSAY} or an attribute
#' dataset named \code{/attrs/SEURAT_ASSAY} for assay name. If not found,
#' defaults to \dQuote{RNA}
#' @param cells Name of dataset in \code{/col_attrs} with cell names
#' @param features Name of dataset in \code{/row_attrs} with feature names
#' @param normalized Name of matrix in \code{/layers} to store normalized data
#' as; pass \dQuote{/matrix} to store \code{/matrix} as normalized data instead
#' of raw counts
#' @param scaled Name of dataset in \code{/layers} to store scaled data as
#' @param filter Keep only selected cells and/or features as specified by
#' \code{/col_attrs/Valid} and \code{/row_attrs/Valid}, respectively
#' @param verbose Show progress updates
#' @param ... Arguments passed to other methods
#'
#' @return A \code{\link[Seurat]{Seurat}} object
#'
#' @name LoadLoom
#' @rdname LoadLoom
#'
#' @inherit LoomLoading details
#'
#' @inheritSection LoomLoading Loom 0.1 Loading
#'
#' @inheritSection LoomLoading Loom 3.0.0 Loading
#'
#' @seealso
#' \href{http://linnarssonlab.org/loompy/conventions/index.html}{Loom
#' file conventions}
#'
#' @export
#'
LoadLoom <- function(
  file,
  assay = NULL,
  cells = 'CellID',
  features = 'Gene',
  normalized = NULL,
  scaled = NULL,
  filter = c('cells', 'features', 'all', 'none'),
  verbose = TRUE,
  ...
) {
  UseMethod(generic = 'LoadLoom', object = file)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname LoadLoom
#' @method LoadLoom character
#' @export
#'
LoadLoom.character <- function(file, ...) {
  loom <- loom$new(filename = file, mode = 'r')
  on.exit(expr = loom$close_all())
  return(LoadLoom(file = loom, ...))
}

#' @rdname LoadLoom
#' @method LoadLoom H5File
#' @export
#'
LoadLoom.H5File <- function(file, ...) {
  return(LoadLoom(
    file = as.loom(x = file),
    ...
  ))
}

#' @rdname LoadLoom
#' @method LoadLoom loom
#' @export
#'
LoadLoom.loom <- function(file, ...) {
  return(as.Seurat(x = file, ...))
}

#' @importFrom stringi stri_count_fixed
#'
#' @aliases as.Seurat
#'
#' @rdname LoadLoom
#' @method as.Seurat loom
#' @export
#'
as.Seurat.loom <- function(
  x,
  assay = NULL,
  cells = 'CellID',
  features = 'Gene',
  normalized = NULL,
  scaled = NULL,
  filter = c('cells', 'features', 'all', 'none'),
  verbose = TRUE,
  ...
) {
  cells <- H5Path('col_attrs', cells)
  features <- H5Path('row_attrs', features)
  if (!Exists(x = x, name = features) || !inherits(x = x[[features]], what = 'H5D')) {
    stop("Cannot find feature names dataset at ", features, call. = FALSE)
  } else if (!Exists(x = x, name = cells) || !inherits(x = x[[cells]], what = 'H5D')) {
    stop("Cannot find cell names dataset at ", cells, call. = FALSE)
  } else if (length(x = Dims(x = x[[cells]])) != 1) {
    stop("The cell names dataset must be one-dimensional", call. = FALSE)
  } else if (length(x = Dims(x = x[[features]])) != 1) {
    stop("The feature names dataset must be one-dimensional", call. = FALSE)
  }
  if (!is.null(x = normalized) && !grepl(pattern = '^[/]?matrix$', x = normalized)) {
    if (!grepl(pattern = '^[/]?layers', x = normalized)) {
      normalized <- H5Path('layers', normalized)
    }
    if (stri_count_fixed(str = normalized, pattern = '/') != 1) {
      stop(
        "'normalized' must be a dataset within the 'layers' group",
        call. = FALSE
      )
    }
    if (!Exists(x = x, name = normalized)) {
      warning(
        "Cannot find ",
        basename(path = normalized),
        " in this loom file, not loading normalized data",
        call. = FALSE,
        immediate. = TRUE
      )
      normalized <- NULL
    }
  }
  if (!is.null(x = scaled)) {
    if (!grepl(pattern = '^[/]?layers', x = scaled)) {
      scaled <- H5Path('layers', scaled)
    }
    if (stri_count_fixed(str = scaled, pattern = '/') != 1) {
      stop(
        "'scaled' must be a dataset within the 'layers' group",
        call. = FALSE
      )
    }
    if (!Exists(x = x, name = scaled)) {
      warning(
        "Cannot find ",
        basename(path = scaled),
        " in this loom file, not loading scaled data",
        call. = FALSE,
        immediate. = TRUE
      )
      scaled <- NULL
    }
  }
  version <- ClosestVersion(query = x$version(), targets = c('0.1.0', '3.0.0'))
  load.fxn <- switch(
    EXPR = version,
    # '0.1.0' = LoadLoom0.1,
    # '3.0.0' = LoadLoom3.0
    LoadLoom3.0
  )
  object <- load.fxn(
      file = x,
      assay = assay,
      cells = cells,
      features = features,
      normalized = normalized,
      scaled = scaled,
      filter = filter,
      verbose = verbose
    )
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Loom-file Loading
#'
#' Version-specific loom-file loading functions
#'
#' @inheritParams LoadLoom
#'
#' @inherit LoadLoom return
#'
#' @name LoomLoading
#' @rdname LoomLoading
#'
#' @importFrom Seurat CreateAssayObject Key<- CreateSeuratObject SetAssayData as.Graph DefaultAssay<- CreateDimReducObject Loadings Loadings<- Stdev
#'
#' @details
#' \code{LoadLoom} will try to automatically fill slots of a \code{Seurat}
#' object based on data presence or absence in a given loom file. This method
#' varies by loom specification version. For version-specific details, see
#' sections below
#'
#' @section Loom 0.1 Loading:
#' Loading data from loom files less than version 3.0.0 is not
#' currently supported
#'
#' @keywords internal
#'
LoadLoom0.1 <- function(
  file,
  assay = NULL,
  cells = 'col_atts/CellID',
  features = 'row_attrs/Gene',
  normalized = NULL,
  scaled = NULL,
  filter = c('cells', 'features', 'all', 'none'),
  verbose = TRUE
) {
  filter <- filter[1]
  filter <- match.arg(arg = filter)
  assay <- assay %||% suppressWarnings(expr = DefaultAssay(object = file)) %||% 'RNA'
  if (isTRUE(x = verbose)) {
    message("Reading in /matrix")
  }
  counts <- t(x = as.matrix(x = file[['matrix']]))
  dnames <- list(
    features = file[[features]][],
    cells = file[[cells]][]
  )
  if (anyDuplicated(x = dnames$features)) {
    warning(
      "Duplicate feature names found, making unique",
      call. = FALSE,
      immediate. = TRUE
    )
    dnames$features <- make.unique(names = dnames$features)
  }
  dimnames(x = counts) <- dnames
}

#' @name LoomLoading
#' @rdname LoomLoading
#'
#' @section Loom 3.0.0 Loading:
#' For loom files version 3.0.0 and higher, \code{LoadLoom3.0} provides
#' comprehensive loading with support for filtering, dimensional reductions,
#' and cell graphs.
#'
LoadLoom3.0 <- function(
  file,
  assay = NULL,
  cells = 'col_attrs/CellID',
  features = 'row_attrs/Gene',
  normalized = NULL,
  scaled = NULL,
  filter = c('cells', 'features', 'all', 'none'),
  verbose = TRUE
) {
  # Implement filtering
  filter <- filter[1]
  filter <- match.arg(arg = filter)
  assay <- assay %||% suppressWarnings(expr = DefaultAssay(object = file)) %||% 'RNA'

  # Check for Valid filters
  cells.valid <- NULL
  features.valid <- NULL
  if (filter != 'none') {
    if (filter %in% c('cells', 'all')) {
      valid.path <- H5Path('col_attrs', 'Valid')
      if (Exists(x = file, name = valid.path)) {
        cells.valid <- as.logical(file[[valid.path]][])
        if (isTRUE(x = verbose)) {
          message("Found cell filter, keeping ", sum(cells.valid), " of ", length(cells.valid), " cells")
        }
      }
    }
    if (filter %in% c('features', 'all')) {
      valid.path <- H5Path('row_attrs', 'Valid')
      if (Exists(x = file, name = valid.path)) {
        features.valid <- as.logical(file[[valid.path]][])
        if (isTRUE(x = verbose)) {
          message("Found feature filter, keeping ", sum(features.valid), " of ", length(features.valid), " features")
        }
      }
    }
  }

  # Read in /matrix
  if (isTRUE(x = verbose)) {
    message("Reading in /matrix")
  }
  counts <- t(x = as.matrix(x = file[['matrix']]))
  dnames <- list(
    features = file[[features]][],
    cells = file[[cells]][]
  )

  # Apply filtering
  if (!is.null(x = features.valid)) {
    counts <- counts[features.valid, , drop = FALSE]
    dnames$features <- dnames$features[features.valid]
  }
  if (!is.null(x = cells.valid)) {
    counts <- counts[, cells.valid, drop = FALSE]
    dnames$cells <- dnames$cells[cells.valid]
  }

  if (anyDuplicated(x = dnames$features)) {
    warning(
      "Duplicate feature names found, making unique",
      call. = FALSE,
      immediate. = TRUE
    )
    dnames$features <- make.unique(names = dnames$features)
  }
  if (anyDuplicated(x = dnames$cells)) {
    warning(
      "Duplicate cell names found, making unique",
      call. = FALSE,
      immediate. = TRUE
    )
    dnames$cells <- make.unique(names = dnames$cells)
  }
  dimnames(x = counts) <- dnames
  # Assemble the initial assay
  if (!is.null(x = normalized) && grepl(pattern = '^[/]?matrix$', x = normalized)) {
    if (isTRUE(x = verbose)) {
      message("Storing /matrix as data")
    }
    assay.obj <- CreateAssayObject(data = counts)
  } else {
    if (isTRUE(x = verbose)) {
      message("Storing /matrix as counts")
    }
    assay.obj <- CreateAssayObject(counts = counts)
  }
  Key(object = assay.obj) <- UpdateKey(key = tolower(x = assay))
  if (isTRUE(x = verbose)) {
    message("Saving /matrix to assay '", assay, "'")
  }
  object <- CreateSeuratObject(counts = assay.obj, assay = assay)
  # Load in normalized data
  if (!is.null(x = normalized) && !grepl(pattern = '^[/]?matrix$', x = normalized)) {
    if (verbose) {
      message("Saving ", basename(path = normalized), " as data")
    }
    norm.data <- t(x = as.matrix(x = file[[normalized]]))
    # Apply filtering to normalized data
    if (!is.null(x = features.valid)) {
      norm.data <- norm.data[features.valid, , drop = FALSE]
    }
    if (!is.null(x = cells.valid)) {
      norm.data <- norm.data[, cells.valid, drop = FALSE]
    }
    dimnames(x = norm.data) <- dnames
    object <- SetAssayData(
      object = object,
      assay = assay,
      slot = 'data',
      new.data = norm.data
    )
  }
  # Load in scaled data
  if (!is.null(x = scaled)) {
    if (verbose) {
      message("Saving ", basename(path = scaled), " as scaled data")
    }
    scaled.data <- t(x = as.matrix(x = file[[scaled]]))
    # Apply filtering to scaled data
    if (!is.null(x = features.valid)) {
      scaled.data <- scaled.data[features.valid, , drop = FALSE]
    }
    if (!is.null(x = cells.valid)) {
      scaled.data <- scaled.data[, cells.valid, drop = FALSE]
    }
    dimnames(x = scaled.data) <- dnames
    object <- SetAssayData(
      object = object,
      assay = assay,
      slot = 'scale.data',
      new.data = scaled.data
    )
  }
  # Load in feature-level metadata
  mf.remove <- c(basename(path = features), 'Valid')
  meta.features <- as.data.frame(x = file[['row_attrs']])
  meta.features <- meta.features[, !colnames(x = meta.features) %in% mf.remove, drop = FALSE]
  # Apply filtering to feature metadata
  if (!is.null(x = features.valid)) {
    meta.features <- meta.features[features.valid, , drop = FALSE]
  }
  rownames(x = meta.features) <- dnames$features
  if (ncol(x = meta.features)) {
    object[[assay]] <- AddMetaData(
      object = object[[assay]],
      metadata = meta.features
    )
  }
  # Load in cell-level metadata
  md.remove <- c(basename(path = cells), 'Valid')
  meta.data <- as.data.frame(x = file[['col_attrs']])
  meta.data <- meta.data[, !colnames(x = meta.data) %in% md.remove, drop = FALSE]
  # Apply filtering to cell metadata
  if (!is.null(x = cells.valid)) {
    meta.data <- meta.data[cells.valid, , drop = FALSE]
  }
  rownames(x = meta.data) <- dnames$cells
  if (ncol(x = meta.data)) {
    object <- AddMetaData(object = object, metadata = meta.data)
  }
  # Load in dimensional reductions
  # Loom files may store dimensional reductions in /col_attrs as 2D datasets
  # We'll look for a '/reductions' group if it exists (non-standard extension)
  if (Exists(x = file, name = 'reductions')) {
    for (reduc.name in names(x = file[['reductions']])) {
      tryCatch({
        if (verbose) {
          message("Loading dimensional reduction ", reduc.name)
        }
        reduc.group <- file[['reductions']][[reduc.name]]
        # Read embeddings
        if (Exists(x = reduc.group, name = 'embeddings')) {
          embeddings <- t(x = as.matrix(x = reduc.group[['embeddings']]))
          # Apply filtering
          if (!is.null(x = cells.valid)) {
            embeddings <- embeddings[cells.valid, , drop = FALSE]
          }
          rownames(x = embeddings) <- colnames(x = object)
          # Create DimReduc object
          reduc.obj <- CreateDimReducObject(
            embeddings = embeddings,
            key = paste0(reduc.name, '_'),
            assay = assay
          )
          # Load additional data if available
          if (Exists(x = reduc.group, name = 'loadings')) {
            loadings <- as.matrix(x = reduc.group[['loadings']])
            if (!is.null(x = features.valid)) {
              # Only filter if loadings match feature count
              if (nrow(loadings) == length(features.valid)) {
                loadings <- loadings[features.valid, , drop = FALSE]
              }
            }
            reduc.obj <- Loadings(object = reduc.obj) <- loadings
          }
          if (Exists(x = reduc.group, name = 'stdev')) {
            Stdev(object = reduc.obj) <- reduc.group[['stdev']][]
          }
          object[[reduc.name]] <- reduc.obj
        }
      }, error = function(e) {
        warning(
          "Failed to load dimensional reduction ", reduc.name, ": ", e$message,
          call. = FALSE,
          immediate. = TRUE
        )
      })
    }
  }
  # Load in cell graphs
  for (i in names(x = file[['col_graphs']])) {
    if (verbose) {
      message("Loading graph ", i)
    }
    graph <- LoadGraph(graph = file[['col_graphs']][[i]])
    # Apply filtering to graphs
    if (!is.null(x = cells.valid)) {
      graph <- graph[cells.valid, cells.valid, drop = FALSE]
    }
    colnames(x = graph) <- rownames(x = graph) <- colnames(x = object)
    graph <- as.Graph(x = graph)
    DefaultAssay(object = graph) <- assay
    object[[i]] <- graph
  }
  return(object)
}

#' @importFrom Matrix sparseMatrix
#'
#' @keywords internal
#'
LoadGraph <- function(graph) {
  if (!inherits(x = graph, what = 'H5Group')) {
    stop(graph$get_obj_name(), " is not a loom graph", call. = FALSE)
  }
  OneD <- function(x) {
    check <- Exists(x = graph, name = x)
    if (isTRUE(x = check)) {
      check <- inherits(x = graph[[x]], what = 'H5D')
    }
    if (isTRUE(x = check)) {
      check <- length(x = Dims(x = graph[[x]])) == 1
    }
    return(check)
  }
  for (ds in c('a', 'b', 'w')) {
    if (!isTRUE(x = OneD(x = ds))) {
      stop("Missing one dimensional dataset '", ds, "'", call. = FALSE)
    }
  }
  return(sparseMatrix(
    i = graph[['a']][] + 1,
    j = graph[['b']][] + 1,
    x = graph[['w']][]
  ))
}
