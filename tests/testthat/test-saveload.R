# Test save/load functionality

test_that("Basic save and load works", {
  skip_if_not_installed("Seurat")

  counts <- matrix(1:12, nrow = 3)
  rownames(counts) <- paste0("Gene", 1:3)
  colnames(counts) <- paste0("Cell", 1:4)
  obj <- Seurat::CreateSeuratObject(counts = counts)

  tmp <- tempfile(fileext = ".h5seurat")
  on.exit(unlink(tmp), add = TRUE)

  SeuratDisk::SaveH5Seurat(obj, tmp, overwrite = TRUE, verbose = FALSE)
  loaded <- SeuratDisk::LoadH5Seurat(tmp, verbose = FALSE)

  expect_s4_class(loaded, "Seurat")
  expect_equal(ncol(loaded), 4)
  expect_equal(nrow(loaded), 3)
})

test_that("Multi-assay objects work correctly", {
  skip_if_not_installed("Seurat")

  rna <- matrix(rpois(100, 1), 10, 10)
  rownames(rna) <- paste0("Gene", seq_len(nrow(rna)))
  colnames(rna) <- paste0("Cell", seq_len(ncol(rna)))
  adt <- matrix(rpois(50, 5), 5, 10)
  rownames(adt) <- paste0("ADT", seq_len(nrow(adt)))
  colnames(adt) <- colnames(rna)

  obj <- Seurat::CreateSeuratObject(counts = rna)
  obj[["ADT"]] <- Seurat::CreateAssayObject(counts = adt)

  tmp <- tempfile(fileext = ".h5seurat")
  on.exit(unlink(tmp), add = TRUE)

  SeuratDisk::SaveH5Seurat(obj, tmp, overwrite = TRUE, verbose = FALSE)
  loaded <- SeuratDisk::LoadH5Seurat(tmp, verbose = FALSE)

  expect_equal(SeuratObject::Assays(loaded), c("RNA", "ADT"))
  expect_equal(nrow(loaded[["RNA"]]), 10)
  expect_equal(nrow(loaded[["ADT"]]), 5)
})

test_that("V5 objects are handled correctly", {
  skip_if_not_installed("Seurat")
  skip_if_not(packageVersion("Seurat") >= "5.0.0")

  counts <- matrix(rpois(100, 1), 10, 10)
  obj <- Seurat::CreateSeuratObject(counts = counts)

  is_v5 <- inherits(obj[["RNA"]], "Assay5")

  if (is_v5) {
    tmp <- tempfile(fileext = ".h5seurat")
    on.exit(unlink(tmp), add = TRUE)

    SeuratDisk::SaveH5Seurat(obj, tmp, overwrite = TRUE, verbose = FALSE)
    loaded <- SeuratDisk::LoadH5Seurat(tmp, verbose = FALSE)

    expect_s4_class(loaded[["RNA"]], "Assay5")
  } else {
    skip("Seurat not running in V5 mode")
  }
})

test_that("Visium spatial coordinates align with metadata", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("hdf5r")

  counts <- matrix(rpois(15, 10), nrow = 5, ncol = 3)
  rownames(counts) <- paste0("Gene", seq_len(nrow(counts)))
  cell_ids <- paste0("Spot", seq_len(ncol(counts)))
  colnames(counts) <- cell_ids

  obj <- Seurat::CreateSeuratObject(counts = counts, assay = "Spatial")

  coords <- data.frame(
    imagerow = c(10, 20, 30),
    imagecol = c(5, 15, 25),
    row.names = cell_ids,
    stringsAsFactors = FALSE
  )

  obj@meta.data$spatial_x <- coords$imagecol
  obj@meta.data$spatial_y <- coords$imagerow

  tmp_h5ad <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp_h5ad), add = TRUE)

  h5 <- hdf5r::H5File$new(tmp_h5ad, mode = "w")
  on.exit(h5$close_all(), add = TRUE)

  obsm <- h5$create_group("obsm")
  spatial_mat <- as.matrix(coords[, c("imagerow", "imagecol")])
  obsm$create_dataset(
    name = "spatial",
    robj = spatial_mat,
    dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE
  )

  uns <- h5$create_group("uns")
  spatial_group <- uns$create_group("spatial")
  lib_group <- spatial_group$create_group("library1")
  sf_group <- lib_group$create_group("scalefactors")
  sf_group$create_dataset(
    name = "tissue_hires_scalef",
    robj = 0.1,
    dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE
  )

  converted <- SeuratDisk:::ConvertH5ADSpatialToSeurat(
    h5ad_file = h5,
    seurat_obj = obj,
    assay_name = "Spatial",
    verbose = FALSE
  )

  expect_equal(converted@meta.data$spatial_x, coords$imagecol)
  expect_equal(converted@meta.data$spatial_y, coords$imagerow)
  expect_identical(converted@misc$spatial_technology, "Visium")
})
