# SeuratDisk V5

[![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)]()
[![R](https://img.shields.io/badge/R-%3E%3D4.0-blue.svg)]()
[![Seurat](https://img.shields.io/badge/Seurat-v5-green.svg)]()

## Overview

SeuratDisk provides interfaces for HDF5-based single cell file formats, enabling efficient storage and conversion between Seurat objects and AnnData/h5ad formats. This fork extends the original [SeuratDisk package](https://github.com/mojaveazure/seurat-disk) with full Seurat V5 compatibility and improved interoperability with Python's scanpy/anndata ecosystem.

## Key Features

- **Seurat V5 Support**: Native handling of Seurat V5 Assay5 objects and layered data structures
- **Format Conversion**: Bidirectional conversion between Seurat, h5Seurat, h5ad, and h5mu formats
- **Multimodal Data Support**: Full h5mu (MuData) format support for CITE-seq, multiome, and other multimodal datasets
- **Direct h5ad Loading**: Load AnnData h5ad files directly into Seurat without intermediate conversion
- **Spatial Visium Support**: Automatically rebuilds Visium spatial images, scalefactors, and coordinates when reading h5ad files
- **Metadata Preservation**: Complete transfer of cell/gene metadata, dimensional reductions, and graphs
- **Modern AnnData Compatibility**: Supports anndata 0.8+ categorical encoding and structure
- **Python Interoperability**: Seamless data exchange with Python's scanpy, muon, and mudata ecosystems

## Installation

```r
# Install dependencies
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Install from GitHub
remotes::install_github("mianaz/seuratdisk-V5")

# Or install locally from source
R CMD INSTALL .
```

## Usage

### Basic Operations

```r
library(SeuratDisk)

# Save Seurat object to h5Seurat format
SaveH5Seurat(seurat_obj, filename = "data.h5Seurat")

# Load h5Seurat file
seurat_obj <- LoadH5Seurat("data.h5Seurat")

# Convert to h5ad for Python/scanpy
Convert("data.h5Seurat", dest = "data.h5ad")
```

### Working with AnnData Files

```r
# Load h5ad directly into Seurat (recommended for V5)
seurat_obj <- LoadH5AD("scanpy_output.h5ad")

# Convert h5ad to h5Seurat
Convert("data.h5ad", dest = "data.h5Seurat")
```

### Working with Multimodal Data (H5MU)

```r
# Load multimodal h5mu file (CITE-seq, multiome, etc.)
seurat_obj <- LoadH5MU("multimodal_data.h5mu")

# Load specific modalities only
seurat_obj <- LoadH5MU("data.h5mu", modalities = c("rna", "prot"))

# Save multimodal Seurat object to h5mu
SaveH5MU(seurat_obj, "output.h5mu")

# Convert between formats
Convert("data.h5mu", dest = "data.h5Seurat")  # h5mu → h5Seurat
Convert("data.h5Seurat", dest = "data.h5mu")  # h5Seurat → h5mu
Convert("data.h5mu", dest = "rna_only.h5ad")  # Extract single modality
```

### Spatial Transcriptomics

```r
# Load a Visium h5ad created with scanpy
brain <- LoadH5AD("visium_brain.h5ad", assay.name = "Spatial")

# Inspect reconstructed spatial images
Images(brain)

# Plot a gene straight away
SpatialFeaturePlot(brain, features = "PTPRC", images = Images(brain)[1])
```

When present, low- and high-resolution Visium images, scalefactors, and spot coordinates are restored directly from the AnnData file, enabling immediate use of standard Seurat spatial workflows.

### Layer Handling

Seurat V5 supports multiple data layers per assay (counts, data, scale.data). These are mapped as follows:

**In h5Seurat files:**
- Stored under `/assays/{assay_name}/layers/{layer_name}`
- Each layer is a separate HDF5 group/dataset

**In h5ad files:**
- Primary layer (typically `data` or `scale.data`) stored as `/X`
- Raw counts stored as `/raw/X` (if available)
- Additional assay layers stored in `/layers/{assay_name}`

**Conversion behavior:**
- h5Seurat → h5ad: `data` layer becomes `/X`, `counts` becomes `/raw/X`
- h5ad → Seurat: `/X` becomes `data` layer, `/raw/X` becomes `counts`
- Missing layers are handled gracefully (e.g., counts-only objects supported)

## System Requirements

- R >= 4.0.0
- Seurat >= 5.0.0
- SeuratObject >= 5.0.0
- hdf5r >= 1.3.0

## Improvements in This Fork

### V5 Compatibility
- Native Assay5 object handling with layer-based data storage
- Patched hdf5r methods for improved multidimensional array access
- Support for counts-only objects (data layer optional)

### AnnData Format Updates
- Modern categorical encoding (`categories`/`codes` instead of `__categories`)
- Proper encoding attributes for anndata 0.8+ compatibility
- Complete metadata group structure (`obsm`, `obsp`, `varm`, `varp`, `layers`, `uns`)
- Correct `shape` attributes for sparse matrices
- Visium spatial libraries are reconstructed into fully functional `VisiumV2` images inside Seurat

### Bug Fixes
- Fixed V5 metadata transfer issues
- Corrected dimension handling for 1D/2D feature datasets
- Improved error handling for missing data layers

## Known Limitations

- Round-trip conversions (h5Seurat ↔ Seurat ↔ h5Seurat) may lose some V5-specific metadata
- H5MU multimodal support requires MuDataSeurat package (automatically installed from GitHub)
- Spatial data in h5mu format is stored per-modality; full Visium reconstruction is in development
- Visium hires imagery is embedded inside the Seurat object; raw TIFF/PNG assets are not exported automatically
- Large datasets (>100K cells) may require substantial memory

## Troubleshooting

### HDF5 Method Patching

If you see warnings about H5D method patching:
```
Error patching H5D methods: no method found for function '[' and signature H5D
```
This is expected if hdf5r methods are not yet loaded. The package uses fallback methods automatically.

### Memory Management

For large datasets:
```r
# Monitor memory usage
gc()

# Windows only: increase memory limit
memory.limit(size = 32000)
```

## Upstream Reference

This fork is based on the original SeuratDisk package by [mojaveazure](https://github.com/mojaveazure/seurat-disk). The original repository may no longer be actively maintained. This fork focuses specifically on Seurat V5 compatibility and modern AnnData format support.

## Citation

If you use SeuratDisk in published research, please cite the relevant Seurat papers:

```
Hao and Hao et al. (2024) Dictionary learning for integrative, multimodal and
scalable single-cell analysis. Nature Biotechnology. [Seurat V5]

Hao et al. (2021) Integrated analysis of multimodal single-cell data.
Cell, 184(13), 3573-3587. [Seurat V4]
```

## Contributing

Contributions are welcome. Key areas for improvement:
- Enhanced multi-modal assay support
- Spatial transcriptomics integration
- Performance optimization for very large datasets
- Additional file format conversions

## License

GPL-3 (inherited from original SeuratDisk)

## Support

For issues related to V5 compatibility or this fork, please open an issue in this repository. For general Seurat questions, visit the [Seurat discussions](https://github.com/satijalab/seurat/discussions).
