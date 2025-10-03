# Spatial Data Support Implementation Plan for SeuratDisk V5

## Overview
This document outlines the implementation plan for adding comprehensive spatial data support to SeuratDisk, enabling seamless conversion between h5Seurat and h5ad formats for spatial transcriptomics data.

## Current State Assessment

### Existing Multimodal Support
✅ **Already Supported:**
- Multiple assays (RNA, ADT, SCT) in h5Seurat format
- Basic h5ad to h5Seurat conversion for expression matrices
- Layer support in V5 compatibility
- Image data structure in h5Seurat
- obsm, obsp, uns data reading from h5ad

❌ **Missing:**
- Spatial coordinate conversion between formats
- Visium/SlideSeq specific metadata handling
- FOV object support for imaging-based spatial data
- Spatial image data conversion
- MuData integration for multimodal spatial data

## Implementation Phases

### Phase 1: Core Spatial Data Structure Support

#### 1.1 Spatial Coordinates Conversion
**Files to modify:** `R/LoadH5AD.R`, `R/Convert.R`

```r
# h5ad → Seurat:
# Read from adata.obsm['spatial'] → Store in Seurat image objects
# Seurat → h5ad:
# Extract from GetTissueCoordinates() → Write to obsm['spatial']
```

**Implementation tasks:**
- [ ] Add spatial coordinate detection in LoadH5AD
- [ ] Create helper function for coordinate extraction from Seurat images
- [ ] Implement coordinate writing in H5SeuratToH5AD
- [ ] Handle coordinate transformations and scaling factors

#### 1.2 Spatial Image Data
**Files to modify:** `R/LoadH5AD.R`, `R/SaveH5Seurat.R`, `R/Convert.R`

```r
# h5ad structure: uns['spatial'][library_id]['images']
# Seurat structure: @images slot with VisiumV1/V2 or SlideSeq objects
```

**Implementation tasks:**
- [ ] Parse uns['spatial'] structure from h5ad
- [ ] Convert image data between formats
- [ ] Handle scale factors and spot diameter metadata
- [ ] Support multiple image resolutions (hires, lowres)

### Phase 2: Technology-Specific Support

#### 2.1 Visium Data Support
**New file:** `R/SpatialVisium.R`

```r
#' Convert Visium-specific data between formats
#' Handles: spot coordinates, tissue positions, scale factors
#' Supports both Visium and Visium HD
```

**Implementation tasks:**
- [ ] Create VisiumV1/V2 objects from h5ad data
- [ ] Handle tissue position lists
- [ ] Convert spot-barcode mappings
- [ ] Support Visium HD high-resolution data

#### 2.2 SlideSeq Support
**New file:** `R/SpatialSlideSeq.R`

```r
#' Convert SlideSeq data between formats
#' Handles: bead locations, puck IDs
```

**Implementation tasks:**
- [ ] Parse SlideSeq coordinate systems
- [ ] Handle bead array positions
- [ ] Convert between continuous and binned coordinates

#### 2.3 FOV/Imaging Data Support
**New file:** `R/SpatialFOV.R`

```r
#' Support for imaging-based spatial data (Xenium, MERFISH, etc.)
#' Handles: FOV objects, molecular positions, cell boundaries
```

**Implementation tasks:**
- [ ] Create FOV object converter
- [ ] Handle single-molecule positions
- [ ] Convert cell segmentation data
- [ ] Support multiple FOVs per dataset

### Phase 3: Multimodal Spatial Integration

#### 3.1 Multiple Assays with Spatial
**Files to modify:** `R/LoadH5AD.R`, `R/SaveH5Seurat.R`

```r
# Support spatial data across multiple assays
# Example: RNA + Protein (ADT) in Visium
```

**Implementation tasks:**
- [ ] Extend assay iteration to preserve spatial links
- [ ] Handle assay-specific spatial metadata
- [ ] Maintain coordinate consistency across assays

#### 3.2 MuData Support
**New file:** `R/MuDataSupport.R`

```r
#' Read and write MuData h5mu files
#' Handle modality-specific spatial information
```

**Implementation tasks:**
- [ ] Detect h5mu file structure
- [ ] Parse .mod dictionary for multiple modalities
- [ ] Preserve global vs modality-specific spatial data
- [ ] Write h5mu files from Seurat objects

### Phase 4: Enhanced Metadata Handling

#### 4.1 Spatial Metadata Preservation
**Files to modify:** `R/Convert.R`, `R/h5Seurat_bindings.R`

```r
# Comprehensive spatial metadata mapping between formats
```

**Mapping table:**
| h5ad | h5Seurat | Description |
|------|----------|-------------|
| obsm['spatial'] | @images[[x]]@coordinates | Spatial coordinates |
| uns['spatial'][lib]['scalefactors'] | @images[[x]]@scale.factors | Scale factors |
| uns['spatial'][lib]['images'] | @images[[x]]@image | Image data |
| obs['in_tissue'] | @images[[x]]@coordinates$tissue | Tissue detection |

**Implementation tasks:**
- [ ] Create metadata mapping functions
- [ ] Handle missing metadata gracefully
- [ ] Preserve custom spatial annotations
- [ ] Support library-specific metadata

### Phase 5: Testing and Validation

#### 5.1 Test Data Preparation
**Directory:** `tests/testthat/test-spatial/`

**Datasets to include:**
```r
# Visium datasets
- SeuratData::stxBrain (Mouse Brain Visium)
- SeuratData::stxKidney (Mouse Kidney Visium)

# SlideSeq datasets
- SeuratData::ssHippo (Mouse Hippocampus SlideSeq)

# Imaging datasets
- Example Xenium data (if available)
```

#### 5.2 Test Suite Implementation
**Files:** `tests/testthat/test-spatial-*.R`

```r
test-spatial-visium.R      # Visium-specific tests
test-spatial-slideseq.R    # SlideSeq tests
test-spatial-fov.R         # FOV/imaging tests
test-spatial-convert.R     # Round-trip conversion tests
test-spatial-mudata.R      # MuData integration tests
```

**Test coverage:**
- [ ] Coordinate preservation
- [ ] Image data integrity
- [ ] Metadata completeness
- [ ] Round-trip conversion accuracy
- [ ] Multi-assay spatial data
- [ ] Edge cases and error handling

### Phase 6: Documentation and Examples

#### 6.1 Vignettes
**New vignettes:**
- `vignettes/spatial_conversion.Rmd` - Spatial data conversion guide
- `vignettes/mudata_integration.Rmd` - MuData integration examples

#### 6.2 Function Documentation
- Update all modified functions with spatial parameters
- Add @examples with spatial data
- Document spatial data structures

## Technical Considerations

### 1. Backward Compatibility
- Maintain compatibility with existing h5Seurat files
- Gracefully handle missing spatial data
- Version detection for format changes

### 2. Performance Optimization
- Lazy loading for large images
- Chunked reading/writing for coordinates
- Memory-efficient image conversion

### 3. Format Standards Compliance
- Follow h5ad 0.10+ spatial standards
- Adhere to Seurat v5 spatial object specifications
- Support scanpy/squidpy conventions

## Implementation Timeline

**Week 1-2: Phase 1** - Core spatial structure support
**Week 3-4: Phase 2** - Technology-specific implementations
**Week 5: Phase 3** - Multimodal integration
**Week 6: Phase 4** - Metadata handling
**Week 7-8: Phase 5** - Testing and validation
**Week 9: Phase 6** - Documentation

## Dependencies and Requirements

### R Package Dependencies
```r
# Existing dependencies (already in DESCRIPTION)
- Seurat (>= 5.0.0)
- SeuratObject (>= 5.0.0)
- hdf5r (>= 1.3.0)

# Consider adding (optional):
- SeuratData  # For test datasets
- imager      # For image processing
- EBImage     # Alternative image processing
```

### Python Interoperability Testing
```python
# Python packages for validation
- scanpy >= 1.9
- squidpy >= 1.3
- mudata >= 0.2
- anndata >= 0.10
```

## Success Criteria

1. ✅ Full round-trip conversion of Visium data without data loss
2. ✅ Successful loading of spatial h5ad files in Seurat
3. ✅ Preservation of spatial coordinates and images
4. ✅ Support for multiple spatial technologies
5. ✅ Integration with existing multimodal workflows
6. ✅ Comprehensive test coverage (>90%)
7. ✅ Clear documentation and examples

## Future Enhancements

- Support for newer spatial technologies (CosMx, STOmics)
- Integration with SpatialData framework
- Streaming support for large spatial datasets
- Cloud storage compatibility
- GUI for spatial data preview during conversion