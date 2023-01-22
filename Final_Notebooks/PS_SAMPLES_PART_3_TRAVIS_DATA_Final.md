PS SAMPLES - PART 3
================

# PSORIASIS SAMPLES - CONTINUED - PART 3

## PSO + PSA + NORMAL SKIN COMBINED ANALYSIS

### Integration with PS skin scRNA data from (dataset 1) Hughes et al.

### (DOI - <https://doi.org/10.1016/j.immuni.2020.09.015>).

| LINK / NUMBER             | DESCRIPTION                                                                                               |
|---------------------------|-----------------------------------------------------------------------------------------------------------|
| [FIGURE 5A](#figure-5a)   | Spatial plots for cell type enrichment using (PS only) data-set 1 (Hughes et al.) - Structural cell types |
| [FIGURE 5B](#figure-5b)   | MIA for (PS only) dataset 1 - Immune cell types                                                           |
| [FIGURE 5C](#figure-5c)   | MIA for TRM (Bulk RNAseq dataset)                                                                         |
| [FIGURE 5E](#figure-5e)   | B CELL - PATHWAYS                                                                                         |
| [FIGURE S9B](#figure-s9b) | UMAP for (PS only) dataset 1                                                                              |
| [FIGURE S9D](#figure-s9d) | MIA for (PS only) dataset 1 - Structural cell types                                                       |

### LOAD ALL PACKAGES

``` r
library(tidyverse)
```

    ## Warning: package 'tidyverse' was built under R version 4.1.2

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
    ## ✔ ggplot2 3.4.0      ✔ purrr   0.3.5 
    ## ✔ tibble  3.1.8      ✔ dplyr   1.0.10
    ## ✔ tidyr   1.2.1      ✔ stringr 1.4.1 
    ## ✔ readr   2.1.3      ✔ forcats 0.5.2

    ## Warning: package 'ggplot2' was built under R version 4.1.2

    ## Warning: package 'tibble' was built under R version 4.1.2

    ## Warning: package 'tidyr' was built under R version 4.1.2

    ## Warning: package 'readr' was built under R version 4.1.2

    ## Warning: package 'purrr' was built under R version 4.1.2

    ## Warning: package 'dplyr' was built under R version 4.1.2

    ## Warning: package 'stringr' was built under R version 4.1.2

    ## Warning: package 'forcats' was built under R version 4.1.2

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

``` r
library(Seurat)
```

    ## Warning: package 'Seurat' was built under R version 4.1.2

    ## Attaching SeuratObject
    ## Attaching sp

``` r
library(cowplot)
library(ggsci)
library(RColorBrewer)
```

    ## Warning: package 'RColorBrewer' was built under R version 4.1.2

``` r
library(pheatmap)
library(SeuratDisk)
```

    ## Registered S3 method overwritten by 'SeuratDisk':
    ##   method            from  
    ##   as.sparse.H5Group Seurat

### LOAD HELPER FUNCTIONS

``` r
source("../SPATIAL_FUNCTIONS.R")
```

    ## Warning: package 'reticulate' was built under R version 4.1.2

    ## Warning: package 'clusterProfiler' was built under R version 4.1.1

    ## 

    ## Registered S3 method overwritten by 'ggtree':
    ##   method      from 
    ##   identify.gg ggfun

    ## clusterProfiler v4.0.5  For help: https://yulab-smu.top/biomedical-knowledge-mining-book/
    ## 
    ## If you use clusterProfiler in published research, please cite:
    ## T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan, X Fu, S Liu, X Bo, and G Yu. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. The Innovation. 2021, 2(3):100141. doi: 10.1016/j.xinn.2021.100141

    ## 
    ## Attaching package: 'clusterProfiler'

    ## The following object is masked from 'package:purrr':
    ## 
    ##     simplify

    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

    ## Loading required package: AnnotationDbi

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## Warning: package 'BiocGenerics' was built under R version 4.1.1

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min

    ## Loading required package: Biobase

    ## Warning: package 'Biobase' was built under R version 4.1.1

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Loading required package: IRanges

    ## Warning: package 'IRanges' was built under R version 4.1.1

    ## Loading required package: S4Vectors

    ## Warning: package 'S4Vectors' was built under R version 4.1.3

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:clusterProfiler':
    ## 
    ##     rename

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     expand

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:clusterProfiler':
    ## 
    ##     slice

    ## The following object is masked from 'package:sp':
    ## 
    ##     %over%

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice

    ## The following object is masked from 'package:purrr':
    ## 
    ##     reduce

    ## 
    ## Attaching package: 'AnnotationDbi'

    ## The following object is masked from 'package:clusterProfiler':
    ## 
    ##     select

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

    ## 

### LOAD COLOR PALETTE USED FOR ALL PLOTS

``` r
col.pal <- RColorBrewer::brewer.pal(9, "OrRd")
```

### Integration with dataset 1 (Hughes et al)

Importing Harmonized (Batch corrected) ST data (produced in PS Samples
part 1 notebook).

``` r
# 1. All Spatial Samples with Harmony Batch Correction
skin_data.hm.sct <- readRDS(file = "/Volumes/Extreme Pro/GITHUB-DATA/ST-DATA/PSORIASIS-DATA/RDS-Files/ALL_SPATIAL_SAMPLES(HM_BATCH_CORRECTED).RDS")

# 2. Single Cell Markers
skin_data.hm.sct.markers <- readRDS("/Volumes/Extreme Pro/GITHUB-DATA/ST-DATA/PSORIASIS-DATA/RDS-Files/ALL_ST_HARMONY_ALIGNED_MARKERS.RDS")
```

``` r
## COLOR FOR LABELS
color.labels <- c("0 Fibroblasts"="#87CEFA",
"1 Macs + fibroblasts"="#4876FF",
"2 Eccrine + melanocyte precursors"="#CD853F",
"3 Epidermis"="#BF96FF",
"4 Epidermis"="#FF0000",
"5 Connective tissue"="#CAF178",
"6 Mixed"="#E0BFB6",
"7 Epidermis"="#68228B",
"8 Hair follicle and sebaceous glands"="#7B0000",
"9 Adipose"="#FFC71A",
"10 Suprabasal keratinocytes"="#C355A0",
"11 Smooth muscle"="#00B923",
"12 Endothelial cells"="#8B5A2B",
"13 Immunoglobulins, fibroblasts"="#838B8B",
"14 Smooth muscle"="#005947",
"15 Mixed"="#C1CDCD",
"16 Adipose, fibroblasts"="#FF7545")
```

Load data-set 1 (Hughes et al.) scRNA data

``` r
# Load Travis Data Seurat Object
travis.sc.seu.obj <- readRDS("/Volumes/Extreme Pro/GITHUB-DATA/SC-RNA-DATA/TRAVIS-DATA/RDS-Files/Travis_Data_Seurat_Object.RDS")
```

``` r
# META DATA FROM TRAVIS Paper
meta.data <- read.csv("/Volumes/Extreme Pro/GITHUB-DATA/SC-RNA-DATA/TRAVIS-DATA/Meta-Data/Cell_Level_Metadata.csv")
```

Subset the data to only include the Psoriasis group samples

``` r
# Gather Travis Meta Data & Subset Psoriasis Data
skin.meta.data <- as.data.frame(travis.sc.seu.obj@meta.data) %>% rownames_to_column("CellID")
travis.sc.seu.obj@meta.data <- inner_join(x=skin.meta.data, y=meta.data,by="CellID", keep=FALSE) %>% column_to_rownames("CellID")

# Subset Data to be only Psoriasis
travis.psoriatic_data <- subset(travis.sc.seu.obj, Condition == c("Psoriasis"))
```

Pre-process the subset object

``` r
## PROCESS SC RNA DATA
travis.psoriatic_data <- NormalizeData(travis.psoriatic_data)
travis.psoriatic_data <- FindVariableFeatures(travis.psoriatic_data, selection.method = "vst", nfeatures = 2000)
travis.psoriatic_data <- ScaleData(travis.psoriatic_data)
travis.psoriatic_data <- RunPCA(travis.psoriatic_data, features = VariableFeatures(object = travis.psoriatic_data))
travis.psoriatic_data <- FindNeighbors(travis.psoriatic_data, dims = 1:40)
travis.psoriatic_data <- FindClusters(travis.psoriatic_data, resolution = 1)
travis.psoriatic_data <- RunUMAP(travis.psoriatic_data, dims = 1:40)
```

``` r
saveRDS(travis.psoriatic_data, file="/Volumes/Extreme Pro/GITHUB-DATA/SC-RNA-DATA/TRAVIS-DATA/RDS-Files/TRAVIS_PSORIASIS_ONLY_DATA.RDS")
```

Import the processed data from the RDS file (PSORIASIS ONLY)

``` r
travis.psoriatic_data <- readRDS("/Volumes/Extreme Pro/GITHUB-DATA/SC-RNA-DATA/TRAVIS-DATA/RDS-Files/TRAVIS_PSORIASIS_ONLY_DATA.RDS")
```

``` r
# Save Markers Genes for Psoriatic Travis Data
Idents(travis.psoriatic_data) <- "Specific_CellType"
Marker_genes.TRAVIS <- FindAllMarkers(travis.psoriatic_data, max.cells.per.ident=1000, min.pct = 0.25, assay = "RNA")
saveRDS(Marker_genes.TRAVIS, file="/Volumes/Extreme Pro/GITHUB-DATA/SC-RNA-DATA/TRAVIS-DATA/RDS-Files/MARKER_GENES_PSORIASIS_ONLY_TRAVIS.RDS")}
```

### FIGURE S9B

### UMAP for (PS only) dataset-1

``` r
#pdf(width = 12,height=8,file = "../TRAVIS_scRNA_DATA/UMAP_TRAVIS_DATA_PS_SAMPLES_ONLY.pdf")
DimPlot(travis.psoriatic_data,group.by = "CellType", pt.size = 3.5) + scale_color_igv()
```

![](PS_SAMPLES_PART_3_TRAVIS_DATA_Final_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
#dev.off()
```

Load Travis (PS Only) Marker Genes

``` r
Marker_genes.TRAVIS <- readRDS("/Volumes/Extreme Pro/GITHUB-DATA/SC-RNA-DATA/TRAVIS-DATA/RDS-Files/MARKER_GENES_PSORIASIS_ONLY_TRAVIS.RDS")
```

``` r
filtered_single_cell.markers <- Marker_genes.TRAVIS %>% filter(p_val_adj <= 0.05) %>% group_by(cluster) %>% top_n(n = 300, wt = avg_log2FC) %>% filter(avg_log2FC>0.25)
filtered_spatial_markers <- skin_data.hm.sct.markers %>% filter(p_val_adj <= 0.05) %>% group_by(cluster) %>% top_n(n = 300, wt = avg_log2FC) %>% filter(avg_log2FC>0.25)

##INTERSECT GENES BETWEEN scRNA and Spatial data
st.genes <- unique(rownames(skin_data.hm.sct@assays$Spatial@counts))
sc.genes <- unique(rownames(travis.psoriatic_data@assays$RNA@counts))
all.genes.scrna_and_spt <- unique(intersect(sc.genes,st.genes))

MIA_results <- MIA(total_genes = length(all.genes.scrna_and_spt), single_cell.markers = filtered_single_cell.markers,spatial.markers = filtered_spatial_markers)

# MIA Plot for All Cell Types
E.data <- MIA_results %>% column_to_rownames("cluster")
E.data <- E.data[,order(colnames(E.data))]
pheatmap(E.data,cluster_cols = FALSE, cluster_rows = FALSE, fontsize=15, color = col.pal)
```

![](PS_SAMPLES_PART_3_TRAVIS_DATA_Final_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

### FIGURE 5B

### MIA ENRICHMENT (IMMUNE CELLS ONLY)

``` r
# 1. MIA Plot for Immune Cell Types Only
immune_only.E.data <- E.data[,c("B cell","T cell-1","T cell-2","T cell-3","Myeloid-1","Myeloid-2","Langerhans","Mast-1","Mast-2")]

#pdf(file = "MIA_regions_IMMUNE_CELLS_ALL_PS_SAMPLES.pdf",width = 10,height = 10)
pheatmap(immune_only.E.data,cluster_cols = FALSE,cluster_rows = FALSE,fontsize=15,color = col.pal)
```

![](PS_SAMPLES_PART_3_TRAVIS_DATA_Final_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
#dev.off()
```

### FIGURE S9D

### MIA ENRICHMENT (STRUCTURAL CELLS ONLY)

``` r
# 2. MIA Plot for Structural Cell Types Only
structure_only.E_data <- E.data[,c("Fibroblast-1","Fibroblast-2","Fibroblast-3","Fibroblast-4","Fibroblast-5","Fibroblast-6","HairFollicle","Keratinocyte-1","Keratinocyte-2","Keratinocyte-3","Keratinocyte-4","Keratinocyte-5","Keratinocyte-6","Keratinocyte-7","Keratinocyte-8","Sebocyte","Venule-1","Venule-2","Venule-3","Venule-4","Venule-5","VSMC-1","Lymphatic")]

#pdf(file = "MIA_regions_STRUCTURAL_CELLS_ALL_PS_SAMPLES.pdf",width = 10,height = 10)
pheatmap(structure_only.E_data,cluster_cols = FALSE,cluster_rows = FALSE,fontsize=15,color = col.pal)
```

![](PS_SAMPLES_PART_3_TRAVIS_DATA_Final_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
#dev.off()
```

### APPROACH \#2 - SEURAT scRNA INTEGRATION

``` r
# DO NOT RUN THIS LOCALLY
skin_reference.3 <- SCTransform(travis.psoriatic_data, ncells = 3000, verbose = FALSE) %>% 
  RunPCA(verbose = FALSE) %>% 
  RunUMAP(dims = 1:40)

anchors <- FindTransferAnchors(reference = skin_reference.3, query = skin_data.hm.sct, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = skin_reference.3$Specific_CellType, prediction.assay = TRUE, weight.reduction = skin_data.hm.sct[["harmony"]],dims = 1:40)
skin_data.hm.sct[["predictions_travis_data_harmony_version"]] <- predictions.assay
```

Save & Reload Seurat scRNA Integration RDS Object

``` r
# Save RDS 
#saveRDS(skin_data.hm.sct, file="/Volumes/Extreme Pro/GITHUB-DATA/SC-RNA-DATA/TRAVIS-DATA/RDS-Files/SKIN_SPATIAL_AND_SC_MERGED_SCORES_HARMONY_VERSION.RDS")

# Reload RDS
skin_data.hm.combined.v2 <- readRDS("/Volumes/Extreme Pro/GITHUB-DATA/SC-RNA-DATA/TRAVIS-DATA/RDS-Files/SKIN_SPATIAL_AND_SC_MERGED_SCORES_HARMONY_VERSION.RDS")
```

Define color scheme for spatial plots

``` r
SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))
cols <- SpatialColors(n = 100)
```

### FIGURE 5A

### CELL TYPE ENRICHMENT PLOTS FOR STRUCTURAL CELLS (AS SHOWN IN FINAL FIGURE)

``` r
Cell_types <- c("Fibroblast-1","Fibroblast-2","HairFollicle","Keratinocyte-1","Keratinocyte-2","Keratinocyte-3","Sebocyte","VSMC-1","Melanocyte","Lymphatic")
DefaultAssay(skin_data.hm.combined.v2) <- "predictions_travis_data_harmony_version"
```

HEALTHY SAMPLE

``` r
# Spatial Plot Healthy Sample
SpatialDimPlot(skin_data.hm.combined.v2, images = c("ST.HF.1.R2"))
```

![](PS_SAMPLES_PART_3_TRAVIS_DATA_Final_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

``` r
# Spatial Feature Plot of All Structural Cells
for (x in Cell_types){
  print(SpatialFeaturePlot(skin_data.hm.combined.v2, features = c(x), pt.size.factor = 2.5, ncol = 2, crop = TRUE, images = "ST.HF.1.R2") + 
          scale_fill_gradientn(colors=cols,limits = c(0,1.01))) 
}
```

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.
    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](PS_SAMPLES_PART_3_TRAVIS_DATA_Final_files/figure-gfm/unnamed-chunk-22-2.png)<!-- -->

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](PS_SAMPLES_PART_3_TRAVIS_DATA_Final_files/figure-gfm/unnamed-chunk-22-3.png)<!-- -->

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](PS_SAMPLES_PART_3_TRAVIS_DATA_Final_files/figure-gfm/unnamed-chunk-22-4.png)<!-- -->

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](PS_SAMPLES_PART_3_TRAVIS_DATA_Final_files/figure-gfm/unnamed-chunk-22-5.png)<!-- -->

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](PS_SAMPLES_PART_3_TRAVIS_DATA_Final_files/figure-gfm/unnamed-chunk-22-6.png)<!-- -->

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](PS_SAMPLES_PART_3_TRAVIS_DATA_Final_files/figure-gfm/unnamed-chunk-22-7.png)<!-- -->

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](PS_SAMPLES_PART_3_TRAVIS_DATA_Final_files/figure-gfm/unnamed-chunk-22-8.png)<!-- -->

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](PS_SAMPLES_PART_3_TRAVIS_DATA_Final_files/figure-gfm/unnamed-chunk-22-9.png)<!-- -->

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](PS_SAMPLES_PART_3_TRAVIS_DATA_Final_files/figure-gfm/unnamed-chunk-22-10.png)<!-- -->![](PS_SAMPLES_PART_3_TRAVIS_DATA_Final_files/figure-gfm/unnamed-chunk-22-11.png)<!-- -->

LESIONAL SKIN SAMPLE

``` r
# Spatial Plot LESIONAL SAMPLE
SpatialDimPlot(skin_data.hm.combined.v2, images = c("ST_21_L_Batch_6"))
```

![](PS_SAMPLES_PART_3_TRAVIS_DATA_Final_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

``` r
# Spatial Feature Plot of All Structural Cells
for (x in Cell_types){
  print(SpatialFeaturePlot(skin_data.hm.combined.v2, features = c(x), pt.size.factor = 2.5, ncol = 2, crop = TRUE,images = c("ST_21_L_Batch_6")) + 
          scale_fill_gradientn(colors=cols,limits = c(0,1.01))) 
}
```

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.
    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](PS_SAMPLES_PART_3_TRAVIS_DATA_Final_files/figure-gfm/unnamed-chunk-23-2.png)<!-- -->

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](PS_SAMPLES_PART_3_TRAVIS_DATA_Final_files/figure-gfm/unnamed-chunk-23-3.png)<!-- -->

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](PS_SAMPLES_PART_3_TRAVIS_DATA_Final_files/figure-gfm/unnamed-chunk-23-4.png)<!-- -->

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](PS_SAMPLES_PART_3_TRAVIS_DATA_Final_files/figure-gfm/unnamed-chunk-23-5.png)<!-- -->

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](PS_SAMPLES_PART_3_TRAVIS_DATA_Final_files/figure-gfm/unnamed-chunk-23-6.png)<!-- -->

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](PS_SAMPLES_PART_3_TRAVIS_DATA_Final_files/figure-gfm/unnamed-chunk-23-7.png)<!-- -->

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](PS_SAMPLES_PART_3_TRAVIS_DATA_Final_files/figure-gfm/unnamed-chunk-23-8.png)<!-- -->

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](PS_SAMPLES_PART_3_TRAVIS_DATA_Final_files/figure-gfm/unnamed-chunk-23-9.png)<!-- -->

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](PS_SAMPLES_PART_3_TRAVIS_DATA_Final_files/figure-gfm/unnamed-chunk-23-10.png)<!-- -->![](PS_SAMPLES_PART_3_TRAVIS_DATA_Final_files/figure-gfm/unnamed-chunk-23-11.png)<!-- -->

## CELL TYPE ENRICHMENT PLOTS FOR IMMUNE CELLS

``` r
Cell_types <- c("T cell-1","T cell-2","T cell-3","B cell","Myeloid-1","Myeloid-2","Langerhans")
DefaultAssay(skin_data.hm.combined.v2) <- "predictions_travis_data_harmony_version"
```

HEALTHY SAMPLE

``` r
## HEALTHY SAMPLE
for (x in Cell_types){
  print(SpatialFeaturePlot(skin_data.hm.combined.v2, features = c(x), pt.size.factor = 2.5, ncol = 2, crop = TRUE,images = "ST.HF.1.R2") + 
          scale_fill_gradientn(colors=cols,limits = c(0,1.01))) 
}
```

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.
    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](PS_SAMPLES_PART_3_TRAVIS_DATA_Final_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](PS_SAMPLES_PART_3_TRAVIS_DATA_Final_files/figure-gfm/unnamed-chunk-25-2.png)<!-- -->

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](PS_SAMPLES_PART_3_TRAVIS_DATA_Final_files/figure-gfm/unnamed-chunk-25-3.png)<!-- -->

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](PS_SAMPLES_PART_3_TRAVIS_DATA_Final_files/figure-gfm/unnamed-chunk-25-4.png)<!-- -->

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](PS_SAMPLES_PART_3_TRAVIS_DATA_Final_files/figure-gfm/unnamed-chunk-25-5.png)<!-- -->

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](PS_SAMPLES_PART_3_TRAVIS_DATA_Final_files/figure-gfm/unnamed-chunk-25-6.png)<!-- -->![](PS_SAMPLES_PART_3_TRAVIS_DATA_Final_files/figure-gfm/unnamed-chunk-25-7.png)<!-- -->

LESIONAL SKIN SAMPLE

``` r
for (x in Cell_types){
  print(SpatialFeaturePlot(skin_data.hm.combined.v2, features = c(x), pt.size.factor = 2.5, ncol = 2, crop = TRUE, images = c("ST_21_L_Batch_6")) + 
          scale_fill_gradientn(colors=cols,limits = c(0,1.01))) 
}
```

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.
    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](PS_SAMPLES_PART_3_TRAVIS_DATA_Final_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](PS_SAMPLES_PART_3_TRAVIS_DATA_Final_files/figure-gfm/unnamed-chunk-26-2.png)<!-- -->

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](PS_SAMPLES_PART_3_TRAVIS_DATA_Final_files/figure-gfm/unnamed-chunk-26-3.png)<!-- -->

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](PS_SAMPLES_PART_3_TRAVIS_DATA_Final_files/figure-gfm/unnamed-chunk-26-4.png)<!-- -->

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](PS_SAMPLES_PART_3_TRAVIS_DATA_Final_files/figure-gfm/unnamed-chunk-26-5.png)<!-- -->

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](PS_SAMPLES_PART_3_TRAVIS_DATA_Final_files/figure-gfm/unnamed-chunk-26-6.png)<!-- -->![](PS_SAMPLES_PART_3_TRAVIS_DATA_Final_files/figure-gfm/unnamed-chunk-26-7.png)<!-- -->

### FIGURE 5E

### B CELL SPECIFIC MARKERS - USED IN ENRICH PATHWAY ANALYSIS (PERFORMED USING ENRICHR)

EnrichR (<https://maayanlab.cloud/Enrichr/>)

``` r
B.CELL.MARKERS <- Marker_genes.TRAVIS %>% 
  filter(cluster == "B cell") %>%
  filter(avg_log2FC > 0.25)

write.csv(B.CELL.MARKERS,file="B_CELL_MARKERS.csv")
```

### FIGURE 5C

### TRM (BULK RNAseq data) - MIA Analysis

``` r
# Read in TRM Data
trm.signature <- read.csv("/Volumes/Extreme Pro/GITHUB-DATA/SC-RNA-DATA/TRM_DATA/TRM_SIGNATURE.csv")
trm.exclusive <- read.csv("/Volumes/Extreme Pro/GITHUB-DATA/SC-RNA-DATA/TRM_DATA/TRM_exclusive_up.csv")
trm.up <- read.csv("/Volumes/Extreme Pro/GITHUB-DATA/SC-RNA-DATA/TRM_DATA/TRM_up.csv")
trm.high <- read.csv("/Volumes/Extreme Pro/GITHUB-DATA/SC-RNA-DATA/TRM_DATA/TRM_high.csv")
```

``` r
# Extract Unique Genes from All TRM data frames
all.trm.genes <- c(trm.exclusive$Gene, trm.high$Gene,trm.up$Gene, trm.signature$Gene) %>% 
  unique()

# Filter Spatial Markers
filtered_spatial_markers <- skin_data.hm.sct.markers %>% 
  filter(p_val_adj <= 0.05) %>% 
  group_by(cluster) %>% 
  top_n(n = 350,wt = avg_log2FC) %>% 
  filter(avg_log2FC > 0.25)

##INTERSECT GENES BETWEEN scRNA and Spatial data
st.genes <- unique(rownames(skin_data.hm.sct@assays$Spatial@counts))

MIA_results <- MIA_bulk(markers = all.trm.genes, spatial.markers = filtered_spatial_markers,total_genes = length(st.genes), name = "TRM-all-genes")

E.data <- MIA_results %>% column_to_rownames("cluster")

#pdf(width = 7,height = 10,file = "TRM_MIA_PLOT.pdf")
pheatmap(E.data, cluster_cols = FALSE, cluster_rows = FALSE, fontsize = 15, color = col.pal)
```

![](PS_SAMPLES_PART_3_TRAVIS_DATA_Final_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

``` r
#dev.off()
```
