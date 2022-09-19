R Notebook
================

# PS SAMPLES PART 4

## PSEUDO-BULK ANALYSIS

### *DIFFERENTIAL TESTING BETWEEN CONDITIONS - LESIONAL, NON-LESIONAL AND HEALTHY*

| FIGURE NO | DESCRIPTION                                                | LINK         |
|-----------|------------------------------------------------------------|--------------|
| 6A        | Hierarchical clustering based heatmap for all samples      | [link](#6a)  |
| 6B        | PCA plot for all samples - grouped by PS group             | [link](#6b)  |
| 6C        | PCA plot for all samples - grouped by Disease severity     | [link](#6c)  |
| 8A        | PCA plot for - Cluster 1 (Dermis macrophages, fibroblasts) | [link](#8a)  |
| 8B        | PCA plot for - Cluster 12 (Dermis lymphatics)              | [link](#8a)  |
| S9        | Dendogram for all samples                                  | [link](#s10) |
| S10       | PCA plots - per cluster basis (PS samples only)            | [link](#s10) |

``` r
library(DESeq2)
```

    ## Loading required package: S4Vectors

    ## Warning: package 'S4Vectors' was built under R version 4.1.2

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## Warning: package 'BiocGenerics' was built under R version 4.1.1

    ## 
    ## Attaching package: 'BiocGenerics'

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

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## Warning: package 'IRanges' was built under R version 4.1.1

    ## Loading required package: GenomicRanges

    ## Warning: package 'GenomicRanges' was built under R version 4.1.2

    ## Loading required package: GenomeInfoDb

    ## Warning: package 'GenomeInfoDb' was built under R version 4.1.2

    ## Loading required package: SummarizedExperiment

    ## Warning: package 'SummarizedExperiment' was built under R version 4.1.1

    ## Loading required package: MatrixGenerics

    ## Warning: package 'MatrixGenerics' was built under R version 4.1.1

    ## Loading required package: matrixStats

    ## Warning: package 'matrixStats' was built under R version 4.1.2

    ## 
    ## Attaching package: 'MatrixGenerics'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
    ##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
    ##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
    ##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
    ##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
    ##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
    ##     colWeightedMeans, colWeightedMedians, colWeightedSds,
    ##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
    ##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
    ##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
    ##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
    ##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
    ##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
    ##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
    ##     rowWeightedSds, rowWeightedVars

    ## Loading required package: Biobase

    ## Warning: package 'Biobase' was built under R version 4.1.1

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## 
    ## Attaching package: 'Biobase'

    ## The following object is masked from 'package:MatrixGenerics':
    ## 
    ##     rowMedians

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     anyMissing, rowMedians

``` r
library(Seurat)
```

    ## Warning: package 'Seurat' was built under R version 4.1.2

    ## Attaching SeuratObject

    ## Attaching sp

    ## 
    ## Attaching package: 'Seurat'

    ## The following object is masked from 'package:SummarizedExperiment':
    ## 
    ##     Assays

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──

    ## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
    ## ✔ tibble  3.1.7     ✔ dplyr   1.0.9
    ## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
    ## ✔ readr   2.1.2     ✔ forcats 0.5.1

    ## Warning: package 'ggplot2' was built under R version 4.1.2

    ## Warning: package 'tibble' was built under R version 4.1.2

    ## Warning: package 'tidyr' was built under R version 4.1.2

    ## Warning: package 'readr' was built under R version 4.1.2

    ## Warning: package 'dplyr' was built under R version 4.1.2

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::collapse()   masks IRanges::collapse()
    ## ✖ dplyr::combine()    masks Biobase::combine(), BiocGenerics::combine()
    ## ✖ dplyr::count()      masks matrixStats::count()
    ## ✖ dplyr::desc()       masks IRanges::desc()
    ## ✖ tidyr::expand()     masks S4Vectors::expand()
    ## ✖ dplyr::filter()     masks stats::filter()
    ## ✖ dplyr::first()      masks S4Vectors::first()
    ## ✖ dplyr::lag()        masks stats::lag()
    ## ✖ ggplot2::Position() masks BiocGenerics::Position(), base::Position()
    ## ✖ purrr::reduce()     masks GenomicRanges::reduce(), IRanges::reduce()
    ## ✖ dplyr::rename()     masks S4Vectors::rename()
    ## ✖ dplyr::slice()      masks IRanges::slice()

``` r
source("SPATIAL_FUNCTIONS.R")
```

    ## Warning: package 'RColorBrewer' was built under R version 4.1.2

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

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     slice

    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     rename

    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

    ## Loading required package: AnnotationDbi

    ## 
    ## Attaching package: 'AnnotationDbi'

    ## The following object is masked from 'package:clusterProfiler':
    ## 
    ##     select

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

    ## 

Import the Seurat object + Meta-data

``` r
skin_data.hm.sct <- readRDS(file="../ALL_SPATIAL_SAMPLES.RDS")

meta.data <- read.csv("all_skin_corrected.meta.data.csv")
skin_data.hm.sct$sample.id <- meta.data$sample.id
```

RE-NORMALIZING TO REGRESS OUT BATCH

``` r
skin_data.hm.sct_re_normalized <- SCTransform(skin_data.hm.sct,assay = "Spatial",new.assay.name ="SCT_BATCH_REGRESSED",vars.to.regress = c("sample.id"))
```

    ## Calculating cell attributes from input UMI matrix: log_umi

    ## Variance stabilizing transformation of count matrix of size 21901 by 16424

    ## Model formula is y ~ log_umi

    ## Get Negative Binomial regression parameters per gene

    ## Using 2000 genes, 5000 cells

    ##   |                                                                              |                                                                      |   0%  |                                                                              |==================                                                    |  25%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================================| 100%

    ## Found 70 outliers - those will be ignored in fitting/regularization step

    ## Second step: Get residuals using fitted parameters for 21901 genes

    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   2%  |                                                                              |===                                                                   |   5%  |                                                                              |=====                                                                 |   7%  |                                                                              |======                                                                |   9%  |                                                                              |========                                                              |  11%  |                                                                              |==========                                                            |  14%  |                                                                              |===========                                                           |  16%  |                                                                              |=============                                                         |  18%  |                                                                              |==============                                                        |  20%  |                                                                              |================                                                      |  23%  |                                                                              |==================                                                    |  25%  |                                                                              |===================                                                   |  27%  |                                                                              |=====================                                                 |  30%  |                                                                              |======================                                                |  32%  |                                                                              |========================                                              |  34%  |                                                                              |=========================                                             |  36%  |                                                                              |===========================                                           |  39%  |                                                                              |=============================                                         |  41%  |                                                                              |==============================                                        |  43%  |                                                                              |================================                                      |  45%  |                                                                              |=================================                                     |  48%  |                                                                              |===================================                                   |  50%  |                                                                              |=====================================                                 |  52%  |                                                                              |======================================                                |  55%  |                                                                              |========================================                              |  57%  |                                                                              |=========================================                             |  59%  |                                                                              |===========================================                           |  61%  |                                                                              |=============================================                         |  64%  |                                                                              |==============================================                        |  66%  |                                                                              |================================================                      |  68%  |                                                                              |=================================================                     |  70%  |                                                                              |===================================================                   |  73%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================                |  77%  |                                                                              |========================================================              |  80%  |                                                                              |=========================================================             |  82%  |                                                                              |===========================================================           |  84%  |                                                                              |============================================================          |  86%  |                                                                              |==============================================================        |  89%  |                                                                              |================================================================      |  91%  |                                                                              |=================================================================     |  93%  |                                                                              |===================================================================   |  95%  |                                                                              |====================================================================  |  98%  |                                                                              |======================================================================| 100%

    ## Computing corrected count matrix for 21901 genes

    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   2%  |                                                                              |===                                                                   |   5%  |                                                                              |=====                                                                 |   7%  |                                                                              |======                                                                |   9%  |                                                                              |========                                                              |  11%  |                                                                              |==========                                                            |  14%  |                                                                              |===========                                                           |  16%  |                                                                              |=============                                                         |  18%  |                                                                              |==============                                                        |  20%  |                                                                              |================                                                      |  23%  |                                                                              |==================                                                    |  25%  |                                                                              |===================                                                   |  27%  |                                                                              |=====================                                                 |  30%  |                                                                              |======================                                                |  32%  |                                                                              |========================                                              |  34%  |                                                                              |=========================                                             |  36%  |                                                                              |===========================                                           |  39%  |                                                                              |=============================                                         |  41%  |                                                                              |==============================                                        |  43%  |                                                                              |================================                                      |  45%  |                                                                              |=================================                                     |  48%  |                                                                              |===================================                                   |  50%  |                                                                              |=====================================                                 |  52%  |                                                                              |======================================                                |  55%  |                                                                              |========================================                              |  57%  |                                                                              |=========================================                             |  59%  |                                                                              |===========================================                           |  61%  |                                                                              |=============================================                         |  64%  |                                                                              |==============================================                        |  66%  |                                                                              |================================================                      |  68%  |                                                                              |=================================================                     |  70%  |                                                                              |===================================================                   |  73%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================                |  77%  |                                                                              |========================================================              |  80%  |                                                                              |=========================================================             |  82%  |                                                                              |===========================================================           |  84%  |                                                                              |============================================================          |  86%  |                                                                              |==============================================================        |  89%  |                                                                              |================================================================      |  91%  |                                                                              |=================================================================     |  93%  |                                                                              |===================================================================   |  95%  |                                                                              |====================================================================  |  98%  |                                                                              |======================================================================| 100%

    ## Calculating gene attributes

    ## Wall clock passed: Time difference of 2.269348 mins

    ## Determine variable features

    ## Place corrected count matrix in counts slot

    ## Regressing out sample.id

    ## Centering data matrix

    ## Warning: Keys should be one or more alphanumeric characters followed by an
    ## underscore, setting key from sct_batch_regressed_ to sctbatchregressed_

    ## Set default assay to SCT_BATCH_REGRESSED

IDENTIFY DE GENES - LESIONAL VS HEALTHY

``` r
LES_vs_HEALTHY_v1 <- FindMarkers(skin_data.hm.sct_re_normalized,group.by = "DISEASE_STATUS",ident.1 = "Lesional",ident.2 = "Healthy skin")
LES_vs_HEALTHY_v2  <- FindMarkers(skin_data.hm.sct_re_normalized,group.by = "DISEASE_STATUS",ident.1 = "Lesional",ident.2 = "Healthy skin")
```

``` r
NON_LES_vs_HEALTHY <- FindMarkers(skin_data.hm.sct_re_normalized,group.by = "DISEASE_STATUS",ident.1 = "Non-Lesional",ident.2 = "Healthy skin")
LES_VS_NON_LES <- FindMarkers(skin_data.hm.sct_re_normalized,group.by = "DISEASE_STATUS",ident.1 = "Lesional",ident.2 = "Non-Lesional")
```

## PSEUDO-BULK APPROACH

Import Pseudo-counts from the Seurat object

``` r
pseudo.counts <- Seurat:::PseudobulkExpression(object = skin_data.hm.sct,assays = "Spatial",group.by = "sample.id",slot = "counts",pb.method = "aggregate") 

summary(colSums(pseudo.counts$Spatial))
```

    ##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    ##   211576   963567  1826699  2826222  4051790 10120020

``` r
colSums(pseudo.counts$Spatial>10)
```

    ##      HV1_S1_R1      HV1_S1_R2      HV2_S1_R1      HV2_S1_R2         HV2_S2 
    ##           9100           9105          10190          12563           7489 
    ##         HV3_S1         HV3_S2     PSA_LES_P1  PSA_LES_P2_R1  PSA_LES_P2_R2 
    ##          11816          10507          13152          10717          13651 
    ##  PSA_LES_P3_R1  PSA_LES_P3_R2     PSA_LES_P4     PSA_LES_P5 PSA_NON_LES_P1 
    ##           8690          13803          12429          12463          10040 
    ## PSA_NON_LES_P2 PSA_NON_LES_P3 PSA_NON_LES_P4     PSO_LES_P1  PSO_LES_P2_R1 
    ##           8428          12828          11984          12697           8031 
    ##  PSO_LES_P2_R2     PSO_LES_P3     PSO_LES_P4     PSO_LES_P5     PSO_LES_P6 
    ##          13411          10764          13216          12779          13118 
    ## PSO_NON_LES_P1 PSO_NON_LES_P2 PSO_NON_LES_P3 PSO_NON_LES_P4 PSO_NON_LES_P5 
    ##           7391           3140           3910          10904          11626

``` r
pseudo.counts.df <- as.data.frame(pseudo.counts$Spatial) %>% rownames_to_column("Gene")

Genes <- pseudo.counts.df$Gene

library(DESeq2)

groups.table <- read.csv(file="PSEUDO-BULK-DATA/groups.table.csv",stringsAsFactors = TRUE) %>% filter(Sample.ID!="") %>% column_to_rownames("Sample.ID") %>% unite("DISEASE_and_SEVERITY","GROUP_I","SEVERITY", sep = "_", remove = FALSE, na.rm = FALSE)
```

``` r
coldata <- groups.table[,c("DISEASE_and_SEVERITY","GROUP_I","SEVERITY","BATCH")]
```

``` r
rownames(groups.table)
```

    ##  [1] "HV1_S1_R1"      "HV2_S1_R1"      "HV3_S1"         "HV2_S2"        
    ##  [5] "HV1_S1_R2"      "HV2_S1_R2"      "HV3_S2"         "PSO_LES_P2_R1" 
    ##  [9] "PSO_NON_LES_P2" "PSO_LES_P2_R2"  "PSA_LES_P5"     "PSA_LES_P2_R1" 
    ## [13] "PSA_NON_LES_P2" "PSA_LES_P2_R2"  "PSO_LES_P3"     "PSO_NON_LES_P3"
    ## [17] "PSA_LES_P4"     "PSA_NON_LES_P4" "PSO_LES_P6"     "PSO_LES_P1"    
    ## [21] "PSO_NON_LES_P1" "PSA_LES_P3_R1"  "PSA_NON_LES_P3" "PSA_LES_P3_R2" 
    ## [25] "PSO_LES_P5"     "PSO_NON_LES_P5" "PSA_LES_P1"     "PSA_NON_LES_P1"
    ## [29] "PSO_LES_P4"     "PSO_NON_LES_P4"

``` r
colnames(pseudo.counts.df)
```

    ##  [1] "Gene"           "HV1_S1_R1"      "HV1_S1_R2"      "HV2_S1_R1"     
    ##  [5] "HV2_S1_R2"      "HV2_S2"         "HV3_S1"         "HV3_S2"        
    ##  [9] "PSA_LES_P1"     "PSA_LES_P2_R1"  "PSA_LES_P2_R2"  "PSA_LES_P3_R1" 
    ## [13] "PSA_LES_P3_R2"  "PSA_LES_P4"     "PSA_LES_P5"     "PSA_NON_LES_P1"
    ## [17] "PSA_NON_LES_P2" "PSA_NON_LES_P3" "PSA_NON_LES_P4" "PSO_LES_P1"    
    ## [21] "PSO_LES_P2_R1"  "PSO_LES_P2_R2"  "PSO_LES_P3"     "PSO_LES_P4"    
    ## [25] "PSO_LES_P5"     "PSO_LES_P6"     "PSO_NON_LES_P1" "PSO_NON_LES_P2"
    ## [29] "PSO_NON_LES_P3" "PSO_NON_LES_P4" "PSO_NON_LES_P5"

``` r
# Ordering the file
counts.file <- pseudo.counts.df[,rownames(groups.table)]
#counts.file <- mutate_all(counts.file, function(x) as.numeric(as.character(x)))

group.name <- colnames(groups.table)[1]
groups.levels <- factor(groups.table[,group.name]) %>% levels()
```

``` r
rownames(counts.file) <- Genes
counts.final <- counts.file
```

``` r
# Run DE here
dds <- DESeqDataSetFromMatrix(countData = counts.final,colData=groups.table,design = ~BATCH + GROUP_I)
dds <- DESeq(dds,quiet = TRUE)
```

``` r
normalized.counts <- counts(dds, normalized=TRUE)
colSums(normalized.counts)
```

    ##      HV1_S1_R1      HV2_S1_R1         HV3_S1         HV2_S2      HV1_S1_R2 
    ##        1794826        1740896        1544520        1374940        1413584 
    ##      HV2_S1_R2         HV3_S2  PSO_LES_P2_R1 PSO_NON_LES_P2  PSO_LES_P2_R2 
    ##        1427695        1596350        1743215        1695625        1682631 
    ##     PSA_LES_P5  PSA_LES_P2_R1 PSA_NON_LES_P2  PSA_LES_P2_R2     PSO_LES_P3 
    ##        1879940        3007601        1637527        2111038        2014346 
    ## PSO_NON_LES_P3     PSA_LES_P4 PSA_NON_LES_P4     PSO_LES_P6     PSO_LES_P1 
    ##        1871857        1663323        1505774        1818696        2911079 
    ## PSO_NON_LES_P1  PSA_LES_P3_R1 PSA_NON_LES_P3  PSA_LES_P3_R2     PSO_LES_P5 
    ##        1583144        1958126        1684964        2165046        2297236 
    ## PSO_NON_LES_P5     PSA_LES_P1 PSA_NON_LES_P1     PSO_LES_P4 PSO_NON_LES_P4 
    ##        1663951        2658807        1609042        1845900        1512737

``` r
vsd <- vst(dds, blind=TRUE)
```

### FIGURE 6B / 6C <a id="6b"></a><a id="6c"></a>

PCA plots showing clustering of samples based on disease group and
severity / PASI score.

``` r
plotPCA(vsd, intgroup=c("GROUP_I","GROUP_II"))+ scale_color_manual(values=c("dimgray","#cc3333","#ffcc33","#ff9999","#ff9933"))
```

![](PS_SAMPLES_PART_4_PSEUDO_BULK_ANALYSIS__files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
pcaData <- plotPCA(vsd, intgroup=c("SEVERITY","PASI_SCORE","DISEASE_and_SEVERITY"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

#pdf("ALL_SAMPLES_PCA_PLOT.pdf",height = 8,width = 10)
ggplot(pcaData, aes(PC1, PC2, fill=PASI_SCORE, shape=DISEASE_and_SEVERITY)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  scale_shape_manual(values=c(21,22,23,24,25))
```

![](PS_SAMPLES_PART_4_PSEUDO_BULK_ANALYSIS__files/figure-gfm/unnamed-chunk-17-2.png)<!-- -->

``` r
#dev.off()
```

``` r
sampleDists <- dist(t(assay(vsd)))
```

### FIGURE 6A <a id="6a"></a>

### Hierarchical Clustering - Samples cluster by disease severity not systemic co-morbidity

``` r
library("RColorBrewer")
library("pheatmap")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$GROUP_I, vsd$SEVERITY,vsd$GROUP_II, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(brewer.pal(9, "RdYlBu")) (255)
#pdf("HC_WITH_HEATMAP.pdf",height = 15,width = 15)
print(pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors))
```

![](PS_SAMPLES_PART_4_PSEUDO_BULK_ANALYSIS__files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

``` r
#dev.off()
```

### FIGURE S9 <a id="s9"></a>

## DENDOGRAM

``` r
hc <- hclust(sampleDists)

print(plot(hc, labels=paste(vsd$GROUP_I,vsd$SEVERITY,vsd$GROUP_II,sep = ":")))
```

![](PS_SAMPLES_PART_4_PSEUDO_BULK_ANALYSIS__files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

    ## NULL

``` r
cluster.labels <- unique(skin_data.hm.sct@meta.data$Spatial.regions) %>% as.vector() %>% sort(decreasing = TRUE)
```

``` r
cluster.labels.filtered <- c("9 Adipose","10 Suprabasal keratinocytes","11 Smooth muscle","12 Endothelial cells","13 Immunoglobulins, fibroblasts","14 Smooth muscle","15 Mixed","16 Adipose, fibroblasts")
```

### Figure S9 <a id="s9"></a>

### CLUSTER SPECIFIC (PSEUDO-BULK) PCA-plots

#### PART 1

``` r
for(x in cluster.labels){
  subset.data <- subset(skin_data.hm.sct,Spatial.regions %in% c(x))
  dds <- pseudo_bulk_out(subset.data,group_label = "sample.id",groups_tbl_path = "PSEUDO-BULK-DATA/groups.table.csv")
  if(!is.null(dds)){
    tryCatch({
      vsd <- vst(dds, blind=TRUE)
      
      # ALL PSORIASIS SAMPLES
      vsd_subset <- vsd[,vsd$GROUP_I %in% c("Lesional","Non-Lesional")]
      
        ## PCA PLOT
      #pdf(file=paste("PSEUDO_BULK_OUTPUT/ALL_PS_SAMPLES/",x,"_PCA_PLOT_PSEUDOBULK_ALL_SAMPLES.pdf"),height = 8,width = 10)
      print(plotPCA(vsd_subset, intgroup=c("GROUP_I", "SEVERITY"))+ ggtitle(paste("PCA Plot for -",x))+ scale_color_manual(values=c("#ff9999","#cc3333","#ffcc33","#ff9933")))
      #dev.off()
      
  }, error=function(e){ skip_to_next <<- TRUE})
  }
  
}
```

![](PS_SAMPLES_PART_4_PSEUDO_BULK_ANALYSIS__files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->![](PS_SAMPLES_PART_4_PSEUDO_BULK_ANALYSIS__files/figure-gfm/unnamed-chunk-23-2.png)<!-- -->![](PS_SAMPLES_PART_4_PSEUDO_BULK_ANALYSIS__files/figure-gfm/unnamed-chunk-23-3.png)<!-- -->![](PS_SAMPLES_PART_4_PSEUDO_BULK_ANALYSIS__files/figure-gfm/unnamed-chunk-23-4.png)<!-- -->![](PS_SAMPLES_PART_4_PSEUDO_BULK_ANALYSIS__files/figure-gfm/unnamed-chunk-23-5.png)<!-- -->![](PS_SAMPLES_PART_4_PSEUDO_BULK_ANALYSIS__files/figure-gfm/unnamed-chunk-23-6.png)<!-- -->

#### PART-2

``` r
#cluster.ids <- skin_data.hm.sct@meta.data$Spatial.regions %>% unique() %>% as.vector()
for(x in cluster.labels.filtered){
  subset.data <- subset(skin_data.hm.sct,Spatial.regions %in% c(x))
  dds <- pseudo_bulk_out(subset.data,group_label = "sample.id",groups_tbl_path = "PSEUDO-BULK-DATA/groups.table.csv")
  if(!is.null(dds)){
    tryCatch({
      vsd <- vst(dds, blind=TRUE,nsub = 100)
      normalized.counts <- counts(dds, normalized=TRUE)
      
      # ONLY PSORIASIS SAMPLES
      vsd_subset <- vsd[,vsd$GROUP_I %in% c("Lesional","Non-Lesional")]
      ## PCA PLOT
      #pdf(file=paste("PSEUDO_BULK_OUTPUT/ALL_SAMPLES(VERSION_2)/",x,"_PCA_PLOT_PSEUDOBULK_ALL_SAMPLES.pdf"),height = 8,width = 10)
      print(plotPCA(vsd_subset, intgroup=c("GROUP_I", "SEVERITY"))+ ggtitle(paste("PCA Plot for -",x))+ scale_color_manual(values=c("#ff9999","#cc3333","#ffcc33","#ff9933"))) 
      #dev.off()
      
  }, error=function(e){ skip_to_next <<- TRUE})
  }
  
}
```

![](PS_SAMPLES_PART_4_PSEUDO_BULK_ANALYSIS__files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->![](PS_SAMPLES_PART_4_PSEUDO_BULK_ANALYSIS__files/figure-gfm/unnamed-chunk-24-2.png)<!-- -->![](PS_SAMPLES_PART_4_PSEUDO_BULK_ANALYSIS__files/figure-gfm/unnamed-chunk-24-3.png)<!-- -->![](PS_SAMPLES_PART_4_PSEUDO_BULK_ANALYSIS__files/figure-gfm/unnamed-chunk-24-4.png)<!-- -->![](PS_SAMPLES_PART_4_PSEUDO_BULK_ANALYSIS__files/figure-gfm/unnamed-chunk-24-5.png)<!-- -->![](PS_SAMPLES_PART_4_PSEUDO_BULK_ANALYSIS__files/figure-gfm/unnamed-chunk-24-6.png)<!-- -->![](PS_SAMPLES_PART_4_PSEUDO_BULK_ANALYSIS__files/figure-gfm/unnamed-chunk-24-7.png)<!-- -->![](PS_SAMPLES_PART_4_PSEUDO_BULK_ANALYSIS__files/figure-gfm/unnamed-chunk-24-8.png)<!-- -->

### FIGURE 8A and 8B<a id="8a"></a>

``` r
subset.data <- subset(skin_data.hm.sct,Spatial.regions %in% c("1 Macs + fibroblasts"))

dds <- pseudo_bulk_out(subset.data,group_label = "sample.id",groups_tbl_path = "PSEUDO-BULK-DATA/groups.table.csv")

vsd <- vst(dds, blind=TRUE,nsub = 100)

# ALL PSORIASIS SAMPLES
vsd_subset_2 <- vsd[,vsd$GROUP_I %in% c("Lesional","Non-Lesional")]
#pdf(file=paste("PSEUDO_BULK_OUTPUT/ALL_PS_SAMPLES/","1 Macs + fibroblasts","_PCA_PLOT_PSEUDOBULK_ALL_SAMPLES(WITH_LABELS)(V2).pdf"),height = 8,width = 10)
plotPCA(vsd_subset_2, intgroup=c("GROUP_I", "SEVERITY")) + geom_text(label = vsd_subset_2$PASI_SCORE,nudge_x = 0.50, nudge_y = 1) + scale_color_manual(values=c("#ff9999","#cc3333","#ffcc33","#ff9933"))
```

![](PS_SAMPLES_PART_4_PSEUDO_BULK_ANALYSIS__files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

``` r
#dev.off()

#12 Endothelial cells
subset.data <- subset(skin_data.hm.sct,Spatial.regions %in% c("12 Endothelial cells"))
dds <- pseudo_bulk_out(subset.data,group_label = "sample.id",groups_tbl_path = "PSEUDO-BULK-DATA/groups.table.csv")
vsd <- vst(dds, blind=TRUE,nsub = 100)

# ALL PSORIASIS SAMPLES
vsd_subset_2 <- vsd[,vsd$GROUP_I %in% c("Lesional","Non-Lesional")]
#pdf(file=paste("PSEUDO_BULK_OUTPUT/ALL_PS_SAMPLES/","12 Endothelial cells","_PCA_PLOT_PSEUDOBULK_ALL_SAMPLES(WITH_LABELS)(V2).pdf"),height = 8,width = 10)
plotPCA(vsd_subset_2, intgroup=c("GROUP_I", "SEVERITY")) + geom_text(label = vsd_subset_2$PASI_SCORE,nudge_x = 0.50, nudge_y = 1)  + scale_color_manual(values=c("#ff9999","#cc3333","#ffcc33","#ff9933"))
```

![](PS_SAMPLES_PART_4_PSEUDO_BULK_ANALYSIS__files/figure-gfm/unnamed-chunk-25-2.png)<!-- -->

``` r
#dev.off()
```
