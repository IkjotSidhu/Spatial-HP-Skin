R Notebook
================

\##PSO + PSA + NORMAL SKIN COMBINED ANALYSIS - PART 2 \## CONTINUED

## Harmony-batch correction

## Table of content

| FIGURE NO | DESCRIPTION                                                             | LINK        |
|-----------|-------------------------------------------------------------------------|-------------|
| S4A       | UMAP WITH ALL SPATIAL SAMPLES (Using Harmony batch correction)          | [link](#3b) |
| 3B        | UMAP WITH ALL SPATIAL SAMPLES (Using Harmony batch correction)          | [link](#3b) |
| 3A        | NON-LESIONAL SKIN SAMPLE (ST 21 NL), and LESIONAL SKIN SAMPLE (ST 22 L) | [link](#3a) |
| S5        | Heatmap (After Harmony batch correction)                                | [link](#s5) |
| 3C        | Percentage composition                                                  | [link](#3c) |

### LOAD ALL PACKAGES

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
library(ggpubr)
```

    ## 
    ## Attaching package: 'ggpubr'

    ## The following object is masked from 'package:cowplot':
    ## 
    ##     get_legend

``` r
load("../all_samples_merge.RData")
```

### Running Harmony batch correction

``` r
DefaultAssay(all_samples.merge) <- "SCT"
all_samples.hm.sct <- RunHarmony(all_samples.merge,assay.use = "SCT",project.dim = FALSE,group.by.vars = "sample.id")
```

``` r
ElbowPlot(skin_data.hm.sct,ndims = 50)
```

``` r
skin_data.hm.sct <- skin_data.hm.sct %>% 
    RunUMAP(reduction = "harmony", dims = 1:40) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:40) %>% 
    FindClusters(resolution = 0.35) %>% 
    identity()
```

``` r
skin_data.hm.sct <- readRDS(file="../ALL_SPATIAL_SAMPLES.RDS")
```

``` r
## ADDING CLUSTER LABELS
cluster.labels <- c("0 Fibroblasts",
"1 Macs + fibroblasts",
"2 Eccrine + melanocyte precursors",
"3 Epidermis",
"4 Epidermis",
"5 Connective tissue",
"6 Mixed",
"7 Epidermis",
"8 Hair follicle and sebaceous glands",
"9 Adipose",
"10 Suprabasal keratinocytes",
"11 Smooth muscle",
"12 Endothelial cells",
"13 Immunoglobulins, fibroblasts",
"14 Smooth muscle",
"15 Mixed",
"16 Adipose, fibroblasts")

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

``` r
## COLOR FOR LABELS
color.bar <- c("lightskyblue",
"royalblue1",
"tan3","darkorchid1",
"violetred",
"paleturquoise1",
"azure2",
"darkorchid4",
"darkorange",
"gold",
"violet",
"tomato3",
"tan4",
"azure4",
"tomato",
"azure3",
"goldenrod1")
```

``` r
names(cluster.labels) <- levels(skin_data.hm.sct)
skin_data.hm.sct <- RenameIdents(skin_data.hm.sct, cluster.labels)
skin_data.hm.sct <- StashIdent(skin_data.hm.sct, save.name = "Spatial.regions")
```

    ## With Seurat 3.X, stashing identity classes can be accomplished with the following:
    ## skin_data.hm.sct[["Spatial.regions"]] <- Idents(object = skin_data.hm.sct)

``` r
Idents(skin_data.hm.sct) <- "Spatial.regions"

DimPlot(skin_data.hm.sct)
```

![](PS_SAMPLES_PART_2_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
skin_data.hm.sct <- PrepSCTFindMarkers(skin_data.hm.sct,assay = "SCT")
```

    ## Only one SCT model is stored - skipping recalculating corrected counts

``` r
skin_data.hm.sct.markers <- FindAllMarkers(skin_data.hm.sct, min.pct = 0.25, assay = "SCT") %>% filter(p_val_adj<0.05)
```

    ## Calculating cluster 0 Fibroblasts

    ## Calculating cluster 1 Macs + fibroblasts

    ## Calculating cluster 2 Eccrine + melanocyte precursors

    ## Calculating cluster 3 Epidermis

    ## Calculating cluster 4 Epidermis

    ## Calculating cluster 5 Connective tissue

    ## Calculating cluster 6 Mixed

    ## Calculating cluster 7 Epidermis

    ## Calculating cluster 8 Hair follicle and sebaceous glands

    ## Calculating cluster 9 Adipose

    ## Calculating cluster 10 Suprabasal keratinocytes

    ## Calculating cluster 11 Smooth muscle

    ## Calculating cluster 12 Endothelial cells

    ## Calculating cluster 13 Immunoglobulins, fibroblasts

    ## Calculating cluster 14 Smooth muscle

    ## Calculating cluster 15 Mixed

    ## Calculating cluster 16 Adipose, fibroblasts

``` r
levels(skin_data.hm.sct@meta.data$Spatial.regions) 
```

    ##  [1] "0 Fibroblasts"                       
    ##  [2] "1 Macs + fibroblasts"                
    ##  [3] "2 Eccrine + melanocyte precursors"   
    ##  [4] "3 Epidermis"                         
    ##  [5] "4 Epidermis"                         
    ##  [6] "5 Connective tissue"                 
    ##  [7] "6 Mixed"                             
    ##  [8] "7 Epidermis"                         
    ##  [9] "8 Hair follicle and sebaceous glands"
    ## [10] "9 Adipose"                           
    ## [11] "10 Suprabasal keratinocytes"         
    ## [12] "11 Smooth muscle"                    
    ## [13] "12 Endothelial cells"                
    ## [14] "13 Immunoglobulins, fibroblasts"     
    ## [15] "14 Smooth muscle"                    
    ## [16] "15 Mixed"                            
    ## [17] "16 Adipose, fibroblasts"

``` r
DimPlot(skin_data.hm.sct,split.by="DISEASE_STATUS",cols=color.labels,pt.size=3.5)
```

![](PS_SAMPLES_PART_2_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

## FIGURE 3B and S4A

\###UMAP WITH ALL SPATIAL SAMPLES (AFTER HARMONY BATCH CORRECTION)

<a id="3b">

``` r
DimPlot(skin_data.hm.sct,cols=color.labels,pt.size=3.5)
```

![](PS_SAMPLES_PART_2_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

</a>

## FIGURE 3A -

### NON-LESIONAL SKIN SAMPLE (ST 21 NL)

<a id="3a">

``` r
SpatialDimPlot(skin_data.hm.sct,images = c("ST_21_NL_Batch_6"),cols=color.labels,pt.size.factor = 2.5)
```

![](PS_SAMPLES_PART_2_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

## FIGURE 3A -

### LESIONAL SKIN SAMPLE (ST 22 L)

``` r
SpatialDimPlot(skin_data.hm.sct,images = c("ST_22L_Batch_8"),cols=color.labels,pt.size.factor = 2.5)
```

![](PS_SAMPLES_PART_2_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

</a>

``` r
DefaultAssay(skin_data.hm.sct) <- "SCT"
top10 <- skin_data.hm.sct.markers %>%
    group_by(cluster) %>%
    filter(gene %in% rownames(skin_data.hm.sct@assays$SCT@scale.data)) %>%
    top_n(n = 10, wt = avg_log2FC)
```

<a id="s5"> \## FIGURE SUPPLEMENTARY S5 \### HEATMAP showing top marker
genes after Harmony batch correction

``` r
DoHeatmap(skin_data.hm.sct, features = top10$gene,assay = "SCT",group.colors = color.labels,angle=90) + NoLegend()
```

![](PS_SAMPLES_PART_2_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

</a>

### Fraction of spots per group (DISEASE STATUS)

``` r
seurat_clusters.df.v2 <- table(skin_data.hm.sct@active.ident,skin_data.hm.sct@meta.data$DISEASE_STATUS,skin_data.hm.sct@meta.data$sample.id) %>% as.data.frame() %>%
  dplyr::rename(Cluster_id=Var1,Group=Var2,Sample_id=Var3) %>% filter(Freq!=0)  %>% group_by(Sample_id)%>% mutate(Fraction=Freq*100/sum(Freq)) %>% mutate(SUM_OF_FRACTIONS=sum(Fraction)) 
```

<a id="3c"> \## FIGURE 3c \### PERCENTAGE COMPOSITION PLOT

``` r
black.bold.16.text <- element_text(face = "bold", color = "black", size = 14,angle = 90, vjust = 0.5, hjust=1)
brks <- c(0, 0.25, 0.5, 0.75, 1)

ggplot(seurat_clusters.df.v2,aes(x=Group,y=Freq,fill=Cluster_id)) + geom_bar(stat="identity",position="fill") + ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(),axis.text.x =black.bold.16.text) + scale_y_continuous(breaks = brks, labels = scales::percent(brks)) + scale_fill_manual(values = color.labels)
```

![](PS_SAMPLES_PART_2_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->
</a>
