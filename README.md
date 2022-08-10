# Spatial-HP-Skin

This repo will cover all of the code linked to **Castillo et al.**

(**Manuscript Title:** *Spatial transcriptomics stratifies health and psoriatic disease severity by emergent cellular ecosystems*)

## GEO DATA LINK

All of the raw files linked to this project can be found under GEO Accession number - ([**GSE202011**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202011))

## CODE AND REPOSITORY NAVIGATION

The code is organised into several R notebooks with each notebook covering different parts of the analysis. The figure numbers and references are found within each notebook and a table content will also be added to link each figure to different parts of code for ease of navigation. Reproducibility is the key here, please contact **XXXX** for any issues with the code.

The notebooks will also require you to import the data from the linked Zenodo / Figshare repository (Will be added soon) in order to reproduce all of the respective figures.

Please make sure to install all of the dependencies prior to running the code by running the script **XXXX**.

Each notebook link will take you to the notebook with necessary code, please select the figure number in the table of contents found on top of the notebook to .

### FIGURE 1

+------------+--------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------+
| FIGURE NO  | NOTEBOOK LINK                                                                                                | DESCRIPTION                                                                      |
+============+==============================================================================================================+==================================================================================+
| 1B         | [link](https://github.com/IkjotSidhu/Spatial-HP-Skin/blob/main/ST_HEALTHY_SAMPLES_FIGURE_1.md)               | Spatial Plot - Healthy sample with cluster labels                                |
+------------+--------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------+
| 1C         | [link](https://github.com/IkjotSidhu/Spatial-HP-Skin/blob/main/ST_HEALTHY_SAMPLES_FIGURE_1.md)               | UMAP (Healthy samples only) - with spatial regions labels                        |
+------------+--------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------+
| 1D         | [link](https://github.com/IkjotSidhu/Spatial-HP-Skin/blob/main/ST_HEALTHY_SAMPLES_FIGURE_1.md)               | Percentage composition (Healthy samples only)                                    |
+------------+--------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------+
| 1E         | [link](https://github.com/IkjotSidhu/Spatial-HP-Skin/blob/main/ST_HEALTHY_SAMPLES_PART_2_TRAVIS_DATA_md)     | Seurat Integeration (scRNA + ST data) showing cell type enrichment               |
+------------+--------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------+
| 1G         | [link](https://github.com/IkjotSidhu/Spatial-HP-Skin/blob/main/ST_HEALTHY_SAMPLES_PART_2_TRAVIS_DATA_.md)    | MIA - healthy skin with dataset 1 (Hughes et al.) - Structural cell types only   |
+------------+--------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------+
| 1H         | [link](https://github.com/IkjotSidhu/Spatial-HP-Skin/blob/main/ST_HEALTHY_SAMPLES_PART_3(REYNOLDS_DATA).md)  | MIA - healthy skin with dataset 2 (Reynolds et al.) - Structural cell types only |
+------------+--------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------+
| 1I         | [link](https://github.com/IkjotSidhu/Spatial-HP-Skin/blob/main/ST_HEALTHY_SAMPLES_PART-4(SULLIVAN_DATA).Rmd) | MIA - healthy skin with dataset 3 (Sullivan et al.) - APCs and Preadiplocytes    |
+------------+--------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------+

### FIGURE 2

+------------+-------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------+
| FIGURE NO  | NOTEBOOK LINK                                                                                               | DESCRIPTION                                                                  |
+============+=============================================================================================================+==============================================================================+
| 2A         | [link](https://github.com/IkjotSidhu/Spatial-HP-Skin/blob/main/ST_HEALTHY_SAMPLES_PART_2_TRAVIS_DATA_.md)   | MIA - healthy skin with dataset 1 (Hughes et al.) - Immune cell types only   |
+------------+-------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------+
| 2B         | [link](https://github.com/IkjotSidhu/Spatial-HP-Skin/blob/main/ST_HEALTHY_SAMPLES_PART_3(REYNOLDS_DATA).md) | MIA - healthy skin with dataset 2 (Reynolds et al.) - Immune cell types only |
+------------+-------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------+

### FIGURE 3

+------------+--------------------------------------------------------------------------------------+--------------------------------------------------------------------------------+
| FIGURE NO  | NOTEBOOK LINK                                                                        | DESCRIPTION                                                                    |
+============+======================================================================================+================================================================================+
| 3A         | [link](https://github.com/IkjotSidhu/Spatial-HP-Skin/blob/main/PS_SAMPLES_PART_2.md) | Spatial Plots with cluster labels - Lesional and Non-Lesional skin             |
+------------+--------------------------------------------------------------------------------------+--------------------------------------------------------------------------------+
| 3B         | [link](https://github.com/IkjotSidhu/Spatial-HP-Skin/blob/main/PS_SAMPLES_PART_2.md) | UMAP - All samples combined (PS + Healthy skin) after Harmony batch correction |
+------------+--------------------------------------------------------------------------------------+--------------------------------------------------------------------------------+
| 3C         | [link](https://github.com/IkjotSidhu/Spatial-HP-Skin/blob/main/PS_SAMPLES_PART_2.md) | Percentage composition plots for different clusters / spatial regions          |
+------------+--------------------------------------------------------------------------------------+--------------------------------------------------------------------------------+

### SUPPLEMENTARY FIGURE 2

+------------+------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------+
| FIGURE NO  | NOTEBOOK LINK                                                                                  | DESCRIPTION                                                                |
+============+================================================================================================+============================================================================+
| S2-B       |                                                                                                |                                                                            |
+------------+------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------+
| S2-C       | [link](https://github.com/IkjotSidhu/Spatial-HP-Skin/blob/main/ST_HEALTHY_SAMPLES_FIGURE_1.md) | Heat-map showing top 8 marker genes per cluster (Healthy sample data only) |
+------------+------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------+
| S2-D       | [link](https://github.com/IkjotSidhu/Spatial-HP-Skin/blob/main/ST_HEALTHY_SAMPLES_FIGURE_1.md) | UMI counts                                                                 |
+------------+------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------+
|            |                                                                                                |                                                                            |
+------------+------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------+

### SUPPLEMENTARY FIGURE 3

+-----------+-------------------------------------------------------------------------------------------------------------+-------------+
| FIGURE NO | NOTEBOOK LINK                                                                                               | DESCRIPTION |
+===========+=============================================================================================================+=============+
| S3-A      | [link](https://github.com/IkjotSidhu/Spatial-HP-Skin/blob/main/ST_HEALTHY_SAMPLES_PART_2_TRAVIS_DATA_.md)   |             |
+-----------+-------------------------------------------------------------------------------------------------------------+-------------+
| S3-B      | [link](https://github.com/IkjotSidhu/Spatial-HP-Skin/blob/main/ST_HEALTHY_SAMPLES_PART_3(REYNOLDS_DATA).md) |             |
+-----------+-------------------------------------------------------------------------------------------------------------+-------------+

### SUPPLEMENTARY FIGURE 4

+------------+----------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------+
| FIGURE NO  | NOTEBOOK LINK                                                                                            | DESCRIPTION                                                        |
+============+==========================================================================================================+====================================================================+
| S4-A       | [link](https://github.com/IkjotSidhu/Spatial-HP-Skin/blob/main/PS_SAMPLES_PART_1.md) (Anchor data UMAP)  | UMAP plots - Harmony batch correction vs Anchor batch correction   |
|            |                                                                                                          |                                                                    |
|            | [link](https://github.com/IkjotSidhu/Spatial-HP-Skin/blob/main/PS_SAMPLES_PART_2.md) (Harmony data UMAP) |                                                                    |
+------------+----------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------+
| S4-C       | [link](https://github.com/IkjotSidhu/Spatial-HP-Skin/blob/main/PS_SAMPLES_PART_1.md) (Anchor data)       | Spatial plots- Harmony batch correction vs Anchor batch correction |
|            |                                                                                                          |                                                                    |
|            | [link](https://github.com/IkjotSidhu/Spatial-HP-Skin/blob/main/PS_SAMPLES_PART_2.md) (Harmony data)      |                                                                    |
+------------+----------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------+

### SUPPLEMENTARY FIGURE 5

+------------+--------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------+
| FIGURE NO  | NOTEBOOK LINK                                                                        | DESCRIPTION                                                                              |
+============+======================================================================================+==========================================================================================+
| S5         | [link](https://github.com/IkjotSidhu/Spatial-HP-Skin/blob/main/PS_SAMPLES_PART_2.md) | Heat-map showing top 10 marker genes per cluster (All samples) - Harmony batch corrected |
+------------+--------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------+

### SUPPLEMENTARY FIGURE 6

+------------+--------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------+
| FIGURE NO  | NOTEBOOK LINK                                                                        | DESCRIPTION                                                                             |
+============+======================================================================================+=========================================================================================+
| S6         | [link](https://github.com/IkjotSidhu/Spatial-HP-Skin/blob/main/PS_SAMPLES_PART_1.md) | Heat-map showing top 10 marker genes per cluster (All samples) - Anchor batch corrected |
+------------+--------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------+

### SUPPLEMENTARY FIGURE 7

+------------+--------------------------------------------------------------------------------------+--------------------------------------------------+
| FIGURE NO  | NOTEBOOK LINK                                                                        | DESCRIPTION                                      |
+============+======================================================================================+==================================================+
| S7-B       | [link](https://github.com/IkjotSidhu/Spatial-HP-Skin/blob/main/PS_SAMPLES_PART_2.md) | UMIs per clusters - Harmony batch corrected data |
+------------+--------------------------------------------------------------------------------------+--------------------------------------------------+
