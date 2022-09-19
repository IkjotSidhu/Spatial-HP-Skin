library(tidyverse)
library(Seurat)
library(cowplot)
library(tidyverse)
library(RColorBrewer)
library(reticulate)


st_filter_by_genes <- function(st.data,x){
  st.data <- subset(st.data, subset = nFeature_Spatial > x)
  pdf(file=paste(unique(st.data@meta.data$orig.ident),"_filtered_by",x,".pdf"),width = 10,height=12)
  print(SpatialDimPlot(st.data,pt.size.factor = 2.3))
  dev.off()
  return(st.data)
}


st_single_cell_integerate <- function(data.list){
  for (i in 1:length(data.list)) {
    data.list[[i]] <- SCTransform(data.list[[i]], verbose = FALSE,assay="Spatial")
  }
  features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = 3000)
  data.list <- PrepSCTIntegration(object.list = data.list, anchor.features = features)
  data.anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT", 
                                         anchor.features = features)
  data.combined.sct <- IntegrateData(anchorset = data.anchors, normalization.method = "SCT")
  return(data.combined.sct)
}


st_plot_genes <- function(genes,st.data){
  for(x in genes){
    pdf(file =paste(st.data,x,"_SPATIAL_EXPRESSION_PLOT.pdf"),width = 5,height=7)
    print(SpatialFeaturePlot(st.data,features = c(x),pt.size.factor = 2.3))
    dev.off()
  }
}

st_plot <- function(st.data){
  pdf(file =paste(unique(st.data@meta.data$orig.ident),"_SPATIAL_PLOT.pdf"),width = 10,height=12)
  print(SpatialPlot(st.data,pt.size.factor =1.9))
  dev.off()
}

st_scatter_QC <- function(st.data){
  DATA.COUNTS <- data.frame(UMI_Count=st.data$nCount_Spatial,FEATURE_Count=st.data$nFeature_Spatial)
  pdf(file=paste(unique(st.data@meta.data$orig.ident),"_scatter_plot_QC",".pdf"),width = 7,height=7)
  print(ggplot2::ggplot(DATA.COUNTS,aes(UMI_Count,FEATURE_Count)) + geom_point() + xlim(0,125000) + ylim(0,8500))
  dev.off()
}

st_plot_QC <- function(st.data){
  pdf(file =paste(unique(st.data@meta.data$orig.ident),"_nCount_VIOLIN_PLOT.pdf"),width = 5,height=7)
  print(VlnPlot(st.data, features = "nCount_Spatial", pt.size = 0.1) + NoLegend())
  dev.off()
  pdf(file =paste(unique(st.data@meta.data$orig.ident),"_nCount_SPATIAL_PLOT.pdf"),width = 6,height=8)
  print(SpatialFeaturePlot(st.data, features = "nCount_Spatial") + theme(legend.position = "right"))
  dev.off()
  pdf(file =paste(unique(st.data@meta.data$orig.ident),"_nFeature_VIOLIN_PLOT.pdf"),width = 5,height=7)
  print(VlnPlot(st.data, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend())
  dev.off()
  pdf(file =paste(unique(st.data@meta.data$orig.ident),"_nFeature_SPATIAL_PLOT.pdf"),width = 6,height=8)
  print(SpatialFeaturePlot(st.data, features = "nFeature_Spatial") + theme(legend.position = "right"))
  dev.off()
}

st_filter_spots <- function(st.data,x){
  remove.spots <- read.csv(file=x)
  subset_spots <- Cells(st.data)[which((!(rownames(st.data@meta.data) %in% remove.spots$Barcode)))]
  st.data <- subset(st.data,cells=subset_spots)
  return(st.data)
}

st_combine <- function(data.list,ndim,res,seed=100){
  set.seed(seed)
  for (i in 1:length(data.list)) {
    data.list[[i]] <- SCTransform(data.list[[i]], verbose = FALSE,assay="Spatial")
  }
  features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = 3000)
  data.list <- PrepSCTIntegration(object.list = data.list, anchor.features = features)
  data.anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT", 
                                         anchor.features = features)
  data.combined.sct <- IntegrateData(anchorset = data.anchors, normalization.method = "SCT")
  data.combined.sct <- RunPCA(data.combined.sct, verbose = FALSE)
  data.combined.sct <- FindNeighbors(data.combined.sct, dims = 1:ndim)
  data.combined.sct <- FindClusters(data.combined.sct, verbose = FALSE,resolution = res)
  data.combined.sct <- RunUMAP(data.combined.sct, dims = 1:ndim)
  return(data.combined.sct)
}

MIA <- function(stlist,sclist,total=14705){
  overlap <- length(intersect(stlist$gene,sclist$gene))
  C <- length(sclist$gene)
  D <- length(stlist$gene)
  
  # ENRICHMENT CALCULATION
  e <- -log10(phyper(overlap,C, total - C, D, lower.tail = FALSE))
  E.SCORES <- append(E.SCORES,e)
  return(e)
}

library(clusterProfiler)
library(org.Hs.eg.db)


GO_PATHWAYS_ENRICH <- function(MARKERS,ONT,OUTPUT){
  ## gene UPREGULATED
  up.gene <- MARKERS %>% filter(avg_log2FC>0.25) %>% arrange(desc(avg_log2FC)) %>% dplyr::select(c("gene"))
  
  ## geneS DOWNREGULATED
  down.gene <- MARKERS %>% filter(avg_log2FC<(-0.25)) %>% arrange(avg_log2FC) %>% dplyr::select(c("gene"))
  
  GO_res_UP <- enrichGO(up.gene$gene,OrgDb = org.Hs.eg.db,keyType = 'SYMBOL',pvalueCutoff = 0.1,ont= ONT) 
  
  GO_res_DOWN <- enrichGO(down.gene$gene,OrgDb = org.Hs.eg.db,keyType = 'SYMBOL',pvalueCutoff = 0.1,ont= ONT)
  if(is.null(GO_res_UP)){
  }
  else{
    GO_res_UP <- as.data.frame(GO_res_UP)
    write.csv(GO_res_UP,file=paste(deparse(substitute(MARKERS)),OUTPUT,"_UP_GO_RESULTS_",ONT,".csv",sep=""))
  }
  if(is.null(GO_res_DOWN)){
  } else {
    GO_res_DOWN <- as.data.frame(GO_res_DOWN)
    write.csv(GO_res_DOWN,file=paste(deparse(substitute(MARKERS)),OUTPUT,"_DOWN_GO_RESULTS_",ONT,".csv",sep=""))
  }
  #results_list <- c(GO_res_UP,GO_res_DOWN)
  return(NULL)
}

KEGG_PATHWAYS <- function(MARKERS,OUTPUT){
  up.gene <- MARKERS %>% filter(avg_log2FC>0.25) %>% arrange(desc(avg_log2FC)) %>% dplyr::select(c("gene"))
  ge_list_2 <- bitr(up.gene$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb =org.Hs.eg.db,drop=TRUE)
  KEGG_res_UP <- enrichKEGG(ge_list_2$ENTREZID,organism= 'hsa',pvalueCutoff = 0.1)
  if (is.null(KEGG_res_UP)){
  } else {
    KEGG_res_UP <- setReadable(KEGG_res_UP, OrgDb = org.Hs.eg.db,keyType = "ENTREZID") %>% as.data.frame()
    write.csv(KEGG_res_UP,file=paste(deparse(substitute(MARKERS)),OUTPUT,"_UP_KEGG_RESULTS_",".csv",sep=""))
  }
  down.gene <- MARKERS %>% filter(avg_log2FC<(-0.25)) %>% arrange(avg_log2FC) %>% dplyr::select(c("gene"))
  ge_list_2 <- bitr(down.gene$gene,fromType = "SYMBOL",toType ="ENTREZID",OrgDb =org.Hs.eg.db,drop=TRUE)
  KEGG_res_DOWN <- enrichKEGG(ge_list_2$ENTREZID,organism= 'hsa',pvalueCutoff = 0.1)
  if (is.null(KEGG_res_DOWN)) {
    return(NULL)
  } else {
    KEGG_res_DOWN <- setReadable(KEGG_res_DOWN, OrgDb = org.Hs.eg.db, keyType="ENTREZID") %>% as.data.frame()
    write.csv(KEGG_res_DOWN,file=paste(deparse(substitute(MARKERS)),OUTPUT,"_DOWN_KEGG_RESULTS_",".csv",sep=""))
  }
}









Cell_types <- c("Fibroblast-1","Fibroblast-2","HairFollicle","Keratinocyte-1","Keratinocyte-2","Keratinocyte-3","Keratinocyte-4","Keratinocyte-5","Sebocyte","VSMC-1","Melanocyte","Lymphatic")
Cell_types <- c("T cell-1","T cell-2","T cell-3","B cell","Myeloid-1","Myeloid-2","Langerhans")

st_enrichment <- function(sc.data,st.data,Cell_types,seed=100){
  set.seed(seed)
  st.data <- st_filter_by_genes(st.data = st.data,x = 200)
  st.data <- SCTransform(st.data,assay = "Spatial")
  st.data <- RunPCA(st.data, assay = "SCT", verbose = FALSE)
  
  data_reference <- SCTransform(sc.data, ncells = 3000, verbose = FALSE) %>% RunPCA(verbose = FALSE) %>% RunUMAP(dims = 1:20)
  anchors <- FindTransferAnchors(reference = data_reference, query = st.data, normalization.method = "SCT")
  predictions.assay <- TransferData(anchorset = anchors, refdata = data_reference$final_clustering, prediction.assay = TRUE, weight.reduction = st.data[["pca"]],dims = 1:20)
  st.data[["predictions_sc_data"]] <- predictions.assay
  
  DefaultAssay(st.data) <- "predictions_sc_data"
  image_name <- Images(st.data)
  Cell_types <- rownames(st.data) %>% as.vector()
  for (x in Cell_types){
    pdf(file = paste(image_name[1],x,".pdf"),width = 10,height = 15)
    print(SpatialFeaturePlot(st.data, features = c(x), pt.size.factor = 2.5, ncol = 2, crop = TRUE) + scale_fill_gradientn(colors=cols,limits = c(0,1.00001))) 
    dev.off()
  }
  
}
st_enrichment_v2 <- function(sc.data,st.data,Cell_types){
  st.data <- st_filter_by_genes(st.data = st.data,x = 200)
  st.data <- SCTransform(st.data,assay = "Spatial")
  st.data <- RunPCA(st.data, assay = "SCT", verbose = FALSE)
  
  data_reference <- SCTransform(sc.data, ncells = 3000, verbose = FALSE) %>% RunPCA(verbose = FALSE) %>% RunUMAP(dims = 1:20)
  anchors <- FindTransferAnchors(reference = data_reference, query = st.data, normalization.method = "SCT")
  predictions.assay <- TransferData(anchorset = anchors, refdata = data_reference$final_clustering, prediction.assay = TRUE, weight.reduction = st.data[["pca"]],dims = 1:20)
  st.data[["predictions_sc_data"]] <- predictions.assay
  
  DefaultAssay(st.data) <- "predictions_sc_data"
  return(st.data)
}

pseudo_bulk_out <- function(Seurat_obj,group_label,groups_tbl_path){
  pseudo.counts <- Seurat:::PseudobulkExpression(object = Seurat_obj,assays = "Spatial",group.by = group_label,slot = "counts",pb.method = "aggregate") 
  pseudo.counts.df <- as.data.frame(pseudo.counts$Spatial) %>% rownames_to_column("Gene")
  groups.table <- read.csv(file=groups_tbl_path,stringsAsFactors = TRUE) %>% filter(Sample.ID!="") %>% filter(Sample.ID %in% colnames(pseudo.counts.df)) %>% column_to_rownames("Sample.ID")
  
  # Ordering the file
  counts.file <- pseudo.counts.df[,rownames(groups.table)]
  
  group.name <- colnames(groups.table)[1]
  groups.levels <- levels(groups.table[,group.name])
  counts.final <- lapply(counts.file, as.integer) %>% as.data.frame(row.names = pseudo.counts.df$Gene)
  dds <- DESeqDataSetFromMatrix(countData = counts.final,colData=groups.table,design = ~BATCH + GROUP_I)
  dds <- DESeq(dds,quiet = TRUE)
  return(dds)
  #pdf(file=paste("PCA",,"PLOT_PSEUDOBULK.pdf"),height = 8,width = 10)
  #print(plotPCA(vsd, intgroup=c("GROUP_I", "SEVERITY")))
  #dev.off()
  
  #sampleDists <- dist(t(assay(vsd)))
  
  #sampleDistMatrix <- as.matrix(sampleDists)
  #rownames(sampleDistMatrix) <- paste(vsd$GROUP_I, vsd$SEVERITY,vsd$GROUP_II, sep="-")
  #colnames(sampleDistMatrix) <- NULL
  #colors <- colorRampPalette(brewer.pal(9, "RdYlBu")) (255)
  #pdf("HC_WITH_HEATMAP.pdf",height = 15,width = 15)
  #print(pheatmap(sampleDistMatrix,
  #               clustering_distance_rows=sampleDists,
  #               clustering_distance_cols=sampleDists,
  #               col=colors))
  #dev.off()
  
  #sampleDists <- dist(t(assay(vsd)))
  
  #hc <- hclust(sampleDists)
  
  #pdf("DENDOGRAM_ALL_SAMPLES.pdf",height = 10,width = 10)
  #print(plot(hc, labels=paste(vsd$GROUP_I,vsd$SEVERITY,vsd$GROUP_II,sep = ":")))
  #dev.off()
}

## MULTI-MODAL INTERSECTION ANALYSIS

MIA_ENRICH <- function(stlist,sclist,total){
  overlap <- length(intersect(stlist$gene,sclist$gene))
  C <- length(sclist$gene)
  D <- length(stlist$gene)
  
  # ENRICHMENT CALCULATION
  e <- -log10(phyper(overlap,C, total - C, D, lower.tail = FALSE))
  return(e)
}
MIA_DEPLETE <- function(stlist,sclist){
  overlap <- length(intersect(stlist$gene,sclist$gene))
  C <- length(sclist$gene)
  D <- length(stlist$gene)
  
  # ENRICHMENT CALCULATION
  d <- -log10(1- (phyper(overlap,C, total - C, D, lower.tail = FALSE)))
  return(d)
}
MIA <- function(total_genes,single_cell.markers,spatial.markers)
{
  #D.SCORES <- c()
  # Perform this operation for every cell type
  #single_cell.markers <- FindAllMarkers(single_cell,assay = assay_use,logfc.threshold = 0.25,return.thresh = p_val_adj < 0.1)
  #spatial.markers <- FindAllMarkers(spatial_data,assay = assay_use,logfc.threshold = 0.25,return.thresh = p_val_adj < 0.1)
  cell.types <- single_cell.markers %>% dplyr::select(cluster) %>% unique() %>% as.list()
  spatial.regions <- spatial.markers %>% dplyr::select(cluster) %>% unique() %>% as.list()
  E.SCORES <- data.frame(spatial.regions)
  for(i in cell.types){
    for (x in i){
      e_list <- c()
      #list.append(e_list,i)
      for(y in spatial.regions){
        for(z in y){
          single_cell <- single_cell.markers %>% filter(cluster==x)
          spatial_data <- spatial.markers %>% filter(cluster==z)
          e <- MIA_ENRICH(single_cell,spatial_data,total = total_genes)
          #d <- MIA_DEPLETE(single_cell,spatial_data,total = total_genes)
          e_list <- c(e_list,e)
        }
        #D.SCORES <- append(D.SCORES,d)
      }
      E.SCORES[paste(x)] <- e_list
      #E.SCORES <- append(E.SCORES,e)
    }
  }
  #e.data <- data.frame("GA"=E.SCORES[1],"ER"=E.SCORES[2],"C0L17A1+"=E.SCORES[3])
  #d.data <- data.frame("GA"=D.SCORES[1],"ER"=D.SCORES[2],"C0L17A1+"=D.SCORES[3])
  #res <- list(E.SCORES,e.data)
  return(E.SCORES)
}

