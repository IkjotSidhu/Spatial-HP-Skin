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

MIA_ENRICH_bullk <- function(stlist,bulklist,total){
  overlap <- length(intersect(stlist$gene,bulklist))
  C <- length(bulklist)
  D <- length(stlist$gene)
  
  # ENRICHMENT CALCULATION
  e <- -log10(phyper(overlap,C, total - C, D, lower.tail = FALSE))
  return(e)
}
MIA_bulk <- function(total_genes,markers,spatial.markers,name)
{
  spatial.regions <- spatial.markers %>% dplyr::select(cluster) %>% unique() %>% as.list()
  E.SCORES <- data.frame(spatial.regions)
  e_list <- c()
      #list.append(e_list,i)
    for(y in spatial.regions){
      for(z in y){
        spatial_data <- spatial.markers %>% filter(cluster==z)
        e <- MIA_ENRICH_bullk(spatial_data,markers,total = total_genes)
        e_list <- c(e_list,e)
      }
    }
   E.SCORES[paste(name)] <- e_list
  return(E.SCORES)
}

n_cells = nrow(combo_tbl)
for (ac in integrated_clusters) {
  ac_cells = combo_tbl %>% filter(cluster_integrated == ac) %>% pull(cell)
  n_ac = length(ac_cells)
  for (bc in independent_clusters) {
    bc_cells = combo_tbl %>% filter(cluster_independent == bc) %>% pull(cell)
    n_bc = length(bc_cells)
    n_common = length(intersect(ac_cells, bc_cells))
    n_either = length(union(ac_cells, bc_cells))
    
    # calculate jaccard index
    jaccard_index = n_common / n_either
    # message(glue("{ac} : {bc} : {length(ac_cells)} : {length(bc_cells)} : {jaccard_index}"))
    jaccard_mat[ac, bc] = round(jaccard_index, 3)
    
    # perform hypergeometric test
    ph = 1.0
    if (n_common > 0) {
      ph = phyper(n_common - 1, n_ac, n_cells - n_ac, n_bc, lower.tail = FALSE, log.p = FALSE)
    }
    # round(-log10(ph), 1)
    ph_mat[ac, bc] = round(ph, 3)
    # perform fisher's exact test
    fe = fisher.test(
      matrix(c(n_common,
               n_ac - n_common,
               n_bc - n_common,
               n_cells - n_ac - n_bc + n_common),
             2, 2),
      alternative="greater"
    )$p.value
    fe_mat[ac, bc] = round(ph, 3)
    
  }
}