#' Run Velocyto analysis on your Seurat2 object
#'
#' This function allows you to un Velocyto analysis on your Seurat2 object and visualise it on your umap embeddings,
#' you need to have fully pre-proccessed Seurat object with QC done, 2d embeddings and clustering pre-calculated,
#'
#'
#' @param Seurat_obj  Seurat object
#' @param loom_path path to the loom file, pre-calculated with python command line tool
#'
#' @return umap 2d plot with velocity
#'
#' @keywords Seurat, single cell sequencing, RNA-seq, RNA velocity
#'
#' @examples
#'
#'
#'
#' @export
#' 

Seurat2_monocle3 <- function(Seurat_obj, 
                             clus_var = c('res.0.6', 'res.0.8', 'res.1', 'celltype'), 
                             multi_samples = T,
                             samples_var = c('condition', 'sample'),
                             cell_cycle = F)
  {
  library(monocle3)

  features_select <- grepl("^MT[-,{RNR}]|^IG[H,L,K][A-Z][0-9].*|^RP[L,S][0-9].*|^RP[0-9].*|^FO[0-9]{2,}|^AP[0-9]{2,}|\\.", 
                           toupper(rownames(Seurat_obj@data))
                           )

  raw.matr <- Seurat_obj@raw.data[,colnames(Seurat_obj@raw.data) %in% rownames(Seurat_obj@meta.data)]
  raw.matr <- raw.matr[!(features_select),]

  gene_meta <- data.frame(gene_short_name = rownames(raw.matr))
  rownames(gene_meta) <- rownames(raw.matr)

  Monocle_obj <- new_cell_data_set(expression_data = raw.matr,
                                  cell_metadata = Seurat_obj@meta.data[colnames(raw.matr),], 
                                  gene_metadata = gene_meta)

  # pre-processing stage 
  Monocle_obj <- preprocess_cds(Monocle_obj, 
                               num_dim = 50) 

  if (cell_cycle == F & multi_samples == T){
  Monocle_obj <- align_cds(Monocle_obj, 
                          alignment_group = samples_var, 
                          residual_model_formula_str = paste0("~ ", samplesvar)
                          )
  } 
  if (cell_cycle == T & multi_samples == T){
    Monocle_obj <- align_cds(Monocle_obj, 
                             alignment_group = samples_var, 
                             residual_model_formula_str = '~ Phase + S.Score + G2M.Score'
    )
  }
  if (cell_cycle == T & multi_samples == F){
    Monocle_obj <- align_cds(Monocle_obj, 
                             alignment_group = 'Phase', 
                             residual_model_formula_str = '~ S.Score + G2M.Score'
    )
  }

  # calculate umap
  Monocle_obj <- reduce_dimension(Monocle_obj)

  monoclePCA <- list()

  monoclePCA$gene.loadings <- Monocle_obj@preprocess_aux$gene_loadings
  monoclePCA$cell.loadings <- Monocle_obj@reducedDims$PCA

  Seurat_obj@misc$monocle <- Monocle_obj@reducedDims$UMAP 
  Seurat_obj@misc$monoclePCA <- monoclePCA

  # monocle clustering
  Monocle_obj <- cluster_cells(Monocle_obj)

  # learn trajectory
  Monocle_obj <- learn_graph(Monocle_obj)
  
  if (cell_cycle == F & multi_samples == T){
    p1 <- plot_cells(Monocle_obj,
                     color_cells_by = clus_var,
                     label_groups_by_cluster=FALSE,
                     label_leaves=FALSE,
                     label_branch_points=FALSE, cell_size = 1, 
                     group_label_size = 6)
    
    
    p2 <- plot_cells(Monocle_obj,
                     color_cells_by = samples_var,
                     label_groups_by_cluster=FALSE,
                     label_leaves=FALSE,
                     label_branch_points=FALSE, cell_size = 1, 
                     group_label_size = 6)
    
    save.plot_grid(p1,p2, 
                   filename = filename)
    
  } 
  if (cell_cycle == T & multi_samples == T){
    p1 <- plot_cells(Monocle_obj,
                     color_cells_by = clus_var,
                     label_groups_by_cluster=FALSE,
                     label_leaves=FALSE,
                     label_branch_points=FALSE, cell_size = 1, 
                     group_label_size = 6)
    
    
    p2 <- plot_cells(Monocle_obj,
                     color_cells_by = samples_var,
                     label_groups_by_cluster=FALSE,
                     label_leaves=FALSE,
                     label_branch_points=FALSE, cell_size = 1, 
                     group_label_size = 6)
    
    p3 <- plot_cells(Monocle_obj,
                     color_cells_by = 'Phase',
                     label_groups_by_cluster=FALSE,
                     label_leaves=FALSE,
                     label_branch_points=FALSE, cell_size = 1, 
                     group_label_size = 6)
    
    save.plot_grid(p1,p2,p3, 
                   filename = filename)
  }
  if (cell_cycle == T & multi_samples == F){
    p1 <- plot_cells(Monocle_obj,
                     color_cells_by = clus_var,
                     label_groups_by_cluster=FALSE,
                     label_leaves=FALSE,
                     label_branch_points=FALSE, cell_size = 1, 
                     group_label_size = 6)
    p3 <- plot_cells(Monocle_obj,
                     color_cells_by = 'Phase',
                     label_groups_by_cluster=FALSE,
                     label_leaves=FALSE,
                     label_branch_points=FALSE, cell_size = 1, 
                     group_label_size = 6)
    
    save.plot_grid(p1,p3, 
                   filename = filename)
  }
  
  return(list(
    Seurat_obj,
    Monocle_obj
    )
      )
    
}