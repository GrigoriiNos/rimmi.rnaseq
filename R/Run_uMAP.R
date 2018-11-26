#' Calculate uMAP for Seurat object
#'
#' This function allows you in addition to classical dimentionality reductions approaches as tSNE and PCA calculate also uMAP and embedd it to your Seurat object in the dr slot.
#'
#' @param Seurat_obj your Seurat object
#'
#' @return Seurat object with new dimentionality reduction slot
#'
#' @keywords Seurat, single cell sequencing, RNA-seq, uMAP, dimentionality reduction
#'
#' @examples
#'
#' Seurat_obj <- Run_uMAP(Seurat_obj)
#' uMAPPlot(Seurat_obj, sample_name = 'b06')
#'
#' @export
#'

Run_uMAP <- function(Seurat_obj, variable_genes = TRUE){
  library(umap)
  library(Seurat)
  ### get expression matrix
  if (variable_genes == T) {
    exp.matrix <- Seurat_obj@data[Seurat_obj@var.genes,]
  } else {
    exp.matrix <- Seurat_obj@data
  }

  print('started to calculate uMAPs, can take up to 1 hour')
  ### calculate uMAPs
  Seurat_obj_map <- umap(
    t(as.matrix(
      exp.matrix)
    )
  )

  ### put umap into the new slot
  Seurat_obj@dr$umap <- Seurat_obj_map

  print('done')
  print('you now have your Seurat object with a new dimentionality rediction slot')
  Seurat_obj
}
