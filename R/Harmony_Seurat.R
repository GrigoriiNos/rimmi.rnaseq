#' use the batch effect correction method harmony on Seurat objects
#'
#' This function allows you to use an additional BE correction method still using Seurat objects
#'
#' @param merged_object pre-merged Seurat object
#'
#' @return Seurat object with BE corrected expression matrix in the data slot
#'
#' @keywords Seurat, single cell sequencing, RNA-seq, gene signature
#'
#' @examples
#'
#' mnn_for_seurat(B04_B05)
#' B05@meta.data[, "sample"] <- "B05"
#' B04@meta.data[, "sample"] <- "B04"
#'
#' add prefix to cellnames
#' B05 <- RenameCells(B05,
#'                   add.cell.id = 'B05')
#'
#' B04 <- RenameCells(B04,
#'                   add.cell.id = 'B04')
#'
#' B04_B05 <- MergeSeurat(B04, B05)
#' B04_B05 <- mnn_for_seurat(B04_B05)
#' @export
#'

harmony_for_seurat <- function(merged_object){
  library(harmony)
  
  samples <- unique(merged_object@meta.data$orig.ident)
  
  pca <- merged_object@dr$pca@cell.embeddings
  
  harmony_emb <- HarmonyMatrix(pca, 
                               merged_object@meta.data$orig.ident, 
                               theta=2, do_pca=FALSE)
  
  merged_object@dr$harmony <- harmony_emb
  
  merged_object
}
