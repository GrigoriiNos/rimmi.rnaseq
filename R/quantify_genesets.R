#'
#' This function allows you to calculate an average expression for the list of genes of your interest.
#'
#' @param Seurat_obj Seurat: your Seurat object
#' @param genesets, list: named list with genesets
#' @param reset_genesets, bool: if you want reset already existing genesets with a new assay object
#'
#' @return Seurat object with geneset assay
#'
#' @keywords geneset, gene signature,
#'
#' @examples
#'
#' @export
#'

quantify_genesets <- function(Seurat_obj, genesets, reset_genesets=FALSE){

  if (is.null((names(genesets)))){
    throw('genesets must be a names list')
  }

  genesets <- lapply(genesets, function(gset) gset[gset %in% rownames(Seurat_obj@assays$RNA@counts)])

  DefaultAssay(Seurat_obj) <- 'RNA'
  signs <- AddModuleScore(Seurat_obj, genesets)
  signs_matrix <- t(sign@meta.data[grepl("Cluster" ,colnames(signs@meta.data))])
  rownames(signs_matrix) <- names(genesets)

  if (reset_genesets){
    Seurat_obj@assays$genesets <- NULL
  }

  if (is.null(Seurat_obj@assays$genesets)){
    Seurat_obj_signs <- CreateAssayObject(data = as.matrix(signs_matrix))
    Seurat_obj@assays$genesets <- Seurat_obj_signs
    DefaultAssay(Seurat_obj) <- 'genesets'
    Seurat_obj <- ScaleData(Seurat_obj)
    return(Seurat_obj)
  } else {
    signs_matrix <- rbind(Seurat_obj@assays$genesets@data, signs_matrix)
    signs_matrix <- signs_matrix[unique(rownames(signs_matrix)),]
    Seurat_obj@assays$genesets <- CreateAssayObject(data = signs_matrix)
    DefaultAssay(Seurat_obj) <- 'genesets'
    Seurat_obj <- ScaleData(Seurat_obj)
    return(Seurat_obj)
  }

}
