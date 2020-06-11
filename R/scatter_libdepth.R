#' Plot clusters in 2 umaps with the point size corresponting to the library size
#'
#' This function allows you to see if the library size influences in your clustering
#'
#' @param Seurat_obj your Seurat object
#'
#' @return scatter plot with the library size correponding to the point size
#'
#' @keywords Seurat, single cell sequencing, RNA-seq, gene signature
#'
#' @examples
#'
#' scatter_libdepth(A07)
#'
#' @export
#'

lib.qc_plot <- function(Seurat_obj){
  library(ggplot2)

  qc.df <- data.frame(row.names = rownames(Seurat_obj@meta.data),
                      umap1 = Seurat_obj@reductions$umap@cell.embeddings[,1],
                      umap2 = Seurat_obj@reductions$umap@cell.embeddings[,2],
                      nUMI = Seurat_obj$nCount_RNA,
                      pMito = Seurat_obj$percent.mt,
                      nGenes = Seurat_obj$nFeature_RNA,
                      cluster = Idents(Seurat_obj))

  p1 <- ggplot(qc.df, aes(x = umap1, y = umap2, colour = cluster)) +
    geom_point(aes(size = nUMI)) +
    theme_bw()


  p2 <- ggplot(qc.df, aes(x = umap1, y = umap2, colour = cluster)) +
    geom_point(aes(size = nGenes)) +
    theme_bw()

  p3 <- ggplot(qc.df, aes(x = umap1, y = umap2, colour = cluster)) +
    geom_point(aes(size = pMito)) +
    theme_bw()

  cowplot::plot_grid(p1, p2, p3)
}
