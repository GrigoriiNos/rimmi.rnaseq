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

scatter_libdepth <- function(Seurat_obj){
  library(ggplot2)
  library(cowplot)

  df <- data.frame(umap1 = Seurat_obj@dr$umap@cell.embeddings[,1 ],
                   umap2 = Seurat_obj@dr$umap@cell.embeddings[,2 ],
                   cells = as.character(Seurat_obj@ident),
                   lsize = Seurat_obj@meta.data$nUMI,
                   ngene = Seurat_obj@meta.data$nGene
                   )

  p1 <- ggplot(df, aes(x = umap1, y = umap2, col = cells))+
          geom_point(aes(size = lsize)) +
          scale_size_continuous(range = c(0.001, 3))
  p2 <- ggplot(df, aes(x = umap1, y = umap2, col = cells))+
          geom_point(aes(size = ngene)) +
          scale_size_continuous(range = c(0.001, 3))
  plot_grid(p1, p2)
}
