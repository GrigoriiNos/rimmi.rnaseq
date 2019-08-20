#' Plot your clusters in 3 dimensions
#'
#' This function allows you to make a 3D scatter plot for your single cell gene expression measurements using UMAP dimensions by default from your Seurat object
#'
#' @param obj your Seurat object
#' @param dims how many of pca/cca dimentions you want to use, 10 is default
#'
#' @return plot
#'
#' @keywords Seurat, single cell sequencing, RNA-seq, uMAP, dimentionality reduction
#'
#' @examples
#'
#'
#' umap_plot3d(tcells)
#'
#' @export
#'

umap_plot3d <- function(obj, dims = 1:10){
  library(plotly)
  library(Seurat)
  # recalculate umaps with 3 dimensions
  if (is.null(obj@dr$cca.aligned)){
    obj <- RunUMAP(obj,
                   reduction.use = "pca",
                   dims.use = dims,
                   max.dim = 3)
  } else {
    obj <- RunUMAP(obj,
                   reduction.use = "cca.aligned",
                   dims.use = dims,
                   max.dim = 3)
  }
  # construct new data frame with umaps and cell identities
  df <- data.frame(umap1 = obj@dr$umap@cell.embeddings[,1],
                   umap2 = obj@dr$umap@cell.embeddings[,2],
                   umap3 = obj@dr$umap@cell.embeddings[,3],
                   cell = obj@ident)

  # plot the fancy 3d scatter plot ; D
  plot_ly(df, x = ~umap1, y = ~umap2, z = ~umap3, color = ~cell) %>%
    add_markers() %>%
    layout(scene = list(xaxis = list(title = 'umap 1'),
                        yaxis = list(title = 'umap 2'),
                        zaxis = list(title = 'umap 3')))
}
