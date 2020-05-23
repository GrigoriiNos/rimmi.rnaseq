#' Plot uMAP in 2 dimentions
#'
#' This function allows you in addition to plot your clustering results in 2 uMAPs dimentions
#' @param Seurat_obj your Seurat object
#'
#' @param title title for your plot
#' @param sample_name sample name
#' @param col should be 'clusters' if you want to colour cells according to the clustering, change on 'samples' for colouring according to the batches.
#'
#' @return ggplot2 object, umap plot
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


uMAPPlot <- function(Seurat_obj,
                     sample_name = 'X',
                     col = 'clusters'){
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  # data frame with 2 dims of dr umap and cell labels
  to_plot <- data_frame(umap1 = Seurat_obj@dr$umap$layout[,1],
                        umap2 = Seurat_obj@dr$umap$layout[,2],
                        )

  if (col == 'clusters')
  {
  to_plot$cols <-  as.character(unname(
                        Seurat_obj@ident)
                        )
  } else {
    to_plot$cols <- Seurat_obj@meta.data$protocol
  }

  # plot it!
  ggplot(to_plot, aes(x = umap1, y = umap2, col = cols))+
    geom_point()+
    ggtitle(paste0('uMAP plot for ', sample_name))+
    guides(colour = guide_legend(override.aes = list(size=8)))+
    theme_bw()+
    theme(legend.position = 'top', legend.text = element_text(size = 8),
          legend.background = element_rect(color = "black",
                                           fill = "grey70", size = 0.01, linetype = "solid"),
          plot.title = element_text(size = 15, face = "bold", hjust=0.5)
    )
}


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

umap_plot3d <- function(Seurat_obj,
                        coloring = c('clusters', 'sample', 'tissue', 'condition'),
                        reduction = c('pca', 'harmony', 'mnn', 'cca.aligned')){
  library(plotly)
  library(Seurat)
  # recalculate umaps with 3 dimensions

  if (reduction == 'pca'){
    dims <- Seurat_obj@commands$RunUMAP.RNA.pca$dims
  } else if ((reduction == 'cca.aligned')){
    dims <- Seurat_obj@commands$RunUMAP.RNA.cca.aligned$dims
  } else if ((reduction == 'harmony')){
    dims <- Seurat_obj@commands$RunUMAP.RNA.harmony$dims
  } else if (reduction == 'mnn'){
    dims <- Seurat_obj@commands$RunUMAP.RNA.mnn$dims
  }

  Seurat_obj <- RunUMAP(Seurat_obj,
                 reduction = reduction,
                 dims = dims,
                 n.components = 3)

  if (coloring == 'clusters'){
    coloring <- Idents(Seurat_obj)
  } else if (coloring == 'sample'){
    coloring <- Seurat_obj$sample
  } else if (coloring == 'tissue'){
    coloring <- Seurat_obj$tissue
  } else if (coloring == 'condition'){
    coloring <- Seurat_obj$condition
  }
  # construct new data frame with umaps and cell identities
  df <- data.frame(umap1 = Seurat_obj@reductions$umap@cell.embeddings[,1],
                   umap2 = Seurat_obj@reductions$umap@cell.embeddings[,2],
                   umap3 = Seurat_obj@reductions$umap@cell.embeddings[,3],
                   cell = coloring)

  # plot the fancy 3d scatter plot ; D
  plot_ly(df, x = ~umap1, y = ~umap2, z = ~umap3, color = ~cell) %>%
    add_markers() %>%
    layout(scene = list(xaxis = list(title = 'umap 1'),
                        yaxis = list(title = 'umap 2'),
                        zaxis = list(title = 'umap 3')))
}



