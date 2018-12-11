#' Plot uMAP in 2 dimentions
#'
#' This function allows you in addition to plot your clustering results in 2 uMAPs dimentions
#' @param Seurat_obj your Seurat object
#'
#' @param title title for your plot
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



