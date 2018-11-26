#' Plot uMAP in 2 dimentions
#'
#' This function allows you in addition to plot your clustering results in 2 uMAPs dimentions
#' @param Seurat_obj your Seurat object
#' @param title title for your plot
#' @keywords Seurat, single cell sequencing, RNA-seq, uMAP, dimentionality reduction
#' @export
#' @examples
#' Seurat_obj <- Run_uMAP(Seurat_obj)
#' uMAPPlot(Seurat_obj, sample_name = 'b06')


uMAPPlot <- function(Seurat_obj, sample_name = 'X'){
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  #get the cells names
  Seurat_obj_labels <- as.character(
    unname(Seurat_obj@ident)
  )

  #select the colours
  cols<-brewer.pal(n=length(
    unique(Seurat_obj_labels)
  ),
  name="Set3")

  names(cols) <- unique(Seurat_obj_labels)
  for (i in 1:length(cols)){
    Seurat_obj_labels[Seurat_obj_labels == names(cols)[i]] <- cols[i]
  }
  # plot it!
  to_plot <- data_frame(umap1 = Seurat_obj@dr$uMAP$layout[,1],
                       umap2 = Seurat_obj@dr$uMAP$layout[,2],
                   cells = as.character(
                     unname(Seurat_obj@ident)
                   ))

  ggplot(to_plot, aes(x = umap1, y = umap2, col = cells))+
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



