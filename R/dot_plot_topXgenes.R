#' Make a dot plot for all your markers
#'
#' This function allows you to plot each top X DEG for each cell cluster using a Seurat::FindAllMarkers function output
#'
#' @param markers_table The output of the Seurat::FndAllMarkers function, don't change anything in the column names within it.
#'
#' @param X How much of the genes you need to analyse in each cluster. 5 is a default.
#'
#' @param Seurat_obj your Seurat object to plot markers from
#'
#' @keywords dot plot, Seurat, markers, single cell sequencing, RNA-seq
#'
#' @examples
#'
#' dot_plot_topXgenes(markers, 10)
#'
#' @export
#'


dot_plot_topXgenes <- function(markers_table, X = 5, Seurat_obj){
  library(Seurat)
  library(dplyr)

  ###
  plot_top_X_genes <- function(markers_table, X = 5){
    top_X <- markers_table %>%
      group_by(cluster) %>%
      top_n(X, avg_logFC)
    top_X$gene
  }
  ###
  top_5_B06 <- unique(
    plot_top_X_genes(markers_table, X)
  )
  DotPlot(
    b06,
    Seurat_obj,
    x.lab.rot = T
  )
}
