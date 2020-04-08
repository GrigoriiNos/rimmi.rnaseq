#' Annotation for the cell markers.
#'
#' This function allows you to get multiple functional annotations for DE genes for each cell cluster using a Seurat::FindAllMarkers function output
#'
#' @param markers_table The output of the Seurat::FndAllMarkers function, don't change anything in the column names within it.
#'
#' @param X How much of the genes you need to analyse in each cluster. 50 is a default.
#'
#' @param organism What species you work with? "mmusculus" and "hsapiens" is available in gProfiler. Human is a default.
#'
#' @return table with annotation for each cluster
#'
#' @keywords GO, KEGG, HPA, single cell, RNA-seq
#'
#' @examples
#'
#' annotate_markers(stromal_markers)
#'
#' @export
#'

annotate_markers <- function(markers_table,
                             X = 50,
                             organism = 'hsapiens',
                             reg = 'up'){
  library(gProfileR)
  library(dplyr)

  ###
  if (reg == 'up'){
    q <- markers_table %>%
      group_by(cluster) %>%
      filter(p_val_adj < 0.05 & avg_logFC > 0) %>%
      top_n(X, abs(avg_logFC)) %>%
      dplyr::select(gene) %>%
      do(gprofiler(., organism = organism))
  }
  if (reg == 'down'){
    q <- markers_table %>%
      group_by(cluster) %>%
      filter(p_val_adj < 0.05 & avg_logFC < 0) %>%
      top_n(X, abs(avg_logFC)) %>%
      dplyr::select(gene) %>%
      do(gprofiler(., organism = organism))
  }
  return(q)
}
