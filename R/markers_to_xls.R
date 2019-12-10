#' getting your markers into the excel shpreadsheet
#'
#' This function allows you write Seurat::FindAllMarkers function output into the excel file having each cluster markers written in the distinct sheet
#'
#' @param markers_table The output of the Seurat::FndAllMarkers function, don't change anything in the column names within it.
#'
#' @param filename
#'
#' @keywords excel, single cell, RNA-seq, Seurat, markers
#'
#' @examples
#'
#' annotate_markers(my_markers, "stromal_markers")
#'
#' @export
#'

markers_to_xls <- function(markers_table, filename = 'my_markers'){
  library(dplyr)
  ###
  clusters <- unique(markers_table$cluster)
  tables <- list()
  for (i in 1:length(clusters)){
    sub_markers_table <- markers_table %>%
      filter(cluster == clusters[i])
    tables[[i]] <- sub_markers_table
  }
  WriteXLS::WriteXLS(x = tables,
                     ExcelFileName = paste0(filename, '.xls'),
                     SheetNames = paste(clusters, 'cluster')
  )
}
