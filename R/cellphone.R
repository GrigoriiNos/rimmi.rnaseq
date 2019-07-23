#' cell-cell interaction analysis using method CellPhoneDB
#'
#' This function allows you to use an CellPhoneDB on Seurat Object
#'
#' @param Seurat_obj Seurat Object
#'
#' @return output tables for cell-cell interction inference
#'
#' @keywords Seurat, single cell sequencing, RNA-seq, cell-cell interactions
#'
#' @examples
#'
#'
#'
#'
#' @export
#'

cellphone_for_seurat <- function(Seurat_obj){

  counts <- as.matrix(
    Seurat_obj@data)

  metadata <- data.frame(Cell = rownames(Seurat_obj@meta.data),
                         cell_type = Seurat_obj@meta.data$res.0.6
                         )

  write.table(couts,
            file = 'counts.txt',
            quote = F,
            col.names = T)

  write.table(metadata,
            file = 'metadata.txt',
            quote = F,
            col.names = T)

  system('cellphonedb method statistical_analysis metadata.txt counts.txt --iterations=10 --threads=2')

  system('cellphonedb plot dot_plot')

  system('cellphonedb plot heatmap_plot metadata.txt')
}
