#' use the batch effect correction method scanorama on Seurat objects
#'
#' This function allows you to use an additional BE correction method still using Seurat objects
#'
#' @param merged_object pre-merged Seurat object
#'
#' @return Seurat object with BE corrected expression matrix in the data slot
#'
#' @keywords Seurat, single cell sequencing, RNA-seq, gene signature
#'
#' @examples
#'
#' mnn_for_seurat(B04_B05)
#' B05@meta.data[, "sample"] <- "B05"
#' B04@meta.data[, "sample"] <- "B04"
#'
#' add prefix to cellnames
#' B05 <- RenameCells(B05,
#'                   add.cell.id = 'B05')
#'
#' B04 <- RenameCells(B04,
#'                   add.cell.id = 'B04')
#'
#' B04_B05 <- MergeSeurat(B04, B05)
#'
#' @export
#'
#'

scanorama_seurat <- function(merged_object){
  samples <- unique(merged_object@meta.data$orig.ident)

  batch1 <- merged_object@data[, merged_object@meta.data$orig.ident == samples[1]]
  cells1 <- colnames(batch1)

  batch2 <- merged_object@data[, merged_object@meta.data$orig.ident == samples[2]]
  cells2 <- colnames(batch2)

  # List of data sets (matrices of cells-by-genes):
  datasets <- list(t(as.matrix(batch1)),
                   t(as.matrix(batch2)))
  # List of gene lists:
  genes_list <- list(rownames(batch1),
                     rownames(batch2))

  library(reticulate)
  scanorama <- import('scanorama')

  # Batch correction.
  corrected.data <- scanorama$correct(datasets,
                                      genes_list,
                                      return_dense=TRUE)

  batch1 <- corrected.data[[1]][[1]]
  rownames(batch1) <- cells1
  colnames(batch1) <- as.character(corrected.data[[2]])

  batch2 <- corrected.data[[1]][[2]]
  rownames(batch2) <- cells2
  colnames(batch2) <- as.character(corrected.data[[2]])

  corrected <- t(
    rbind(batch1,
          batch2)
  )

  colnames(corrected) <- c(cells1, cells2)

  merged_object@data <- as(corrected, 'dgCMatrix')
  merged_object
}
