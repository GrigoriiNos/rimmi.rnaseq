#' use the batch effect correction method mnn correct on Seurat objects
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
#' B04_B05 <- mnn_for_seurat(B04_B05)
#' @export
#'

mnn_for_seurat <- function(merged_object){
  samples <- unique(merged_object@meta.data$orig.ident)

  batch1 <- merged_object@data[, merged_object@meta.data$orig.ident == samples[1]]
  batch2 <- merged_object@data[, merged_object@meta.data$orig.ident == samples[2]]

  be_corrected <- scran::mnnCorrect(batch1,
                                    batch2,
                                    k = 20,
                                    sigma = 0.1,
                                    cos.norm.in = TRUE,
                                    svd.dim = 2)

  batch11 <- be_corrected$corrected[[1]]
  batch22 <- be_corrected$corrected[[2]]

  colnames(batch11) <- colnames(batch1)
  colnames(batch22) <- colnames(batch2)

  corrected <- cbind(batch11,
                     batch22
                     )

  merged_object@raw.data <- corrected

  merged_object
}
