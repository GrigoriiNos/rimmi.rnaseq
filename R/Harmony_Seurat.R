#' use the batch effect correction method harmony on Seurat objects
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

harmony_for_seurat <- function(merged_object,
                               pc.genes = 'default',
                               pcs = 50,
                               theta = 2,
                               lambda = 1,
                               sigma = 0.1){
  library(harmony)

  DefaultAssay(merged_object) <- "RNA"

  merged_object <- Seurat::FindVariableFeatures(object = merged_object,
                                                selection.method = "vst",
                                                nfeatures = 2000,
                                                verbose = F)

  if (pc.genes == 'default'){
    pc.gen <- VariableFeatures(merged_object)
  }
  if (pc.genes == 'no.ig'){
    pc.gen <- VariableFeatures(merged_object)[!grepl(pattern = "^IG[H,L,K][A-Z][0-9].*", x = VariableFeatures(merged_object))]
  }
  if (pc.genes == 'no.tcr'){
    pc.gen <- VariableFeatures(merged_object)[!grepl(pattern = "^TR[A,B][D,J,V].*", x = VariableFeatures(merged_object))]
  }

  merged_object <- merged_object %>%
    ScaleData() %>%
    RunPCA(pc.genes = pc.gen,
           pcs.compute = pcs) %>%
    RunTSNE()

  samples <- unique(merged_object$orig.ident)

  pca <- merged_object@reductions$pca@cell.embeddings

  harmony_emb <- HarmonyMatrix(pca,
                               merged_object@meta.data$orig.ident,
                               theta=theta,
                               lambda = lambda,
                               sigma = sigma,
                               do_pca=FALSE)

  merged_object@reductions$tsne@cell.embeddings <- harmony_emb

  merged_object
}


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

scanorama_seurat <- function(merged_object, raw.data = T){
  samples <- unique(merged_object@meta.data$orig.ident)

  if (raw.data == F){
    batch1 <- merged_object@data[, merged_object@meta.data$orig.ident == samples[1]]
    cells1 <- colnames(batch1)

    batch2 <- merged_object@data[, merged_object@meta.data$orig.ident == samples[2]]
    cells2 <- colnames(batch2)
  } else {
    batch1 <- merged_object@raw.data[, merged_object@meta.data$orig.ident == samples[1]]
    cells1 <- colnames(batch1)

    batch2 <- merged_object@raw.data[, merged_object@meta.data$orig.ident == samples[2]]
    cells2 <- colnames(batch2)
  }

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
