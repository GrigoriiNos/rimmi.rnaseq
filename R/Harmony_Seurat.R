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

harmony_for_seurat <- function(merged_object, pc.genes = 'default'){
  library(harmony)

  merged_object <- Seurat::FindVariableGenes(object = merged_object,
                           mean.function = ExpMean,
                           dispersion.function = LogVMR,
                           x.low.cutoff = 0.0125,
                           x.high.cutoff = 3,
                           y.cutoff = 0.5)

  if (pc.genes == 'default'){
    pc.gen <-  merged_object@var.genes
  }
  if (pc.genes == 'no.ig'){
    pc.gen <- merged_object@var.genes[-(grep(pattern = "^IG[H,L,K][A-Z][0-9].*", x = rownames(x = B03_A07@data)))]
  }
  if (pc.genes == 'no.tcr'){
    pc.gen <- merged_object@var.genes[-(grep(pattern = "^TR[A,B][D,J,V].*", x = rownames(x = B03_A07@data), value = TRUE))]
  }

  merged_object <- RunPCA(object = merged_object,
                pc.genes = pc.gen,
                pcs.compute = 50,
                do.print = F, pcs.print = 1:5,
                genes.print = 5)

  samples <- unique(merged_object@meta.data$orig.ident)

  pca <- merged_object@dr$pca@cell.embeddings

  harmony_emb <- HarmonyMatrix(pca,
                               merged_object@meta.data$orig.ident,
                               theta=2, do_pca=FALSE)

  merged_object@dr$harmony <- harmony_emb

  merged_object
}
