#' Merge clonotype and isotype data to GEX Seurat object
#'
#' This function allows you to fuse the output of VDJ library to Seurat object
#'
#' @param vdj_location location of cel ranger vdj output
#' @param Seurat_obj  Seurat object
#'
#' @return umap 2d plot with velocity
#'
#' @keywords Seurat, single cell sequencing, RNA-seq, RNA velocity
#'
#' @examples
#'
#' A07 <- add_vdj('path_to_vdj', A07)
#'
#'
#' @export

add_vdj <- function(vdj_location, Seurat_obj){
  vdj <- read.csv(paste(vdj_location,"filtered_contig_annotations.csv", sep=""))

  # Remove the -1 at the end of each barcode.
  # Subsets so only the first line of each barcode is kept,
  # as each entry for given barcode will have same clonotype.
  vdj$barcode <- gsub("-1", "", vdj$barcode)
  vdj <- vdj[!duplicated(vdj$barcode), ]

  # Only keep the barcode, isotype and clonotype columns.
  # We'll get additional clonotype info from the clonotype table.
  vdj <- vdj[,c("barcode", "raw_clonotype_id", "c_gene")]
  names(vdj)[names(vdj) == "raw_clonotype_id"] <- "clonotype_id"

  # Clonotype/isotype-centric info
  clono <- read.csv(paste(vdj_location,"clonotypes.csv", sep=""))

  # Slap the AA sequences onto our original table by clonotype_id.
  vdj <- merge(vdj, clono[, c("clonotype_id", "cdr3s_aa")])

  # Reorder so barcodes are first column and set them as rownames.
  vdj <- vdj[, c('barcode', 'clonotype_id', 'c_gene', 'cdr3s_aa')]
  rownames(vdj) <- vdj[,'barcode']
  vdj[,'bacrode'] <- NULL

  # Add to the Seurat object's metadata.
  clono_seurat <- AddMetaData(object=Seurat_obj, metadata=vdj)
  return(clono_seurat)
}
