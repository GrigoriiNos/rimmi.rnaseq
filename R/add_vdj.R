#' Merge clonotype and isotype data to GEX Seurat object
#'
#' This function allows you to fuse the output of VDJ library to Seurat object
#'
#' @param vdj_location location of cel ranger vdj output
#'
#' @return table with vdj info
#'
#' @keywords Seurat, single cell sequencing, RNA-seq, RNA velocity
#'
#' @examples
#'
#' A07 <- add_vdj('path_to_vdj', A07)
#'
#'
#' @export

pull_vdj <- function(vdj_location){
  vdj <- read.csv(paste(vdj_location,"filtered_contig_annotations.csv", sep=""))

  # Remove the -1 at the end of each barcode.
  # Subsets so only the first line of each barcode is kept,
  # as each entry for given barcode will have same clonotype.
  vdj$barcode <- gsub("-1", "", vdj$barcode)
  vdj <- vdj[!duplicated(vdj$barcode), ]

  # Only keep the barcode, isotype and clonotype columns.
  # We'll get additional clonotype info from the clonotype table.
  vdj <- vdj[,c("barcode", "raw_clone_id", "ccall", "vcall", "dcall", "jcall")]
  names(vdj)[names(vdj) == "raw_clone_id"] <- "clone_id"

  # Clonotype/isotype-centric info
  clono <- read.csv(paste(vdj_location,"clonotypes.csv", sep=""))

  # Slap the AA sequences onto our original table by clone_id.
  vdj <- merge(vdj, clono[, c("clone_id", "cdr3s_aa", "cdr3s_nt")])

  # Reorder so barcodes are first column and set them as rownames.
  vdj <- vdj[, c("barcode", "clone_id", "ccall", "vcall", "dcall", "jcall", "cdr3s_aa", "cdr3s_nt")]
  rownames(vdj) <- vdj[,'barcode']
  vdj$barcode <- NULL
  vdj
}


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
  vdj <- pull_vdj(vdj_location)

  # Add to the Seurat object's metadata.
  Seurat_obj@meta.data <- merge(Seurat_obj@meta.data,
                                vdj,
                                by = 0,
                                all.x = T)
  rownames(Seurat_obj@meta.data) <- Seurat_obj@meta.data$Row.names
  Seurat_obj@meta.data$Row.names <- NULL

  return(Seurat_obj)
}


#' This function allows you pull a summary information about your clonotypes
#'
#' @param Seurat_obj Seurat object with integrated data from VDJ analysis with Isotype assigment
#' @param n how much clones
#'
#' @keywords DEG, single cell, VDJ analysis, B cells, T cells
#'
#' @examples
#'
#' clonotypes_summary(A07, n = 10)
#'
#' @export
#'

clonotypes_summary <- function(Seurat_obj, n = 10, separate_cdr3 = T){

  # take top n clonotypes
  clonotypes <- table(Seurat_obj@meta.data$clone_id)

  topnclones <- clonotypes[
    order(clonotypes, decreasing = T)
    ][1:n]

  topnclones <- names(
    topnclones
  )

  # pull unique info about clonotypes from meta data
  summary <- Seurat_obj@meta.data %>%
    filter(clone_id %in% topnclones) %>%
    group_by(clone_id, c_call, v_call, j_call, junction, junction_aa) %>%
    summarise(n = n()) %>%
    #tidyr::separate_rows(cdr3s_aa, cdr3s_nt, sep = ';', convert = T) %>%
    arrange(as.numeric(gsub('[a-z]','', clone_id)))

  write.csv(summary,
            file = 'clonotypes_summary.csv',
            quote = F,
            row.names = F)
  print('clonotypes summary table is written')

  print(summary)
}
