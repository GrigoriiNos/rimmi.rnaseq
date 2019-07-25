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

  counts <- as.data.frame(
    as.matrix(
      Seurat_obj@data)
    )

  #colnames(counts) <- paste('d-pos_', colnames(counts), sep = '')

  library("biomaRt")

  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  genes  <-  getBM(filters='hgnc_symbol',
            attributes = c('ensembl_gene_id','hgnc_symbol'),
            values = rownames(counts),
            mart = ensembl)

  counts <- counts[rownames(counts) %in% genes$hgnc_symbol,]

  counts <- tibble::rownames_to_column(
    as.data.frame(counts), var = 'hgnc_symbol')

  counts <- plyr::join(counts, genes)

  counts$hgnc_symbol <- NULL

  counts <- cbind(counts[,which(colnames(counts) == 'ensembl_gene_id')], counts)

  colnames(counts)[1] <- 'Gene'
  counts$ensembl_gene_id <- NULL

  metadata <- data.frame(Cell = rownames(Seurat_obj@meta.data),
                         cell_type = Seurat_obj@ident
                         )

  #metadata$Cell <- paste('d-pos_', metadata$Cell, sep = '')

  write.table(counts,
            file = 'counts.txt',
            quote = F,
            col.names = T,
            row.names = F,
            sep = '\t')

  write.table(metadata,
            file = 'metadata.txt',
            quote = F,
            col.names = T,
            row.names = F,
            sep = '\t')

  system('cellphonedb method statistical_analysis metadata.txt counts.txt --iterations=10 --threads=2')

  system('cellphonedb plot dot_plot')

  system('cellphonedb plot heatmap_plot metadata.txt')
}












