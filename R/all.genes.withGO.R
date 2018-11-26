#' Get all genes with partiular GO
#'
#' This function allows you to get a list of genes from the biomaRt database that have a GO of your choice.
#' @param GO a GO term that you want to get all genes with
#' @keywords GO, all genes with GO
#' @import biomaRt
#' @export
#' @examples
#' all.genes.withGO("GO:0008217")

all.genes.withGO <- function(GO){
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations
  #gets gene symbol, transcript_id and go_id for all genes annotated with GO:0007507
  gene.data <- getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id', 'go_id'),
                     filters = 'go', values = GO, mart = ensembl)
  # list of genes assosiated with BP term
  unique(gene.data$hgnc_symbol)
}
