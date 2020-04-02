#' heatmap genes with certain function
#'
#' This function allows you to plot a heatmap of some genes with the certain function
#'
#' @param query string containing a regular expression to find a certain gene function
#' @param Seurat_obj Seurat object
#'
#' @return heatmap
#'
#' @keywords heatmap, genes, gene function
#'
#' @examples
#'
#' query_heatmap('chemokine', pp.comp)
#'
#' @export
#'

query_heatmap <- function(query, Seurat_obj, species = 'Mouse'){
  if (species == 'Mouse'){
    library(org.Mm.eg.db)
    db <- org.Mm.eg.db
  } else {
    library(org.Hs.eg.db)
    db <- org.Hs.eg.db
  }

  allkeys <- keys(db, keytype="GENENAME")

  get_genes <- select(db,
                      keys=allkeys[grep(query, allkeys)],
                      columns="SYMBOL",
                      keytype="GENENAME")

  genes <- unique(
      (get_genes$SYMBOL)
      )

  genes <- genes[genes %in% rownames(Seurat_obj@assays$RNA@scale.data)]

  DoHeatmap(Seurat_obj,
            features = genes,
            size = 3.5, )
}
