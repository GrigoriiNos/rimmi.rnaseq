#' Calculate the contribution of each gene for the given gene signature
#'
#' This function allows you to see how did each gene from the given transcriptional program contribute to the pre-calculated gene signature
#'
#' @param Seurat_obj your Seurat object
#' @param ident number of cluster you want to look at
#' @param signature the gene set of your interest
#'
#' @return heatmap with the gene contribution
#'
#' @keywords Seurat, single cell sequencing, RNA-seq, gene signature
#'
#' @examples
#'
#' Gene_contribution(A07, 7, blooddata[[1]])
#'
#' @export
#'


Gene_contribution <- function(Seurat_obj, ident, signature, signature_name){
  library(gplots)
  library(Seurat)

  divisors <- function(x){
    #  Vector of numberes to test against
    y <- seq_len(x)
    #  Modulo division. If remainder is 0 that number is a divisor of x so return it
    y[ x%%y == 0 ]
  }


  # make an object subset with the given cluster
  Seurat_obj <- SubsetData(Seurat_obj,
                       ident.use = as.character(ident)
  )
  signature <- unique(
    intersect(signature, rownames(Seurat_obj@data))
    )

  sign_index <- which(rownames(Seurat_obj@data) == signature_name)

  genes_percentage <- c()
  # take a percentage for each gene into a vector
  for (i in 1:length(signature)){

    gene_index <- which(rownames(Seurat_obj@data) == signature[i])

    gene_contrib <- Seurat_obj@data[gene_index,] / Seurat_obj@data[sign_index,]

    q <- sum(gene_contrib) / length(gene_contrib)

    genes_percentage <- c(genes_percentage, q)

  }
  names(genes_percentage) <- signature
  # scale it
  genes_percentage <- scale(
    sort(genes_percentage,
         decreasing = T)
  )

  df <- data.frame(
    genes = rownames(genes_percentage),
    data = unname(genes_percentage[,1])
  )

  # get a divisor for the length of signature so the heatmap is simmetrical
  q <- divisors(length(signature))
  q <- q[q > 2][1]

  # plot it
  heatmap.2(matrix(df$data, ncol = q),
            trace="n",
            Colv = FALSE,
            Rowv=FALSE,
            cellnote = matrix(df$genes, ncol = q),
            notecol = 1,
            scale = "none",
            dendrogram = "none",
            labCol = "",
            labRow = "",
            cexRow = 0.75)

  print(genes_percentage)
}
