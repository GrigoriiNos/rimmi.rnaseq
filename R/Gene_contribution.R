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


Gene_contribution <- function(Seurat_obj, ident, signature){
  library(gplots)
  library(Seurat)
  
  print('you store the signature after calculation in the 1st row of the matrix, it should have a name that you gave it')
  print('so make sure that 1st row of the matrix is signature')
  print(rownames(Seurat_obj@data)[1])
  print('if that was a gene name, then it was wrong, you didnt calculate the avarage expression')
  
  Seurat_obj <- SubsetData(Seurat_obj, 
                       ident.use = as.character(ident)
  )
  
  genes_percentage <- c()
  for (i in 1:length(signature)){
    
    gene_index <- which(rownames(Seurat_obj@data) == signature[i])
    
    gene_contrib <- Seurat_obj@data[gene_index,] / Seurat_obj@data[1,]
    
    q <- sum(gene_contrib) / length(gene_contrib)
    
    genes_percentage <- c(genes_percentage, q)
    
  }
  names(genes_percentage) <- signature
  genes_percentage <- scale(
    sort(genes_percentage, 
         decreasing = T)
  )
  
  df <- data.frame(
    genes = rownames(genes_percentage),
    data = unname(genes_percentage[,1])
  )
  
  heatmap.2(matrix(df$data, ncol = 3), 
            trace="n", 
            Colv = FALSE, 
            Rowv=FALSE,
            cellnote = matrix(df$genes, ncol = 3), 
            notecol = 1, 
            scale = "none",
            dendrogram = "none", 
            labCol = "", 
            labRow = "", 
            cexRow = 0.75)
  
  print(genes_percentage)
}
