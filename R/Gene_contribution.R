#' Calculate the contribution of each gene for the given gene signature
#'
#' This function allows you to see how did each gene from the given transcriptional program contribute to the pre-calculated gene signature
#'
#' @param Seurat_obj your Seurat object
#' @param resolution what resolution you wanna check, either 0.6 or 1 or 2
#' @param signature the gene set of your interest
#' @param signature_name the name you gave to the signature
#'
#' @return heatmap with the gene contribution
#'
#' @keywords Seurat, single cell sequencing, RNA-seq, gene signature
#'
#' @examples
#'
#' Gene_contribution(A07, 0.6, signature, 'sLZ')
#'
#' @export
#'


Gene_contribution <- function(Seurat_obj,
                              resolution,
                              signature,
                              signature_name){

  library(Seurat)

  sign_index <- which(rownames(Seurat_obj@data) == signature_name)

  if (resolution == 0.6){
    ident <- sort(
      unique(
        Seurat_obj@meta.data$res.0.6))
  }
  if (resolution == 1){
    ident <- sort(
      unique(
        Seurat_obj@meta.data$res.1))
  }
  if (resolution == 2){
    ident <- sort(
      unique(
        Seurat_obj@meta.data$res.2))
  }
  if (!resolution %in% c(0.6, 1, 2)){
    print('nono please use only resolution 0.6 1 or 2')
    break
  }

  signature <- signature[signature %in% rownames(Seurat_obj@data)]

  output1 <- data.frame(genes = signature)
  output2 <- data.frame(genes = signature)
  for (i in 1:length(ident)){
  # make an object subset with the given cluster
  sub_Seurat_obj <- SubsetData(Seurat_obj,
                       ident.use = ident[i]
  )

  genes_percentage <- c()
  # take a percentage for each gene into a vector
    for (j in 1:length(signature)){

      # contibution of genes from the signature
      gene_index <- which(rownames(sub_Seurat_obj@data) == signature[j])

      scaler <- max(
        abs(
          sub_Seurat_obj@data[sign_index,])
        )

      gene_contrib <- sub_Seurat_obj@data[gene_index,] + scaler / sub_Seurat_obj@data[sign_index,] + scaler

      q <- sum(gene_contrib) / length(gene_contrib)

      genes_percentage <- c(genes_percentage, q)
      }

  names(genes_percentage) <- signature
  # scale it
  genes_percentage <- scale(
    sort(genes_percentage,
         decreasing = T)
  )

  df1 <- data.frame(
    genes = rownames(genes_percentage),
    data = unname(genes_percentage[,1])
  )
  colnames(df1)[2] <- paste('cluster', ident[i])
  output1 <- merge(output1, df1, 'genes')

  #the average expression per clsuter for a given gene
  genes_average <- rowSums(
    sub_Seurat_obj@data[rownames(sub_Seurat_obj@data) %in% signature, ]
  ) / ncol(sub_Seurat_obj@data)

  df2 <- data.frame(
    genes = names(genes_average),
    data = unname(genes_average)
  )

  colnames(df2)[2] <- paste('cluster', ident[i])
  output2 <- merge(output2, df2, 'genes')

  }

  rownames(output1) <- output1$genes
  rownames(output2) <- output2$genes

  write.csv(output1,
            quote = F,
            file = paste0('genes from the ',signature_name, ' signature, resolution', resolution, '.csv')
            )

  write.csv(output2,
            quote = F,
            file = paste0('average gene exp from the ',signature_name, ' signature, resolution', resolution, '.csv')
  )

  heatmap(
    as.matrix(output1[,-1])
    )
}


