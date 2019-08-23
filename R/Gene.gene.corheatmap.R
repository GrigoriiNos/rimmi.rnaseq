#' heatmap for gene-gene correlation matrix
#'
#' This function allows you to plot a heatmap for gene-gene correlation
#'
#' @param Seurat_obj Seurat object
#'
#' @return a heaatmap
#'
#' @keywords single cell, correlation analysis, heatmap, data visualisation
#'
#' @examples
#'
#' Gene.gene.corheatmap(Seurat_obj, do.cluster = T)
#'
#' @export
#'
  
Gene.gene.corheatmap <- function(Seurat_obj, 
                                 cor.method = 'spearman', 
                                 cor.cut.off = 0.3){  
  
  print('preprocessing expression matrix')
  matr <- as.matrix(Seurat_obj@data)
  
  #get rid of ribo, mito and ig genes
  features_select <- grep("^MT[-,{RNR}]|^IG[H,L,K][A-Z][0-9].*|^RP[L,S][0-9].*|^RP[0-9].*|^FO[0-9]{2,}|^AP[0-9]{2,}|\\.", rownames(Seurat_obj@data), value = T)
  
  # take only higly variable genes into account
  features <- Seurat_obj@var.genes[!(Seurat_obj@var.genes %in% features_select)]
  
  # get a matrix with only variable and functional genes
  matr <- matr[rownames(matr) %in% features,]
  
  # get a correlation matrix
  cor.mat <- cor(t(matr), method = cor.method)
  
  print('calculating cor matrix')
  # clean in case of na values
  cor.mat[is.na(cor.mat)] <- 0
  
  # plot the heatmap using ComplexHeatmap package
  library(ComplexHeatmap)
  library(circlize)
  
  print('subsetting cor matrix')
  high.corr <- c()
  for (i in 1:nrow(cor.mat)){
    a <- sum(cor.mat[,i] > cor.cut.off) > 1
    high.corr <- c(high.corr, a)
  }
  
  sub_matr <- matr[rownames(matr) %in% rownames(cor.mat)[high.corr], ]
  sub_corr.mat <- cor(t(sub_matr), method = cor.method)
  
  col_fun = colorRamp2(c(-1, 0, 1),
                       c("blue", "white", "red"),
                       transparency = 0.5)
  
  print('started to plot heatmap, it can take quite a lot of time, grab the coffee or take a break')
  Heatmap(sub_corr.mat,
          name = "corr",
          col = col_fun,
          cluster_rows = T,
          cluster_columns = T,
          show_row_names = T,
          show_column_names = T)
}
