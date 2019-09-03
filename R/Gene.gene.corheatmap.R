#' heatmap for gene-gene correlation matrix
#'
#' This function allows you to plot a heatmap for gene-gene correlation
#'
#' @param Seurat_obj Seurat object
#' @param cor.method method to par ro cor function for correlation calculation, spearman is by default, bayesian (package psycho), pearson and kendall can also be used
#' @param coef.cut.off what monimum correlation coeffitient to choose to cut off the noise
#' @param gene.cut.off how much genes should have this correlation coefficient
#' @param imputed should MAGIC imputation be used on an expression matrix primarily to correlation analysis, FALSE by default
#' @param gene.names do you want do show gene names? Set to FALSE is you gonna have a large matrix
#' @param impute_after do you want to impute cor matrix that was generated after analysing unimputed data?
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
                                 coef.cut.off = 0.3,
                                 gene.cut.off = 1,
                                 use.hvg = F,
                                 imputed = F,
                                 impute_after = T,
                                 gene.names = T){

  print('preprocessing expression matrix')
  if (imputed == F){
    matr <- as.matrix(Seurat_obj@data)
  } else {
    matr <- Rmagic::magic(
      t(as.matrix(Seurat_obj@data))
      )
    matr <- t(matr[[1]])
  }

  #Seurat_obj@var.genes <- toupper(Seurat_obj@var.genes)
  rownames(Seurat_obj@data) <- toupper(rownames(Seurat_obj@data))
  Seurat_obj@var.genes <- toupper(Seurat_obj@var.genes)

  #get rid of ribo, mito and ig genes
  features_select <- grep("^MT[-,{RNR}]|^IG[H,L,K][A-Z][0-9].*|^RP[L,S][0-9].*|^RP[0-9].*|^FO[0-9]{2,}|^AP[0-9]{2,}|\\.|[0-9].*RIK", rownames(Seurat_obj@data), value = T)

  library(dplyr)
  if (use.hvg == T){
    var.genes <- Seurat_obj@hvg.info %>%
      tibble::rownames_to_column() %>%
      mutate(genes = toupper(rowname)) %>%
      top_n(n = 1500, wt = gene.dispersion) %>%
      pull(genes)

  # take only higly variable genes into account
    features <- var.genes[!(var.genes %in% features_select)]
  } else {
    features <- Seurat_obj@var.genes[!(Seurat_obj@var.genes %in% features_select)]
  }
  # get a matrix with only variable and functional genes
  matr <- matr[rownames(matr) %in% features,]

  # get a correlation matrix
  if (cor.method == 'bayesian'){
    cor.mat <- psycho::bayes_cor(t(matr))
  } else {
    cor.mat <- cor(t(matr),
                   method = cor.method)
  }

  print('calculating cor matrix')
  # clean in case of na values
  cor.mat[is.na(cor.mat)] <- 0

  # plot the heatmap using ComplexHeatmap package
  library(ComplexHeatmap)
  library(circlize)

  print('subsetting cor matrix')
  high.corr <- c()
  for (i in 1:nrow(cor.mat)){
    a <- sum(abs(cor.mat[,i]) > coef.cut.off) > gene.cut.off
    high.corr <- c(high.corr, a)
  }

  sub_matr <- matr[rownames(matr) %in% rownames(cor.mat)[high.corr], ]
  nrow(sub_matr)

  if (impute_after == T){
    sub_matr <- Rmagic::magic(
      t(as.matrix(sub_matr))
    )
    sub_matr <- t(sub_matr[[1]])
  }

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
          show_row_names = gene.names,
          show_column_names = gene.names)
}
