#' heatmap for cell-cell correlation matrix
#'
#' This function allows you to plot a heatmap for cell-cell correlation
#'
#' @param Seurat_obj Seurat object
#'
#' @param do.cluster set to T if you want a hierarchical clustering on the top of the heatmap, defulat is FALSE
#'
#' @param colours set the colours for clusters that you want, random colours by default
#'
#' @return a heaatmap
#'
#' @keywords single cell, correlation analysis, heatmap, data visualisation
#'
#' @examples
#'
#' Cell.cell.corheatmap(Seurat_obj, do.cluster = T)
#'
#' @export
#'

Cell.cell.corheatmap <- function(Seurat_obj,
                                 do.cluster = F,
                                 colours = 'random',
                                 cor.method = 'spearman',
                                 title = "Cell-Cell corellation heatmap"){

  print('converting matrix')
  matr <- as.matrix(Seurat_obj@data)

  print('selecting genes')
  #get rid of ribo, mito and ig genes
  features_select <- grep("^MT[-,{RNR}]|^IG[H,L,K][A-Z][0-9].*|^RP[L,S][0-9].*|^RP[0-9].*|^FO[0-9]{2,}|^AP[0-9]{2,}|\\.", rownames(Seurat_obj@data), value = T)

  # take only higly variable genes into account
  features <- Seurat_obj@var.genes[!(Seurat_obj@var.genes %in% features_select)]

  # get a matrix with only variable and functional genes
  matr <- matr[rownames(matr) %in% features,]

  print('calculating correlation matrix')
  # get a correlation matrix
  cor.mat <- cor(matr, method = cor.method)

  # clean in case of na values
  cor.mat[is.na(cor.mat)] <- 0

  # take a meta data for heat map
  meta.data <- data.frame(Phase = Seurat_obj@meta.data$Phase,
                          Cluster = Seurat_obj@meta.data$res.0.6,
                          Cellname = rownames(Seurat_obj@meta.data)
  )

  print('preparing metadata for the plot')
  # take groups
  group_size <- table(Seurat_obj@ident)

  # assign each cell barcode to cluster
  gl  <-  lapply(1:length(group_size), function(i) {
    rownames(cor.mat)[sum(group_size[seq_len(i-1)]) + 1:group_size[i]]
  })

  # name them by the cluster
  names(gl) <-  paste0("Cluster ", 0:(length(group_size)-1))

  gd <-  structure(rep(names(gl),
                 times = sapply(gl, length)),
                 names = unlist(gl)
                 )

  # give a random colour to each cluster got a heatmap
  if (colours == 'random'){
    group_color <-  structure(circlize::rand_color(length(group_size)), names = names(gl))
  } else {
    group_color <- colours
  }

  n_group = length(gl)

  # plot the heatmap using ComplexHeatmap package
  library(ComplexHeatmap)
  library(circlize)

  col_fun = colorRamp2(c(-1, 0, 1),
                       c("darkblue", "white", "red"),
                       transparency = 0.5)

  print('started to plot heatmap, it can take quite a lot of time, grab the coffee or take a break')
  Heatmap(cor.mat,
          name = "corr",
          col = col_fun,
          cluster_rows = do.cluster,
          cluster_columns = do.cluster,
          show_row_names = FALSE,
          show_column_names = FALSE,
          column_title = title,
          top_annotation = HeatmapAnnotation(group = gd, col = list(group = group_color), show_legend = FALSE)) +
    rowAnnotation(group = gd, col = list(group = group_color), width = unit(0.5, "cm"))
}















