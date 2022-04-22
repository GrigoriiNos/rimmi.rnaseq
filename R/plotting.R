#' heatmap for cell-cell correlation matrix
#'
#' This function allows you to plot a heatmap for cell-cell correlation
#'
#' @param Seurat_obj Seurat object
#' @param do.cluster set to T if you want a hierarchical clustering on the top of the heatmap, defulat is FALSE
#' @param colours set the colours for clusters that you want, random colours by default
#' @param cor.method what method parameter should be passed to cor function, spearman is by default, bayesian (package psycho), pearson and kendall can also be used
#' @param @param imputed should MAGIC imputetion be used on an expression matrix primarily to correlation analysis, FALSE by default
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
                                 imputed = F,
                                 title = "Cell-Cell corellation heatmap"){

  print('converting matrix')
  if (imputed == F){
    matr <- as.matrix(Seurat_obj@data)
  } else {
    matr <- Rmagic::magic(t(as.matrix(Seurat_obj@data)))
    matr <- t(matr[[1]])
  }

  print('selecting genes')
  #get rid of ribo, mito and ig genes
  features_select <- grep("^MT[-,{RNR}]|^IG[H,L,K][A-Z][0-9].*|^RP[L,S][0-9].*|^RP[0-9].*|^FO[0-9]{2,}|^AP[0-9]{2,}|\\.", rownames(Seurat_obj@data), value = T)

  # take only higly variable genes into account
  features <- Seurat_obj@var.genes[!(Seurat_obj@var.genes %in% features_select)]

  # get a matrix with only variable and functional genes
  matr <- matr[rownames(matr) %in% features,]

  print('calculating correlation matrix')
  # get a correlation matrix
  if (cor.method == 'bayesian'){
    cor.mat <- psycho::bayes_cor(matr)
  } else {
    cor.mat <- cor(matr,
                   method = cor.method)
  }


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



#' ggplot extension for split violin plot
#'
#' @return ggplot extension
#'
#' @keywords violin plot, ggplot
#'
#' @examples
#'
#' ggplot(data, aes(x=x, y=y, fill = factor))+
#' geom_split_violin
#'
#' @export

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ...,
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE,
                              show.legend = NA, inherit.aes = TRUE) {
  GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin,
                             draw_group = function(self, data, ..., draw_quantiles = NULL) {
                               data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                               grp <- data[1, "group"]
                               newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                               newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                               newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])

                               if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                                 stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                           1))
                                 quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                                 aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                                 aesthetics$alpha <- rep(1, nrow(quantiles))
                                 both <- cbind(quantiles, aesthetics)
                                 quantile_grob <- GeomPath$draw_panel(both, ...)
                                 ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                               }
                               else {
                                 ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                               }
                             })
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin,
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}


#' save multiple gg plots in pdf
#'
#' This function allows you to save multiple gg plots in pdf
#'
#'
#' @param ...  ggplot objects
#' @param filename filename
#'
#' @keywords ggplot, save
#'
#' @examples
#' p1 <- qplot(1:10)
#' p2 <- qplot(1:20)
#' p3 <- qplot(1:30)
#'
#' save.plot.grid(p1,p2,p3, 'whatever')
#'
#'
#' @export

save.plot_grid <- function(..., filename){
  library(ggplot2)
  library(gridExtra)

  grid.arrange(...) #arranges plots within grid

  #save
  g <- arrangeGrob(...) #generates g

  ggsave(file = paste0(filename, ".pdf"), g) #saves
}


#' Make a dot plot for all your markers
#'
#' This function allows you to plot each top X DEG for each cell cluster using a Seurat::FindAllMarkers function output
#'
#' @param markers_table The output of the Seurat::FndAllMarkers function, don't change anything in the column names within it.
#'
#' @param X How much of the genes you need to analyse in each cluster. 5 is a default.
#'
#' @param Seurat_obj your Seurat object to plot markers from
#'
#' @keywords dot plot, Seurat, markers, single cell sequencing, RNA-seq
#'
#' @examples
#'
#' dot_plot_topXgenes(markers, 10)
#'
#' @export
#'


dot_plot_topXgenes <- function(markers_table, X = 5, Seurat_obj){
  library(Seurat)
  library(dplyr)

  ###
  plot_top_X_genes <- function(markers_table, X = 5){
    top_X <- markers_table %>%
      group_by(cluster) %>%
      top_n(X, avg_logFC)
    top_X$gene
  }
  ###
  top_5_B06 <- unique(
    plot_top_X_genes(markers_table, X)
  )
  DotPlot(
    b06,
    Seurat_obj,
    x.lab.rot = T
  )
}


set.seed(100)
rnorm(n = 100, mean = 0, sd = 2)




#' Plot clusters in 2 umaps with the point size corresponting to the library size
#'
#' This function allows you to see if the library size influences in your clustering
#'
#' @param Seurat_obj your Seurat object
#'
#' @return scatter plot with the library size correponding to the point size
#'
#' @keywords Seurat, single cell sequencing, RNA-seq, gene signature
#'
#' @examples
#'
#' scatter_libdepth(A07)
#'
#' @export
#'

lib.qc_plot <- function(Seurat_obj){
  library(ggplot2)

  qc.df <- data.frame(row.names = rownames(Seurat_obj@meta.data),
                      umap1 = Seurat_obj@reductions$umap@cell.embeddings[,1],
                      umap2 = Seurat_obj@reductions$umap@cell.embeddings[,2],
                      nUMI = Seurat_obj$nCount_RNA,
                      pMito = Seurat_obj$percent.mt,
                      nGenes = Seurat_obj$nFeature_RNA,
                      cluster = Idents(Seurat_obj))

  p1 <- ggplot(qc.df, aes(x = umap1, y = umap2, colour = cluster)) +
    geom_point(aes(size = nUMI)) +
    theme_bw()


  p2 <- ggplot(qc.df, aes(x = umap1, y = umap2, colour = cluster)) +
    geom_point(aes(size = nGenes)) +
    theme_bw()

  p3 <- ggplot(qc.df, aes(x = umap1, y = umap2, colour = cluster)) +
    geom_point(aes(size = pMito)) +
    theme_bw()

  cowplot::plot_grid(p1, p2, p3)
}

