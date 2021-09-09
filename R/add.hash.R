#' Merge cell hashcoding demultiplexed data with GEX Seurat 3 object with VDJ metadata
#'
#' This function allows you to fuse the output of hashcode demultiplexing with GEX+VDJ object
#'
#' @param hash.matrix a cellular bacrode vs hashtags matrix
#' @param Seurat_obj  Seurat object
#'
#' @return Seurat 3 object merged with hashtag demiltiplexed identities in the metadata slot
#'
#' @keywords Seurat, single cell sequencing, RNA-seq, RNA velocity
#'
#' @examples
#'
#' A07 <- add_hash(Seurat_objhash, Seurat_obj, 'A14')
#'
#'
#' @export

add_hash <- function(hash.matrix, Seurat_obj, sample_name, quantile_theshold = 0.99, method=c('HTODemux', 'MULTIhash')){

  #cells0tags <- which(apply(hash.matrix, 2, function(x) all(x == 0)))
  #tags0everywhere <- which(apply(hash.matrix, 1, function(x) all(x == 0)))

  Seurat_obj.umis <- Seurat_obj@assays$RNA@counts

  # Select cell barcodes detected by both RNA and HTO In the example datasets we have already
  # filtered the cells for you, but perform this step for clarity.
  #hash.matrix <- hash.matrix[,!apply(hash.matrix[c(3,5),], 2, function(x) all(x == 0))]

  joint.bcs <- intersect(colnames(Seurat_obj.umis), colnames(hash.matrix))

  # Subset RNA and HTO counts by joint cell barcodes
  Seurat_obj.umis <- Seurat_obj.umis[, joint.bcs]
  hash.matrix <- as.matrix(hash.matrix[, joint.bcs])

  # Confirm that the HTO have the correct names
  rownames(hash.matrix)
  #hash.matrix <- hash.matrix + 1

  # Setup Seurat object
  Seurat_obj.hashtag <- CreateSeuratObject(counts = Seurat_obj.umis)
  rm(Seurat_obj.umis)
  # Normalize RNA data with log normalization
  Seurat_obj.hashtag <- NormalizeData(Seurat_obj.hashtag)
  # Find and scale variable features
  Seurat_obj.hashtag <- FindVariableFeatures(Seurat_obj.hashtag, selection.method = "mean.var.plot")
  Seurat_obj.hashtag <- ScaleData(Seurat_obj.hashtag, features = VariableFeatures(Seurat_obj.hashtag))

  #Add HTO data as a new assay independent from RNA
  Seurat_obj.hashtag[["HTO"]] <- CreateAssayObject(counts = hash.matrix)
  # Normalize HTO data, here we use centered log-ratio (CLR) transformation
  Seurat_obj.hashtag <- NormalizeData(Seurat_obj.hashtag, assay = "HTO", normalization.method = "CLR")

  # If you have a very large dataset we suggest using k_function = 'clara'. This is a k-medoid
  # clustering function for large applications You can also play with additional parameters (see
  # documentation for HTODemux()) to adjust the threshold for classification Here we are using the
  # default settings
  if (method == 'HTODemux'){
    Seurat_obj.hashtag <- HTODemux(Seurat_obj.hashtag, init=2, assay = "HTO", positive.quantile = quantile_theshold)
    print(
      table(Seurat_obj.hashtag$HTO_classification)
    )
  } else if (method == 'MULTIhash'){
    Seurat_obj.hashtag <- MULTIseqDemux(Seurat_obj.hashtag, assay = "HTO")
    print(
      table(Seurat_obj.hashtag$MULTI_classification)
    )
  }

  # Group cells based on the max HTO signal
  Idents(Seurat_obj.hashtag) <- "HTO_maxID"

  DefaultAssay(Seurat_obj.hashtag) <- 'HTO'
  try({
    Ridge <- RidgePlot(Seurat_obj.hashtag,
                       assay = "HTO",
                       features = rownames(Seurat_obj.hashtag[["HTO"]])[1:2], ncol = 2)
    ggsave(paste0(sample_name, '_RidgePlot_QC.jpg'), Ridge, width = 15, height = 12)

    Feat <- FeatureScatter(Seurat_obj.hashtag,
                           feature1 = rownames(Seurat_obj.hashtag[["HTO"]])[1],
                           feature2 = rownames(Seurat_obj.hashtag[["HTO"]])[2])

    ggsave(paste0(sample_name, '_FeatureScatter_QC12.jpg'), Feat, width = 15, height = 12)

    Feat <- FeatureScatter(Seurat_obj.hashtag,
                           feature1 = rownames(Seurat_obj.hashtag[["HTO"]])[3],
                           feature2 = rownames(Seurat_obj.hashtag[["HTO"]])[4])

    ggsave(paste0(sample_name, '_FeatureScatter_QC34.jpg'), Feat, width = 15, height = 12)

    Feat <- FeatureScatter(Seurat_obj.hashtag,
                           feature1 = rownames(Seurat_obj.hashtag[["HTO"]])[2],
                           feature2 = rownames(Seurat_obj.hashtag[["HTO"]])[5])

    ggsave(paste0(sample_name, '_FeatureScatter_QC25.jpg'), Feat, width = 15, height = 12)

    Feat <- FeatureScatter(Seurat_obj.hashtag,
                           feature1 = rownames(Seurat_obj.hashtag[["HTO"]])[1],
                           feature2 = rownames(Seurat_obj.hashtag[["HTO"]])[5])

    ggsave(paste0(sample_name, '_FeatureScatter_QC15.jpg'), Feat, width = 15, height = 12)

    Idents(Seurat_obj.hashtag) <- Seurat_obj.hashtag@meta.data$HTO_classification
    Seurat_obj.hashtag$HTO_classification_short <- gsub('-|[A-Z]+','', Seurat_obj.hashtag$HTO_classification)
    hto.heatmap(Seurat_obj.hashtag)

    Vld <- VlnPlot(Seurat_obj.hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE)+
      theme(legend.position = "none")
    ggsave(paste0(sample_name, '_VlnPlot_QC.jpg'), Vld, width = 20, height = 12)
  })

  meta <- merge(Seurat_obj@meta.data[, !(colnames(Seurat_obj@meta.data) %in% c('nCount_HTO', 'nFeature_HTO', 'HTO_maxID',
                                                                               'HTO_secondID', 'HTO_margin', 'HTO_classification',
                                                                               'HTO_classification.global', 'hash.ID'))],
                Seurat_obj.hashtag@meta.data[, c('nCount_HTO', 'nFeature_HTO', 'HTO_maxID',
                                                 'HTO_secondID', 'HTO_margin', 'HTO_classification',
                                                 'HTO_classification.global', 'hash.ID')],
                by = 'row.names',
                all.x = T)

  rownames(meta) <- meta$Row.names
  meta$Row.names <- NULL

  Seurat_obj@meta.data <- meta[Cells(Seurat_obj),]

  Seurat_obj
}

#' mark tags that you think are useless as such
#'
#' This function allows you to  mark tags that you think are useless as such
#'
#' @param Seurat_obj  Seurat object
#' @param useless_tags indexes of tags that you think are useless
#'
#' @return Seurat 3 object
#'
#' @keywords Seurat, single cell sequencing, RNA-seq, RNA velocity
#'
#' @examples
#'
#' throw_away_tags(Seurat_obj,
#' useless_tags =  c(1,2,3))
#'
#'
#' @export

throw_away_tags <- function(Seurat_obj,
                            useless_tags){
  # variable with all the tags
  tags <- unique(Seurat_obj$HTO_classification)
  # get actual tags here instead of indexes
  useless_tags <- tags[useless_tags]
  # save old tags for a back uo
  Seurat_obj@meta.data$HTO_classification_back_up <- Seurat_obj@meta.data$HTO_classification
  # change uselesss as such
  Seurat_obj@meta.data$HTO_classification[Seurat_obj@meta.data$HTO_classification == useless_tags] <- 'useless_tags'

  Seurat_obj
}

#' hto heatmap
#'
#'
#' @param Seurat_obj.hashtag  Seurat object
#'
#' @export

hto.heatmap <- function(Seurat_obj.hashtag, scale=T){
  library(ComplexHeatmap)

  hto.mat <- Seurat_obj.hashtag@assays$HTO@data
  hto.mat <- hto.mat[rownames(hto.mat)!='unmapped',]
  rownames(hto.mat) <- gsub('[A-Z]|-', '', rownames(hto.mat))
  tagshort <- Seurat_obj.hashtag$HTO_classification_short
  tagglobal <- Seurat_obj.hashtag$HTO_classification.global

  if (scale){
    hto.mat <- t(scale(t(hto.mat)))
  }

  tagsann <- HeatmapAnnotation(hto.class = tagshort,
                               global = tagglobal,
                               which = 'col',
                               annotation_width = unit(c(1, 4), 'cm'),
                               gap = unit(1, 'mm'))

  hmap <- Heatmap(
    hto.mat[rownames(hto.mat) != 'tag4',],
    column_order = order(tagshort),
    show_row_names = T,
    show_column_names = FALSE,
    cluster_rows = F,
    cluster_columns = F,
    clustering_method_rows = "ward.D2",
    clustering_method_columns = "ward.D2",
    top_annotation=tagsann)

  draw(hmap, heatmap_legend_side="left", annotation_legend_side="right")
}




