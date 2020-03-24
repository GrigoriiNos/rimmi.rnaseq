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

add_hash <- function(hash.matrix, Seurat_obj, sample_name){
  Seurat_obj.umis <- Seurat_obj@assays$RNA@counts

  # Select cell barcodes detected by both RNA and HTO In the example datasets we have already
  # filtered the cells for you, but perform this step for clarity.
  joint.bcs <- intersect(colnames(Seurat_obj.umis), colnames(hash.matrix))
  
  # Subset RNA and HTO counts by joint cell barcodes
  Seurat_obj.umis <- Seurat_obj.umis[, joint.bcs]
  hash.matrix <- as.matrix(hash.matrix[, joint.bcs])
  
  # Confirm that the HTO have the correct names
  rownames(hash.matrix)
  
  # Setup Seurat object
  Seurat_obj.hashtag <- CreateSeuratObject(counts = Seurat_obj.umis)
  
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
  colnames(Seurat_obj.hashtag@assays$HTO@counts)[colSums(Seurat_obj.hashtag@assays$HTO@counts) < 100]
  rownames(Seurat_obj.hashtag@assays$HTO@counts)[rowSums(Seurat_obj.hashtag@assays$HTO@counts) < 10000]
  
  #Seurat_obj.hashtag <- HTODemux(Seurat_obj.hashtag, assay = "HTO", positive.quantile = 0.99)
  Seurat_obj.hashtag <- MULTIseqDemux(Seurat_obj.hashtag, assay = "HTO")
  
  # Global classification results
  print(
    table(Seurat_obj.hashtag$MULTI_classification)
    )
  
  # Group cells based on the max HTO signal
  Idents(Seurat_obj.hashtag) <- "HTO_maxID"

  Ridge <- RidgePlot(Seurat_obj.hashtag, assay = "HTO", features = rownames(Seurat_obj.hashtag[["HTO"]])[1:2], ncol = 2)
  ggsave(paste0(sample_name, '_RidgePlot_QC.jpg'), Ridge)  
  
  Feat <- FeatureScatter(Seurat_obj.hashtag, feature1 = "tag1-GTCAACTCTTTAGCG", feature2 = "tag2-TGATGGCCTATTGGG")
  ggsave(paste0(sample_name, '_FeatureScatter_QC.jpg'), Feat)
  
  Idents(Seurat_obj.hashtag) <- Seurat_obj.hashtag@meta.data$HTO_classification
  #Idents(Seurat_obj.hashtag) <- Seurat_obj.hashtag@meta.data$MULTI_classification
  
  Vld <- VlnPlot(Seurat_obj.hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE)+
    theme(legend.position = "none")
  ggsave(paste0(sample_name, '_VlnPlot_QC.jpg'), Vld)
  
  meta <- merge(Seurat_obj@meta.data,
                Seurat_obj.hashtag@meta.data[, c('nCount_HTO', 'nFeature_HTO', 'HTO_maxID', 'HTO_secondID', 'HTO_margin', 'HTO_classification', 'HTO_classification.global', 'hash.ID')],
                by = 0, 
                all.x = T)
  
  rownames(meta) <- meta$Row.names
  meta$Row.names <- NULL
  
  Seurat_obj@meta.data <- meta
  
  Seurat_obj
}







