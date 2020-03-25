#' This function allows you to calculate DEG for the IG hashtags
#'
#' @param Seurat_obj Seurat object with integrated data from VDJ analysis with Ihashtag assigment
#'
#' @keywords DEG, single cell, VDJ analysis, B cells, T cells
#'
#' @examples
#'
#' hashtags_DEG(A07)
#'
#' @export
#'

hashtags_DEG <- function(Seurat_obj, demult.method = 'HTODemux'){
  library(Seurat)
  # align cells in the same order
  Seurat_obj@meta.data <- Seurat_obj@meta.data[names(Idents(Seurat_obj)),]
  
  # replace cell identity from cluster to ihashtag
  if (demult.method == 'HTODemux'){
    Idents(Seurat_obj) <- factor(Seurat_obj@meta.data$HTO_classification)
  } else {
    Idents(Seurat_obj) <- factor(Seurat_obj@meta.data$MULTI_ID)
  }
  names(Idents(Seurat_obj)) <- rownames(Seurat_obj@meta.data)
  
  Seurat_obj_sub <- subset(Seurat_obj, 
                           idents = Idents(Seurat_obj)[!is.na(Idents(Seurat_obj)) & Idents(Seurat_obj) != 'Negative'])
  
  print('calculating markers for short hashtags')
  hash_markers <- FindAllMarkers(object = Seurat_obj_sub,
                                      only.pos = F,
                                      min.pct = 0.25,
                                      thresh.use = 0.25)
  
  rimmi.rnaseq::markers_to_xls(hash_markers, filename = 'hashtags_DEGs')
  print('saved hashtags DEG  as an excel spreadsheet')
  
  print(
    hash_markers %>%
      group_by(cluster) %>%
      top_n(3, avg_logFC)
  )
}

#' This function allows you to calculate DEG for some specific IG hashtags
#'
#' @param Seurat_obj Seurat object with integrated data from VDJ analysis with Ihashtag assigment
#'
#' @keywords DEG, single cell, VDJ analysis, B cells, T cells
#'
#' @examples
#'
#' hashtags_DEG(A07)
#'
#' @export
#'

Specific_hashtags_DEG <- function(Seurat_obj,
                                  hash1,
                                  hash2,
                                  filename,
                                  demult.method = 'HTODemux',
                                  clusters = NULL){
  library(Seurat)
  # resort cells in the right order
  Seurat_obj@meta.data <- Seurat_obj@meta.data[names(Idents(Seurat_obj)),]
  
  if (!is.null(clusters)){
    Seurat_obj <- subset(Seurat_obj,
                         idents = clusters)
  }
  # replace cell identity from cluster to ihashtag
  # replace cell identity from cluster to ihashtag
  if (demult.method == 'HTODemux'){
    Idents(Seurat_obj) <- factor(Seurat_obj@meta.data$HTO_classification)
  } else {
    Idents(Seurat_obj) <- factor(Seurat_obj@meta.data$MULTI_ID)
  }
  names(Idents(Seurat_obj)) <- rownames(Seurat_obj@meta.data)
  
  Seurat_obj_sub <- subset(Seurat_obj, 
                           idents = Idents(Seurat_obj)[!is.na(Idents(Seurat_obj)) & Idents(Seurat_obj) != 'Negative'])
  
  print('calculating markers for specific hashtags')
  
  hash_markers <- FindMarkers(object = Seurat_obj_sub,
                             ident.1 = hash1,
                             ident.2 = hash2,
                             only.pos = F,
                             min.pct = 0.25,
                             thresh.use = 0.25)
  
  write.csv(hash_markers,
            file = filename,
            quote = F,
            row.names = T)
  print('saved hashtags DEG as a csv file')
}


#' This function allows you to pull specific barcodes from the object with specific IG hashtags from specific clusters
#'
#'
#' @param Seurat_obj Seurat object with integrated data from VDJ analysis with Ihashtag assigment
#' @param hashtags Ihashtag names
#' @param clusters cluster names
#'
#' @keywords VDJ analysis, B cells, T cells
#'
#' @examples
#'
#' pull_ihashtag(A07, 'IGHA1')
#'
#' @export
#'

pull_hashtag <- function(Seurat_obj,
                         hashtags, 
                         demult.method = 'HTODemux',
                         clusters = NULL){
  
  Seurat_obj@meta.data <- Seurat_obj@meta.data[names(Idents(Seurat_obj)),]
  if (demult.method == 'HTODemux'){
    if (is.null(clusters)){
      hash.cells <- rownames(
        Seurat_obj@meta.data)[Seurat_obj@meta.data$HTO_classification %in% hashtags]
    } else{
      hash.cells <- rownames(
        Seurat_obj@meta.data)[(Seurat_obj@meta.data$HTO_classification %in% hashtags) & (Idents(Seurat_obj) %in% clusters)]
    }
  } else {
    if (is.null(clusters)){
      hash.cells <- rownames(
        Seurat_obj@meta.data)[Seurat_obj@meta.data$MULTI_ID %in% hashtags]
    } else{
      hash.cells <- rownames(
        Seurat_obj@meta.data)[(Seurat_obj@meta.data$MULTI_ID %in% hashtags) & (Idents(Seurat_obj) %in% clusters)]
    }
  }
  hash.cells
}
