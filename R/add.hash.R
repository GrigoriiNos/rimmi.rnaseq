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
    Seurat_obj.hashtag <- HTODemux(Seurat_obj.hashtag, init=NULL, assay = "HTO", positive.quantile = quantile_theshold)
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
    Seurat_obj.hashtag$HTO_classification_short <- gsub('-[A-Z]+','', Seurat_obj.hashtag$HTO_classification)
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

#' HTODemux
#' in house modified Seurat::HTODemux for identifiyng multuplets
#'
#' @param object  Seurat object
#'
#' @export

HTODemux <- function(
  object,
  assay = "HTO",
  positive.quantile = 0.99,
  init = NULL,
  nstarts = 100,
  kfunc = "clara",
  nsamples = 100,
  seed = 42,
  verbose = TRUE,
  patient = 'patient1'
) {
  if (!is.null(x = seed)) {
    set.seed(seed = seed)
  }
  #initial clustering
  assay <- assay %||% DefaultAssay(object = object)
  data <- GetAssayData(object = object, assay = assay)
  data <- data[rownames(data) != 'unmapped',]
  counts <- GetAssayData(
    object = object,
    assay = assay,
    slot = 'counts'
  )[, colnames(x = object)]
  counts <- as.matrix(x = counts)
  counts <- counts[rownames(counts)!='unmapped',]
  ncenters <- init %||% (nrow(x = data) + 1)
  switch(
    EXPR = kfunc,
    'kmeans' = {
      init.clusters <- kmeans(
        x = t(x = GetAssayData(object = object, assay = assay)),
        centers = ncenters,
        nstart = nstarts
      )
      #identify positive and negative signals for all HTO
      Idents(object = object, cells = names(x = init.clusters$cluster)) <- init.clusters$cluster
    },
    'clara' = {
      #use fast k-medoid clustering
      init.clusters <- clara(
        x = t(x = GetAssayData(object = object, assay = assay)),
        k = ncenters,
        samples = nsamples
      )
      #identify positive and negative signals for all HTO
      Idents(object = object, cells = names(x = init.clusters$clustering), drop = TRUE) <- init.clusters$clustering
    },
    stop("Unknown k-means function ", kfunc, ", please choose from 'kmeans' or 'clara'")
  )
  #average hto signals per cluster
  #work around so we don't average all the RNA levels which takes time
  average.expression <- AverageExpression(
    object = object,
    assays = assay,
    verbose = FALSE
  )[[assay]]
  #checking for any cluster with all zero counts for any barcode
  if (sum(average.expression == 0) > 0) {
    stop("Cells with zero counts exist as a cluster.")
  }
  #create a matrix to store classification result
  discrete <- GetAssayData(object = object, assay = assay)
  discrete[discrete > 0] <- 0
  # for each HTO, we will use the minimum cluster for fitting
  ndata <- object@assays$HTO@data
  raw.data <- object@assays$HTO@counts

  rownames(ndata) <- gsub('-','_', rownames(ndata))
  rownames(raw.data) <- gsub('-','_', rownames(raw.data))

  for (iter in rownames(x = data)) {
    values <- counts[iter, colnames(object)]
    #commented out if we take all but the top cluster as background
    #values_negative=values[setdiff(object@cell.names,WhichCells(object,which.max(average.expression[iter,])))]
    values.use <- values[WhichCells(
      object = object,
      idents = levels(x = Idents(object = object))[[which.min(x = average.expression[iter, ])]]
    )]
    fit <- suppressWarnings(expr = fitdist(data = values.use, distr = "nbinom"))
    cutoff <- as.numeric(x = quantile(x = fit, probs = positive.quantile)$quantiles[1])
    ##### ARTIFICIAL DEMUXING
    #cutoff <- cut_offs_pat1[iter]
    ##### ARTIFICIAL DEMUXING
    discrete[iter, names(x = which(x = values > cutoff))] <- 1
    if (verbose) {
      message(paste0("Cutoff for ", iter, " : ", cutoff, " reads"))
    }

    tag <- gsub('-','_', iter)
    p <- ggplot(as_tibble(t(as.matrix(raw.data))), aes_string(x=tag)) +
      geom_histogram(aes(y=..density..), colour="black", fill="white")+
      geom_density(alpha=.2, fill="blue") + xlim(c(-.1, cutoff+200)) + ggtitle(paste0('raw counts dist ', gsub('_.+','', tag))) +
      geom_vline(aes(xintercept=cutoff), color="red", linetype="dashed", size=1) + theme_bw()

    cells <- names(raw.data[tag,][raw.data[tag,] > cutoff-3 & raw.data[tag,] < cutoff+3])
    cutoffs.norm <- mean(ndata[, cells][tag,])

    p2 <- ggplot(as_tibble(t(ndata)), aes_string(x=tag)) +
      geom_histogram(aes(y=..density..), colour="black", fill="white")+
      geom_density(alpha=.2, fill="blue") + ggtitle(paste0('normalised counts dist ', gsub('_.+','', tag))) +
      geom_vline(aes(xintercept=cutoff.norm), color="red", linetype="dashed", size=1) + theme_bw()

    setwd('/Volumes/Elements/RIMMI/Adam_germinal_cenrte/T_cells_together/hto_panT/')
    ggsave(plot=p, paste0('hists/', patient, '/raw_dist_',tag, 'cutoff=', positive.quantile, '.jpg'))
    ggsave(plot=p2, paste0('hists/', patient, '/norm_dist_', tag,'cutoff=' ,positive.quantile, '.jpg'))

  }
  # now assign cells to HTO based on discretized values
  npositive <- colSums(x = discrete)
  classification.global <- npositive
  classification.global[npositive == 0] <- "Negative"
  classification.global[npositive == 1] <- "Singlet"
  classification.global[npositive > 1] <- "Doublet"
  classification.global[npositive > 2] <- "Triplet"
  classification.global[npositive > 3] <- "Quadruplet"
  classification.global[npositive > 4] <- "Pentaplet"

  print(table(classification.global))

  donor.id = rownames(x = data)
  hash.max <- apply(X = data, MARGIN = 2, FUN = max)
  hash.maxID <- apply(X = data, MARGIN = 2, FUN = which.max)
  hash.second <- apply(X = data, MARGIN = 2, FUN = MaxN, N = 2)
  hash.third <- apply(X = data, MARGIN = 2, FUN = MaxN, N = 3)
  hash.maxID <- as.character(x = donor.id[sapply(
    X = 1:ncol(x = data),
    FUN = function(x) {
      return(which(x = data[, x] == hash.max[x])[1])
    }
  )])
  hash.secondID <- as.character(x = donor.id[sapply(
    X = 1:ncol(x = data),
    FUN = function(x) {
      return(which(x = data[, x] == hash.second[x])[1])
    }
  )])
  hash.thirdID <- as.character(x = donor.id[sapply(
    X = 1:ncol(x = data),
    FUN = function(x) {
      return(which(x = data[, x] == hash.third[x])[1])
    }
  )])
  hash.margin <- hash.max - hash.second
  doublet_id <- sapply(
    X = 1:length(x = hash.maxID),
    FUN = function(x) {
      return(paste(sort(x = c(hash.maxID[x], hash.secondID[x])), collapse = "_"))
    }
  )
  triplet_id <- sapply(
    X = 1:length(x = hash.maxID),
    FUN = function(x) {
      return(paste(sort(x = c(hash.maxID[x], hash.secondID[x], hash.thirdID[x])), collapse = "_"))
    }
  )
  # doublet_names <- names(x = table(doublet_id))[-1] # Not used
  classification <- classification.global
  classification[classification.global == "Negative"] <- "Negative"
  classification[classification.global == "Singlet"] <- hash.maxID[which(x = classification.global == "Singlet")]
  classification[classification.global == "Doublet"] <- doublet_id[which(x = classification.global == "Doublet")]
  classification[classification.global == "Triplet"] <- triplet_id[which(x = classification.global == "Triplet")]
  classification.metadata <- data.frame(
    hash.maxID,
    hash.secondID,
    hash.margin,
    classification,
    classification.global
  )
  colnames(x = classification.metadata) <- paste(
    assay,
    c('maxID', 'secondID', 'margin', 'classification', 'classification.global'),
    sep = '_'
  )
  object <- AddMetaData(object = object, metadata = classification.metadata)
  Idents(object) <- paste0(assay, '_classification')
  # Idents(object, cells = rownames(object@meta.data[object@meta.data$classification.global == "Doublet", ])) <- "Doublet"
  doublets <- rownames(x = object[[]])[which(object[[paste0(assay, "_classification.global")]] == "Doublet")]
  Idents(object = object, cells = doublets) <- 'Doublet'
  # object@meta.data$hash.ID <- Idents(object)
  object$hash.ID <- Idents(object = object)
  return(object)
}


`%||%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(lhs)
  } else {
    return(rhs)
  }
}

MaxN <- function(x, N = 2){
  len <- length(x)
  if (N > len) {
    warning('N greater than length(x).  Setting N=length(x)')
    N <- length(x)
  }
  sort(x, partial = len - N + 1)[len - N + 1]
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




