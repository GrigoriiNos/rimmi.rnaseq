#' Merge clonotype and isotype data to GEX Seurat object
#'
#' This function allows you to fuse the output of VDJ library to Seurat object
#'
#' @param vdj_location location of cel ranger vdj output
#'
#' @return table with vdj info
#'
#' @keywords Seurat, single cell sequencing, RNA-seq, RNA velocity
#'
#' @examples
#'
#' A07 <- add_vdj('path_to_vdj', A07)
#'
#'
#' @export

pull_vdj <- function(vdj_location){
  vdj <- read.csv(paste(vdj_location,"filtered_contig_annotations.csv", sep=""))

  # Remove the -1 at the end of each barcode.
  # Subsets so only the first line of each barcode is kept,
  # as each entry for given barcode will have same clonotype.
  vdj$barcode <- gsub("-1", "", vdj$barcode)
  vdj <- vdj[!duplicated(vdj$barcode), ]

  # Only keep the barcode, isotype and clonotype columns.
  # We'll get additional clonotype info from the clonotype table.
  vdj <- vdj[,c("barcode", "raw_clone_id", "ccall", "vcall", "dcall", "jcall")]
  names(vdj)[names(vdj) == "raw_clone_id"] <- "clone_id"

  # Clonotype/isotype-centric info
  clono <- read.csv(paste(vdj_location,"clonotypes.csv", sep=""))

  # Slap the AA sequences onto our original table by clone_id.
  vdj <- merge(vdj, clono[, c("clone_id", "cdr3s_aa", "cdr3s_nt")])

  # Reorder so barcodes are first column and set them as rownames.
  vdj <- vdj[, c("barcode", "clone_id", "ccall", "vcall", "dcall", "jcall", "cdr3s_aa", "cdr3s_nt")]
  rownames(vdj) <- vdj[,'barcode']
  vdj$barcode <- NULL
  vdj
}


#' Merge clonotype and isotype data to GEX Seurat object
#'
#' This function allows you to fuse the output of VDJ library to Seurat object
#'
#' @param vdj_location location of cel ranger vdj output
#' @param Seurat_obj  Seurat object
#'
#' @return umap 2d plot with velocity
#'
#' @keywords Seurat, single cell sequencing, RNA-seq, RNA velocity
#'
#' @examples
#'
#' A07 <- add_vdj('path_to_vdj', A07)
#'
#'
#' @export

add_vdj <- function(vdj_location, Seurat_obj){
  vdj <- pull_vdj(vdj_location)

  # Add to the Seurat object's metadata.
  Seurat_obj@meta.data <- merge(Seurat_obj@meta.data,
                                vdj,
                                by = 0,
                                all.x = T)
  rownames(Seurat_obj@meta.data) <- Seurat_obj@meta.data$Row.names
  Seurat_obj@meta.data$Row.names <- NULL

  return(Seurat_obj)
}


#' Gives you a summary of your clonotypes
#'
#' @param Seuraat_obj Seurat object
#' @param clusters_column column in the @meta.data that stores your clusters
#'
#' @return data.frame with clonotypes summary
#'
#' @export

clonotypes_summary <- function(Seurat_obj, clusters_column){
  clonometa <- Seurat_obj@meta.data[!is.na(Seurat_obj@meta.data$cdr3),]

  colnames(clonometa)[colnames(clonometa) == clusters_column] <- 'cell.types'

  clonometa$c_call <- factor(clonometa$c_call)
  clonometa$GC <- factor(clonometa$GC)
  clonometa$cell.types <- factor(clonometa$cell.types)

  clonometa <- clonometa[,c('clone_id',
                            'GC',
                            'cell.types',
                            'c_call'
  )]

  summarise_clone <- function(df){
    cbind(
      df$clone_id %>% table %>% as.data.frame %>% column_to_rownames('.'),
      df$c_call %>% table %>% as_tibble() %>% column_to_rownames('.') %>% t,
      df$GC %>% table %>% as_tibble() %>% column_to_rownames('.') %>% t,
      df$cell.types %>% table %>% as_tibble() %>% column_to_rownames('.') %>% t
    )
  }

  clonotypes <- data.frame()
  for (clone in unique(clonometa$clone_id)){
    df <- clonometa[clonometa$clone_id == clone,]

    df <- summarise_clone(df)

    clonotypes <- rbind(clonotypes, df[,colnames(df) != "V1"])
  }

  clonotypes <- clonotypes[order(clonotypes$Freq, decreasing = T),]
  clonotypes
}

#' Rank clonotypes based of the spread in GEX space
#'
#' @param Seuraat_obj Seurat object
#' @param cells_th minimum cells in the clonotype to analyse
#'
#' @return data.frame ranked clonotypes
#'
#' @export


rank_clonotypes <- function(Seurat_obj, cells_th=10){

  freq.cl <- table(Seurat_obj@meta.data$clone_id)
  freq.cl <- freq.cl[freq.cl > 2]
  big.clones <- freq.cl[freq.cl > cells_th]

  mean.dist <- c()
  max.dist <- c()

  Seurat_obj <- SubsetData(Seurat_obj, cells = complete.cases(Seurat_obj$clone_id))

  for (clone_id in names(big.clones)){
    cells <- Cells(Seurat_obj)[Seurat_obj$clone_id == clone_id]

    distances <- dist(Seurat_obj@reductions$pca@cell.embeddings[cells,])

    mean.d <- mean(distances)
    names(mean.d) <- clone_id

    max.d <- max(distances)
    names(max.d) <- clone_id

    mean.dist <- c(mean.dist, mean.d)
    max.dist <- c(max.dist, max.d)
  }

  distances <- data.frame(max.dist = max.dist,
                          mean.dist = mean.dist,
                          row.names = names(max.dist))


  distances[order(-distances$max.dist, -distances$mean.dist), ]
}
