#' Splits object into a list of subsetted objects.
#'
#' Splits object based on a single attribute into a list of subsetted objects,
#' one for each level of the attribute. For example, useful for taking an object
#' that contains cells from many patients, and subdividing it into
#' patient-specific objects.
#'
#' @param object Seurat object
#' @param split.by Attribute for splitting. Default is "ident". Currently
#' only supported for class-level (i.e. non-quantitative) attributes.
#'
#' @return A named list of Seurat objects, each containing a subset of cells
#' from the original object.
#'
#' @export
#'
#' @examples
#' # Assign the test object a three level attribute
#' groups <- sample(c("group1", "group2", "group3"), size = 80, replace = TRUE)
#' names(groups) <- colnames(pbmc_small)
#' pbmc_small <- AddMetaData(object = pbmc_small, metadata = groups, col.name = "group")
#' obj.list <- SplitObject(pbmc_small, split.by = "group")
#'
SplitObject <- function(object, split.by = "ident") {
  if (split.by == 'ident') {
    groupings <- Idents(object = object)
  } else {
    groupings <- FetchData(object = object, vars = split.by)[, 1]
  }
  groupings <- unique(x = as.character(x = groupings))
  obj.list <- list()
  for (i in groupings) {
    if (split.by == "ident") {
      obj.list[[i]] <- SubsetData(object, idents = i)
    }
    else {
      cells <- which(x = object[[split.by, drop = TRUE]] == i)
      cells <- colnames(x = object)[cells]
      obj.list[[i]] <- SubsetData(object, cells = cells)
    }
  }
  return(obj.list)
}


#' Prepare 10x v3 out files files for Seurat 2
#'
#'
#' @param path path to filtered gene bc matrix
#'
#'
#' @keywords Seurat, single cell sequencing, RNA-seq, 10x
#'
#' unarchive('"~/.../filtered_feature_bc_matrix/"')
#'
#' @export

unarchive_10x <- function(path){
  setwd(paste0(path, 'filtered_feature_bc_matrix/'))
  {
    sink('unarchive.sh')
    cat('#!/bin/bash \n')
    cat(paste0('cd ', path, 'filtered_feature_bc_matrix/ \n'))
    cat('ls | grep -e .gz | xargs gunzip \n')
    cat('mv features.tsv genes.tsv')
    sink()
  }

  system('chmod +x unarchive.sh')
  system('./unarchive.sh')

  system('rm unarchive.sh')
}



#' Merge two loom files into one and produce rvel.cd
#'
#' This function allows you to merge emat and nmat matrices into one and then calculate rvel.cd for velocyto analysis
#'
#'
#' @param loom1 path to loom file 1
#' @param loom2 path to loom file 2
#' @param sample1 sample name 1, for the prefix in cellular barcodes
#' @param sample2 sample name 2, for the prefix in cellular barcodes
#' @param emb 2d embedding
#'
#' @return umap 2d plot with velocity
#'
#' @keywords Seurat, single cell sequencing, RNA-seq, RNA velocity
#'
#' @examples
#'
#'
#'
#' @export

merge_loom <- function(loom1,
                       loom2,
                       Seurat_obj,
                       sample1,
                       sample2,
                       emb){
  library(velocyto.R)
  library(Matrix.utils)

  ldat1 <- read.loom.matrices(loom1)
  # exonic read (spliced) expression matrix
  emat1 <- ldat1$spliced
  # intronic read (unspliced) expression matrix
  nmat1 <- ldat1$unspliced

  ldat2 <- read.loom.matrices(loom2)
  # exonic read (spliced) expression matrix
  emat2 <- ldat2$spliced
  # intronic read (unspliced) expression matrix
  nmat2 <- ldat2$unspliced

  colnames(emat1) <- paste0(sample1, '_', colnames(emat1))
  colnames(nmat1) <- paste0(sample1, '_', colnames(nmat1))
  colnames(emat2) <- paste0(sample2, '_', colnames(emat2))
  colnames(nmat2) <- paste0(sample2, '_', colnames(nmat2))

  if (all(rownames(emat1) %in% rownames(emat2)) & all(rownames(emat1) %in% rownames(emat2))){
    print('starting merging')

    emat <- merge.Matrix(emat1, emat2, by.x = rownames(emat1), by.y = rownames(emat2))
    nmat <- merge.Matrix(nmat1, nmat2, by.x = rownames(nmat1), by.y = rownames(nmat2))
  } else {
    stop('features are not matching')
  }
  rm(ldat1, ldat2, emat1, emat2, nmat1, nmat2)

  emat <- clean_spmat(emat, Seurat_obj)
  nmat <- clean_spmat(nmat, Seurat_obj)

  cell.dist <- as.dist(1-armaCor(t(emb)))
  fit.quantile <- 0.02
  # Main velocity estimation
  rvel.cd <- gene.relative.velocity.estimates(emat,nmat,
                                              deltaT=2,
                                              kCells=10,
                                              cell.dist=cell.dist,
                                              fit.quantile=fit.quantile,
                                              n.cores=24)
}

#' Run Velocyto analysis on your Seurat2 object
#'
#' This function allows you to get 2d embeddings from your Seurat object that you want to use for velocyto analysis
#'
#'
#' @param mat matrix
#' @param emb Seurat_obj Seurat object
#'
#' @return cleaned matrix
#'
#' @keywords Seurat, single cell sequencing, RNA-seq, RNA velocity
#'
#' @examples
#'
#'
#'
#' @export

clean_spmat <- function(mat, Seurat_obj, emb, sample){
  print('trimming cellular barcodes')

  colnames(mat) <- gsub('possorted_genome_bam_.*:|x|[A-Z][0-9]{2}:', '', colnames((mat)))

  colnames(mat) <- paste0('A',gsub('.+-|\\.loom','', sample), '_', colnames(mat))

  print(paste0('filtering ', sum(duplicated(rownames(mat))), ' dublicated genes'))
  mat <- mat[!duplicated(rownames(mat)),]

  notnas <- colSums(is.na(mat))==0
  mat <- mat[, notnas]

  #colnames(mat) <- gsub('SJ-2365-Adam-gex-', 'A', colnames(mat))
  #colnames(mat) <- gsub(':', '_', colnames(mat))
  #colnames(mat) <- gsub('x', '', colnames(mat))
  # subset dead and bad quality cells
  mat <- mat[,colnames(mat) %in% rownames(emb)]

  mat
}

