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



