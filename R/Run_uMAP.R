#' Calculate uMAP for Seurat object
#'
#' This function allows you in addition to canonical dimentionality reductions approache as tSNE calculate also uMAP and embedd it to your Seurat object in the dr slot, can be used for visualisation of the cell clusters.
#'
#' @param Seurat_obj your Seurat object
#' @param dims how many of pca dimentions you want to use, 10 is default
#' @param method method for umap manifold calculation, can be, "naive" or "umap-learn". Naive method is default. Further reading is in ??umap::umap
#' @param CCA set TRUE of you want to build uMAPs on the CCAs dimentionality rediction
#'
#' @return Seurat object with new dimentionality reduction slot
#'
#' @keywords Seurat, single cell sequencing, RNA-seq, uMAP, dimentionality reduction
#'
#' @examples
#'
#' Seurat_obj <- Run_uMAP(Seurat_obj)
#' uMAPPlot(Seurat_obj, sample_name = 'b06')
#'
#' @export
#'

Run_uMAP <- function(Seurat_obj,
                     dims = 10,
                     method = 'naive',
                     CCA = F)
  {
  library(umap)
  library(Seurat)

  ### get pca/cca of expression matrix
  if (CCA == F)
    {
    print('calculating PCA')
    dr <- GetDimReduction(object = Seurat_obj,
                         reduction.type = 'pca',
                         slot = 'cell.embeddings')
    } else {
    print('calculating CCA')
    dr <- GetDimReduction(object = Seurat_obj,
                          reduction.type = 'cca',
                          slot = 'cell.embeddings')
    }

  print('started to calculate uMAPs, can take up to 1 hour')
  ### calculate uMAPs
  Seurat_obj_map <- umap(
    dr[,1:dims],
    method = method
    )

  ### put umap into the new slot
  Seurat_obj@dr$umap <- Seurat_obj_map

  print('done')
  print('you now have your Seurat object with a new dimentionality rediction slot')
  Seurat_obj
  }
