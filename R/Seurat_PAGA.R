#' Run PAGA analysis on your Seurat2 object
#'
#' This function allows you to run PAGA analysis on your Seurat2 object, visualise and then store in your Seurat object
#'
#' @param Seurat_obj Seurat object
#'
#' @return umap plot, Seurat object with PAGA information
#'
#' @keywords Seurat, single cell sequencing, RNA-seq, pseudotime analysis, PAGA
#'
#' @examples
#'
#'
#' @export


Seurat2_PAGA <- function(Seurat_obj){
  library(Seurat)
  library(tidyverse)
  library(reticulate)

  # convert the object
  seurat_ad <- Convert(
    from = seurat,
    to = "anndata",
    filename = "seurat.h5ad"
  )

  # load scanpy
  sc <- import("scanpy", convert = FALSE)

  # load the exported anndata file (sc.read())
  adata = sc$read("seurat.h5ad")

  # find neighbors (sc.pp.neighbors())
  sc$pp$neighbors(adata)

  # cluster cells using the Leiden algorithm (sc.tl.leiden())
  sc$tl$leiden(adata, resolution = 1.0)

  # run PAGA analysis (sc.tl.paga()) and give it a visual check (sc.pl.paga())
  sc$tl$paga(adata, groups = "leiden")
  sc$pl$paga(adata)

  # generate a UMAP using the PAGA results as initial positions (sc.tl.umap())
  sc$tl$umap(adata, init_pos = "paga")

  # plot UMAP and PAGA side-by-side (sc.pl.paga_compare())
  sc$pl$paga_compare(
    adata,
    basis = "umap",
    threshold = 0.15,
    edge_width_scale = 0.5,
    save = TRUE
  )

  # put a Leiden clustering info to the seurat object
  seurat@meta.data$leiden <- py_to_r(adata$obs$leiden)

  # create a list to be stored on your seurat object with PAGA information  
  paga <- list(
    connectivities = py_to_r(adata$uns$paga$connectivities$todense()),
    connectivities_tree = py_to_r(adata$uns$paga$connectivities$todense()),
    group_name = py_to_r(adata$uns$paga$groups),
    groups = levels(py_to_r(adata$obs$leiden)),
    group_colors = setNames(py_to_r(adata$uns$leiden_colors), c(0:(nrow(py_to_r(adata$uns$paga$pos))-1))),
    position = tibble(
      group = levels(py_to_r(adata$obs$leiden)),
      x = as.data.frame(py_to_r(adata$uns$paga$pos))$V1,
      y = as.data.frame(py_to_r(adata$uns$paga$pos))$V2
    ),
    umap = tibble(
      UMAP_1 = as.data.frame(py_to_r(adata$obsm$X_umap))$V1,
      UMAP_2 = as.data.frame(py_to_r(adata$obsm$X_umap))$V2
    )
  )

  rownames(paga$connectivities) <- c(1:nrow(paga$pos))
  colnames(paga$connectivities) <- c(1:nrow(paga$pos))

  # store it in the object
  seurat@misc$paga <- paga

  seurat
}