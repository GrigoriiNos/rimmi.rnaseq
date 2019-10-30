#' Run PAGA analysis on your Seurat_obj2 object
#'
#' This function allows you to run PAGA analysis on your Seurat_obj2 object, visualise and then store in your Seurat_obj object
#'
#' @param Seurat_obj_obj Seurat_obj object
#'
#' @return umap plot, Seurat_obj object with PAGA information
#'
#' @keywords Seurat_obj, single cell sequencing, RNA-seq, pseudotime analysis, PAGA
#'
#' @examples
#'
#'
#' @export


Seurat_obj2_PAGA <- function(Seurat_obj_obj,
                            python_path = "/usr/local/bin/anaconda3/bin/python3.6"){

  library(Seurat)
  library(tidyverse)
  library(reticulate)

  use_python("/usr/local/bin/anaconda3/bin/python3.6")

  an <- import("anndata2ri", convert = FALSE)
  # Activate the anndata2ri conversion between SingleCellExperiment and AnnData
  an$activate

  #Loading the rpy2 extension enables cell magic to be used
  #This runs R code in jupyter notebook cells

  # convert the object
  Seurat_obj_ad <- Convert(
    from = A07,
    to = "anndata",
    filename = "Seurat_obj.h5ad"
  )

  # load scanpy
  sc <- import("scanpy", convert = FALSE)

  sc$logging$print_versions()
  # load the exported anndata file (sc.read())
  adata = sc$read("Seurat_obj.h5ad")

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

  # put a Leiden clustering info to the Seurat_obj object
  Seurat_obj@meta.data$leiden <- py_to_r(adata$obs$leiden)

  # create a list to be stored on your Seurat_obj object with PAGA information
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
  Seurat_obj@misc$paga <- paga

  Seurat_obj
}
