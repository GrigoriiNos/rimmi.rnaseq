#' Run Velocyto analysis on your Seurat2 object
#'
#' This function allows you to get 2d embeddings from your Seurat object that you want to use for velocyto analysis
#'
#'
#' @param Seurat_obj  Seurat object
#' @param emb embedding type of your choice
#'
#' @return umap 2d plot with velocity
#'
#' @keywords Seurat, single cell sequencing, RNA-seq, RNA velocity
#'
#' @examples
#'
#' emb <- get_emb(B06, emb = 'umap')
#'
#' rvel.cd <- preprocess_loom('~/Desktop/RIMMI/Adam_germinal_cenrte/B06/velocyto/B06.loom', B06, emb)
#'
#' Seurat2_velocyto(rvel.cd, B06, emb)
#'
#'
#' @export

get_emb <- function(Seurat_obj, emb = c('umap', 'monocle', 'custom')){
  # take embedding from the Seurat data object
  # NOTE: This assumes you have a seurat data object loaded
  # into R memory prior to using this script. STOP and rerun seurat
  # pipeline if you do not have this loaded. In my case, my seurat object is simply myData
  if (ncol(Seurat_obj@data) == ncol(Seurat_obj@raw.data)){
    stop('Perform cell filtering primarily to velocity analysis')
  }
  if (is.null(Seurat_obj@dr$umap)){
    stop('umap embeddings was npt precalculated for this Seurat object')
  }
  if (is.null(Seurat_obj@ident)){
    stop('clustering was not precalculated for this Seurat object')
  }

  if (emb == 'umap'){
    emb <- Seurat_obj@dr$umap@cell.embeddings

  }
  if (emb == 'monocle'){
    emb <- Seurat_obj@misc$monocle

  }
  emb
}

#' Run Velocyto analysis on your Seurat2 object
#'
#' This function allows you to preprocess loom file into the object that fits to plotting velocyto arrows on 2d embeddings
#'
#'
#' @param Seurat_obj  Seurat object
#' @param loom_path path to the loom file, pre-calculated with python command line tool
#' @param emb the matrix with 2d embeddings that you got from get_emb() function
#'
#' @return umap 2d plot with velocity
#'
#' @keywords Seurat, single cell sequencing, RNA-seq, RNA velocity
#'
#' @examples
#'
#' emb <- get_emb(B06, emb = 'umap')
#'
#' rvel.cd <- preprocess_loom('~/Desktop/RIMMI/Adam_germinal_cenrte/B06/velocyto/B06.loom', B06, emb)
#'
#' Seurat2_velocyto(rvel.cd, B06, emb)
#'
#'
#' @export


preprocess_loom <- function(loom_path, Seurat_obj, emb){
  library(velocyto.R)
  # This is generated from the Velocyto python command line tool.
  # You need a loom file before you can proceed
  ldat <- read.loom.matrices(loom_path)

  # exonic read (spliced) expression matrix
  emat <- ldat$spliced
  # intronic read (unspliced) expression matrix
  nmat <- ldat$unspliced

  emat <- clean_spmat(emat)
  nmat <- clean_spmat(nmat)

  # I'm not sure what this parameter does to be honest. 0.02 default
  # perform gamma fit on a top/bottom quantiles of expression magnitudes
  fit.quantile <- 0.02

  # Estimate the cell-cell distances
  cell.dist <- as.dist(1-armaCor(t(emb)))
  # Main velocity estimation
  rvel.cd <- gene.relative.velocity.estimates(emat,nmat,
                                              deltaT=2,
                                              kCells=10,
                                              cell.dist=cell.dist,
                                              fit.quantile=fit.quantile,
                                              n.cores=24)
  rvel.cd
}

#' Run Velocyto analysis on your Seurat2 object
#'
#' This function allows you to plot velocyto on your 2d embeddings
#'
#'
#' @param rvel.cd the output of preprocess_loom() function
#' @param Seurat_obj  Seurat object
#' @param emb a matrix with 2d embeddings, the oitput of get_emb() function
#'
#' @return umap 2d plot with velocity
#'
#' @keywords Seurat, single cell sequencing, RNA-seq, RNA velocity
#'
#' @examples
#'
#' emb <- get_emb(B06, emb = 'umap')
#'
#' rvel.cd <- preprocess_loom('~/Desktop/RIMMI/Adam_germinal_cenrte/B06/velocyto/B06.loom', B06, emb)
#'
#' Seurat2_velocyto(rvel.cd, B06, emb)
#'
#'
#' @export


Seurat2_velocyto <- function(rvel.cd, Seurat_obj, emb){
  library(velocyto.R)

  emb <- as.data.frame(emb)
  emb <- emb[rownames(Seurat_obj@meta.data),]
  colnames(emb) <- c('UMAP1', 'UMAP2')
  emb$cluster <- Seurat_obj@meta.data$res.0.6
  Seurat_obj@dr$monocle <- NULL

  gg <- ggplot(data = emb,
               aes(x = UMAP1, y = UMAP2, colour = cluster))+
    geom_point()
  emb$cluster <- NULL
  emb <- as.matrix(emb)

# This section gets the colors out of the seurat tSNE object so that my seurat and velocyto plots use the same color scheme.
 ggplot_build(gg)$data
 colors <- as.list(ggplot_build(gg)$data[[1]]$colour)
 names(colors) <- rownames(emb)

 {
 png('My_velocity.png')
 p1 <- show.velocity.on.embedding.cor(emb,
                                      rvel.cd,
                                      n=30,
                                      scale='sqrt',
                                      cell.colors=ac(colors,alpha=0.5),
                                      cex=0.8,
                                      arrow.scale=2,
                                      show.grid.flow=T,
                                      min.grid.cell.mass=1.0,
                                      grid.n=50,
                                      arrow.lwd=1,
                                      do.par=F,
                                      cell.border.alpha = 0.1,
                                      n.cores=24,
                                      main="Cell Velocity")
  dev.off()
  }

}

