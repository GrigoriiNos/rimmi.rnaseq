#' Run Velocyto analysis on your Seurat2 object
#'
#' This function allows you to un Velocyto analysis on your Seurat2 object and visualise it on your umap embeddings
#'
#' @param Seurat_obj  Seurat object
#' @param loom_path path to the loom file, pre-calculated with python command line tool
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


Seurat2_velocyto <- function(loom_path, Seurat_obj){
  library(velocyto.R)
 # This is generated from the Velocyto python command line tool.
 # You need a loom file before you can proceed
 ldat <- read.loom.matrices(loom_path)

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

 emb <- Seurat_obj@dr$umap@cell.embeddings

 # Estimate the cell-cell distances
 cell.dist <- as.dist(1-armaCor(t(emb)))

 # exonic read (spliced) expression matrix
 emat <- ldat$spliced
 # intronic read (unspliced) expression matrix
 nmat <- ldat$unspliced

 colnames(emat) <- gsub('possorted_genome_bam_[A-Z].*:|x|[A-Z][0-9]{2}:', '', colnames((emat)))
 colnames(nmat) <- gsub('possorted_genome_bam_[A-Z].*:|x|[A-Z][0-9]{2}:', '', colnames((nmat)))

 # subset dead and bad quality cells
 emat <- emat[,colnames(emat) %in% colnames(Seurat_obj@data)]
 nmat <- nmat[,colnames(nmat) %in% colnames(Seurat_obj@data)]

 emat <- emat[,names(Seurat_obj@ident)]
 nmat <- nmat[,names(Seurat_obj@ident)]

  # I'm not sure what this parameter does to be honest. 0.02 default
 # perform gamma fit on a top/bottom quantiles of expression magnitudes
 fit.quantile <- 0.02

 # Main velocity estimation
 rvel.cd <- gene.relative.velocity.estimates(emat,nmat,
                                             deltaT=2,
                                             kCells=10,
                                             cell.dist=cell.dist,
                                             fit.quantile=fit.quantile,
                                             n.cores=24)

 # This section gets the colors out of the seurat tSNE object so that my seurat and velocyto plots use the same color scheme.
 gg <- DimPlot(Seurat_obj,
               reduction.use = "umap",
               do.label = T)

 ggplot_build(gg)$data
 colors <- as.list(ggplot_build(gg)$data[[1]]$colour)
 names(colors) <- rownames(emb)


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
}



