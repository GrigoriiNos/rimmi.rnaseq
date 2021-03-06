% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat2_monole3.R
\name{Seurat2_monocle3}
\alias{Seurat2_monocle3}
\title{Run Velocyto analysis on your Seurat2 object}
\usage{
Seurat2_monocle3(
  Seurat_obj,
  clus_var = c("res.0.6", "res.0.8", "res.1", "celltype"),
  multi_samples = T,
  samples_var = c("condition", "sample"),
  cell_cycle = F
)
}
\arguments{
\item{Seurat_obj}{Seurat object}

\item{loom_path}{path to the loom file, pre-calculated with python command line tool}
}
\value{
umap 2d plot with velocity
}
\description{
This function allows you to un Velocyto analysis on your Seurat2 object and visualise it on your umap embeddings,
you need to have fully pre-proccessed Seurat object with QC done, 2d embeddings and clustering pre-calculated,
}
\examples{



}
\keyword{RNA}
\keyword{RNA-seq,}
\keyword{Seurat,}
\keyword{cell}
\keyword{sequencing,}
\keyword{single}
\keyword{velocity}
