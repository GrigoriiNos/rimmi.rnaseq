% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Harmony_Seurat.R
\name{scanorama_seurat}
\alias{scanorama_seurat}
\title{use the batch effect correction method scanorama on Seurat objects}
\usage{
scanorama_seurat(merged_object, raw.data = T)
}
\arguments{
\item{merged_object}{pre-merged Seurat object}
}
\value{
Seurat object with BE corrected expression matrix in the data slot
}
\description{
This function allows you to use an additional BE correction method still using Seurat objects
}
\examples{

mnn_for_seurat(B04_B05)
B05@meta.data[, "sample"] <- "B05"
B04@meta.data[, "sample"] <- "B04"

add prefix to cellnames
B05 <- RenameCells(B05,
                  add.cell.id = 'B05')

B04 <- RenameCells(B04,
                  add.cell.id = 'B04')

B04_B05 <- MergeSeurat(B04, B05)

}
\keyword{RNA-seq,}
\keyword{Seurat,}
\keyword{cell}
\keyword{gene}
\keyword{sequencing,}
\keyword{signature}
\keyword{single}
