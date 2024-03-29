% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{all.genes.withGO}
\alias{all.genes.withGO}
\title{Get all genes with partiular GO}
\usage{
\method{all}{genes.withGO}(GO)
}
\arguments{
\item{GO}{a GO term that you want to get all genes with}
}
\value{
vector of genes that are involved in GO that was put as an argument
}
\description{
This function allows you to get a list of genes from the biomaRt database that have a GO of your choice.
}
\examples{

all.genes.withGO("GO:0008217")

}
\keyword{GO}
\keyword{GO,}
\keyword{all}
\keyword{genes}
\keyword{with}
