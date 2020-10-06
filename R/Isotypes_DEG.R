#' This function allows you to calculate DEG for the IG Isotypes
#'
#' @param Seurat_obj Seurat object with integrated data from VDJ analysis with Isotype assigment
#'
#' @keywords DEG, single cell, VDJ analysis, B cells, T cells
#'
#' @examples
#'
#' Isotypes_DEG(A07)
#'
#' @export
#'

Isotypes_DEG <- function(Seurat_obj){
  library(Seurat)

  Seurat_obj@meta.data <- Seurat_obj@meta.data[names(Idents(Seurat_obj)),]

  Seurat_obj$c_call_short <- gsub('[0-9]', '', Seurat_obj$c_call)
  # replace cell identity from cluster to isotype
  Idents(Seurat_obj) <- factor(Seurat_obj@meta.data$c_call_short)
  names(Idents(Seurat_obj)) <- rownames(Seurat_obj@meta.data)

  short_isotypes <- grep('IGH[A-Z]', levels(Idents(Seurat_obj)), value = T)

  Seurat_obj_sub <- subset(Seurat_obj,
                           idents = short_isotypes)

  print('calculating markers for short isotypes')
  print(grep('IGH[A-Z].*', levels(Idents(Seurat_obj)), value = T))
  Short_iso_markers <- FindAllMarkers(object = Seurat_obj_sub,
                                      test.use = 'poisson', # another method
                                      latent.vars = 'orig.ident', # adjust for sample
                                      only.pos = F,
                                      min.pct = 0.25,
                                      thresh.use = 0.25)

  rimmi.rnaseq::markers_to_xls(Short_iso_markers, filename = 'Short_isotypes_DEGs')
  print('saved short isotypes as an excel spreadsheet')

  # replace cell identity from cluster to isotype
  Idents(Seurat_obj) <- factor(Seurat_obj@meta.data$c_call)
  names(Idents(Seurat_obj)) <- rownames(Seurat_obj@meta.data)

  full_isotypes <- grep('IGH[A-Z].*', levels(Idents(Seurat_obj)), value = T)

  Seurat_obj_sub <- subset(Seurat_obj,
                           idents = full_isotypes)

  print('calculating markers for full isotypes')
  print(grep('IGH[A-Z].*', levels(Idents(Seurat_obj)), value = T))
  Full_iso_markers <- FindAllMarkers(object = Seurat_obj_sub,
                                     test.use = 'poisson', # another method
                                     latent.vars = 'orig.ident', # adjust for sample
                                     only.pos = F,
                                     min.pct = 0.25,
                                     thresh.use = 0.25)

  rimmi.rnaseq::markers_to_xls(Full_iso_markers, filename = 'Full_isotypes_DEGs')
  print('saved full isotypes as an excel spreadsheet')

  print(
    Full_iso_markers %>%
      group_by(cluster) %>%
      top_n(5, avg_logFC)
  )
}

#' This function allows you to calculate DEG for some specific IG Isotypes
#'
#' @param Seurat_obj Seurat object with integrated data from VDJ analysis with Isotype assigment
#'
#' @keywords DEG, single cell, VDJ analysis, B cells, T cells
#'
#' @examples
#'
#' Isotypes_DEG(A07)
#'
#' @export
#'

Specific_isotypes_DEG <- function(Seurat_obj,
                                  iso1,
                                  iso2,
                                  filename,
                                  clusters = NULL){
  library(Seurat)

  Seurat_obj@meta.data <- Seurat_obj@meta.data[names(Idents(Seurat_obj)),]

  if (!is.null(clusters)){
    Seurat_obj <- subset(Seurat_obj,
                         idents = clusters)
  }
  # replace cell identity from cluster to isotype
  Idents(Seurat_obj) <- factor(Seurat_obj@meta.data$c_call)
  names(Idents(Seurat_obj)) <- rownames(Seurat_obj@meta.data)

  isotypes <- grep('IGH[A-Z].*', levels(Idents(Seurat_obj)), value = T)

  Seurat_obj_sub <- SubsetData(object = Seurat_obj,
                               ident.use = isotypes
  )
  print('calculating markers for specific isotypes')
  print(grep('IGH[A-Z].*', levels(Idents(Seurat_obj)), value = T))
  iso_markers <- FindMarkers(object = Seurat_obj_sub,
                             test.use = 'poisson', # another method
                             latent.vars = 'orig.ident', # adjust for sample
                             ident.1 = iso1,
                             ident.2 = iso2,
                             only.pos = F,
                             min.pct = 0.25,
                             thresh.use = 0.25)

  write.csv(iso_markers,
            file = filename,
            quote = F,
            row.names = T)
  print('saved isotypes as a csv file')
}


#' This function allows you to pull specific barcodes from the object with specific IG Isotypes from specific clusters
#'
#'
#' @param Seurat_obj Seurat object with integrated data from VDJ analysis with Isotype assigment
#' @param isotypes Isotype names
#' @param clusters cluster names
#'
#' @keywords VDJ analysis, B cells, T cells
#'
#' @examples
#'
#' pull_isotype(A07, 'IGHA1')
#'
#' @export
#'

pull_isotype <- function(Seurat_obj,
                         isotypes,
                         clusters = NULL){

  Seurat_obj@meta.data <- Seurat_obj@meta.data[names(Idents(Seurat_obj)),]

  if (is.null(clusters)){
  iso.cells <- rownames(
    Seurat_obj@meta.data)[Seurat_obj@meta.data$c_call %in% isotypes]
  } else{
    iso.cells <- rownames(
      Seurat_obj@meta.data)[(Seurat_obj@meta.data$c_call %in% isotypes) & (Idents(Seurat_obj) %in% clusters)]
  }
  iso.cells
}


#' This function allows you to pull specific barcodes from the object with specific clonotypes from specific clusters
#'
#' @param Seurat_obj Seurat object with integrated data from VDJ analysis with Isotype assigment
#' @param clonotypes clonotype names
#' @param clusters cluster names
#'
#' @keywords VDJ analysis, B cells, T cells
#'
#' @examples
#'
#' pull_clonotype(A07, 'IGHA1')
#'
#' @export
#'

pull_clonotype <- function(Seurat_obj,
                           clonotypes,
                           clusters = NULL){

  Seurat_obj@meta.data <- Seurat_obj@meta.data[names(Idents(Seurat_obj)),]

  if (is.null(clusters)){
    clono.cells <- rownames(
      Seurat_obj@meta.data)[Seurat_obj@meta.data$clone_id %in% clonotypes]
  } else{
    clono.cells <- rownames(
      Seurat_obj@meta.data)[(Seurat_obj@meta.data$clone_id %in% clonotypes) & (Idents(Seurat_obj) %in% clusters)]
  }
  clono.cells
}

