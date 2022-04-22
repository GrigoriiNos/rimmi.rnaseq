#' Annotation for the cell markers.
#'
#' This function allows you to get multiple functional annotations for DE genes for each cell cluster using a Seurat::FindAllMarkers function output
#'
#' @param markers_table The output of the Seurat::FndAllMarkers function, don't change anything in the column names within it.
#'
#' @param X How much of the genes you need to analyse in each cluster. 50 is a default.
#'
#' @param organism What species you work with? "mmusculus" and "hsapiens" is available in gProfiler. Human is a default.
#'
#' @return table with annotation for each cluster
#'
#' @keywords GO, KEGG, HPA, single cell, RNA-seq
#'
#' @examples
#'
#' annotate_markers(stromal_markers)
#'
#' @export
#'

annotate_markers <- function(markers_table,
                             X = 50,
                             method = 'gprofiler',
                             organism = 'hsapiens',
                             reg = 'up'){

  library(gProfileR)
  library(enrichR)

  enrich <- function(df, method=method){
    gprofiler1 <- function(df){
      x <- gprofiler(df$gene, organism = organism)
      x$cluster <- df$cluster[1]
      x
    }

    enrichr1 <- function(df){
      dbs <- c("GO_Biological_Process_2018",
               "GO_Molecular_Function_2018", "Transcription_Factor_PPIs",
               "Human_Gene_Atlas", "Ligand_Perturbations_from_GEO_down",
               "KEGG_2019_Human","Reactome_2016", "WikiPathways_2016",
               "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO",
               "TRRUST_Transcription_Factors_2019")

      geneset <- df$gene

      enriched <- enrichr(toupper(geneset), dbs)

      enriched <- lapply(1:length(enriched), function(i){
        enriched[[i]]$database <- dbs[i]

        enriched[[i]]
      }
      )

      enriched <- bind_rows(enriched)
      enriched$cluster <- df$cluster[1]

      enriched %>% filter(str_count(Genes, ';') != 0)
    }
    if (method == 'gprofiler'){
      gprofiler1(df)
    } else if (method == 'enrichr'){
      enrichr1(df)
    }
  }

  ###
  if (reg == 'up'){
    q <- markers_table %>%
      group_by(cluster) %>%
      filter(p_val_adj < 0.05 & avg_logFC > 0) %>%
      top_n(X, abs(avg_logFC)) %>%
      dplyr::select(gene)
  } else if (reg == 'down'){
    q <- markers_table %>%
      group_by(cluster) %>%
      filter(p_val_adj < 0.05 & avg_logFC < 0) %>%
      top_n(X, abs(avg_logFC)) %>%
      dplyr::select(gene)
  }

  report <- tibble()
  for (cl in unique(q$cluster)){
    try(
      report <- rbind(report, enrich(df = q %>% filter(cluster == cl) #%>% pull(gene)
                                     ,
                                     method=method))
    )
  }
  return(report)
}



#' getting your markers into the excel shpreadsheet
#'
#' This function allows you write Seurat::FindAllMarkers function output into the excel file having each cluster markers written in the distinct sheet
#'
#' @param markers_table The output of the Seurat::FndAllMarkers function, don't change anything in the column names within it.
#'
#' @param filename
#'
#' @keywords excel, single cell, RNA-seq, Seurat, markers
#'
#' @examples
#'
#' annotate_markers(my_markers, "stromal_markers")
#'
#' @export
#'

markers_to_xls <- function(markers_table, filename = 'my_markers'){
  library(dplyr)
  ###
  clusters <- unique(markers_table$cluster)
  tables <- list()
  for (i in 1:length(clusters)){
    sub_markers_table <- markers_table %>%
      filter(cluster == clusters[i])
    tables[[i]] <- sub_markers_table
  }
  WriteXLS::WriteXLS(x = tables,
                     ExcelFileName = paste0(filename, '.xls'),
                     SheetNames = paste(clusters, 'cluster')
  )
}
