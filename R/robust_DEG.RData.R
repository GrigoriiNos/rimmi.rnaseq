#' apply DEG with multiple glms on Seruat object with case control meta data
#'
#' This function allows you to calculate an average expression for the list of genes of your interest.
#'
#' @param Seurat_obj your Seurat object
#' @param model type of the model to choose
#'
#' @return DEG table
#'
#' @keywords GO, KEGG, HPA, single cell, RNA-seq
#'
#' @export


Hglm.deg <- function(Seurat_obj,
                     condition = 'condition',
                     cell.types = 'RNA_snn_res.0.5',
                     model = c('glm0', 'glm', 'glmer')){

  library(tidyverse)
  library(lme4)
  library(lmtest)
  library(ggrepel)
  library(ggplot2)

  metadata <- Seurat_obj@meta.data[,c('nCount_RNA',
                                      'orig.ident',
                                      condition,
                                      cell.types)]

  gene_data <- t(Seurat_obj@assays$RNA@counts)
  gene_data <- gene_data[, colnames(gene_data) %in% VariableFeatures(Seurat_obj)]

  genes <- VariableFeatures(Seurat_obj)[1:50]
  genes <- union(genes, c('CCL19','CXCL13','CR2', 'CCL21A', 'INMT', 'STMN2', 'GAS6',
                              'POSTN', 'RBP1', 'LRAT', 'ACTA2', 'MADCAM1', 'MAF', 'IL7',
                              'TAGLN', 'CD74','CCL2','GREM1','BGN','RPB1', 'SERPING1', 'C3', 'MAF'))
  # taking gene set so far to get it done quickly
  gene_data <- gene_data[,colnames(gene_data) %in% genes]

  metadata <- metadata[rownames(metadata) %in% rownames(gene_data),]

  # Function for fitting the two models to each gene and doing an LR test.

  # multilevel glm
  fit_lme <- function(expr_vec) {
    tmp_df <- cbind(metadata, data.frame(expr = unname(expr_vec)))

    mres1 <- glmer(f1, data = tmp_df, family = poisson())
    mres2 <- glmer(f2, data = tmp_df, family = poisson())
    anova_res <- anova(mres1, mres2)
    pval <- anova_res$`Pr(>Chisq)`[2]
    effect_size <- -(ranef(mres2)$condition[,1][1] - ranef(mres2)$condition[,1][2])

    c(pval = pval, effect_size = effect_size)
  }

  # commom glm function
  fit_glm <- function(expr_vec) {
    tmp_df <- cbind(metadata, data.frame(expr = unname(expr_vec)))

    mres1 <- glm(f1, data = tmp_df, family = poisson())
    mres2 <- glm(f2, data = tmp_df, family = poisson())

    lrt_res <- lrtest(mres1, mres2)
    pval <- lrt_res$`Pr(>Chisq)`[2]
    effect_size <- unname(coef(mres2)[2])

    c(pval = pval, effect_size = effect_size)
  }


  if (model == 'glm'){
    f1 <- 'expr ~ offset(log(nCount_RNA)) + sample'
    f2 <- 'expr ~ offset(log(nCount_RNA)) + condition + sample'
    e_fit_glm <- function(expr_vec) {
      tryCatch(fit_glm(expr_vec), error = function(e) c(pval = NA, effect_size = NA))
    }
    fit_results <- apply(gene_data, 2, e_fit_glm)
  } else if (model == 'glm0'){
    f1 <- 'expr ~ offset(log(nCount_RNA))'
    f2 <- 'expr ~ offset(log(nCount_RNA)) + condition'
    e_fit_glm <- function(expr_vec) {
      tryCatch(fit_glm(expr_vec), error = function(e) c(pval = NA, effect_size = NA))
    }
    fit_results <- apply(gene_data, 2, e_fit_glm)
  } else if (model == 'glmer'){
    f1 <- 'expr ~ offset(log(nCount_RNA)) + (1 | sample)'
    f2 <- 'expr ~ offset(log(nCount_RNA)) + (1 | condition/sample)'
    e_fit_lme <- function(expr_vec) {
      tryCatch(fit_lme(expr_vec), error = function(e) c(pval = NA, effect_size = NA))
    }
    fit_results <- apply(gene_data, 2, e_fit_lme)
    }
  # return results
  return(
    as.tibble(cbind(t(fit_results),
                    data.frame(gene = rownames(t(fit_results)))))
    )
}

#' apply DESeq for DEG pseudo bulk collapsed samples
#'
#' This function allows you to calculate an average expression for the list of genes of your interest.
#'
#' @param Seurat_obj your Seurat object
#' @param thresh.pv display genes with at least this p value
#' @param show.genes additional genes to display
#'
#' @return DEG table
#'
#' @keywords GO, KEGG, HPA, single cell, RNA-seq
#'
#' @export

pseudo.bulk.DEG <- function(Seurat_sub.obj,
                            sample = 'orig.ident',
                            condition = 'condition',
                            fitType = 'local',
                            show.genes = c('CCL19', 'CCL21A')){

  collapsed.list <- lapply(Seurat::SplitObject(Seurat_sub.obj, sample),
                           function(sample)
                             {rowSums(as.matrix(sample@assays$RNA@counts))
                             }
                           )

  design <- unique(
    Seurat_sub.obj@meta.data[,c(sample, condition)]
  )

  counts <- as.data.frame(collapsed.list)

  counts <- counts[,design[sample][,1]] %>%
    rownames_to_column() %>%
    filter(rowname %in% VariableFeatures(Seurat_sub.obj))

  rownames(design) <- 1:nrow(design)


  dds <- DESeqDataSetFromMatrix(countData=counts,
                                colData=design,
                                design=model.matrix(~design[condition][,1]),
                                tidy = TRUE)


  dds <- DESeq(dds,
               fitType = fitType
               )

  res <- as.data.frame(results(dds)) %>%
    rownames_to_column('gene')

  res <- res[,c('gene', 'log2FoldChange','pvalue',  'padj')]
  colnames(res) <- c('gene', 'effect_size','pvalue',  'pval')

  return(res)
}

#' unite list of deg table
#'
#' better version volcano plot
#'
#' @param deg.table deg table
#'
#' @return deg table united
#'
#' @keywords DEG, differential expression analysis
#'
#' @export

list.to.table <- function(deg.table){

  if (is.null(names(deg.table))){
    deg.table <- lapply(deg.table, function(df) {
      df <- df[order(abs(df$effect_size), decreasing = T),]
      df$pval[is.na(df$pval)] <- 1
      df$sign <- ifelse(df$pval < 0.05, '**','')
      df
    })

    deg.tables <- do.call('bind_rows', deg.table)
    deg.tables$cluster <- rep(0:(length(deg.table)-1), each = nrow(deg.table[[1]]))
    deg.tables
  } else {
    deg.table <- lapply(deg.table, function(df) {
      df <- df[order(abs(df$effect_size), decreasing = T),]
      df$pval[is.na(df$pval)] <- 1
      df$sign <- ifelse(df$pval < 0.05, '**','')
      df
    })

    deg.tables <- do.call('bind_rows', deg.table)
    deg.tables$cluster <- rep(names(deg.table), each = nrow(deg.table[[1]]))
    deg.tables
  }
}

#' plot volcano plot 2
#'
#' better version volcano plot
#'
#' @param response table with DEG
#' @param  ES_th ES
#' @param title the name of your plot
#'
#' @return volcano plot
#'
#' @keywords DEG, differential expression analysis, volcano plot
#'
#' @examples
#'
#' volcano_pl(responce)
#'
#' @export

volcano_pl2 <- function(response, ES_th = 0.5, title){
  table <- response %>%
    mutate(`-log10p` = -log10(pval)) %>%
    mutate(conclusion = ifelse(response$pval < 0.05 & response$effect_size > ES_th,
                               'sign_up',
                               ifelse(response$pval < 0.05 & response$effect_size < -ES_th,
                                      'sign_down',
                                      'not_sign'))) %>%
    filter(!is.na(pval))

  if (length(unique(table$conclusion)) == 3){
    colours <- c("black", "blue", "red")
  }
  if (all(table$conclusion %in% c('sign_up', 'not_sign'))){
    colours <- c("black", "red")
  }
  if (all(table$conclusion %in% c('sign_down', 'not_sign'))){
    colours <- c("black", "blue")
  }
  if (length(unique(table$conclusion)) == 1){
    colours <- c("black")
  }

  ggplot(table, aes(x = effect_size, y = `-log10p`, color = conclusion)) +
    geom_point() +
    scale_color_manual(values=colours) +
    geom_label_repel(data = table %>%
                       filter(conclusion %in% c('sign_up', 'sign_down')),
                     aes(label = gene),
                     box.padding   = 0.35,
                     point.padding = 0.5,
                     segment.color = 'grey50') +
    ggtitle(title)
}



#' plot volcano plot 2
#'
#' better version volcano plot
#'
#' @param response table with DEG
#' @param avg_logFC_th avg_logFC threshold
#' @param log10pval log10pval threshold to plot
#' @param y1 y1 limut
#' @param y2 y2 limit
#' @param x1 x1 limit
#' @param x2 x2 limit
#' @param title the name of your plot
#' @param genes_plot force these genes into the plot
#'
#'
#' @return volcano plot
#'
#' @keywords DEG, differential expression analysis, volcano plot
#'
#' @examples
#'
#' volcano_pl(responce)
#'
#' @export

# volcano plots for markers
volcano_pl <- function (response, ES_th = .5, ES_th_plot = .5, log10pval = -log10(.05), y1=0, y2=20, x1=-1, x2=1, title, genes_plot=NA)
{
  response <- response[!is.na(response$p_val_adj),]

  table <- response %>% mutate(`-log10p` = -log10(p_val_adj)) %>%
    mutate(conclusion = ifelse(p_val_adj < 0.05 &
                                 avg_logFC > ES_th, "upregulated", ifelse(p_val_adj <
                                                                            0.05 & avg_logFC < -ES_th, "downregulated", "not_sign"))) %>%
    mutate(to_plot = ifelse(`-log10p` > log10pval &
                              avg_logFC > ES_th_plot, "upregulated", ifelse(`-log10p` >
                                                                              log10pval & avg_logFC < -ES_th_plot, "downregulated", "not_sign")))

  colours <- c('not_sign'='black', 'upregulated'='red', 'downregulated'='blue')

  if (!is.na(genes_plot)){
    table$to_plot[table$gene %in% genes_plot] <- table$conclusion[table$gene %in% genes_plot]
  }

  ggplot(table, aes(x = avg_logFC, y = `-log10p`, color = conclusion)) +
    geom_point() + scale_color_manual(values = colours) +
    geom_label_repel(data = table[table$to_plot != 'not_sign',], aes(label = gene),  size = 3,
                     point.padding = 0.5, segment.color = "grey50") +
    theme_classic() +
    ylim(y1,y2) +
    xlim(x1,x2) +
    ggtitle(title)
}
