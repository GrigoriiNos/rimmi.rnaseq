#' plor volcano plot
#'
#' This function allows you to vizualise differential expression analisys for bulk stidy by volcano plot
#'
#' @param edgeR_output your edgeR output
#'
#' @param title the name of your plot
#'
#' @return volcano plot
#'
#' @keywords DEG, differential expression analysis, volcano plot
#'
#' @examples
#'
#' volcano_plot(edgeR_output)
#'
#' @export

volcano_plot <- function(edgeR_output, title = 'volcano plot'){
  edgeR_output
  colour <- ifelse(
    abs(edgeR_output$avg_logFC) > 1.5 & edgeR_output$FDR < 0.05,
    "red", "black")
  {
    plot(edgeR_output$avg_logFC,
         -log10(edgeR_output$FDR),
         col = colour,
         pch = 16,
         xlab = "avg_logFC",
         ylab = "-10log(q-value)",
         main = title)
    abline(h = -1.5, v = -1.5, col = "gray60")
    abline(h = 1.3, v = -1.5, col = "gray60")
    abline(h = -1.5, v = 1.5, col = "gray60")
  }
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
