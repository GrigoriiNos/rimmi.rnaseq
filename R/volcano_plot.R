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
    abs(edgeR_output$logFC) > 1.5 & edgeR_output$FDR < 0.05,
    "red", "black")
  {
    plot(edgeR_output$logFC,
         -log10(edgeR_output$FDR),
         col = colour,
         pch = 16,
         xlab = "LogFC",
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
#' @param logFC_th LogFC threshold
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

# volcano plots for markers
volcano_pl <- function(response, ES_th = 0.5, title){
  table <- response %>%
    mutate(`-log10p` = -log10(adj.P.Val)) %>%
    mutate(conclusion = ifelse(response$adj.P.Val < 0.05 & response$logFC > ES_th,
                               'sign_up',
                               ifelse(response$adj.P.Val < 0.05 & response$logFC < -ES_th,
                                      'sign_down',
                                      'not_sign')))

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

  ggplot(table, aes(x = logFC, y = `-log10p`, color = conclusion)) +
    geom_point() +
    scale_color_manual(values=colours) +
    geom_label_repel(aes(label = gene),
                     box.padding   = 0.35,
                     point.padding = 0.5,
                     segment.color = 'grey50') +
    ggtitle(title)
}
