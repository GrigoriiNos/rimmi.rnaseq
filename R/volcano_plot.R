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

volcano_pl <- function(response, logFC_th = 0.5, title){
  table <- response %>%
    mutate(`-log10p` = -log(adj.P.Val, 10)) %>%
    mutate(col = ifelse((adj.P.Val < 0.05 & abs(logFC) > 0.5), 'red', 'grey'))

  ggplot(table, aes(x = logFC, y = `-log10p`, color = col)) +
    geom_point() +
    scale_color_manual(breaks = c("8", "6"), values=c("blue", "red")) +
    geom_label_repel(aes(label = gene),
                     box.padding   = 0.35,
                     point.padding = 0.5,
                     segment.color = 'grey50') +
    ggtitle(title)
}
