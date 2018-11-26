#' plor volcano plot
#'
#' This function allows you to vizualise differential expression analisys for bulk stidy by volcano plot
#' @param edgeR_output your edgeR output
#' @param title the name of your plot
#' @keywords DEG, differential expression analysis, volcano plot
#' @export
#' @examples
#' volcano_plot(edgeR_output)

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
