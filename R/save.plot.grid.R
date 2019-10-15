#' save multiple gg plots in pdf
#'
#' This function allows you to save multiple gg plots in pdf
#'
#'
#' @param ...  ggplot objects
#' @param filename filename
#'
#' @keywords ggplot, save 
#'
#' @examples
#' p1 <- qplot(1:10)
#' p2 <- qplot(1:20)
#' p3 <- qplot(1:30)
#' 
#' save.plot.grid(p1,p2,p3, 'whatever')
#'
#'
#' @export

save.plot_grid <- function(..., filename){
  library(ggplot2)
  library(gridExtra)
  
  grid.arrange(...) #arranges plots within grid
  
  #save
  g <- arrangeGrob(...) #generates g
  
  ggsave(file = paste0(filename, ".pdf"), g) #saves   
}



