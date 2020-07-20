#' Splits object into a list of subsetted objects.
#'
#' Splits object based on a single attribute into a list of subsetted objects,
#' one for each level of the attribute. For example, useful for taking an object
#' that contains cells from many patients, and subdividing it into
#' patient-specific objects.
#'
#' @param object Seurat object
#' @param split.by Attribute for splitting. Default is "ident". Currently
#' only supported for class-level (i.e. non-quantitative) attributes.
#'
#' @return A named list of Seurat objects, each containing a subset of cells
#' from the original object.
#'
#' @export
#'
#' @examples
#' # Assign the test object a three level attribute
#' groups <- sample(c("group1", "group2", "group3"), size = 80, replace = TRUE)
#' names(groups) <- colnames(pbmc_small)
#' pbmc_small <- AddMetaData(object = pbmc_small, metadata = groups, col.name = "group")
#' obj.list <- SplitObject(pbmc_small, split.by = "group")
#'
SplitObject <- function(object, split.by = "ident") {
  if (split.by == 'ident') {
    groupings <- Idents(object = object)
  } else {
    groupings <- FetchData(object = object, vars = split.by)[, 1]
  }
  groupings <- unique(x = as.character(x = groupings))
  obj.list <- list()
  for (i in groupings) {
    if (split.by == "ident") {
      obj.list[[i]] <- SubsetData(object, idents = i)
    }
    else {
      cells <- which(x = object[[split.by, drop = TRUE]] == i)
      cells <- colnames(x = object)[cells]
      obj.list[[i]] <- SubsetData(object, cells = cells)
    }
  }
  return(obj.list)
}
