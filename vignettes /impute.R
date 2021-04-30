Sys.setenv(RETICULATE_PYTHON = "/Users/kbl105/anaconda3/bin/python")

library(Seurat)
library(Rmagic)
E20[['imputed']] <- CreateAssayObject(data = as('dgCMatrix', t(
  magic(data = t(E20@assays$RNA@data), t = 2))))


PP.combined <- merge(
  lapply(SplitObject(PP.combined, split.by = 'sample'),
                      function(Seurat_obj){
  Seurat_obj[['imputed']] <- CreateAssayObject(data = as('dgCMatrix', t(
    magic(data = t(Seurat_obj@assays$RNA@data), t = 2)
          )
        )
      )
    Seurat_obj
  })
)



q <- magic(data = subset(MLNb, cells = Cells(MLNb)[1:10]), t = 2)

??Rmagic::magic

