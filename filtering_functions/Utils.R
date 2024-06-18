
import_into_sobj <- function(data_dir){
  sc_data <- Seurat::Read10X(
    data.dir = data_dir, gene.column = 1
  )
  sobj <- Seurat::CreateSeuratObject(counts = sc_data)
  rm(sc_data)
  return(sobj)
}
