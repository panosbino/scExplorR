
import_into_sobj <- function(data_dir){
  sc_data <- Seurat::Read10X(
    data.dir = data_dir, gene.column = 1
  )
  sobj <- Seurat::CreateSeuratObject(counts = sc_data)
  rm(sc_data)
  return(sobj)
}


normalize_all <- function(sc_input_object,method = "cpm") {
  if(class(sc_input_object) == "Seurat")
  {sce <- Seurat::as.SingleCellExperiment(sc_input_object)}
  else if (class(sc_input_object) == "SingleCellExperiment")
  {sce <- sc_input_object}

  if (method == "cpm") {
    assay(sce, "normed") <- scuttle::calculateCPM(sce)
  }
  else if(method == "log"){
    assay(sce, "normed") <- scuttle::normalizeCounts(sce, log=FALSE, transform = "log",size.factors=librarySizeFactors(sce), pseudo.count=1)
  }
  else if (method == "scran"){
    clusters <- scran::quickCluster(sce)
    size_factors <- scuttle::computePooledFactors(sce, clusters=clusters)
    assay(sce, "normed") <- scuttle::normalizeCounts(sce, log=FALSE, transform = "log",size.factors=librarySizeFactors(size_factors) , pseudo.count=1)
  }
  else if(method == "asinh"){
    assay(sce, "normed") <- scuttle::normalizeCounts(sce, log=FALSE, transform = "asihn",size.factors=librarySizeFactors(sce), pseudo.count=1)
  }
  return(Seurat::as.Seurat(sce, counts = "counts", data = "normed"))
}


