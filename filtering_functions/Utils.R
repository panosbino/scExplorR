
import_into_sobj <- function(data_dir, project_name){

  sc_data <- Seurat::Read10X(
    data.dir = data_dir, gene.column = 1
  )
  sobj <- Seurat::CreateSeuratObject(counts = sc_data, project = project_name)
  rm(sc_data)
  return(sobj)
}

get_mito_genes <- function(organism, ensembl_version){
  ensembl <- useEnsembl(biomart = 'genes',
                        dataset = paste0(organism, '_gene_ensembl'),
                        version = ensembl_version)

  mito_genes <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name'),
                      filters='chromosome_name',
                      values='MT',
                      mart=ensembl)
  return(mito_genes)

}


read_and_merge_sobjs <- function(data_rep, sample_names){
  sobj_list = list()
  for (sample_name in sample_names) {
    print(sample_name)
    ## here would go things that are necessary to do in a sample independent manner, i.e. noise correction and doublet detection, which we are skipping for now
    # Import data into seurat object
    sobj_temp <- import_into_sobj(data_dir = file.path(data_rep, sample_name), project_name = sample_name)

    # Rename cell barcode names
    sobj_temp <- Seurat::RenameCells(sobj_temp, new.names = paste0(Cells(sobj_temp), ".", sample_name))

    # Append to list
    sobj_list[[sample_name]] <- sobj_temp
  }
  # Merge list
  sobj <- merge(sobj_list[[1]], sobj_list[2:length(sobj_list)])
  sobj <- JoinLayers(sobj)

  return(sobj)
}


load_and_annotate_recipe_1 <- function(data_rep, organism, ensembl_version){
  sample_names <- list.files(data_rep)
  mito_genes <- get_mito_genes(organism, ensembl_version)
  sobj <- read_and_merge_sobjs(data_rep, sample_names)
  sobj[["percent.mt"]] <- Seurat::PercentageFeatureSet(sobj, features = mito_genes$ensembl_gene_id)

  return(sobj)
}


#c("cpm","log","scran","asinh")
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


