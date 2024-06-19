get_possible_ensembl_versions <- function(){
  require(biomaRt)
  return(biomaRt::listEnsemblArchives()[['version']])
}

get_possible_organisms <- function(ensembl_version){
  require(biomaRt)
  ensembl <- biomaRt::useEnsembl(biomart = "genes", version = ensembl_version)
  return(biomaRt::listDatasets(ensembl)[["dataset"]])
}

import_into_sobj <- function(data_dir, project_name){
  require(Seurat)
  sc_data <- Seurat::Read10X(
    data.dir = data_dir, gene.column = 1
  )
  sobj <- SeuratObject::CreateSeuratObject(counts = sc_data, project = project_name)
  rm(sc_data)
  return(sobj)
}

get_mito_genes <- function(organism, ensembl_version){
  require(biomaRt)
  ensembl <- biomaRt::useEnsembl(biomart = 'genes',
                        dataset = organism,
                        version = ensembl_version)

  mito_genes <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name'),
                      filters='chromosome_name',
                      values='MT',
                      mart=ensembl)
  return(mito_genes)

}


read_and_merge_sobjs <- function(data_rep, sample_names){
  require(Seurat)
  sobj_list = list()
  for (sample_name in sample_names) {
    print(sample_name)
    ## here would go things that are necessary to do in a sample independent manner, i.e. noise correction and doublet detection, which we are skipping for now
    # Import data into seurat object
    sobj_temp <- import_into_sobj(data_dir = file.path(data_rep, sample_name), project_name = sample_name)

    # Rename cell barcode names
    sobj_temp <- Seurat::RenameCells(sobj_temp, new.names = paste0(SeuratObject::Cells(sobj_temp), ".", sample_name))

    # Append to list
    sobj_list[[sample_name]] <- sobj_temp
  }
  # Merge list
  sobj <- merge(sobj_list[[1]], sobj_list[2:length(sobj_list)])
  sobj <- SeuratObject::JoinLayers(sobj)

  return(sobj)
}


load_and_annotate_recipe_1 <- function(data_rep, organism, ensembl_version){
  mito_genes <- get_mito_genes(organism, ensembl_version)
  files_in_data_rep <- list.files(data_rep)

  if ("matrix.mtx.gz" %in% files_in_data_rep){
    # if we are only selecting one folder with data
    sample_name <- c(basename(data_rep))
    data_rep <- dirname(data_rep)
    sobj <- import_into_sobj(data_dir = file.path(data_rep, sample_name), project_name = sample_name)
  } else {
    # if we are loading in multiple datasets
    sample_names <- files_in_data_rep
    sobj <- read_and_merge_sobjs(data_rep, sample_names)
  }
  sobj[["percent.mt"]] <- Seurat::PercentageFeatureSet(sobj, features = mito_genes$ensembl_gene_id)

  return(sobj)
}

show_QC_plots <- function(sobj){
  require(patchwork)
  top <- VlnPlot(sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

  bottom_1 <- FeatureScatter(sobj, feature1 = "nCount_RNA", feature2 = "percent.mt")
  bottom_2 <- FeatureScatter(sobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

  return(top / (bottom_1 + bottom_2))

}

filter_sobj <- function(sobj, min_genes = NULL, max_genes = NULL,
                        min_UMIs = NULL, max_UMIs = NULL,
                        max_pct_mito = NULL){
  require(Seurat)
  if (!is.null(min_genes)) {
    message("filtering min genes")
    sobj <- subset(sobj, subset = nFeature_RNA > min_genes)
  }
  if (!is.null(max_genes)) {
    message("filtering max genes")
    sobj <- subset(sobj, subset = nFeature_RNA < max_genes)
  }
  if (!is.null(min_UMIs)) {
    sobj <- subset(sobj, subset = nCount_RNA > min_UMIs)
  }
  if (!is.null(max_UMIs)) {
    sobj <- subset(sobj, subset = nCount_RNA > max_UMIs)
  }
  if (!is.null(max_pct_mito)) {
    message("filtering high mitos")
    sobj <- subset(sobj, subset = percent.mt < max_pct_mito)
  }
  return (sobj)
}

preview_filtered_sobj <- function(sobj, min_genes = NULL, max_genes = NULL,
                                  min_UMIs = NULL, max_UMIs = NULL,
                                  max_pct_mito = NULL) {
  # returns a plot with showing QC plots of filtered data
  sobj <- filter_sobj(sobj, min_genes = min_genes, max_genes = max_genes,
                      min_UMIs = min_UMIs, max_UMIs = max_UMIs,
                      max_pct_mito = max_pct_mito)
  preview_plot <- show_QC_plots(sobj)
  return(preview_plot)
}

commit_filtered_sobj <- filter_sobj

export_sobj <- function(sobj, path, save_as = 'rds'){
  if (save_as == 'rds'){
    if (endsWith(path, '.rds')){saveRDS(sobj, file = path)
    }
    else{warning('path should end with .rds!')}
  }
  else {message('only rds is supported for now')}
}
get_max_ngenes <- function(sobj) {return(max(sobj@meta.data$nFeature_RNA))}
get_max_nUMIs <- function(sobj){return(max(sobj@meta.data$nCount_RNA))}
get_max_pct_mito <- function(sobj){return(max(sobj@meta.data$percent.mt))}


#c("cpm","log","scran","asinh")
normalize_all <- function(sc_input_object,method = "cpm") {
  if(class(sc_input_object) == "Seurat")
  {sce <- Seurat::as.SingleCellExperiment(sc_input_object)}
  else if (class(sc_input_object) == "SingleCellExperiment")
  {sce <- sc_input_object}

  if (method == "cpm") {
    assay(sce, "normed") <- scuttle::normalizeCounts(sce, log=FALSE, transform = "none", size.factors=librarySizeFactors(sce), pseudo.count=0)
  }
  else if(method == "log"){
    assay(sce, "normed") <- scuttle::normalizeCounts(sce, log=TRUE, transform = "log",size.factors=librarySizeFactors(sce), pseudo.count=1)
  }
  else if (method == "scran"){
    clusters <- scran::quickCluster(sce)
    sce <- computeSumFactors(sce, clusters=clusters)
    assay(sce, "normed") <- scuttle::normalizeCounts(sce, log=TRUE, transform = "log",size.factors= sizeFactors(sce) , pseudo.count=1)
  }
  else if(method == "log_geom"){
    assay(sce, "normed") <- scuttle::normalizeCounts(sce, log=TRUE, transform = "log",size.factors=geometricSizeFactors(sce), pseudo.count=1)
  }
  else if(method == "asihn"){
    assay(sce, "normed") <- scuttle::normalizeCounts(sce, log=FALSE, transform = "asihn",size.factors=librarySizeFactors(sce), pseudo.count=1)
  }
  #tmp <- assay(sce, "normed")
  return(Seurat::as.Seurat(sce, counts = "counts", data = "normed"))
}


# nfeatures, number of PCs and dims can all be user defined input variables if we want to do that.
dim_reduction <- function(normalized_seurat_obj){
  normalized_seurat_obj <- FindVariableFeatures(normalized_seurat_obj, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(normalized_seurat_obj)
  normalized_seurat_obj <- ScaleData(normalized_seurat_obj, features = all.genes)
  dim_reduced_seurat_obj <- RunPCA(normalized_seurat_obj, features = VariableFeatures(object = normalized_seurat_obj))
  dim_reduced_seurat_obj <- RunUMAP(dim_reduced_seurat_obj, dims = 1:10)
  return(dim_reduced_seurat_obj)
}
