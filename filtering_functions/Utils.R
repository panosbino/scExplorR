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
  top <- plot_QC_violin(sobj)
  bottom_left <- plot_QC_scatter_left(sobj)
  bottom_right <- plot_QC_scatter_right(sobj)
  return(top / (bottom_left + bottom_right))

}

filter_sobj <- function(sobj, min_genes = NULL, max_genes = NULL,
                        min_UMIs = NULL, max_UMIs = NULL,
                        min_pct_mito = NULL,
                        max_pct_mito = NULL){
  require(Seurat)
  if (!is.null(min_genes)) {
    sobj <- subset(sobj, subset = nFeature_RNA >= min_genes)
  }
  if (!is.null(max_genes)) {
    sobj <- subset(sobj, subset = nFeature_RNA <= max_genes)
  }
  if (!is.null(min_UMIs)) {
    sobj <- subset(sobj, subset = nCount_RNA >= min_UMIs)
  }
  if (!is.null(max_UMIs)) {
    sobj <- subset(sobj, subset = nCount_RNA >= max_UMIs)
  }
  if (!is.null(max_pct_mito)) {
    sobj <- subset(sobj, subset = percent.mt <= max_pct_mito)
  }
  if (!is.null(min_pct_mito)) {
    sobj <- subset(sobj, subset = percent.mt >= min_pct_mito)
  }
  return (sobj)
}

preview_filtered_sobj <- function(sobj, min_genes = NULL, max_genes = NULL,
                                  min_UMIs = NULL, max_UMIs = NULL,
                                  min_pct_mito = NULL,
                                  max_pct_mito = NULL) {
  # returns a plot with showing QC plots of filtered data
  sobj <- filter_sobj(sobj, min_genes = min_genes, max_genes = max_genes,
                      min_UMIs = min_UMIs, max_UMIs = max_UMIs, min_pct_mito = min_pct_mito,
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

get_bomb_palette <- function(n_colors){
  bomb_colors <- c('#cb353d', '#f9b64e', '#ed6240', '#563d43')
  bomb_palette <- colorRampPalette(bomb_colors)
  return(bomb_palette)
}

plot_QC_scatter_left <- function(sobj){
  plot_data <- tidyd_qcp_point(sobj)
  bomb_palette <- get_bomb_palette()
  n_samples <- plot_data$orig.ident |> unique() |> length()

  scatter_size = 1200/dim(plot_data)[1]

  pl <- plot_data |>
    ggplot(aes(x = nGenes, y = pct_mito, color = orig.ident),) +
    geom_point(size = scatter_size) + theme_classic() +
    scale_color_manual(values = bomb_palette(n_samples)) +
    scale_y_continuous(expand = expansion(mult = c(0.01, 0))) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0))) +
    theme(legend.position = "none")

  return(pl)
}

plot_QC_scatter_right<- function(sobj){
  plot_data <- tidyd_qcp_point(sobj)
  bomb_palette <- get_bomb_palette()
  n_samples <- plot_data$orig.ident |> unique() |> length()

  scatter_size = 1200/dim(plot_data)[1]

  pl <- plot_data |>
    ggplot(aes(x = nGenes, y = nUMIs, color = orig.ident)) +
    geom_point(size = scatter_size) + theme_classic() +
    scale_color_manual(values = bomb_palette(n_samples)) +
    scale_y_continuous(expand = expansion(mult = c(0.01, 0))) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0))) +
    theme(legend.position = "none")
  return(pl)
}

tidyd_qcp_point <- function(sobj){
  plot_data <- sobj@meta.data |>
    dplyr::rename(nUMIs = 'nCount_RNA',
                  nGenes = 'nFeature_RNA',
                  pct_mito = 'percent.mt')
  return(plot_data)
}

tidyd_qcp_violin <- function(sobj) {
  plot_data <- sobj@meta.data |>
    dplyr::rename(nUMIs = 'nCount_RNA',
           nGenes = 'nFeature_RNA',
           pct_mito = 'percent.mt') |>
    tidyr::gather('key', 'value',-orig.ident)

  return(plot_data)
}

plot_QC_violin <- function(sobj){
  require(patchwork)
  require(ggplot2)
  bomb_palette <- get_bomb_palette()

  plot_data <- tidyd_qcp_violin(sobj)
  n_samples <- plot_data$orig.ident |> unique() |> length()

  nUMI_violin <- plot_data |>
    dplyr::filter(key == 'nUMIs') |>
    ggplot(aes(x = orig.ident, y = value, fill = orig.ident)) +
    geom_violin() + scale_fill_manual(values =
                                        bomb_palette(n_samples)) + facet_grid( ~ key) + theme_classic() +
    theme(legend.position = "none", axis.title = element_blank())
  nGenes_violin <- plot_data |>
    dplyr::filter(key == 'nGenes') |>
    ggplot(aes(x = orig.ident, y = value, fill = orig.ident)) +
    geom_violin() + scale_fill_manual(values =
                                        bomb_palette(n_samples)) + facet_grid( ~ key) + theme_classic() +
    theme(legend.position = "none", axis.title = element_blank())
  pct_mito_violin <- plot_data |>
    dplyr::filter(key == 'pct_mito') |>
    ggplot(aes(x = orig.ident, y = value, fill = orig.ident)) +
    geom_violin() + scale_fill_manual(values =
                                        bomb_palette(n_samples)) + facet_grid( ~ key) + theme_classic() +
    theme(legend.position = "none", axis.title = element_blank())

  return(nUMI_violin + nGenes_violin + pct_mito_violin)
}

#c("cpm","log","scran","asinh")
normalize_all <- function(sc_input_object,method = "cpm") {
  require(scuttle)
  require(scRNAseq)
  require(SingleCellExperiment)
  require(tidyverse)
  require(Seurat)
  #print("0")
  if(class(sc_input_object) == "Seurat")
  {sce <- Seurat::as.SingleCellExperiment(sc_input_object)}
  else if (class(sc_input_object) == "SingleCellExperiment")
  {sce <- sc_input_object}
  #print("1")
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
  print("2")
  seurat_norm_obj <- Seurat::as.Seurat(sce, counts = "counts", data = "normed")
  return(seurat_norm_obj)
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


# Needs to be generalized for different number of samples
plot_UMAP <- function(dim_reduced_seurat_obj) {
  UMAPs <- dim_reduced_seurat_obj@reductions$umap@cell.embeddings %>% cbind(dim_reduced_seurat_obj$orig.ident) %>% as.data.frame()
  UMAPs$umap_1 <- UMAPs$umap_1 %>% as.numeric()
  UMAPs$umap_2 <- UMAPs$umap_2 %>% as.numeric()
  colnames(UMAPs)[3] <- "Sample"
  UMAPs$Sample <- ifelse(str_detect(UMAPs$Sample,pattern = "_"),
                         UMAPs$Sample %>% str_replace(pattern = "_", replacement = " "),
                         UMAPs$Sample)
  ggplot() + geom_point(data = UMAPs, mapping = aes(x = umap_1, y = umap_2, color = Sample), size =0.5) +
    theme_minimal() +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_blank(),
          axis.title = element_text(size = 20),
          legend.position="top",
          legend.text = element_text(size = 18),
          legend.title = element_blank(),
          strip.text = element_blank()) +
    facet_wrap(~Sample, nrow = 2,) +
    guides(colour = guide_legend(override.aes = list(size=10))) +
    scale_color_manual(values = c("#cb353d", "#ed6240", "#f9b64e", "#6a4a57"))
}

# Needs further beautifying, will do later
plot_PCA <- function(dim_reduced_seurat_obj){
  PCAs <- dim_reduced_seurat_obj@reductions$pca@cell.embeddings[,1:3] %>%
    as.data.frame() %>%
    cbind(dim_reduced_seurat_obj$orig.ident)
  colnames(PCAs)[4] <- "Sample"
  PCAs$Sample <- ifelse(str_detect(PCAs$Sample,pattern = "_"),
                        PCAs$Sample %>% str_replace(pattern = "_", replacement = " "),
                        PCAs$Sample)
  plotly::plot_ly(data = PCAs, x=~PC_1, y=~PC_2, z=~PC_3, type="scatter3d", mode = "markers" ,colors= c("#6a4a57"), size = 0.5)
}

# Main function that starts from normalization and produces PCA and UMAP plots.
# It can take more variables, but this is good for now
# This returns a list of 3 objects
# The first is the 3D PCA plot
# The second is the UMAP plot
# And the third is the seurat obeject that has the normalized values and the dimensionality reduction
normalize_and_plot_main <- function(path_to_filtered_seurat_obj, normalization_method = "cpm"){
  sc_filtered_object <- readRDS(path_to_filtered_seurat_obj)
  normalized_seurat_object <- normalize_all(sc_input_object = sc_filtered_object, method = normalization_method)
  dim_reduced_seurat_obj <- dim_reduction(normalized_seurat_object)
  PCA_plot <- plot_PCA(dim_reduced_seurat_obj = dim_reduced_seurat_obj)
  UMAP_plot <- plot_UMAP(dim_reduced_seurat_obj = dim_reduced_seurat_obj)
  return(list(PCA_plot,UMAP_plot,dim_reduced_seurat_obj))
}

# This returns a list of 2 objects
# The first is the clustered UMAP plot
# The second is the the Seurat object with the clustering info
# Would need to check if the number of clusters shown in the plot is correct when you change resolution
cluster <- function(dim_reduced_seurat_obj, clustering_resolution = 0.3){
  clustered_seurat_obj <- FindNeighbors(dim_reduced_seurat_obj, dims = 1:10)
  clustered_seurat_obj <- FindClusters(clustered_seurat_obj, resolution = clustering_resolution)
  clustering_plot <- DimPlot(clustered_seurat_obj, reduction = "umap", group.by = paste0("RNA_snn_res.",clustering_resolution)) +
    theme_minimal() +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    ggtitle(paste0("Clustering resolution = ", clustering_resolution,"  |  ",clustered_seurat_obj@meta.data$seurat_clusters %>% unique() %>% length()," clusters found")) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_blank(),
          axis.title = element_text(size = 20),
          #legend.text = element_text(size = 18),
          legend.position = "none",
          plot.title = element_text(size = 20, hjust = 0.5),
    )
  return(list(clustering_plot,clustered_seurat_obj))
}


