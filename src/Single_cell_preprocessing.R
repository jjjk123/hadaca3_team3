install.packages('Seurat')
library(Seurat)

# Read single reference data
reference_data <- readRDS(file ="/Users/solene/Workspace/HADACA/starting_kit_phase2-3/data/reference_pdac.rds")
sc_ref <- reference_data$ref_scRNA

#create seurat object for each sc dataset

create_and_process_seurat <- function(sc_data, mincell, minfeature){
  temp.obj <- (CreateSeuratObject(counts = sc_data$counts, meta.data = sc_data$metadata, min.cells = mincell, min.features = minfeature))
  temp.obj <- NormalizeData(temp.obj, normalization.method = "LogNormalize", scale.factor = 10000)
  temp.obj <- FindVariableFeatures(temp.obj, selection.method = "vst", nfeatures = 2000)
  temp.obj <- ScaleData(temp.obj)
  return (temp.obj)
}

compute_differential_expression <- function(seurat_processed_obj){
  
  cell_types <- unique(seurat_processed_obj@meta.data[["cell_type"]])
  
  # Set identity classes to an existing column in meta data
  Idents(seurat_processed_obj) <- "cell_type"
  
  # Initialize a list to store results
  de_results <- list()
  de_results_genes <- list()
  
  for (cell_type in cell_types) {
    de_result <- FindMarkers(seurat_processed_obj, ident.1 = cell_type)
    
    # Store results in the list
    de_results[[cell_type]] <- de_result
    de_results_genes[[cell_type]] <- rownames(subset(de_results[[cell_type]], avg_log2FC>1 & p_val_adj<0.05 & pct.1>0.75))
  }
  return(de_results_genes)
}

merge_markers_list <- function(l1, l2, l3){
  inner_lists <- list()
  for (list_name in names(l1)) {
    inner_list <- unique(l1[[list_name]], l2[[list_name]], l3[[list_name]])
    inner_lists[[list_name]] <- inner_list
  }
  return (inner_lists)
}

sc_peng <- create_and_process_seurat(sc_ref$ref_sc_peng, 1, 20)
sc_baron <- create_and_process_seurat(sc_ref$ref_sc_baron, 1, 20)
sc_raghavan <- create_and_process_seurat(sc_ref$ref_sc_raghavan, 1, 20)

sc_peng.marker.genes <- compute_differential_expression(sc_peng)
sc_baron.marker.genes  <- compute_differential_expression(sc_baron)
sc_raghavan.marker.genes <- compute_differential_expression(sc_raghavan)

immune.markers <- unique(c(sc_peng.marker.genes[["immune"]], sc_baron.marker.genes[["immune"]], sc_raghavan.marker.genes[["immune"]]))
endo.markers <- unique(c(sc_peng.marker.genes[["endo"]], sc_baron.marker.genes[["endo"]], sc_raghavan.marker.genes[["endo"]]))
fibro.markers <- unique(c(sc_peng.marker.genes[["fibro"]], sc_baron.marker.genes[["fibro"]]))
classic.markers <- unique(c(sc_peng.marker.genes[["classic"]], sc_raghavan.marker.genes[["classic"]]))
basal.markers <- unique(c(sc_peng.marker.genes[["basal"]], sc_raghavan.marker.genes[["basal"]]))

final_gene_list <- unique(c(immune.markers, endo.markers, fibro.markers, classic.markers, basal.markers))
