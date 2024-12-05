if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19",force = TRUE)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(stringr)
library(dplyr)
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annoIllunimaMeth450<-dplyr::select(as.data.frame(ann), Name, chr, pos, Relation_to_Island, UCSC_RefGene_Name)

make_unique_values <- function(x) {
  value <- sapply(strsplit(x,";"), function(x) paste(unique(x), collapse = ","))
  return(value)
}

# Map intersected genes to their probes
get_probes_for_genes <- function(genes, methylation_data) {
  gene_to_probe_map <- setNames(methylation_data$Name, methylation_data$UCSC_RefGene_Name)
  probe_list <- gene_to_probe_map[genes]
  return(probe_list[!is.na(probe_list)]) # Remove NA values
}

# Extract probes corresponding to the union of indices
get_probes_for_indices <- function(indices, data) {
  return(data[indices, "Name"])
}

# Extract genes corresponding to the probes
get_genes_for_probes <- function(probes, methylation_data) {
  # Map probes back to genes using the UCSC_RefGene_Name column
  probe_to_gene_map <- setNames(methylation_data$UCSC_RefGene_Name, methylation_data$Name)
  
  # Get genes corresponding to the input probes
  genes <- probe_to_gene_map[probes]
  
  # Clean up duplicate or multi-mapped genes (e.g., "GENE1;GENE2")
  genes <- sapply(genes, function(x) {
    if (!is.na(x)) make_unique_values(x) else NA
  })
  
  # Remove NA values and return unique genes
  return(unique(genes[!is.na(genes)]))
}

# Compute explicit MAD for each gene
compute_mad <- function(gene_values) {
  # Calculate the median of the gene values
  median_gene <- median(gene_values, na.rm = TRUE)
  # Compute the absolute deviations
  absolute_deviations <- abs(gene_values - median_gene)
  # Return the median of the absolute deviations
  return(median(absolute_deviations, na.rm = TRUE))
}

ref_met = reference_data$ref_met
ref_bulkrna = reference_data$ref_bulkRNA

# filter data to keep only CpG island from the ref_met and reduce computation time
df_new <- annoIllunimaMeth450[annoIllunimaMeth450$Name %in% intersect(rownames(ref_met), annoIllunimaMeth450$Name),]
methylation.subset <- subset(df_new, df_new$UCSC_RefGene_Name!="")

methylation.subset.gene.list <- lapply(methylation.subset[,c("UCSC_RefGene_Name")], FUN=make_unique_values)
methylation.subset.probe.list <- lapply(methylation.subset[,c("Name")], FUN=make_unique_values)

# filter data to keep only CpG island matching a gene name from the bulk reference
bulk.meth.gene.intersection <- intersect(rownames(ref_bulkrna), unlist(strsplit(unlist(methylation.subset.gene.list), ",")))
#bulk.meth.probe.intersection <- intersect(rownames(ref_bulkrna), unlist(strsplit(unlist(methylation.subset.probe.list), ",")))

# Get probes for the intersected genes
bulk_meth_probe_intersection <- get_probes_for_genes(bulk.meth.gene.intersection, methylation.subset)

# Extract methylation values for the intersected probes
probes_to_keep <- bulk_meth_probe_intersection
filtered_ref_met <- ref_met[probes_to_keep, ]


# Define genes_to_keep based on the probes_to_keep
genes_to_keep <- get_genes_for_probes(probes_to_keep, annoIllunimaMeth450)

filtered_ref_bulkrna <- ref_bulkrna[genes_to_keep, ]

# Map probes to genes (assuming many-to-one mapping)
probes_to_genes <- setNames(methylation.subset$UCSC_RefGene_Name, methylation.subset$Name)
mapped_genes <- probes_to_genes[rownames(filtered_ref_met)]


# Apply the function to each row of filtered_ref_met
mad_values_explicit_met <- apply(filtered_ref_met, 1, compute_mad)
mad_values_explicit_bulkrna <- apply(filtered_ref_bulkrna, 1, compute_mad)

# Sort the genes by MAD values in descending order
sorted_mad_indices_met <- order(mad_values_explicit_met, decreasing = TRUE)
sorted_mad_indices_bulkrna <- order(mad_values_explicit_bulkrna, decreasing = TRUE)

n=300

top_n_indices_met=sorted_mad_indices_met[0:n]
top_n_indices_bulkrna=sorted_mad_indices_bulkrna[0:n]

# Get the union of the indices
top_n_indices_union <- union(top_n_indices_met, top_n_indices_bulkrna)

union_genes=bulk.meth.gene.intersection[top_n_indices_union]
union_probes=get_probes_for_indices(top_n_indices_union, annoIllunimaMeth450)

# Combine cleaned genes and corresponding probes
top_chosen_features <- list(
  genes = union_genes,
  probes = union_probes
)

saveRDS(top_chosen_features, file = "top_chosen_features_v4.rds")
