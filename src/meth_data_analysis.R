if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(stringr)
library(dplyr)
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annoIllunimaMeth450<-dplyr::select(as.data.frame(ann), Name, chr, pos, Relation_to_Island, UCSC_RefGene_Name)

## read reference methylation data
reference_data <- readRDS(file ="/Users/solene/Workspace/HADACA/starting_kit_phase1/data/reference_pdac.rds")
ref_met = reference_data$ref_met
ref_bulkrna = reference_data$ref_bulkRNA

make_unique_values <- function(x) {
  value <- sapply(strsplit(x,";"), function(x) paste(unique(x), collapse = ","))
  return(value)
}

# filter data to keep only CpG island from the ref_met and reduce computation time
df_new <- annoIllunimaMeth450[annoIllunimaMeth450$Name %in% intersect(rownames(ref_met), annoIllunimaMeth450$Name),]
methylation.subset <- subset(df_new, df_new$UCSC_RefGene_Name!="")

methylation.subset.gene.list <- lapply(methylation.subset[,c("UCSC_RefGene_Name")], FUN=make_unique_values)

# filter data to keep only CpG island matching a gene name from the bulk reference
bulk.meth.feature.intersection <- intersect(rownames(ref_bulkrna), unlist(strsplit(unlist(methylation.subset.gene.list), ",")))
