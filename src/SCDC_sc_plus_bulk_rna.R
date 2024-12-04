
references = readRDS("/Users/guillaume/Documents/HADACA/starting_kit_phase1/data/reference_pdac.rds")
mixes = readRDS("/Users/guillaume/Documents/HADACA/starting_kit_phase1/data/mixes1_SDE5_pdac.rds")

if (!require("devtools")) {
    install.packages("devtools")
  }
  devtools::install_github("anijudy/SCDC", ref = "anijudy-patch-1")
  library(SCDC)
  
  eset_peng <- getESET(ref_scRNA$ref_sc_peng$counts, row.names(ref_scRNA$ref_sc_peng$counts),
                  ref_scRNA$ref_sc_peng$metadata)
  eset_baron <- getESET(ref_scRNA$ref_sc_baron$counts, row.names(ref_scRNA$ref_sc_baron$counts),
                      ref_scRNA$ref_sc_baron$metadata)
  eset_raghavan <- getESET(ref_scRNA$ref_sc_raghavan$counts, row.names(ref_scRNA$ref_sc_raghavan$counts),
                        ref_scRNA$ref_sc_raghavan$metadata)
  colnames(mix) <- 1:30

  bulk_rna_mixes_eset <- getESET(mix, row.names(mix), colnames(mix))


  bulk.baron <- SCDC_prop(bulk.eset = bulk_rna_mixes_eset,
                                sc.eset = eset_baron,
                                ct.varname="cell_type",
                                sample = "sample",
                                ct.sub=c("basal", "classic", "endo", "fibro", "immune"))

  bulk.peng <- SCDC_prop(bulk.eset = bulk_rna_mixes_eset,
                          sc.eset = eset_peng,
                          ct.varname="cell_type",
                          sample = "sample",
                          ct.sub=c("basal", "classic", "endo", "fibro", "immune"))

  global_nnls = function(mix=NULL, ref=NULL, ...) {
    idx_feat = intersect(rownames(mix), rownames(ref))
    prop = apply(mix[idx_feat,], 2, function(b, A) {
      tmp_prop = nnls::nnls(b=b, A=A)$x
      tmp_prop = tmp_prop / sum(tmp_prop) # Sum To One
      return(tmp_prop)
    }, A=ref[idx_feat,])  
    rownames(prop) = colnames(ref)
    
    return(prop)
  }
  nnls_bulk_rna <- global_nnls(mix, ref)
  bulk_rna_scdc_prop <- CreateSCDCpropObj(t(nnls_bulk_rna), ref)


  peng_baron_plus_bulk_rna <- SCDC_ENSEMBLE(bulk.eset = bulk_rna_mixes_eset,
                                    ct.varname="cell_type",
                                    sample = "sample",
                                    ct.sub=c("basal", "classic", "endo", "fibro", "immune"),
                                    prop.input = list(bulk.peng, bulk.baron, bulk_rna_scdc_prop))
  
  min_arg <- which.min(as.data.frame(peng_baron_plus_bulk_rna$w_table)$RMSD_Y)

  weights <- peng_baron_plus_bulk_rna$w_table[min_arg, c(1,2,3)]

  weighted_matrices <- lapply(1:length(weights), function(x) peng_baron_plus_bulk_rna$prop.only[[x]]*weights[[x]])
  weighted_matrices

  prop <- weighted_matrices[[1]] + weighted_matrices[[2]] + weighted_matrices[[3]]
  
  return(t(prop))