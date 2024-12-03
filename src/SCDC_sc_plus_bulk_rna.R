library(MuSiC)
library(SCDC)

references = readRDS("/Users/guillaume/Documents/HADACA/starting_kit_phase1/data/reference_pdac.rds")
mixes = readRDS("/Users/guillaume/Documents/HADACA/starting_kit_phase1/data/mixes1_SDE5_pdac.rds")

typeof(references$ref_scRNA$ref_sc_peng)
references$ref_scRNA$ref_sc_raghavan

# Using getESET to create eset object from SCDC/R/Basic_Function.R
eset_peng <- getESET(references$ref_scRNA$ref_sc_peng$counts, row.names(references$ref_scRNA$ref_sc_peng$counts),
                  references$ref_scRNA$ref_sc_peng$metadata)
eset_baron <- getESET(references$ref_scRNA$ref_sc_baron$counts, row.names(references$ref_scRNA$ref_sc_baron$counts),
                     references$ref_scRNA$ref_sc_baron$metadata)
eset_raghavan <- getESET(references$ref_scRNA$ref_sc_raghavan$counts, row.names(references$ref_scRNA$ref_sc_raghavan$counts),
                      references$ref_scRNA$ref_sc_raghavan$metadata)
colnames(mixes$mix_rna) <- 1:30

bulk_mixes_eset <- getESET(mixes$mix_rna, row.names(mixes$mix_rna), colnames(mixes$mix_rna))


bulk.baron <- SCDC_prop(bulk.eset = bulk_mixes_eset,
                              sc.eset = eset_baron,
                              ct.varname="cell_type",
                              sample = "sample",
                              ct.sub=c("basal", "classic", "endo", "fibro", "immune"))

bulk.peng <- SCDC_prop(bulk.eset = bulk_mixes_eset,
                        sc.eset = eset_peng,
                        ct.varname="cell_type",
                        sample = "sample",
                        ct.sub=c("basal", "classic", "endo", "fibro", "immune"))


global_peng_baron <- SCDC_ENSEMBLE(bulk.eset = bulk_mixes_eset,
                       sc.eset = c(eset_peng, eset_baron),
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

ref_bulk_rna = references$ref_bulkRNA
mix_bulk_rna = mixes$mix_rna

nnls_bulk_rna <- global_nnls(mix_bulk_rna, ref_bulk_rna)
bulk_rna_scdc_prop <- CreateSCDCpropObj(t(nnls_bulk_rna), ref_bulk_rna)


peng_baron_plus_bulk_rna <- SCDC_ENSEMBLE(bulk.eset = bulk_mixes_eset,
                                   ct.varname="cell_type",
                                   sample = "sample",
                                   ct.sub=c("basal", "classic", "endo", "fibro", "immune"),
                                   prop.input = list(bulk.baron, bulk.peng, bulk_rna_scdc_prop))

