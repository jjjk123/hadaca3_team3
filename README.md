> Hadaca3 Data Challenge

# Challenge: "Multimodal data integration to quantify tumor heterogeneity in cancer"

The Data Challenge took place on December 2-6, 2024 in Aussois in France.

Description of the challenge: https://hadaca3.sciencesconf.org/

# Team B

- 
-
-
-

# Abstract

**The aim** of the project was to design and develop a bioinformatic deconvolution workflow for pancreatic cancer decomposition.

This documentation provides an overview of the process for performing cell type decomposition using the `SCDC` R libraries. The code provided utilizes reference single-cell RNA-seq datasets and bulk RNA-seq mixtures to infer cell-type proportions.

# Worklow

# Methods

### Dependencies

```R
library(MuSiC) ## Currently not in use, but will be tried
library(SCDC)
```

Ensure these packages are installed before executing the script.

### Data Preparation

### Load Reference and Mixture Datasets

```R
references = readRDS("/path/to/reference_pdac.rds")
mixes = readRDS("/path/to/mixes1_SDE5_pdac.rds")
```
- `references` contains single-cell RNA-seq reference data.
- `mixes` contains bulk RNA-seq mixture data.

### Inspect Reference Data

```R
typeof(references$ref_scRNA$ref_sc_peng)
references$ref_scRNA$ref_sc_raghavan
```
This code verifies the structure of the reference datasets.

## Create ExpressionSet Objects

The `getESET` function transforms the single-cell RNA-seq and bulk RNA-seq data into `ExpressionSet` objects required by the `SCDC` package.

```R
eset_peng <- getESET(references$ref_scRNA$ref_sc_peng$counts, row.names(references$ref_scRNA$ref_sc_peng$counts),
                  references$ref_scRNA$ref_sc_peng$metadata)
eset_baron <- getESET(references$ref_scRNA$ref_sc_baron$counts, row.names(references$ref_scRNA$ref_sc_baron$counts),
                     references$ref_scRNA$ref_sc_baron$metadata)
eset_raghavan <- getESET(references$ref_scRNA$ref_sc_raghavan$counts, row.names(references$ref_scRNA$ref_sc_raghavan$counts),
                      references$ref_scRNA$ref_sc_raghavan$metadata)

colnames(mixes$mix_rna) <- 1:30
bulk_mixes_eset <- getESET(mixes$mix_rna, row.names(mixes$mix_rna), colnames(mixes$mix_rna))
```

## Estimate Cell Type Proportions

### Using SCDC_prop

The `SCDC_prop` function estimates cell-type proportions using single-cell references and bulk mixtures.

```R
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
```

### Using SCDC_ENSEMBLE

The `SCDC_ENSEMBLE` function combines results from multiple single-cell references.

```R
global_peng_baron <- SCDC_ENSEMBLE(bulk.eset = bulk_mixes_eset,
                                   sc.eset = c(eset_peng, eset_baron),
                                   ct.varname="cell_type",
                                   sample = "sample",
                                   ct.sub=c("basal", "classic", "endo", "fibro", "immune"))
```

## Non-Negative Least Squares (NNLS) for Proportions

A custom function `global_nnls` calculates cell-type proportions using NNLS.

```R
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
```

**Note: Ensure proper preprocessing of input data for accurate results.**

# Codabench Platform

Competition website: https://www.codabench.org/competitions/4714/



# Results

The `SCDC_ENSEMBLE` function integrates cell-type proportions estimated from multiple methods.

```R
peng_baron_plus_bulk_rna <- SCDC_ENSEMBLE(bulk.eset = bulk_mixes_eset,
                                          ct.varname="cell_type",
                                          sample = "sample",
                                          ct.sub=c("basal", "classic", "endo", "fibro", "immune"),
                                          prop.input = list(bulk.baron, bulk.peng, bulk_rna_scdc_prop))
```

Codabench score (best method): 0.66

# Conclusions

This script performs cell-type decomposition using the `SCDC` library. The Music library is yet to be tested. It combines single-cell references with bulk RNA-seq mixtures. Key steps include:
1. Loading reference and mixture data.
2. Creating `ExpressionSet` objects.
3. Estimating cell-type proportions using `SCDC_prop` and custom methods.
4. Combining results using `SCDC_ENSEMBLE`.

# Future ideas

# References

## Special thank you to the Data Challenge Organizers!
