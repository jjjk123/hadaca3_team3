> Hadaca3 Data Challenge

# Challenge: "Multimodal data integration to quantify tumor heterogeneity in cancer"

The Data Challenge took place on December 2-6, 2024 in Aussois in France.

The detailed description of the challenge can be found here: https://hadaca3.sciencesconf.org/

# Contributors

### Team B

- Solène Weill (solene@epigenelabs.com)
- Jędrzej Kubica (jedrzej.kubica@univ-grenoble-alpes.fr)
- Vesna Lukic (vesna.lukic@centralesupelec.fr)
- Guillaume Appe (guillaume@epigenlabs.com)

# Table of contents

# Overview

**The aim** of the project was to design and develop a bioinformatic workflow to quantify pancreatic tumor heterogeneity using supervised deconvolution methods and multi-omics data. There have been previous studies that introduced various deconvolution methods[^1], however there are a number of challenges that still persist in the field. The first challenge was the integration of multi-omics data (RNA-seq, single cell RNA-seq, and DNA methylation) for a reference in the deconvolution process, and the second challenge was the selection and combination of the best deconvolution software packages. The project results were measured and compared to other approaches on the Codabench platform[^2].

<img width="535" alt="Screenshot 2024-12-06 at 14 09 53" src="https://github.com/user-attachments/assets/f55e096e-43ae-46e5-b7f8-0f0139d3e486">

This documentation provides comprehensive details of our contribution, which focused on performing cell-type deconvolution using bulk RNA and methylation data, and trying both uni- and multi-modal predictions.

# Workflow

<img width="590" alt="Screenshot 2024-12-06 at 09 47 59" src="https://github.com/user-attachments/assets/37987e5c-9105-4c04-9be6-c96107a04eae">

# Data

The challenge was split into three phases as follows:

- Phase 1: Discovery of the data and the Codabench platform
- Phase 2: Estimation of cell type heterogeneity and submissions of methods/results into the platform
- Phase 3: Migration from phase 2 of the best methods and evalution of them

Data was provided in the first two phases. Phase 1 data consisted of one simulated multi-omic dataset having an in-silico mixture of 5 cell types with explicit dependence between genes/CpG probes. Phase 2 data consisted of xxxxx.

# Methods

## Feature selection 

The feature selection is done on both methylation and single-cell RNA data

For methylation data:
- Begin with full bulk methylation data
- Map CpG islands to genes from bulk RNA
- Find cell-type specific methylation sites by a) Ranking CpG islands by difference between max methylation across all cell types and mean of remaining cell types and b) Select top CpG islands based on threshold (minimum of means methylation across all cell types)

For single-cell RNA:
- Identify marker genes for cell types based on differential expression analysis
- Create pseudo-bulk data

## Cell type deconvolution

Deconvolution: Running SCDC + NNLS and an Ensemble method (this was abandoned however, as it achieved poor results)

Unimodal predictions: 
- NNLS on bulk RNA, bulk methylation and pseudo-bulk created from single-cell RNA-seq separately
- Trying to mix bulk RNA and methylation changing reference to (1-bulk_methylation)*bulk_RNA element wise

Multimodal predictions:
- All unimodal from a. + Ensemble method. This allowed using the intersection of CpG sites to gene, bulkRNA and scRNA genes

<img width="466" alt="Screenshot 2024-12-05 at 19 23 27" src="https://github.com/user-attachments/assets/83f2ab94-dca8-4517-a758-6e89bb2177c2">

# Exploratory data analysis

The script cellType_specific_CpGmet.ipynb shows some preliminary analysis of the data. It shows how genes are clustered based on their expression.

<img width="510" alt="Screenshot 2024-12-06 at 14 21 48" src="https://github.com/user-attachments/assets/2c00674d-378f-4a7c-89c6-b6b412c9ece1">

We get the genes that have the most distinct methylation across the 5 cell types.

<img width="619" alt="Screenshot 2024-12-06 at 14 25 26" src="https://github.com/user-attachments/assets/fdc5360c-603a-44af-a450-a47837729e02">

# Pre-processing

The script single_cell_preprocessing.R reads the single cell reference data, then creates and processes a Seurat object for each single cell dataset. The differential expression is computed, and the markers for each cell type are obtained, as well as the final gene list.

### Dependencies



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

# Codabench Platform[^2]

Competition website: https://www.codabench.org/competitions/4714/

# Results

The best deconvolution results: Mmethylation

Approach NNLS + feature selection (CpG to gene mapping, cell-type specific methylation)
(Python)

Codabench score = 0.66

The decomposition of scores for 9 validation datasets:

![image](https://github.com/user-attachments/assets/df4bdc3b-7647-4e6c-91ea-c2e38c59e7a6)


# Conclusions

# Future ideas

## Special thank you to the Data Challenge Organizers!

# References

[^1]: Epigenomic Deconvolution of Breast Tumors Reveals Metabolic Coupling between Constituent Cell Types: 10.1016/j.celrep.2016.10.057
[^2]: Codabench: Flexible, easy-to-use, and reproducible meta-benchmark platform: 10.1016/j.patter.2022.100543
[^3]: SCDC: bulk gene expression deconvolution by multiple single-cell RNA sequencing references: doi.org/10.1093/bib/bbz166
