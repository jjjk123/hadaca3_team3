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

# Abstract

**The aim** of the project was to design and develop a bioinformatic workflow to quantify pancreatic tumor heterogeneity using supervised deconvolution methods and multi-omics data. There have been previous studies that introduced various deconvolution methods[^1], however there are a number of challenges that still persist in the field. The first challenge was the integration of multi-omics data (RNA-seq, single cell RNA-seq, and DNA methylation) for a reference in the deconvolution process, and the second challenge was the selection and combination of the best deconvolution software packages. The project results were measured and compared to other approaches on the Codabench platform[^2].

This documentation provides an overview of the our project, which focused on performing cell type decomposition using the `SCDC` R libraries[^3]. The code provided utilizes reference single-cell RNA-seq datasets and bulk RNA-seq mixtures to infer cell-type proportions.

# Workflow

<img width="590" alt="Screenshot 2024-12-06 at 09 47 59" src="https://github.com/user-attachments/assets/37987e5c-9105-4c04-9be6-c96107a04eae">


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

This script performs cell-type decomposition using the `SCDC` library. The Music library is yet to be tested. It combines single-cell references with bulk RNA-seq mixtures. Key steps include:
1. Loading reference and mixture data.
2. Creating `ExpressionSet` objects.
3. Estimating cell-type proportions using `SCDC_prop` and custom methods.
4. Combining results using `SCDC_ENSEMBLE`.

# Future ideas

## Special thank you to the Data Challenge Organizers!

# References

[^1]: Epigenomic Deconvolution of Breast Tumors Reveals Metabolic Coupling between Constituent Cell Types: 10.1016/j.celrep.2016.10.057
[^2]: Codabench: Flexible, easy-to-use, and reproducible meta-benchmark platform: 10.1016/j.patter.2022.100543
[^3]: SCDC: bulk gene expression deconvolution by multiple single-cell RNA sequencing references: doi.org/10.1093/bib/bbz166
