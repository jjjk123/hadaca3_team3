> Hadaca3 Data Challenge

# Challenge: Multimodal data integration to quantify tumor heterogeneity in cancer

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

<img width="648" alt="Screenshot 2024-12-07 at 09 08 11" src="https://github.com/user-attachments/assets/da7cf24f-4c4f-4951-971c-8074071ebfe2">


# Data

The challenge was split into three phases as follows:

- Phase 1: Discovery of the data and the Codabench platform
- Phase 2: Estimation of cell type heterogeneity and submissions of methods/results into the platform
- Phase 3: Migration from phase 2 of the best methods and evalution of them

Data was provided in the first two phases. Phase 1 data consisted of one simulated multi-omic dataset having an in-silico mixture of 5 cell types with explicit dependence between genes/CpG probes. Phase 2 data consisted of xxxxx.

# Methods

Our approach consisted firstly of feature selection and cell type deconvolution.

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
- a) NNLS on bulk RNA, bulk methylation and pseudo-bulk created from single-cell RNA-seq separately
- b) Trying to mix bulk RNA and methylation changing reference to (1-bulk_methylation)*bulk_RNA element wise

Multimodal predictions:
- All unimodal from a. + Ensemble method. This allowed using the intersection of CpG sites to gene, bulkRNA and scRNA genes

<img width="466" alt="Screenshot 2024-12-05 at 19 23 27" src="https://github.com/user-attachments/assets/83f2ab94-dca8-4517-a758-6e89bb2177c2">

# Exploratory data analysis

The script `cellType_specific_CpGmet.ipynb` shows some preliminary analysis of the data. It shows how genes are clustered based on their expression.

<img width="510" alt="Screenshot 2024-12-06 at 14 21 48" src="https://github.com/user-attachments/assets/2c00674d-378f-4a7c-89c6-b6b412c9ece1">

We get the genes that have the most distinct methylation across the 5 cell types.

<img width="619" alt="Screenshot 2024-12-06 at 14 25 26" src="https://github.com/user-attachments/assets/fdc5360c-603a-44af-a450-a47837729e02">

# Pre-processing

The script `single_cell_preprocessing.R` reads the single cell reference data, then creates and processes a Seurat object for each single cell dataset. The differential expression is computed, and the markers for each cell type are obtained, as well as the final gene list.

`meth_data_analysis.R` installs the annotation file `IlluminaHumanMethylation450kanno.ilmn12.hg19` for Illuminas 450K methylation arrays. After loading this file:

- Read reference methylation data
- Filter data to keep only CpG islands from `ref_met` and reduce computation time
- Further filtering to keep only CpG islands matching a gene name from a bulk reference

`meth_rna_mapping.ipynb` produces a mapping between CpG sites and UCSC refgene names. The mapping is saved as `mapping_meth_rna.csv`

# Model

The main model that was submitted to the Codabench platform is `submission_script.py` which provides a Python-based pipeline for estimating proportions of different components in a biological mixture using RNA and methylation data. The program integrates multiple data modalities and combines them for accurate proportion estimation. The main steps are listed below:

- **Data Alignment and Filtering**:
  - Aligns RNA and methylation datasets based on shared features.
  - Filters features for variability 

- **Proportion Estimation**:
  - Uses **Non-Negative Least Squares (NNLS)** to estimate component proportions in mixtures.
  - Supports RNA, methylation, and pseudo-bulk RNA datasets.

- **Optimal Weighting**:
  - Combines results from RNA and methylation datasets.
  - Finds optimal weights to minimize RMSE using `additionnal_script.py` 

## Usage

Prepare Input Data:
- RNA and methylation mixture data (mix_rna, mix_met).
- RNA and methylation reference data (ref_rna, ref_met).
- Mapping file linking RNA and methylation features (mapping_meth_rna.csv).
- Pseudo-bulk RNA data (peng_pseudo_bulk_sum.csv).

Run the program:
python `submission_script.py`

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

- It is difficult to make conclusions about which is the best-performing deconvolution algorithm given the limited available time in the competition.
- Feature selection is the key to improve performance
- Given our teams expertise in python, we preferred the usage of python rather than R (the InMoose python package, developed by Epigene Labs is a good example)

# Future ideas

- Test different pre-processing methods (reduce noise, batch correction and integrate single cell data, TMM or DESeq2 normalisation)
- Improve Ensemble methods for late integration
- Use M-values for methylation instead of beta-values
- Test feature selection with biological a priori
    - Gene set enrichment analysis using msgdib and hallmark of epigenetic
    - Use genes signatures allowing molecular classification of the cancer (and predict based on subtypes)
- Generalize pipeline to other cancer types


## Special thank you to the Data Challenge Organizers!

# References

[^1]: Epigenomic Deconvolution of Breast Tumors Reveals Metabolic Coupling between Constituent Cell Types: 10.1016/j.celrep.2016.10.057
[^2]: Codabench: Flexible, easy-to-use, and reproducible meta-benchmark platform: 10.1016/j.patter.2022.100543
[^3]: SCDC: bulk gene expression deconvolution by multiple single-cell RNA sequencing references: doi.org/10.1093/bib/bbz166
