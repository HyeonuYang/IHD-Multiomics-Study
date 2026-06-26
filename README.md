# Gut Microbiome and Metabolome Signatures Across the Ischemic Heart Disease Severity Spectrum: A Multi-Omics Study

![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)

## 📌 Project Overview
This repository contains the complete analytical pipeline for investigating the association between the gut microbiome, blood metabolome, and the severity of Ischemic Heart Disease (IHD). The study categorizes IHD severity into Stable Angina (SA), Acute Coronary Syndrome (ACS), and Heart Failure with reduced Ejection Fraction (HFrEF), employing a comprehensive multi-omics and ordinal differential analysis approach.

> **Note:** This code is associated with the paper *"Gut Microbiome and Metabolome Signatures Across the Ischemic Heart Disease Severity Spectrum: A Multi-Omics Study"* (Under Review).
> 
> **Co-first Authors:** Chihoon Park & Hyun-Woo Yang

## 📂 Repository Structure & Workflow

### Part 1: Microbiome Analysis (Authored by Chihoon Park)
The microbiome preprocessing and profiling pipeline relies on shotgun metagenomic sequencing data, encompassing quality control, taxonomic/functional profiling, and microbial risk score calculations.

* **`01_Microbiome/`**: 
  * `01_quality_filtering.sh`: Adapter trimming and quality filtering using **Trimmomatic**.
  * `02_removal_of_human_reads.sh`: Depletion of host reads using **Bowtie2** and **Samtools**.
  * `03_taxonomic_and_functional_profiling.sh` & `03_GMM_table.R`: Profiling via **HUMAnN3/MetaPhlAn**, and Curated Gut Metabolic Modules (GMMs) utilizing `omixerRpm`.
  * `04_statistical_analysis.R`: Alpha/Beta diversity (`vegan`), pathway analysis, and Spearman correlations.
  * `05_calculation_of_MRS.R` & `05_calculation_of_PRS.sh`: Calculation of Abundance/Presence-based Microbial Risk Scores (aMRS, paMRS) and Polygenic Risk Scores (PRS) via **PLINK**.

### Part 2: Metabolomics & Multi-Omics Integration (Authored by Hyun-Woo Yang)
This section covers blood metabolome quality control, ordinal differential analysis, complex multi-omics visualizations, and the development of cumulative ordinal risk models.

* **`02_Metabolomics/`**:
  * `01_metabolomics_QC.ipynb`: Robust quality control of metabolomics data, including coefficient of variation (CV) filtering, Limit of Detection (LOD) cutoffs, and feature-wise minimum imputation.
  * `02_run_maaslin3.R`: Ordinal differential abundance analysis across IHD severity stages (SA -> ACS -> HFrEF) adjusting for clinical covariates (Age, Sex, BMI, etc.) using **MaAsLin3**.
* **`03_Circular_heatmap/`**:
  * `01_data_for_circular_heatmap.ipynb`: Aggregation of multi-omics significant features (Metabolites, Taxa, GMM, MetaCyc) and calculation of group-wise mean abundances.
  * `02_circular_heatmap.R`: Construction of an advanced multi-track circular heatmap using `circlize` and `ComplexHeatmap` to visualize severity-associated features.
* **`04_Correlation/`**:
  * `01_residual_correlation.R`: Computation of stage-residualized Spearman correlations to identify intrinsic interactions between gut microbiome features and circulating metabolites.
* **`05_Ordinal_model/`**:
  * `01_calc_MetRS.ipynb`: Derivation of the Metabolomic Risk Score (MetRS) using a weighted sum of significant metabolites and Z-score normalization.
  * `02_OrdinalModel.ipynb`: Development of cumulative Multi-omics Risk Prediction Models utilizing Ordinal Logistic Regression (`statsmodels`). Evaluates Clinical Risk Factors (CRFs), PRS, MRS, and MetRS via Stratified K-Fold CV and Bootstrap AUC.
  * `03_quantile_plots.ipynb`: Evaluation and visualization of Odds Ratios (OR) across risk score quantiles to assess clinical impact.

## 🛠 Prerequisites & Dependencies

### Software & Bioinformatics Tools
* Trimmomatic (v0.39), Bowtie2, Samtools (v1.20), HUMAnN3, MetaPhlAn, PLINK

### R Packages
* `vegan`, `omixerRpm`, `maaslin3`, `circlize`, `ComplexHeatmap`, `dplyr`, `tidyr`, `tibble`, `readr`, `stringr`

### Python Libraries
* `pandas`, `numpy`, `statsmodels`, `scikit-learn`, `matplotlib`, `seaborn`, `openpyxl`

## 🚀 How to Use
1. Clone the repository:
   ```bash
   git clone [https://github.com/HyeonuYang/IHD-Multiomics-Study.git](https://github.com/HyeonuYang/IHD-Multiomics-Study.git)
