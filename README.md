
# TCGA-BRCA: Breast Cancer Genomics Analysis

## Project Overview

This project focuses on the analysis of breast cancer RNA-seq data from The Cancer Genome Atlas (TCGA-BRCA).

The main goal is to explore gene expression patterns, perform transcriptomic data processing, and investigate cancer-related genes, with an initial focus on **BRCA1**.

The project is being developed using R and Bioconductor.

## Objectives

- Retrieve RNA-seq data from the TCGA-BRCA project.
- Explore and preprocess gene expression count data.
- Perform gene annotation and filtering.
- Separate different RNA types, including mRNA, lncRNA, and miRNA.
- Extract and investigate BRCA1 expression data.
- Prepare sample metadata for downstream analysis.
- Build a DESeq2 dataset for differential gene expression analysis.
- Apply rlog transformation for data exploration and visualization.
- Explore gene expression patterns and their potential clinical relevance.

## Tools & Technologies

- R
- Bioconductor
- TCGAbiolinks
- DESeq2
- ggplot2
- SummarizedExperiment
- Git & GitHub

## Data Source

**The Cancer Genome Atlas (TCGA-BRCA)**

RNA-seq data were retrieved using TCGAbiolinks.

The current working dataset includes:
- 4 Primary Tumor samples
- 4 Solid Tissue Normal samples

This small sample subset is primarily used for learning and developing the analysis workflow.

## Analysis Workflow

### 1. Data Acquisition
- Explore TCGA projects.
- Retrieve TCGA-BRCA RNA-seq data.
- Select tumor and normal samples.
- Download and prepare expression data.

### 2. Data Preprocessing
- Extract raw gene expression counts.
- Filter low-expression genes.
- Annotate genes using Ensembl identifiers.
- Separate protein-coding genes and other RNA types.

### 3. BRCA1 Gene Exploration
- Extract BRCA1 expression data.
- Identify its Ensembl gene ID and gene biotype.
- Prepare BRCA1-specific expression data.

### 4. Sample Metadata
- Identify tumor and normal samples using TCGA barcodes.
- Prepare sample metadata.
- Match metadata with expression matrix columns.

### 5. DESeq2 Data Preparation
- Construct DESeqDataSet objects.
- Define the experimental design using SampleType.
- Apply regularized log (rlog) transformation to the mRNA dataset.

## Project Structure

```text
TCGA-BRCA/
│
├── data/
│   ├── raw/
│   └── processed/
│
├── metadata/
│
├── scripts/
│
├── results/
│
├── figures/
│
└── README.md
```

## Current Project Status

**In Progress**

Completed:
- TCGA-BRCA data retrieval
- Sample selection
- Raw count extraction
- Gene annotation and filtering
- RNA type classification
- BRCA1 expression extraction
- Sample metadata preparation
- DESeq2 dataset construction
- rlog transformation

Upcoming:
- Differential gene expression analysis (DEG)
- Expression data visualization
- PCA and sample clustering
- Survival analysis
- Integration with clinical data
- Further investigation of breast cancer-related genes

## Notes

The current dataset contains only 8 samples and is intended for workflow development and educational purposes. Larger sample cohorts will be required for statistically robust biological conclusions.

