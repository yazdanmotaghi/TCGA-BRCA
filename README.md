# TCGA-BRCA: Breast Cancer Genomics Analysis

## Project Overview

This project explores breast cancer RNA-seq data from **The Cancer Genome Atlas (TCGA-BRCA)** using R and Bioconductor.

The analysis focuses initially on **BRCA1** and follows a transcriptomic workflow, from data retrieval and preprocessing to exploratory differential gene expression analysis using DESeq2.

## Objectives

- Retrieve RNA-seq data from the TCGA-BRCA project.
- Explore and preprocess gene expression count data.
- Annotate genes using Ensembl identifiers and filter low-expression genes.
- Separate different RNA types, including mRNA, lncRNA, and miRNA.
- Extract and investigate BRCA1 expression data.
- Prepare sample metadata for downstream analysis.
- Construct DESeq2 datasets and perform differential gene expression analysis.
- Annotate differential expression results with gene symbols.
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

RNA-seq data were retrieved from TCGA using TCGAbiolinks.

The current working subset includes:

- 4 Primary Tumor samples
- 4 Solid Tissue Normal samples

This small subset is being used to develop and demonstrate the analysis workflow. It is not intended to support definitive biological or clinical conclusions.

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
- Separate protein-coding genes and other RNA types, including lncRNA and miRNA.

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
- Define the experimental design using `SampleType`.
- Estimate size factors for normalization.
- Apply regularized log (rlog) transformation to the mRNA dataset.

### 6. Differential Gene Expression Analysis

- Run the DESeq2 differential expression workflow on the mRNA dataset.
- Extract differential expression results using the false discovery rate (FDR) adjustment method.
- Convert the DESeq2 results object into a data frame.
- Annotate Ensembl gene identifiers with gene symbols.
- Sort the annotated results by adjusted p-value.
- Export the differential expression table for downstream review.

## Source Code

The R analysis script is available here:

The script uses the current working directory as the project directory and creates a local `GDCdata` folder for downloaded data. Open the project in RStudio (or set the working directory to the repository folder) before running it.

[TCGA_BRCA_BRCA1_analysis.R](TCGA_BRCA_BRCA1_analysis.R)

## Current Project Status

**In Progress**

Completed:

- TCGA-BRCA data retrieval
- Tumor and normal sample selection
- Raw count extraction and preprocessing
- Gene annotation and filtering
- RNA type classification
- BRCA1 expression extraction
- Sample metadata preparation
- DESeq2 dataset construction
- rlog transformation of the mRNA dataset
- Initial DESeq2 differential expression workflow
- Gene-symbol annotation of differential expression results
- Sorting and exporting the annotated results table

Upcoming:

- Review and interpret differential expression results
- Additional expression data visualization
- PCA and sample clustering
- Survival analysis
- Integration with clinical data
- Further investigation of breast cancer-related genes

## Limitations

The current working subset contains only 8 samples (4 tumor and 4 normal) and is intended for workflow development and educational purposes. The small sample size—and the presence of multiple aliquots/samples from some patients in this working subset—limits the reliability and generalizability of statistical findings. The current `~SampleType` design does not model patient-level pairing or repeated samples. Therefore, the differential expression output is exploratory and must not be interpreted as a statistically robust result. A larger cohort with a patient-aware design, followed by appropriate validation, is required before drawing biological or clinical conclusions.

