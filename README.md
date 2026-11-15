# TCGA-BRCA: Breast Cancer Genomics Analysis

## Project Overview

This project explores breast cancer RNA-seq data from **The Cancer Genome Atlas (TCGA-BRCA)** using R and Bioconductor.

The analysis focuses initially on **BRCA1** and follows a transcriptomic workflow, from data retrieval and preprocessing to exploratory differential gene expression analysis and visualization.

This is an independent bioinformatics learning and portfolio project. The current analysis is exploratory and is not intended to establish clinical or biological conclusions.

## Objectives

- Retrieve RNA-seq data from the TCGA-BRCA project.
- Explore and preprocess gene expression count data.
- Annotate genes using Ensembl identifiers and filter low-expression genes.
- Separate different RNA types, including mRNA, lncRNA, and miRNA.
- Extract and investigate BRCA1 expression data.
- Prepare sample metadata and match it with the expression matrix.
- Construct DESeq2 datasets and perform exploratory differential gene expression analysis.
- Annotate differential expression results with gene symbols.
- Apply regularized log (rlog) transformation for data exploration.
- Visualize expression data using bar plots, a volcano plot, a heatmap, and PCA.

## Tools & Technologies

- R
- Bioconductor
- TCGAbiolinks
- SummarizedExperiment
- DESeq2
- ggplot2
- EnhancedVolcano
- pheatmap
- gplots
- Git & GitHub

## Data Source

**The Cancer Genome Atlas (TCGA-BRCA)**

RNA-seq data are retrieved from the Genomic Data Commons (GDC) using TCGAbiolinks. The workflow uses the **STAR - Counts** gene expression data.

The current working subset contains:

- 4 Primary Tumor samples
- 4 Solid Tissue Normal samples

The sample metadata are classified using TCGA barcode sample-type codes:

- `01`: Primary Tumor
- `11`: Solid Tissue Normal

The subset is intentionally small and is being used to develop and demonstrate the analysis workflow.

## Analysis Workflow

### 1. Data Acquisition

- Explore available GDC projects and the TCGA-BRCA project summary.
- Query TCGA-BRCA RNA-seq data.
- Select four tumor and four normal samples for the current workflow.
- Download and prepare the selected expression data.

### 2. Data Preprocessing and Annotation

- Extract raw gene expression counts.
- Create a gene annotation table using gene symbols and gene biotypes.
- Filter genes with low total counts.
- Separate protein-coding genes, lncRNA, and miRNA.
- Extract BRCA1 expression data.

### 3. Sample Metadata

- Derive tumor and normal sample labels from TCGA barcode sample-type codes.
- Match metadata to the expression matrix sample order.
- Check sample alignment before constructing DESeq2 objects.

### 4. DESeq2 Analysis

- Construct DESeqDataSet objects for the mRNA expression matrix and BRCA1.
- Use `SampleType` as the design variable.
- Estimate size factors.
- Apply rlog transformation to the mRNA dataset.
- Run the DESeq2 differential expression workflow.
- Adjust p-values using the false discovery rate (FDR) method.
- Annotate Ensembl gene identifiers with gene symbols.
- Sort and export the annotated differential expression results.

### 5. Data Visualization

The current R script includes the following exploratory visualizations:

- **RNA type distribution bar plot:** summarizes the number of genes assigned to each RNA biotype.
- **rlog column-sum bar plot:** displays the sum of rlog-transformed expression values for each sample. This is a descriptive visualization, not a library-size normalization or sample-quality metric.
- **Volcano plot:** displays log2 fold change against p-value and highlights selected gene labels.
- **Heatmap:** displays rlog-transformed expression for up to the top 100 genes from the ranked differential expression results.
- **PCA plot:** uses rlog-transformed expression data to explore sample-level expression patterns.

The figures are generated locally when the script is run. They are not represented as validated biological findings.

## Repository Contents

- `README.md`: project overview, workflow, tools, and limitations.
- `TCGA_BRCA_BRCA1_analysis.R`: main R script for data acquisition, preprocessing, DESeq2 analysis, annotation, and visualization.
- `.gitignore` and `.gitattributes`: repository configuration files.

The script also generates intermediate data tables and visualization files in the working project directory when executed.

## How to Run

1. Clone or download this repository.
2. Open the project folder in RStudio.
3. Install the required R and Bioconductor packages if they are not already installed.
4. Open and run `TCGA_BRCA_BRCA1_analysis.R`.

The script uses the current working directory as the project directory and creates a local `GDCdata` folder for downloaded data. Running the data acquisition steps requires an internet connection and access to the GDC.

The analysis script is available here:

[TCGA_BRCA_BRCA1_analysis.R](TCGA_BRCA_BRCA1_analysis.R)

## Current Project Status

**In Progress — exploratory RNA-seq workflow and visualization implemented.**

Completed in the current workflow:

- TCGA-BRCA RNA-seq query and data retrieval
- Selection of four tumor and four normal samples
- Raw count extraction and preprocessing
- Gene annotation and low-count filtering
- RNA biotype classification
- BRCA1 expression extraction
- Barcode-based sample metadata preparation
- DESeq2 dataset construction and size-factor estimation
- rlog transformation
- Exploratory differential expression analysis and gene-symbol annotation
- RNA-type distribution visualization
- rlog column-sum bar plot
- Volcano plot
- Top-gene heatmap
- PCA plot

Further work:

- Review the analysis design and sample independence.
- Repeat the analysis using a larger, appropriately selected TCGA-BRCA cohort.
- Use a patient-aware design if paired or repeated samples are included.
- Perform more detailed BRCA1-focused analyses.
- Explore associations with clinical variables after appropriate data preparation.
- Consider survival analysis and independent validation in a future stage.

## Limitations and Interpretation

The current working subset contains only **8 samples (4 tumor and 4 normal)**. In addition, the subset includes multiple samples/aliquots from at least one patient. Therefore, the eight columns should not be assumed to represent eight independent patients.

The current DESeq2 design uses `~SampleType` and does not model patient-level pairing or repeated samples. With this small and non-independent subset, differential expression results, volcano plot labels, heatmap patterns, and PCA separation are **exploratory only**. They should not be interpreted as statistically robust, clinically meaningful, or generalizable findings.

A larger cohort, careful sample selection, a patient-aware statistical design where appropriate, and independent validation are needed before drawing biological or clinical conclusions.

---

*This repository documents an ongoing learning project in cancer genomics, transcriptomics, and bioinformatics.*
