################################################################################
# Project: BRCA1 Gene Expression Analysis Using TCGA-BRCA Data
# Research Area: Cancer Genomics, Transcriptomics & Bioinformatics
# Target Gene: BRCA1
# Researcher: Yazdan Mottaghi
# Data Source: Genomic Data Commons (GDC)
# Project ID: TCGA-BRCA
################################################################################


# 1. Load Required Packages
library(TCGAbiolinks)
library(SummarizedExperiment)
library(ggplot2)


# 2. Set Working Directory
# Open this script from the project folder in RStudio, or set the
# working directory to the project folder before running the script.
project_dir <- getwd()
gdc_dir <- file.path(project_dir, "GDCdata")
dir.create(gdc_dir, showWarnings = FALSE, recursive = TRUE)


################################################################################
# 3. Explore Available GDC Projects
################################################################################

All_GDC_Projects <- getGDCprojects()

"TCGA-BRCA" %in% All_GDC_Projects$project_id

View(All_GDC_Projects)
All_GDC_Projects[, 1]


################################################################################
# 4. Explore TCGA-BRCA Project Summary
################################################################################

All_TCGA_BRCA_Data <- getProjectSummary("TCGA-BRCA")

All_TCGA_BRCA_Data$data_categories


################################################################################
# 5. Query TCGA-BRCA Gene Expression Data
################################################################################

query_TCGA_BRCA <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)


################################################################################
# 6. Retrieve Query Results
################################################################################

results <- getResults(query_TCGA_BRCA)

View(results)


################################################################################
# 7. Identify Normal Tissue Samples
################################################################################

Normal_sample_list <- subset(
  results,
  sample_type == "Solid Tissue Normal" &
    access == "open" &
    experimental_strategy == "RNA-Seq"
)

View(Normal_sample_list)

write.table(
  Normal_sample_list,
  file = "all normal samples of BRCA.txt",
  quote = FALSE,
  sep = "\t"
)

nrow(Normal_sample_list)


################################################################################
# 8. Identify Primary Tumor Samples
################################################################################

Tumor_sample_list <- subset(
  results,
  sample_type == "Primary Tumor" &
    access == "open" &
    experimental_strategy == "RNA-Seq"
)

View(Tumor_sample_list)

write.table(
  Tumor_sample_list,
  file = "all tumor samples of BRCA.txt",
  quote = FALSE,
  sep = "\t"
)

nrow(Tumor_sample_list)

nrow(Tumor_sample_list)
nrow(Normal_sample_list)


################################################################################
# 9. Select Four Tumor and Four Normal Samples
################################################################################

sample_T_4 <- Tumor_sample_list[1:4, ]
sample_N_4 <- Normal_sample_list[1:4, ]

# Keep the same sample order used in the original selection.
Selected_samples <- rbind(sample_N_4, sample_T_4)

View(Selected_samples)


################################################################################
# 10. Export Selected Samples
################################################################################

write.table(
  Selected_samples,
  file = "selected samples.txt",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)


################################################################################
# 11. Create a GDC Query for Selected Samples
################################################################################

query_TCGA_2 <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = Selected_samples$sample.submitter_id
)


################################################################################
# 12. Download RNA-seq Data
################################################################################

GDCdownload(
  query_TCGA_2,
  method = "api",
  directory = gdc_dir,
  files.per.chunk = 8
)


################################################################################
# 13. Prepare RNA-seq Data
################################################################################

data_prepare <- GDCprepare(
  query_TCGA_2,
  directory = gdc_dir,
  save = TRUE,
  save.filename = "data.rda"
)


################################################################################
# 14. Inspect the Prepared Dataset
################################################################################

class(data_prepare)
dim(data_prepare)


################################################################################
# 15. Extract Raw Gene Expression Counts
################################################################################

# Use unstranded raw counts when that assay is available.
if ("unstranded" %in% assayNames(data_prepare)) {
  data <- assay(data_prepare, "unstranded")
} else {
  data <- assay(data_prepare)
}

View(data)
dim(data)


################################################################################
# 16. Export Raw Count Matrix
################################################################################

write.table(
  data,
  file = "Raw data.txt",
  quote = FALSE,
  sep = "\t",
  col.names = NA
)


################################################################################
# 17. Create Gene Annotation File
################################################################################

# Extract gene annotation directly from the prepared SummarizedExperiment.
# This avoids relying on an external 'files' object or machine-specific paths.
annotation <- as.data.frame(rowData(data_prepare))

required_annotation_columns <- c("gene_name", "gene_type")
if (!all(required_annotation_columns %in% colnames(annotation))) {
  stop(
    "Expected gene_name and gene_type columns were not found in rowData(data_prepare). ",
    "Please inspect rowData(data_prepare) and update the annotation step."
  )
}

# Match annotation rows to the expression matrix.
annotation_index <- match(rownames(data), rownames(data_prepare))
if (anyNA(annotation_index)) {
  stop("Some expression-matrix genes could not be matched to the annotation.")
}
annotation <- annotation[annotation_index, , drop = FALSE]

# Keep gene symbol and biotype.
annot_df <- annotation[, c("gene_name", "gene_type"), drop = FALSE]
colnames(annot_df) <- c("SYMBOL", "BioType")
rownames(annot_df) <- rownames(data)

# Save annotation for the downstream annotation steps.
write.table(
  annot_df,
  file = "annot_df.txt",
  quote = FALSE,
  sep = "\t",
  col.names = NA
)


################################################################################
# 18. Explore and Filter Raw Counts
################################################################################

class(data)
summary(as.vector(data))

colnames(data)

boxplot(data)
hist(data)
hist(data, breaks = 50, xlim = c(0, 50))

data_filtered <- data[which(rowSums(data) > 10), ]

nrow(data)
nrow(data_filtered)

write.table(
  data_filtered,
  file = "data filtered.txt",
  quote = FALSE,
  sep = "\t",
  col.names = NA
)


################################################################################
# 19. Read Filtered Data and Annotation
################################################################################

annot_df <- read.table(
  "annot_df.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

data_filtered <- read.table(
  "data filtered.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

View(annot_df)
View(data_filtered)


################################################################################
# 20. Add Gene Annotation to Filtered Data
################################################################################

gene_name <- annot_df[rownames(data_filtered), "SYMBOL"]
gene_type <- annot_df[rownames(data_filtered), "BioType"]

data_filtered2 <- data_filtered

data_filtered2 <- cbind(data_filtered2, gene_name)
data_filtered2 <- cbind(data_filtered2, gene_type)

data_filtered2 <- as.data.frame(data_filtered2)

class(data_filtered2)
table(data_filtered2$gene_type)

View(data_filtered2)


################################################################################
# 21. Separate RNA Types
################################################################################

mRNA_df <- subset(
  data_filtered2,
  gene_type == "protein_coding"
)

View(mRNA_df)
nrow(mRNA_df)

lncRNA_df <- subset(
  data_filtered2,
  gene_type == "lncRNA"
)

View(lncRNA_df)
nrow(lncRNA_df)

miRNA_df <- subset(
  data_filtered2,
  gene_type == "miRNA"
)

View(miRNA_df)
nrow(miRNA_df)


################################################################################
# 22. Extract BRCA1 Expression Data
################################################################################

"BRCA1" %in% data_filtered2$gene_name

annot_df["ENSG00000012048.23", ]

BRCA1_mRNA <- subset(
  data_filtered2,
  gene_type == "protein_coding" & gene_name == "BRCA1"
)

View(BRCA1_mRNA)
nrow(BRCA1_mRNA)
dim(BRCA1_mRNA)
head(BRCA1_mRNA)


################################################################################
# 23. Export Expression Data
################################################################################

write.table(
  BRCA1_mRNA,
  file = "BRCA1_mRNA.txt",
  quote = FALSE,
  sep = "\t"
)

write.table(
  lncRNA_df,
  file = "lncRNA_df.txt",
  quote = FALSE,
  sep = "\t"
)

write.table(
  miRNA_df,
  file = "miRNA_df.txt",
  quote = FALSE,
  sep = "\t"
)

write.table(
  mRNA_df,
  file = "mRNA_df.txt",
  quote = FALSE,
  sep = "\t"
)


################################################################################
# 24. Plot RNA Type Distribution
################################################################################

RNA_table <- table(data_filtered2$gene_type)
RNA_table <- as.data.frame(RNA_table)

colnames(RNA_table) <- c("Type", "Number")

View(RNA_table)

Plot_RNA_table <- ggplot(
  RNA_table,
  aes(x = Type, y = Number)
) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(
    title = "Barplot",
    x = "RNA Type",
    y = "Number"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

ggsave(
  "RNA table info.png",
  plot = Plot_RNA_table,
  width = 15,
  height = 7
)


################################################################################
# 25. Create Sample Metadata
################################################################################
# Build metadata directly from the count-matrix column names.
# TCGA barcode sample-type codes: 01 = Primary Tumor; 11 = Solid Tissue Normal.
sample_barcodes <- colnames(data)
sample_type_codes <- substr(sample_barcodes, 14, 15)

Samples_type <- ifelse(
  sample_type_codes == "01", "Tumor",
  ifelse(sample_type_codes == "11", "Normal", NA_character_)
)

if (anyNA(Samples_type)) {
  stop(
    "One or more sample barcodes could not be classified as tumor (01) or normal (11). ",
    "Please inspect colnames(data) and the TCGA sample barcodes."
  )
}

Meta_data <- data.frame(
  SampleName = sample_barcodes,
  SampleType = factor(Samples_type, levels = c("Normal", "Tumor")),
  row.names = sample_barcodes
)

# Confirm metadata order matches the expression matrix.
stopifnot(identical(rownames(Meta_data), colnames(data)))

View(Meta_data)
head(Meta_data)

write.table(
  Meta_data,
  "Metadata.txt",
  quote = FALSE,
  sep = "\t",
  row.names = TRUE,
  col.names = NA
)

#######################################################################################
# DEG

# Read mRNA data
mRNA_df <- read.table(
  "mRNA_df.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

View(mRNA_df)

# Extract raw counts in the exact metadata order.
mRNA_Raw <- as.matrix(mRNA_df[, Meta_data$SampleName, drop = FALSE])
storage.mode(mRNA_Raw) <- "integer"
View(mRNA_Raw)

# Read BRCA1 expression data
BRCA1_mRNA <- read.table(
  "BRCA1_mRNA.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

View(BRCA1_mRNA)

# Extract BRCA1 raw counts in the exact metadata order.
BRCA1_Raw <- as.matrix(BRCA1_mRNA[, Meta_data$SampleName, drop = FALSE])
storage.mode(BRCA1_Raw) <- "integer"
View(BRCA1_Raw)


# Read sample metadata
Meta_data <- read.table(
  "Metadata.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# DESeq2 requires categorical design variables to be factors.
Meta_data$SampleType <- factor(
  Meta_data$SampleType,
  levels = c("Normal", "Tumor")
)

# Confirm metadata and count matrix contain the same samples in the same order.
stopifnot(identical(rownames(Meta_data), colnames(mRNA_Raw)))
stopifnot(!anyDuplicated(rownames(Meta_data)))

View(Meta_data)


# DESeq2
library(DESeq2)

# Confirm count and metadata dimensions and sample order.
stopifnot(identical(colnames(mRNA_Raw), rownames(Meta_data)))
stopifnot(all(is.finite(mRNA_Raw)), all(mRNA_Raw >= 0))

# Create DESeqDataSet for all mRNA genes
d1 <- DESeqDataSetFromMatrix(
  countData = mRNA_Raw,
  colData = Meta_data,
  design = ~SampleType
)

dim(d1)
colData(d1)
design(d1)


# Create DESeqDataSet for BRCA1
db <- DESeqDataSetFromMatrix(
  countData = BRCA1_Raw,
  colData = Meta_data,
  design = ~SampleType
)

dim(db)
colData(db)
design(db)


# Estimate size factors before transformation.
d2 <- estimateSizeFactors(d1)

# rlog Transformation
data_rlogd1 <- rlogTransformation(d2, blind = TRUE)

library(SummarizedExperiment)

data_rlogd1_df <- data.frame(assay(data_rlogd1))

# Export rlog-transformed data
write.table(
  data_rlogd1_df,
  "d1 rlog.txt",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

# Differential Gene Expression Analysis (DEG)
# IMPORTANT: This 8-sample subset is for workflow practice only.
# The small sample size and repeated samples/aliquots from some patients
# mean this design is not suitable for definitive statistical conclusions.
# A larger cohort and an appropriately specified patient-aware design are
# required for a robust analysis.

library(DESeq2)

d3 <- DESeq(d2)

deseq_result <- results(d3, pAdjustMethod = "fdr")

class(deseq_result)

final_deseq_results <- data.frame(deseq_result)

View(final_deseq_results)

write.table(
  final_deseq_results,
  "DEG mRNA.txt",
  quote = F,
  sep = "\t"
)

##########################################################################

# Gene Annotation

DEG <- read.table(
  "DEG mRNA.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)
View(DEG)

annot_df <- read.table(
  "annot_df.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
colnames(annot_df) <- c("Symbol", "Biotype")
View(annot_df)

# Annotation Function

annot_fx <- function(x) {
  k <- annot_df[x, 1]
  return(k)
}

annot_fx("ENSG00000002834.1") # Example

gene_names <- annot_fx(rownames(DEG))

length(gene_names) == nrow(DEG)

final_result_annot <- cbind(gene_names)

final_result_annot <- cbind(
  final_result_annot,
  DEG[1:6]
)

rownames(final_result_annot) <- rownames(DEG)

colnames(final_result_annot)[1] <- "Symbol"

View(final_result_annot)

# Sort by padj value

final_result_annot <- final_result_annot[
  order(final_result_annot$padj, decreasing = F),
]

View(final_result_annot)

# Export Final Results

write.table(
  final_result_annot,
  "Final result.txt",
  quote = F,
  sep = "\t"
)
