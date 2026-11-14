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
setwd("my project directory")


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

length(Tumor_sample_list)
length(Normal_sample_list)


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
  directory = "D:/TCGA_BRCA",
  files.per.chunk = 8
)


################################################################################
# 13. Prepare RNA-seq Data
################################################################################

data_prepare <- GDCprepare(
  query_TCGA_2,
  directory = "D:/TCGA_BRCA",
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

data <- assay(data_prepare)

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

# Read one downloaded STAR Counts file.
# The object 'files' must contain the downloaded STAR Counts file paths.
star_annotation <- read.delim(
  files[1],
  skip = 1,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# Keep only actual gene records.
star_annotation <- star_annotation[
  grepl("^ENSG", star_annotation$gene_id),
]

# Extract gene annotation columns.
annotation <- star_annotation[, c(
  "gene_id",
  "gene_name",
  "gene_type"
)]

# Set Ensembl IDs as row names.
rownames(annotation) <- annotation$gene_id

# Match annotation order with the expression matrix.
annotation <- annotation[rownames(data), ]

# Keep gene symbol and biotype.
annot_df <- annotation[, c("gene_name", "gene_type")]

# Rename annotation columns.
colnames(annot_df) <- c("SYMBOL", "BioType")

# Inspect annotation.
dim(annot_df)
head(annot_df)

# Check missing values.
sum(is.na(annot_df$SYMBOL))
sum(is.na(annot_df$BioType))

# Export annotation.
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
hist(data, breaks = 1000000, xlim = c(0, 50))

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
#metadata

Selected_samples

read.table("selected samples.txt", header = T, sep = "\t")

Samples_type <- c(rep("Normal",4), rep("Tumor",4))

Meta_data <- data.frame(Selected_samples,Samples_type)

View(Meta_data)
head(Meta_data)

Meta_data <- data.frame(
  Sample_ID = Selected_samples$V30,
  Samples_type = Samples_type
)

View(Meta_data)
head(Meta_data)

colnames(Meta_data) <- c("SampleName","SampleType")

# Match sample order with data

Sample_ID <- substr(colnames(data), 1, 16)

Meta_data <- Meta_data[
  match(Sample_ID, Meta_data$SampleName),
]

all(Sample_ID == Meta_data$SampleName)
write.table(Meta_data,"Metadata.txt", quote = F,sep = "\t",row.names = F)

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

# Extract raw counts
mRNA_Raw <- mRNA_df[, -9:-10]
View(mRNA_Raw)

# Convert counts to numeric
mRNA_Raw[,1]
mRNA_Raw[,1] <- as.numeric(mRNA_Raw[,1])

for (i in 1:ncol(mRNA_Raw)) {
  mRNA_Raw[,i] <- as.numeric(mRNA_Raw[,i])
}


# Read BRCA1 expression data
BRCA1_mRNA <- read.table(
  "BRCA1_mRNA.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

View(BRCA1_mRNA)

# Extract BRCA1 raw counts
BRCA1_Raw <- BRCA1_mRNA[, -9:-10]
View(BRCA1_Raw)

# Convert counts to numeric
BRCA1_Raw[,1]
BRCA1_Raw[,1] <- as.numeric(BRCA1_Raw[,1])

for (i in 1:ncol(BRCA1_Raw)) {
  BRCA1_Raw[,i] <- as.numeric(BRCA1_Raw[,i])
}


# Read sample metadata
Meta_data <- read.table(
  "Metadata.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

View(Meta_data)


# DESeq2
library(DESeq2)

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


# rlog Transformation
data_rlogd1 <- rlogTransformation(d1)

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

library(DESeq2)
