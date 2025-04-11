library(biomaRt)
library(DESeq2)
library(ggplot2)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ComplexHeatmap)
library(dplyr)
################################################################################

time_point <- '24h'
sub_dir <- 'UK4'
workdir <- paste0('/home/pospim/Desktop/work/bioinformatics/datasets/SARS/GSE213759_calu3/IC19/',sub_dir, '/',  time_point, "/")
gse_part <- sub(".*/(GSE[^_/]+).*", "\\1", workdir)
res_table <- paste0(gse_part, "-",sub_dir, '-', time_point, ".csv")

## Load data
setwd(workdir)

## Read data
data <- read.table("counts_named.csv", header = TRUE, row.names = 1, comment.char = "", sep="\t")
coldata <- read.table("samplefile.csv", header = T)
contr <- read.table("contrasts.csv", header = T)
rownames(coldata) <- coldata$sample
coldata

control <- "UK4_Calu3_Mock"
virus <- "UK4_Calu3_IC19"

################################################################################

## Process data
# Remove version numbers 
rowNames <- rownames(data)
rowNames <- sub("\\..*", "", rowNames)
rowNames <- make.unique(rowNames)
rownames(data) <- rowNames

all(rownames(coldata) == colnames(data))
data <- round(data)

## Create DESeq2 data object
ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = data,
  colData = coldata,
  design = ~condition
)
ddsFullCountTable

## Differential gene expression analysis
dds <- DESeq(ddsFullCountTable)
dds
################################################################################

## Change virus names
contrast <- c("condition",virus, control)
res.deseq2 <- results(dds, contrast)
res.deseq2 <- as.data.frame(res.deseq2)

## Annotation
res.deseq2$gene_symbol <- mapIds(org.Hs.eg.db,
                                 keys = row.names(res.deseq2),
                                 column = "SYMBOL",
                                 keytype = "ENSEMBL",
                                 multiVals = "first"
)

## clean results
gsea_input <- res.deseq2[!is.na(res.deseq2$log2FoldChange), ]
gsea_input <- gsea_input[!is.na(gsea_input$pvalue), ]

colnames(gsea_input)[colnames(gsea_input) == "log2FoldChange"] ="log2fold"

data_sets_dir <- "/home/pospim/Desktop/work/bioinformatics/proteomics/GSEA/datasets/"
write.table(gsea_input, file = paste0(data_sets_dir, res_table), sep = "\t", quote = FALSE, row.names = FALSE)

