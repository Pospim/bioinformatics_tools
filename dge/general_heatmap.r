## LOAD PACKAGES ##
BiocManager::install(version = "3.20", update = TRUE)
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  library(BiocManager)
}
if (!require(DESeq2)) {
  BiocManager::install("DESeq2")
}
if (!require(biomaRt)) {
  BiocManager::install("biomaRt", force = TRUE)
}
if (!require(ggplot2)) {
  BiocManager::install("ggplot2")
}
if (!require(AnnotationDbi)) {
  BiocManager::install("AnnotationDbi")
}
if (!require(org.Hs.eg.db)) {
  BiocManager::install("org.Hs.eg.db")
}
if (!require(ComplexHeatmap)) {
  BiocManager::install("ComplexHeatmap")
}
if (!require(dplyr)) {
  BiocManager::install("dplyr")
}

library(biomaRt)
library(DESeq2)
library(ggplot2)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ComplexHeatmap)
library(dplyr)

###############################LOAD DATA#################################################
time_point <- ""
sub_dir <- ""
workdir <- paste0("/home/pospim/Desktop/work/bioinformatics/datasets/SARS/sandra_calu3/", sub_dir, '/', time_point, "/")
#gse_part <- sub(".*/(GSE[^_/]+).*", "\\1", workdir)
gse_part <- "sandra_calu3"

setwd(workdir)

## Read data
data <- read.table("counts_named.csv", header = TRUE, row.names = 1, comment.char = "", sep=",")
coldata <- read.table("samplefile.csv", header = T, sep=",")
contr <- read.table("contrasts.csv", header = T, sep=",")
rownames(coldata) <- coldata$sample

coldata

control <- "C3_24_ctrl"
virus <- "C3_24_MOI"
###############################PREPROCESS#################################################
# Remove version numbers 
#rowNames <- rownames(data)
#rowNames <- sub("\\..*", "", rowNames)
#rowNames <- make.unique(rowNames)
#rownames(data) <- rowNames
###############################VERIFY#################################################
all(rownames(coldata) == colnames(data))
data <- round(data)

# Check if any sum differs by more than 20% from the reference sum
sums <- colSums(data)
reference_sum <- max(sums)
if (any(abs(sums - reference_sum) / reference_sum > 0.50)) {
  warning("One or more column sums deviate by more than 50% from the highest sum.")
}
print(sums)
###############################DESEQ2#################################################
# Create DESeq2 data object
ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = data,
  colData = coldata,
  design = ~condition
)
ddsFullCountTable

## Differential gene expression analysis
dds <- DESeq(ddsFullCountTable)
dds

## Change virus names
contrast <- c("condition",virus, control)
res.deseq2 <- results(dds, contrast)
res.deseq2 <- as.data.frame(res.deseq2)

###############################ANNOTATION#################################################
res.deseq2$gene_symbol <- mapIds(org.Hs.eg.db,
                                 keys = row.names(res.deseq2),
                                 column = "SYMBOL",
                                 keytype = "ENSEMBL",
                                 multiVals = "first"
)
###############################COMBINED VIRAL ANNOTATION#################################################
Geneid <- rownames(data)
non_ensg_rows <- subset(data, !grepl("ENSG", Geneid))
rownames(non_ensg_rows)

sarscov2_gene_map <- c(
  "ENSSASG00005000002" = "ORF1ab",
  #"ENSSASG00005000003" = "ORF1ab",
  "ENSSASG00005000004" = "S",
  "ENSSASG00005000006" = "ORF3a",
  "ENSSASG00005000010" = "E",
  "ENSSASG00005000007" = "M",
  "ENSSASG00005000011" = "ORF6",
  "ENSSASG00005000009" = "ORF7a",
  "ENSSASG00005000012" = "ORF7b",
  "ENSSASG00005000008" = "ORF8",
  "ENSSASG00005000005" = "N",
  "ENSSASG00005000013" = "ORF10"
)
res.deseq2$gene_symbol <- ifelse(
  row.names(res.deseq2) %in% names(sarscov2_gene_map),
  sarscov2_gene_map[row.names(res.deseq2)],
  res.deseq2$gene_symbol
)
mapped_rows <- res.deseq2[row.names(res.deseq2) %in% names(sarscov2_gene_map), ]
mapped_rows

###############################aGPCRs#################################################
res.deseq2 <- res.deseq2[, c(7, 1, 2, 3, 4, 5, 6)]
# Remove rows with NA gene_symbol
res.deseq2 <- res.deseq2[!is.na(res.deseq2$gene_symbol), ]

irs_genes <- c(
  # Sensing proteins
  "TLR1", "NOD", "NOD2", "NLRC4", "NLRC5", "LGR4", "POLR3C",
  
  # Regulation of signaling pathway
  "MYD88", "UCHL1", "MALTI", "TBK1",
  
  # Proteins involved in defense response to virus
  "DDX60", "MX1"
)

redox_genes <- c(
  # Mitochondrial functions
  "ALKBH1", "AIFM1", "CYB5R3", "CYB5A",  "DEGS1", 
  "HSD17B10", "MDH2", "PCBD2", "SDHD", "FDX1",
  
  # Mitochondrial respiratory chain complex I
  "NDUFA4", "NDUFS3", "NDUFS7", "NDUFA13", "NDUFA9", "NDUFS6",
  
  # Metabolic enzymes
  "KDSR", "IFI30", "DEGS2", "BLVRA", "DHRSX", "GRIF", 
  "GMPR2", "MOXD1", "PIPOX", "RDH8", "SCD5"
)

#res.irs <- filter(
#  res.deseq2,
#  gene_symbol %in% irs_genes
#)
res.irs <- res.deseq2[!is.na(res.deseq2$gene_symbol) & !is.na(res.deseq2$log2FoldChange), ] 
#up
res.irs <- res.irs[res.irs$log2FoldChange > 2 & res.irs$pvalue < 0.05, ]

res.irs <- res.irs[order(res.irs$gene_symbol), ]
res_table <- paste0(gse_part, "-", time_point, "-irs_analysis_results.csv")
write.table(res.irs, file = res_table, sep = "\t", quote = FALSE, row.names = FALSE)

###############################SAVE FOR GSEA#################################################
#gsea_input <- res.deseq2[!is.na(res.deseq2$log2FoldChange), ]
#gsea_input <- gsea_input[!is.na(gsea_input$pvalue), ]
#colnames(gsea_input)[colnames(gsea_input) == "pvalue"] ="pval"
#colnames(gsea_input)[colnames(gsea_input) == "log2FoldChange"] ="log2fold"

#data_sets_dir <- "/home/pospim/Desktop/work/bioinformatics/GSEA/datasets/"
#write.table(gsea_input, file = paste0(data_sets_dir, gse_part, "_", sub_dir,'_', time_point, ".csv"), sep = "\t", quote = FALSE, row.names = FALSE)

###############################HEATMAP ANALYSIS#################################################
dds2 <- estimateSizeFactors(dds)
norm.counts <- counts(dds2, normalized = TRUE)
norm.irs <- norm.counts[row.names(res.irs), ]

## Add gene symbols to the normalized counts
headers <- as.data.frame(res.irs[, 1])
rownames(headers) <- rownames(res.irs)
colnames(headers) <- "gene_symbol"

norm.irs <- merge(headers, norm.irs, by = "row.names", all = TRUE)
norm.irs <- norm.irs[!duplicated(norm.irs$gene_symbol), ]

rownames(norm.irs) <- norm.irs$gene_symbol
norm.irs <- norm.irs[, -c(1:2)]

## Calculate Z-scores
norm.irs.z <- t(apply(norm.irs, 1, scale))
colnames(norm.irs.z) <- colnames(norm.irs)

norm.irs.z
norm.irs.z <- as.data.frame(norm.irs.z) %>%
  mutate_all(~replace(., is.nan(.), 0))

norm.irs.z <- as.matrix(norm.irs.z)

#sort norm.irs.z by gene_symbol
norm.irs.z <- norm.irs.z[order(rownames(norm.irs.z)), ]
norm.irs.z

###############################HEATMAPS#################################################
## Create heatmap without ordering
title <- "Innate imune response genes"

ht0 <- Heatmap(norm.irs.z,
               name = "Z-score", cluster_rows = F,
               column_title = "conditions",
               row_title = title,
               cluster_columns = F
)
ht0 <- draw(ht0)

## Create pdf
file_name_unorg=paste("heatmap-not_organized-irs_", time_point, ".pdf", sep = "")
pdf(file=file_name_unorg, height = 8)
ht0 <- draw(ht0) 
dev.off()

## Define column splitting variables
column_split <- rep("group", ncol(norm.irs.z)) # total range - all columns
column_split[c(1:5)] <- "control_24h"
column_split[c(6:10)] <- "infected_24h"
column_split <- factor(column_split, levels = c("control_24h", "infected_24h"))

## Exclude 0 rows
norm.irs.z <- norm.irs.z[rowSums(norm.irs.z != 0) > 0, ]

file_name=paste("heatmap-irs_", time_point, ".pdf", sep = "")
  
ht <- Heatmap(norm.irs.z,
                cluster_rows = F,
                cluster_columns = F,
                column_labels = colnames(norm.irs.z),
                name = "Z-score",
                column_split = column_split,
                row_names_gp = grid::gpar(fontsize = 8), 
                column_names_gp = grid::gpar(fontsize = 9),
                row_title = title                
)
ht <- draw(ht)

## create pdf
pdf(file = file_name, height = 8)
ht <- draw(ht)
dev.off()



