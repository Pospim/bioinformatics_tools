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
sub_dir <- "chlor_fluor"
workdir <- paste0("/home/pospim/Desktop/work/bioinformatics/datasets/other/GSE209582/", sub_dir, '/', time_point)
gse_part <- sub(".*/(GSE[^_/]+).*", "\\1", workdir)

setwd(workdir)

## Read data
data <- read.table("counts_named.csv", header = TRUE, row.names = 1, comment.char = "", sep=",")
coldata <- read.table("samplefile.csv", header = T, sep=",")
contr <- read.table("contrasts.csv", header = T, sep=",")
rownames(coldata) <- coldata$sample

coldata$condition

control <- "TE11_control"
virus <- "TE11_Chlor_Fluor"
###############################PREPROCESS#################################################
# Remove version numbers 
#rowNames <- rownames(data)
#rowNames <- sub("\\..*", "", rowNames)
#rowNames <- make.unique(rowNames)
#rownames(data) <- rowNames
###############################VERIFY#################################################
all(rownames(coldata) == colnames(data))
data <- round(data)

# Check if any sum differs by more than 50% from the max sum
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
                                 keytype = "ENTREZID",
                                 multiVals = "first"
)
###############################VIRAL COMBINED ANNOTATION#################################################
#Geneid <- rownames(data)
#non_ensg_rows <- subset(data, !grepl("ENSG", Geneid))
#rownames(non_ensg_rows)

#res.deseq2$gene_symbol <- ifelse(
#  row.names(res.deseq2) %in% names(sarscov2_gene_map),
#  sarscov2_gene_map[row.names(res.deseq2)],
#  res.deseq2$gene_symbol
#)
#mapped_rows <- res.deseq2[row.names(res.deseq2) %in% names(sarscov2_gene_map), ]
#mapped_rows

###############################aGPCRs#################################################
res.deseq2 <- res.deseq2[, c(7, 1, 2, 3, 4, 5, 6)]

gpcrs = c(
  "ADGRA1",
  "ADGRA2",
  "ADGRA3",
  "ADGRB1",
  "ADGRB2",
  "ADGRB3",
  "ADGRD1",
  "ADGRD2",
  "ADGRE1",
  "ADGRE2",
  "ADGRE3",
  "ADGRE5",
  "ADGRF1",
  "ADGRF2",
  "ADGRF3",
  "ADGRF4",
  "ADGRF5",
  "ADGRG1",
  "ADGRG2",
  "ADGRG3",
  "ADGRG4",
  "ADGRG5",
  "ADGRG6",
  "ADGRG7",
  "ADGRL1",
  "ADGRL2",
  "ADGRL3",
  "ADGRL4",
  "ADGRV1",
  "CELSR1",
  "CELSR2",
  "CELSR3"
  )

res.gpcrs <- filter(
  res.deseq2,
  gene_symbol %in% gpcrs
)
res.gpcrs <- res.gpcrs[order(res.gpcrs$gene_symbol), ]

res_table <- paste0(gse_part, "-", time_point, "-analysis_results.csv")
write.table(res.gpcrs, file = res_table, sep = "\t", quote = FALSE, row.names = TRUE)

###############################SAVE FOR GSEA#################################################
gsea_input <- res.deseq2[!is.na(res.deseq2$log2FoldChange), ]
gsea_input <- gsea_input[!is.na(gsea_input$pvalue), ]
colnames(gsea_input)[colnames(gsea_input) == "pvalue"] ="pval"
colnames(gsea_input)[colnames(gsea_input) == "log2FoldChange"] ="log2fold"

data_sets_dir <- "/home/pospim/Desktop/work/bioinformatics/GSEA/datasets/"
write.table(gsea_input, file = paste0(data_sets_dir, gse_part, "_", sub_dir,'_', time_point, ".csv"), sep = "\t", quote = FALSE, row.names = FALSE)

###############################HEATMAP ANALYSIS#################################################
dds2 <- estimateSizeFactors(dds)
norm.counts <- counts(dds2, normalized = TRUE)

raw_counts <- data[row.names(res.gpcrs), ]
rownames(raw_counts) <- res.gpcrs$gene_symbol
raw_counts$gene_symbol <- rownames(raw_counts)
raw_counts <- raw_counts[, c(ncol(raw_counts), 1:(ncol(raw_counts)-1))]
write.table(raw_counts, file = "raw_gpcrs.csv", sep = '\t', quote = FALSE, row.names = F)

norm.counts <- norm.counts[row.names(res.deseq2), ]
rownames(norm.counts) <- res.deseq2$gene_symbol
norm.counts <- as.data.frame(norm.counts)

norm.counts$gene_symbol <- res.deseq2$gene_symbol

norm.gpcrs <- norm.counts[res.gpcrs$gene_symbol, ]
norm.gpcrs <- norm.gpcrs[, c(ncol(norm.gpcrs), 1:(ncol(norm.gpcrs)-1))]
write.table(norm.gpcrs, file = "norm_gpcrs.csv", sep = '\t', quote = FALSE, row.names = F)

## Calculate Z-scores
norm.gpcrs.z <- t(apply(norm.gpcrs, 1, scale))
colnames(norm.gpcrs.z) <- colnames(norm.gpcrs)

norm.gpcrs.z
norm.gpcrs.z <- as.data.frame(norm.gpcrs.z) %>%
  mutate_all(~replace(., is.nan(.), 0))

norm.gpcrs.z <- as.matrix(norm.gpcrs.z)

#sort norm.gpcrs.z by gene_symbol
norm.gpcrs.z <- norm.gpcrs.z[order(rownames(norm.gpcrs.z)), ]
norm.gpcrs.z

###############################HEATMAPS#################################################
## Create heatmap without ordering
ht0 <- Heatmap(norm.gpcrs.z,
               name = "Z-score", km = 4, cluster_rows = T,
               column_title = "conditions",
               row_title = "adhesion GPCRs",
               cluster_columns = F
)
ht0 <- draw(ht0)

## Create pdf
pdf("heatmap-not_organized.pdf", height = 8)
ht0 <- draw(ht0) 
dev.off()

## Define column splitting variables
#column_split <- factor(rep(c("ctrl_24h", "infected_24h"), each = ncol(norm.sig.z)/2))

column_split <- rep(NA, ncol(norm.gpcrs.z)) # total range - all columns
column_split[c(1:4)] <- "Control"
column_split[c(5:7)] <- "Chlor+Fluoroquine"
column_split <- factor(column_split, levels = c("Control", "Chlor+Fluoroquine"))

## Exclude 0 rows
norm.gpcrs.z <- norm.gpcrs.z[rowSums(norm.gpcrs.z != 0) > 0, ]

########### viral genes ##########
row_split <- ifelse(row.names(norm.gpcrs.z) %in% unname(sarscov2_gene_map),
                    "SARS", 
                    "Calu3")
row_split <- factor(row_split, levels = c("SARS", "Calu3"))
norm.gpcrs.z <- norm.gpcrs.z[order(row_split), ]
##################################

## Create heatmap with n clusters
for (cluster in c(2, 3, 4)) {
  file_name= paste("heatmap-organized_", cluster, ".pdf", sep="")
  
  ht <- Heatmap(norm.gpcrs.z,
                km = cluster, 
                cluster_rows = T,
                cluster_columns = F,
                column_labels = colnames(norm.gpcrs.z),
                name = "Z-score",
                column_split = column_split,
                #row_split = row_split,
                row_names_gp = grid::gpar(fontsize = 8), 
                column_names_gp = grid::gpar(fontsize = 9),
                
                
  )
  ht <- draw(ht)
  ## create pdf
  pdf(file = file_name, height = 8)
  ht <- draw(ht)
  dev.off()
}


