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
library(tibble)

###############################LOAD DATA#################################################
time_point <- ""
sub_dir <- "chloroquine"
workdir <- paste0("/home/pospim/Desktop/work/bioinformatics/datasets/other/GSE209582/", sub_dir, '/', time_point)
gse_part <- sub(".*/(GSE[^_/]+).*", "\\1", workdir)

setwd(workdir)


## Read data
data <- read.table("counts_named.csv", header = TRUE, row.names = 1, comment.char = "", sep=",")
coldata <- read.table("samplefile.csv", header = T, sep=",")
contr <- read.table("contrasts.csv", header = T, sep=",")
rownames(coldata) <- coldata$sample

coldata
control <- "TE11_control"
virus <- "TE11_Chlor"

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

## Change sample names
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

###############################aGPCRs#################################################
res.deseq2 <- res.deseq2[, c(7, 1, 2, 3, 4, 5, 6)]

res.sig <- dplyr::filter(res.deseq2, log2FoldChange < -2 & pvalue < 0.05)
res.sig <- res.sig[!is.na(res.sig$gene_symbol), ]

sig_ids <- rownames(res.sig)

rownames(res.sig) <- res.sig$gene_symbol
res.sig <- rownames_to_column(res.sig, var = "gene_id")

res.sig <- res.sig %>% 
  group_by(gene_symbol) %>%
  summarise(gene_id = first(gene_id), across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")

res.sig <- as.data.frame(res.sig)
rownames(res.sig) <- res.sig$gene_id
res.sig$gene_id <- NULL
res.sig <- res.sig[order(-res.sig$log2FoldChange), ]

res_table <- paste0(gse_part, "-", time_point, "-down.csv")
write.table(res.sig, file = res_table, sep = "\t", quote = FALSE, row.names = FALSE)

###############################HEATMAP ANALYSIS#################################################
dds2 <- estimateSizeFactors(dds)
norm.counts <- counts(dds2, normalized = TRUE)

norm.counts <- as.data.frame(norm.counts)

#write.table(norm.counts, file = "counts_normalized.csv", sep = "\t", quote = FALSE, row.names = FALSE)
norm.sig <- norm.counts[sig_ids, ]
norm.sig <- as.data.frame(norm.sig)
norm.sig$gene_symbol <- mapIds(org.Hs.eg.db,
                                  keys = row.names(norm.sig),
                                  column = "SYMBOL",
                                  keytype = "ENTREZID",
                                  multiVals = "first"
)

norm.sig <- norm.sig[!is.na(norm.sig$gene_symbol), ]
rownames(norm.sig) <- norm.sig$gene_symbol
norm.sig$gene_symbol <- NULL

## Calculate Z-scores
norm.sig.z <- t(apply(norm.sig, 1, scale))
colnames(norm.sig.z) <- colnames(norm.sig)
norm.sig.z <- as.data.frame(norm.sig.z) %>%
  mutate_all(~replace(., is.nan(.), 0))
norm.sig.z <- as.matrix(norm.sig.z)

###############################HEATMAP#################################################
group_labels <- factor(rep(c("control", "chlor+fluoroquine"), each = ncol(norm.sig.z)/2))

column_split <- rep(NA, ncol(norm.sig.z)) # total range - all columns
column_split[c(1:4)] <- "Control"
column_split[c(5:7)] <- "Chloroquine"
column_split <- factor(column_split, levels = c("Control", "Chloroquine"))


num_genes <- nrow(norm.sig.z)
batch <- 50
num_batches <- ceiling(num_genes / batch)

for (i in 1:num_batches) {
  start <- (i-1) * batch + 1
  end <- min(i * batch, num_genes)
  norm.batch <- norm.sig.z[start:end, ]
  ht <- Heatmap(norm.batch, name = "Z-score",
                cluster_rows = TRUE,
                cluster_columns = FALSE,
                column_split = column_split
                )
  pdf(paste0("heatmap_down_genes_", i, ".pdf"), height = 8)
  draw(ht)
  dev.off()
}



