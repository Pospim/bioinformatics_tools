if (!require(ggplot2)) {
  BiocManager::install("ggplot2")
}
if (!require(tidyr)) {
  BiocManager::install("tidyr")
}
if (!require(dplyr)) {
  BiocManager::install("dplyr")
}
if (!require(dplyr)) {
  BiocManager::install("scales")
}
# Load libraries
library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)
library(fmsb)
library(RColorBrewer)

# Set working directory
workdir <- "/home/pospim/Desktop/work/bioinformatics/datasets/SARS/DGE/calu3/24+h/"
setwd(workdir)

file_name <- "deseq-Calu3.csv"
basename <- sub("\\.csv$", "", file_name)
sample_file <- sub("deseq", "samples", file_name)

# Read data
data <- read.table(file_name, header = TRUE, sep = "\t")
rownames(data) <- data$gene_symbol
data <- subset(data, select = -c(gene_symbol))
data <- data[!apply(data, 1, function(row) all(is.na(row))),]

metadata <- read.table(sample_file, header = TRUE, sep = "\t")$Sample
metadata <- make.names(metadata, unique=TRUE)

pvalues <- data[, grepl("Pvalue_", colnames(data))]

# Append *, **, *** to gene symbols based on p-value thresholds
annotate_genes <- function(pvalue_row, gene_symbol) {
  min_pvalue <- min(pvalue_row, na.rm = TRUE)
  if (min_pvalue <= 0.001) {
    return(paste0(gene_symbol, "***"))
  } else if (min_pvalue <= 0.01) {
    return(paste0(gene_symbol, "**"))
  } else if (min_pvalue <= 0.05) {
    return(paste0(gene_symbol, "*"))
  } else {
    return(gene_symbol)
  }
}
rownames(data) <- mapply(annotate_genes, split(pvalues, rownames(pvalues)), rownames(data))
rownames(pvalues) <- mapply(annotate_genes, split(pvalues, rownames(pvalues)), rownames(pvalues))

# Extract Log2FoldChange
log2fold <- data[, grepl("Log2FoldChange", colnames(data))]
log2fold[is.na(log2fold)] <- 0

# Rename the columns 
for (i in 1:length(metadata)) {
  old_name <- paste0("Log2FoldChange_", i)  # Create the old column name
  new_name <- metadata[i]               # Get the new name from sample_names vector
  colnames(log2fold)[colnames(log2fold) == old_name] <- new_name  # Replace the old name with the new name
}

# ------------------------------------------------
# STEP 1: P-VALUE FILTER
# ------------------------------------------------
mask <- !is.na(pvalues) & pvalues <= 0.05
sig_log2fold <- log2fold
sig_log2fold[!mask] <- 0

# remove rows with all 0
sig_log2fold <- sig_log2fold[!apply(sig_log2fold, 1, function(row) all(row == 0)),]
sig_log2fold <- sig_log2fold[, !apply(sig_log2fold, 2, function(col) all(col == 0))]

# ------------------------------------------------
# STEP 2: LOG2FC > 2 FILTER
# ------------------------------------------------
mask2 <- sig_log2fold > 2 | sig_log2fold < -2
sig_log2fold[!mask2] <- 0

# remove rows with all 0
sig_log2fold <- sig_log2fold[!apply(sig_log2fold, 1, function(row) all(row == 0)),]
sig_log2fold <- sig_log2fold[, !apply(sig_log2fold, 2, function(col) all(col == 0))]

# ------------------------------------------------
# STEP 3: PREPARE & PLOT THE RADAR CHART
# ------------------------------------------------
# Transpose
confident <- t(sig_log2fold) %>%
  as.data.frame()

# Vector of colors
colors_set <- c(
  "cornflowerblue", "skyblue", "cyan", "aquamarine", "deepskyblue",
  "lightblue1", "slateblue1", "darkorchid", "mediumpurple1", "pink",
  "orchid", "darkgoldenrod", "burlywood", "bisque3", "brown3",
  "lightcoral", "lightseagreen", "khaki", "lavenderblush", "palegreen",
  "lightgoldenrod", "lightsteelblue", "plum", "rosybrown", "thistle",
  "palevioletred", "peachpuff", "lightpink", "lightyellow", "honeydew"
)

# Prepare for plotting
min_vals_confident <- rep(min(confident), ncol(confident))
max_vals_confident <- rep(max(confident), ncol(confident))
#baseline_confident <- rep(0, ncol(confident))

axis_labels_confident <- seq(min(confident),  max(confident), length.out = 5)
#radar_with_baseline <- rbind(max_vals_confident, min_vals_confident, baseline_confident, confident)
radar_data <- rbind(max_vals_confident, min_vals_confident, confident)

colors_confident <- rainbow(nrow(confident))

# Plot radar chart
plot_name <- sub("deseq", "DGE ", basename)
plot_name <- paste0(plot_name, "_radar.png")

png(plot_name, width = 1000, height = 800, bg="white")

radarchart(
  radar_data,
  axistype = 1, 
  pcol = colors_confident,, 
  plwd = 2, 
  plty = 1,
  cglcol = "gray", 
  cglty = 1, 
  cglwd = 0.8,
  axislabcol = "black",
  #title = sub("deseq-", "DGE: ", "Calu3"),
  title = sub("deseq-", "DGE: ", "Calu3 > 2"),
  caxislabels = axis_labels_confident
)
# Add legend with associated colors
legend("topleft",
       legend = colnames(sig_log2fold),       
       bty = "n",               
       pch = 20,                
       col = colors_confident,           
       cex = 0.8,
       pt.cex = 3,
       inset(4,-1)
)
dev.off()

