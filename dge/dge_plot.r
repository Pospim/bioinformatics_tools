## Load packages
if (! require(ggplot2)) {
  install.packages("ggplot2")
}
if (! require(dplyr)) {
  install.packages("dplyr")
}
if (! require(tidyr)) {
  install.packages("tidyr")
}
library(ggplot2)
library(dplyr)
library(tidyr)

# Function to append *, **, *** to gene symbols based on p-value thresholds
annotate_genes <- function(gene_symbol, pvalues) {
  min_pvalue <- min(pvalues, na.rm = TRUE)
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
gene_order <- c(
  "ADGRA1", "ADGRA2", "ADGRA3",
  "ADGRB1", "ADGRB2", "ADGRB3",
  "CELSR1", "CELSR2", "CELSR3",
  "ADGRD1", "ADGRD2",
  "ADGRE1", "ADGRE2", "ADGRE3", "ADGRE5",
  "ADGRF1", "ADGRF2", "ADGRF3", "ADGRF4", "ADGRF5",
  "ADGRG1", "ADGRG2", "ADGRG3", "ADGRG4", "ADGRG5", "ADGRG6", "ADGRG7",
  "ADGRL1", "ADGRL2", "ADGRL3", "ADGRL4",
  "ADGRV1"
)
gene_order <- rev(gene_order)

## Change variables
workdir <- "/home/pospim/Desktop/work/bioinformatics/datasets/SARS/sandra_calu3/cleaned/"
gse_part <- sub(".*/(GSE[^_/]+).*", "\\1", workdir)

setwd(workdir)
file_name <- paste(gse_part, "6h-analysis_results.tsv", sep="-")
plot_name <- "DGE_plot.png"

# Define sample names dynamically
sample_names <- c( "Calu3_Sars_6h","Calu3_Sars_12h", "Calu3_Sars_24h")

data <- read.table(file_name, header = TRUE, sep = "\t", stringsAsFactors=FALSE)
data[] <- lapply(data, function(x) gsub(",",".",x))

## Clean data
for (i in seq_along(sample_names)) {
  log2fc_col <- paste0("Log2FoldChange_", i)
  pvalue_col <- paste0("Pvalue_", i)
  
  # Convert to numeric
  data[[log2fc_col]] <- as.numeric(as.character(data[[log2fc_col]]))
  data[[pvalue_col]] <- as.numeric(as.character(data[[pvalue_col]]))
  
  # Handle NAs
  data[[log2fc_col]][is.na(data[[log2fc_col]])] <- 0
  data[[pvalue_col]][is.na(data[[pvalue_col]])] <- 1
  
  # Calculate log10 p-value
  data[[paste0("log10_pvalue_", i)]] <- -log10(data[[pvalue_col]])
}

data$gene_symbol <- factor(data$gene_symbol, levels = gene_order)

# Extract p-value columns
pvalue_cols <- grep("Pvalue_", colnames(data), value = TRUE)
pvalues <- data[, pvalue_cols]
pvalues <- as.data.frame(pvalues)
rownames(pvalues) <- data$gene_symbol

# Annotate gene symbols and set row names
annotated_gene_symbols <- mapply(annotate_genes, rownames(pvalues), split(pvalues, seq(nrow(pvalues))))
data$gene_symbol <- annotated_gene_symbols

# Convert data to long format for plotting
log2fc_cols <- grep("Log2FoldChange_", colnames(data), value = TRUE)
log10_pvalue_cols <- grep("log10_pvalue_", colnames(data), value = TRUE)

# Match Log2FoldChange and Pvalue columns properly for each sample
data_lng <- data %>%
  pivot_longer(cols = c(all_of(log2fc_cols), all_of(log10_pvalue_cols)),
               names_to = c(".value", "sample_index"),
               names_pattern = "(Log2FoldChange|log10_pvalue)_(\\d+)") %>%
  mutate(sample = sample_names[as.numeric(sample_index)]) %>%
  dplyr::select(gene_symbol, sample, Log2FoldChange, log10_pvalue)

# Ensure the samples are ordered based on the sample_names vector
data_lng$sample <- factor(data_lng$sample, levels = sample_names)
data_lng$gene_symbol <- factor(data_lng$gene_symbol, levels = rev(annotated_gene_symbols))

## Plot DGE
colors <- c("blue","cadetblue3","cadetblue1", "lightblue", "lightblue4", "azure3")  # Add more colors for additional samples

p <- ggplot(data_lng, aes(x = Log2FoldChange, y = gene_symbol, fill = sample)) +
  geom_col(position = "dodge") +
  scale_fill_manual(name="Sample", values=setNames(colors[1:length(sample_names)], sample_names)) +
  scale_alpha_continuous(range = c(0.3, 1), guide = 'none') +
  labs(title="Differential Gene Expression",
       x = "relative expression (log2FoldChange)", y = "") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)) +
  theme(axis.text.y = element_text(size=10, margin=margin(t=100, r=0,b=20,l=0))) 

ggsave(plot_name, plot=p, width=10, height=10, bg="white", dpi=300)
print(p)

