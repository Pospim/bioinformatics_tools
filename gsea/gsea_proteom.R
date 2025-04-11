# https://biostatsquid.com/fgsea-tutorial-gsea/
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage

# Loading relevant libraries
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  library(BiocManager)
}
if (!require(fgsea)) {
  BiocManager::install("fgsea")
}
if (!require(tidyverse)) {
  BiocManager::install("tidyverse")
}
if (!require(RColorBrewer)) {
  BiocManager::install("RColorBrewer")
}

library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(fgsea)
library(data.table)

# Function: Adjacency matrix to list -------------------------
matrix_to_list <- function(pws){
  pws.l <- list()
  for (pw in colnames(pws)) {
    pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws.l)
}
## Function: prepare_gmt --------------------------------------
prepare_gmt <- function(gmt_file, genes_in_data, savefile = FALSE){
  # for debug
  #file <- gmt_files[1]
  #genes_in_data <- df$gene_symbol
  
  # Read in gmt file
  gmt <- gmtPathways(gmt_file)
  hidden <- unique(unlist(gmt))
  
  # Convert gmt file to a matrix with the genes as rows and for each go annotation (columns) the values are 0 or 1
  mat <- matrix(NA, dimnames = list(hidden, names(gmt)),
                nrow = length(hidden), ncol = length(gmt))
  for (i in 1:dim(mat)[2]){
    mat[,i] <- as.numeric(hidden %in% gmt[[i]])
  }
  
  #Subset to the genes that are present in our data to avoid bias
  hidden1 <- intersect(genes_in_data, hidden)
  mat <- mat[hidden1, colnames(mat)[which(colSums(mat[hidden1,])>5)]] # filter for gene sets with more than 5 genes annotated
  # And get the list again
  final_list <- matrix_to_list(mat) # for this we use the function we previously defined
  
  if(savefile){
    saveRDS(final_list, file = paste0(gsub('.gmt', '', gmt_file), '_subset_', format(Sys.time(), '%d%m'), '.RData'))
  }
  
  print('Wohoo! .gmt conversion successfull!:)')
  return(final_list)
}
# Function: Clean and prepare data -------------------------
clean_data <- function(df) {
  df <- dplyr::filter(df, grepl("GN=", Group))
  df$gene_symbol <- sub(".*GN=([A-Za-z0-9]+).*", "\\1", df$Group)
  df <- df %>% 
    dplyr::select(-c(Group, Peak.Name)) %>% 
    dplyr::filter(!is.na(gene_symbol)) %>% 
    dplyr::group_by(gene_symbol) %>% 
    dplyr::filter(abs(log2fold) == max(abs(log2fold))) %>% 
    dplyr::ungroup()
  return(df)
}

# Configuration -----------------------------------------------------------
workdir <- "/home/pospim/Desktop/work/bioinformatics/proteomics/GSEA"
in_path <- file.path(workdir, "datasets")
out_path <- file.path(workdir, "results")
bg_path <- file.path(workdir, "bg_genes")
setwd(workdir)

# Analysis ====================================================
## 1. Read in data -----------------------------------------------------------
list.files(in_path)
raw_file <- "GSE213759_24h.csv"
df <- read.csv(file.path(in_path, raw_file), sep = "\t")
df <- as.data.frame(df)
df <- clean_data(df)

## 2. Prepare background genes -----------------------------------------------
my_genes <- df$gene_symbol
list.files(bg_path)
gmt_files <- list.files(path = bg_path, pattern = '.gmt', full.names = TRUE)
gmt_files

# Change [x] for different pathway
comparison_number <- 2
bg_genes <- prepare_gmt(gmt_files[comparison_number], my_genes, savefile = FALSE)

if (comparison_number == 1){
  comparison_name <- "KEGG"
} else if (comparison_number == 2){
  comparison_name <- "Reactome"
} else if (comparison_number == 3){
  comparison_name <- "GO"
} else {
  comparison_name <- "Unknown"
}

## 3. Rank genes -------------------------------------------------------------
rankings <- sign(df$log2fold)*(-log10(df$pval)) # we will use the signed p values from spatial DGE as ranking
names(rankings) <- df$gene_symbol # genes as names
head(rankings)

# Some genes have such low p values that the signed pval is +- inf, we need to change it to the maximum * constant to avoid problems with fgsea
max_ranking <- max(rankings[is.finite(rankings)])
min_ranking <- min(rankings[is.finite(rankings)])
rankings <- replace(rankings, rankings > max_ranking, max_ranking * 10)
rankings <- replace(rankings, rankings < min_ranking, min_ranking * 10)
rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
head(rankings)

ggplot(data.frame(gene_symbol = names(rankings)[1:100], ranks = rankings[1:100]), aes(gene_symbol, ranks)) + 
  geom_point() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

## 4. Run GSEA ---------------------------------------------------------------
GSEAres <- fgsea(pathways = bg_genes, # List of gene sets to check
                 stats = rankings,
                 scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                 minSize = 10,
                 maxSize = 500,
                 nproc = 1) # for parallelisation

## 5. Plot results -----------------------------------------------------------
head(GSEAres[order(padj), ])
sum(GSEAres[, padj < 0.01])
sum(GSEAres[, pval < 0.01])

num_top_upregulated <- 25
num_top_downregulated <- 25

if (comparison_number == 1){
  prefix <- "KEGG_MEDICUS_"
} else if (comparison_number == 2){
  prefix <- "REACTOME_"
} else if (comparison_number == 3){
  prefix <- "GO"
} else {
  prefix <- "Unknown"
}
# Remove the prefix in pathway names in GSEAres
GSEAres[, pathway := gsub(prefix, "", pathway)]

# Sort GSEA results by absolute NES (most changed pathways)
GSEAres <- GSEAres[order(-abs(NES))]

# Save results
results_file <- paste0(out_path,"/", comparison_name, "_GSEA_results.csv")
fwrite(GSEAres, results_file, sep = "\t")

top_up_pathways <- GSEAres[ES > 0][head(order(pval), num_top_upregulated), pathway]
top_down_pathways <- GSEAres[ES < 0][head(order(pval), num_top_downregulated), pathway]
top_pathways <- c(top_up_pathways, rev(top_down_pathways))


plotGseaTable(bg_genes[top_pathways], stats=rankings, fgseaRes = GSEAres, gseaParam = 0,5)

# Function: Plot mostly affected pathways -------------------------
plot_top_pathways <- function(GSEAres, rankings, bg_genes) {
  collapsed_pathways <- collapsePathways(GSEAres[order(pval)][pval < 0.05], bg_genes, rankings, nperm = 1000)
  main_pathways <- GSEAres[pathway %in% collapsed_pathways$mainPathways][order(-NES), pathway]
  plotGseaTable(bg_genes[main_pathways], rankings, GSEAres, gseaParam = 0)
}

plot_top_pathways(GSEAres, rankings, bg_genes)
