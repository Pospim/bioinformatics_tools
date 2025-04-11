
rm(list = ls(all.names = TRUE))
gc() # free up memory and report the memory usage

library(dplyr)
library(VennDiagram)
library(pheatmap)
library(ggplot2)
library(stringr)

# Configuration ================================================================
in_path <- "/home/pospim/Desktop/work/bioinformatics/GSEA/results/calu3"
setwd(in_path)
out_path <- file.path(in_path, 'trends')

files <- list.files(path = in_path, pattern = "*.csv")
print(files)

# Read all GSEA results into a list of data frames
gsea_results <- lapply(files, read.delim, sep = "\t", header = TRUE)
names(gsea_results) <- gsub(".*-(\\d+h).*", "\\1", basename(files)) # Extract time points (e.g., "24h")

# Data preparation =============================================================

# Extract significant pathways and NES for each time point
significant_pathways <- lapply(names(gsea_results), function(tp) {
  gsea_results[[tp]] %>%
    filter(padj < 0.05) %>%   
    filter(NES > 0) %>%           # UP / DOWN
    select(Pathway = pathway, NES) %>% 
    mutate(TimePoint = tp)
})

# Combine all datasets into a single data frame
pathway_data <- do.call(rbind, significant_pathways)

# Count occurrences of each pathway across all time points
pathway_summary <- pathway_data %>%
  group_by(Pathway) %>%
  summarize(
    Count = n(),                     # Number of time points the pathway is enriched in
    TimePoints = paste(TimePoint, collapse = ", ") # Time points where enriched
  ) %>%
  arrange(desc(Count))

# Analysis =====================================================================
trend_threshold <- floor(length(files) * 0.5)  # Adjust threshold
common_pathways <- pathway_summary %>%
  filter(Count >= trend_threshold)

common_pathway_data <- pathway_data %>%
  filter(Pathway %in% common_pathways$Pathway)

common_pathway_data <- common_pathway_data %>%
  mutate(TimePoint_num = as.numeric(gsub("h", "", TimePoint))) %>%
  mutate(TimePoint = factor(TimePoint, 
                            levels = unique(TimePoint[order(TimePoint_num)]), 
                            ordered = TRUE))

write.csv(common_pathways, file.path(out_path, "common_pathways_up.csv"), row.names = FALSE)

minNES <- min(common_pathway_data$NES, na.rm = TRUE)
maxNES <- max(common_pathway_data$NES, na.rm = TRUE)

padding <- 0.1 * (maxNES - minNES)
ymin <- minNES - padding
ymax <- maxNES + padding

# Visualization ================================================================
pdf(file.path(out_path, "pathway_line_plots2.pdf"), width = 10, height = 6)

# Loop through each unique pathway and create a plot
for (pathway in unique(common_pathway_data$Pathway)) {
  # Filter data for the current pathway
  pathway_data <- common_pathway_data %>% filter(Pathway == pathway)
  wrapped_title <- str_wrap(gsub("_", " ", pathway), width = 60)
  
  min_val <- min(pathway_data$NES, na.rm = TRUE)
  max_val <- max(pathway_data$NES, na.rm = TRUE)
  rng <- if ((max_val - min_val) > 0.7) max_val - min_val else 0.7
  # Generate the line plot
  p <- ggplot(pathway_data, aes(x = TimePoint, y = NES, group = 1)) +
    geom_line(color = "blue", size = 1) +
    geom_point(size = 3, color = "red") +
    labs(title = wrapped_title,
         x = "Time Point",
         y = "NES") +
    theme_minimal(base_size = 15) +
    theme(
      plot.title = element_text(hjust = 0.1, face = "bold"),
      plot.title.position = "plot",
      plot.margin = margin(t=40,),
     ) +
    scale_y_continuous(limits = c(min_val,min_val+rng))
  
  # Print the plot to the PDF
  print(p)
}


# Close the PDF
dev.off()

