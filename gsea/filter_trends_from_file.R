# Load necessary libraries
library(dplyr)

# Paths configuration
in_path <- "/home/pospim/Desktop/work/bioinformatics/GSEA/results/calu3/raw"  # Replace with your input folder path
common_pathways_file <- "/home/pospim/Desktop/work/bioinformatics/GSEA/results/calu3/trends_down/common_pathways_down.csv"  # Replace with the path to your common pathways file
out_path <- '/home/pospim/Desktop/work/bioinformatics/GSEA/results/calu3/filtered_down'

setwd(in_path)
# Create output directory if it doesn't exist
if (!dir.exists(out_path)) {
  dir.create(out_path)
}

# Load common pathways
data_common <- read.csv(common_pathways_file)
common_pathways <- data_common$Pathway

# List all GSEA result files
files <- list.files(path = in_path, pattern = "*.csv", full.names = TRUE)

# Process each file
for (file in files) {
  # Extract the time point from the file name
  base_name <- basename(file)
  time_point <- gsub(".*-(\\d+h).*", "\1", base_name)
  
  # Read GSEA results
  gsea_data <- read.csv(file, sep = "\t", header = TRUE)
  
  # Filter for common pathways
  filtered_data <- gsea_data %>% filter(pathway %in% common_pathways) %>%
    filter(padj<0.05, NES < 0)
  
  # Save filtered results
  out_file <- file.path(out_path, paste0(base_name, 'filtered', ".csv"))
  write.csv(filtered_data, out_file, row.names = FALSE)
  
  cat("Filtered data for", time_point, "saved to", out_file, "\n")
}

cat("All files processed. Filtered results are saved in", out_path, "\n")

