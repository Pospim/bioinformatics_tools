# Load necessary libraries
library(dplyr)
library(ggplot2)
library(stringr)

# Reset workspace
rm(list = ls(all.names = TRUE))
gc()

# Configuration ================================================================
workdir <- "/home/pospim/Desktop/work/bioinformatics/GSEA/results/calu3"
in_path <- file.path(workdir, 'raw')
setwd(workdir)
out_path <- file.path(workdir, 'trends_down')

# Create output directory if it doesn't exist
dir.create(out_path, showWarnings = FALSE)

files <- list.files(path = in_path, pattern = "*.csv")
print(files)

# Read all GSEA results into a single data frame ================================
pathway_data <- list()

for (file in files) {
  # Extract time point from the file name
  time_point <- gsub(".*-(\\d+h).*", "\\1", file)
  
  # Read file and filter significant pathways
  data <- read.delim(file.path(in_path, file), sep = "\t", header = TRUE) %>%
    filter(padj < 0.05, NES < 0) %>% # Filter significant and positive NES
    mutate(
      TimePoint = time_point,
      SourceFile = gsub("^(GSE\\d+).*", "\\1", file) # Extract GSE[number]
    ) %>%
    select(Pathway = pathway, NES, TimePoint, SourceFile) # Select relevant columns
  
  # Append to list
  pathway_data[[file]] <- data
}

# Combine all data frames into a single data frame
pathway_data <- bind_rows(pathway_data)

# Count occurrences of each pathway across all time points
pathway_summary <- pathway_data %>%
  group_by(Pathway) %>%
  summarize(
    Count = n(),                     # Number of time points the pathway is enriched in
    TimePoints = paste(unique(TimePoint), collapse = ", ") # Time points where enriched
  ) %>%
  arrange(desc(Count))

# Analysis =====================================================================
trend_threshold <- floor(length(files) * 0.5)  # Adjust threshold

common_pathways <- pathway_summary %>%
  filter(Count >= trend_threshold) %>%
  pull(Pathway) # Extract only pathway names

# Filter pathway_data for common pathways
filtered_pathway_data <- pathway_data %>%
  filter(Pathway %in% common_pathways)

# Ensure TimePoint is ordered correctly
filtered_pathway_data <- filtered_pathway_data %>%
  mutate(TimePoint_num = as.numeric(gsub("h", "", TimePoint))) %>%
  mutate(TimePoint = factor(TimePoint, 
                            levels = unique(TimePoint[order(TimePoint_num)]), 
                            ordered = TRUE))
#save common pathways
out_file <- file.path(out_path, 'common_pathways_down.csv')
write.csv(common_pathways, out_file, row.names = FALSE)
# Visualization ================================================================
pdf(file.path(out_path, "pathway_line_plots_all_points.pdf"), width = 10, height = 6)

# Loop through each unique pathway and create a plot
for (pathway in unique(filtered_pathway_data$Pathway)) {
  # Filter data for the current pathway
  pathway_subset <- filtered_pathway_data %>% filter(Pathway == pathway)
  wrapped_title <- str_wrap(gsub("_", " ", pathway), width = 60)
  
  # Mean for lineplot
  mean_data <- pathway_subset %>%
    group_by(TimePoint) %>%
    summarize(mean_NES = mean(NES, na.rm = TRUE))
  
  # Generate the plot
  p <- ggplot(pathway_subset, aes(x = TimePoint, y = NES)) +
    geom_point(aes(color = SourceFile), size = 3) + # Individual points
    geom_line(data = mean_data, aes(x = TimePoint, y = mean_NES, group = 1), 
              color = "blue", size = 1, linetype = "dashed") + # Mean line
    labs(
      title = wrapped_title,
      x = "Time Point",
      y = "NES",
      color = "Source File"
    ) +
    theme_minimal(base_size = 15) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.margin = margin(t = 40)
    )
  
  # Print the plot to the PDF
  print(p)
}

# Close the PDF
dev.off()

