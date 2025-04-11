# Load necessary libraries
library(readr)
library(openxlsx)

# Set the directory containing your .csv files
csv_dir <- "/home/pospim/Desktop/work/bioinformatics/GSEA/results/calu3/filtered_down"
setwd(csv_dir)
output_excel <- "combined_data.xlsx"

# Get the list of .csv files
csv_files <- list.files(path = csv_dir, pattern = "\\.csv$", full.names = TRUE)

# Create a new Excel workbook
workbook <- createWorkbook()

# Loop through each .csv file
for (file_path in csv_files) {
  # Extract the file name without extension
  file_name <- tools::file_path_sans_ext(basename(file_path))
  sheet_name <- sub("_.*$", "", file_name)
  # Read the .csv file
  data <- read_delim(file_path, delim=',')
  
  # Add a worksheet to the workbook with the file name as the sheet name
  addWorksheet(workbook, sheet_name)
  
  # Write the data to the worksheet
  writeData(workbook, sheet_name, data)
}

#Save the Excel file
saveWorkbook(workbook, output_excel, overwrite = TRUE)

cat("Excel file created at:", output_excel, "\n")
