## LOAD PACKAGES ##
# library(dplyr)
library(MSstatsConvert)
# library(MSstats)
library(ggplot2)
library(BiocManager)
library(MSstatsTMT)

rm(list=ls())
#-------------------------------------------------------------------------------
## LOAD DATA ##
workdir <- "/home/pospim/Desktop/work/bioinformatics/proteomics/jitka/PXD030395/MQ_search_results_TXT/"
setwd(workdir)

#file_name <- "20210302_115158_InfectionUK_ab_v2_nonorm_Report_PRIDE.csv"
## CONVERT TO MSSTATS ##
#-------------------------------------------------------------------------------
evidence <- read.table("evidence.csv", sep = "\t", header = TRUE)
protein_groups <- read.table("proteinGroups.csv", sep = "\t", header = TRUE)
annotation <- read.csv("annotation.csv", sep = ",")

quant <- MSstatsConvert::MaxQtoMSstatsFormat(evidence = evidence, annotation = annotation,
                            proteinGroups = protein_groups,
                            removeProtein_with1Peptide = TRUE,
                            useUniquePeptide = TRUE,
                            removeFewMeasurements = TRUE,
                            censoredInt = "NA",
                            summaryMethod = "TMP"
                            )

#-------------------------------------------------------------------------------
#raw_data<-read.csv(file_name, header = TRUE, sep = ",")
#head(raw_data)
#str(raw_data)

# Add the PeptideModifiedSequence column (use PeptideSequence as a placeholder)
#raw_data$PeptideModifiedSequence <- raw_data$PeptideSequence
#raw_data$IsotopeLabelType <- "L"

## PREPROCESSING ##
#length(unique(raw_data$ProteinName)) # Check protein names, runs etc.
#unique(raw_data[, c('Run', 'Condition', 'BioReplicate')])

processed_quantile <- dataProcess(raw = quant,
                               logTrans = 2,
                               normalization = 'equalizeMedians',
                               summaryMethod = 'TMP',
                               MBimpute = TRUE,
                               censoredInt = 'NA',
                               maxQuantileforCensored = 0.999,
                               featureSubset = "topN",
                               n_top_feature = 3
                              )

head(processed_quantile$FeatureLevelData)
write.csv(processed_quantile$ProteinLevelData, "protein_level.csv", row.names=FALSE)

# Verify the data
dataProcessPlots(data = processed_quantile, type = "QCPlot",
                 ylimUp=35, width=5, height = 5, originalPlot = FALSE)
dataProcessPlots(data = processed_quantile, type="Profileplot", ylimUp=35,
                 featureName="NA", width=5, height=5, originalPlot = FALSE)
#-------------------------------------------------------------------------------
## STATISTICAL ANALYSIS ##
levels(processed_quantile$ProteinLevelData$GROUP)

group_levels <- c("MERS_12h", "MERS_24h", "MERS_4h", 
                  "Mock_12h", "Mock_24h", "Mock_4h")
# Define the comparison matrix
comparison <- matrix(
  c(0,  0,  1,  0,  0,  -1,    # MERS_4h vs Mock_4h
    1,  0,  0,  -1,  0, 0,     # MERS_12h vs Mock_12h
    0,  1,  0,  0,  -1,  0),   # MERS_24h vs Mock_24h
  nrow = 3,
  byrow = TRUE
)
# Set row and column names for clarity
rownames(comparison) <- c(
  "MERS_4h-Mock_4h",
  "MERS_12h-Mock_12h",
  "MERS_24h-Mock_24h"
)
colnames(comparison) <- group_levels
comparison
# Perform group comparison
comparison_result <- groupComparison(contrast.matrix = comparison, data = processed_quantile)
#-------------------------------------------------------------------------------
## DGE ##
dge_result <- comparison_result$ComparisonResult

dge_result <- dge_result[is.na(dge_result$issue),]
dge_result <- dge_result[!is.na(dge_result$pvalue),]
dge_result <- dge_result[dge_result$pvalue < 0.05,]

#sort on pvalue
dge_result <- dge_result[order(dge_result$pvalue),]
write.csv(dge_result, "dge_result.csv", row.names=FALSE)

top_proteins <- dge_result[order(dge_result$log2FC),][1:20,]
# Plot the results
msstats_data <- list(ComparisonResult = top_proteins)
groupComparisonPlots(
  data = msstats_data$ComparisonResult,
  type = "Heatmap"
)

