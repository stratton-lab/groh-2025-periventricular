#====================================================
# Outline of the Script - Plotting ATAC-seq Data
#====================================================
# 1. Loading Packages
#    - Load necessary libraries for visualization and motif analysis.
#
# 2. QC and Feature Distribution Plots
#    - Generate scatter and violin plots for quality control metrics.
#
# 3. Footprinting Plots
#    - Visualize footprinting for STAT3, BACH2, STAT1, RFX2, and RFX3.
# 
# 4. Motif Visualization
#    - Generate motif sequence plots for key transcription factors.
# 
# 5. Chromatin Accessibility Analysis
#    - Compute and visualize motif accessibility differences between MS and control.
#

#====================================================
# Loading Packages
#====================================================

# Load libraries (Signac should be loaded before Seurat)
library(Signac)
library(Seurat)
library(ggplot2)
library(dplyr)
library(EnsDb.Hsapiens.v86)
library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
library(chromVAR)

# Set working directory
setwd("path_to_your_wd")
set.seed(1)

#====================================================
# QC and Feature Distribution Plots
#====================================================

# Scatter plots for QC metrics
a1 <- DensityScatter(ATAC_ependymal, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
a2 <- DensityScatter(ATAC_ependymal, x = 'nucleosome_signal', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

# Combine scatter plots
a1 | a2

# Violin plots for feature distributions
VlnPlot(object = ATAC_ependymal, 
        features = c('nCount_ATAC', 'nFeature_ATAC', 'TSS.enrichment'),
        pt.size = 0.1,
        ncol = 6)

# UMAP plot with cell labels
DimPlot(object = ATAC_ependymal, label = TRUE) + NoLegend()

#====================================================
# Footprinting Plots
#====================================================

###################### Supp Fig 3A ######################


# Individual TF footprint plots
footprint_tfs <- c("BACH2", "STAT3", "STAT1", "RFX2", "RFX3")

for (tf in footprint_tfs) {
  p2 <- PlotFootprint(ATAC_ependymal, features = tf, group.by = "disease")
  ggsave(paste0("footprinting_", tf, ".png"), plot = p2, width = 8, height = 5, dpi = 300)
}

#====================================================
# Motif Visualization
#====================================================

###################### Supp Fig 3A ######################

# Generate and save motif plots for key TFs
for (tf in footprint_tfs) {
  motif_plot <- MotifPlot(object = ATAC_ependymal, motifs = tf)
  ggsave(filename = paste0(tf, "_Motif.png"), plot = motif_plot, width = 8, height = 3)
}


#====================================================
# Chromatin Accessibility Analysis
#====================================================

###################### Supp Fig 3B ######################


library(ggplot2)
library(reshape2)

# Subset the ATAC object by disease status
ms_cells <- WhichCells(ATAC_ependymal, expression = disease == "MS")
ctrl_cells <- WhichCells(ATAC_ependymal, expression = disease == "CTRL")

# Extract chromVAR data for MS and CTRL
chromvar_ms <- chromvar_data[, ms_cells]
chromvar_ctrl <- chromvar_data[, ctrl_cells]

# Transpose and convert to data frame
chromvar_ms_df <- as.data.frame(t(chromvar_ms))
chromvar_ctrl_df <- as.data.frame(t(chromvar_ctrl))

# Add metadata (disease type)
chromvar_ms_df$disease <- "MS"
chromvar_ctrl_df$disease <- "CTRL"

# Combine data
chromvar_combined <- rbind(chromvar_ms_df, chromvar_ctrl_df)

# Define motifs of interest
selected_motifs <- c("MA1101.2" = "BACH2", "MA0137.3" = "STAT1", "MA0144.2" = "STAT3", 
                     "MA0600.2" = "RFX2", "MA0798.2" = "RFX3") 

# Extract selected motifs and melt data for plotting
chromvar_selected <- chromvar_combined[, c(names(selected_motifs), "disease")]
chromvar_melted <- melt(chromvar_selected, id.vars = "disease", variable.name = "TF", value.name = "Accessibility")

# Violin plots for accessibility differences
for (tf_id in names(selected_motifs)) {
  tf_name <- selected_motifs[tf_id]
  
  # Subset data
  subset_data <- subset(chromvar_melted, TF == tf_id)
  
  # Create violin plot
  plot <- ggplot(subset_data, aes(x = disease, y = Accessibility, fill = disease)) +
    geom_violin(trim = FALSE, scale = "width", alpha = 0.7) + 
    stat_summary(fun = mean, geom = "crossbar", width = 0.2, color = "black", fatten = 2) +
    labs(title = paste("Motif Accessibility for", tf_name), x = "Disease", y = "Accessibility Score") +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 14)
    ) +
    scale_fill_manual(values = c("MS" = "blue", "CTRL" = "red"))
  
  # Save plot
  ggsave(filename = paste0("ViolinPlot_", tf_name, ".png"), plot = plot, width = 4, height = 5)
}
