#====================================================
# Making a Heatmap Plot to Visualize Ligands in MS CSF 
# Driving Changes in Ependymal Cells
#====================================================

# Outline:
# 1. Loading Packages and Setting Seed
# 2. Loading Ependymal Cells and Whole Object
# 3. Heatmap of Ligands Upregulated in MS CSF Most Correlated with Periventricular Damage
# 4. Creating Heatmap of Ligand-Target Interactions
# 5. Creating Heatmap of Ligands Most Predicted to Drive Target Genes

#====================================================
# 1. Loading Packages and Setting Seed
#====================================================

library(nichenetr)
library(tidyverse)
library(ggplot2)
library(Seurat)
library(data.table)
library(patchwork)

set.seed(1)
setwd("your_wd")

#====================================================
# 2. Loading Matrices 
#====================================================

ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))

#====================================================
# 3. Heatmap of Ligands Upregulated in MS CSF 
# Most Correlated with Periventricular Damage
#====================================================

###################### Fig 5E ######################

# Top 102 genes with log2FC above a threshold
top_upregulated_genes <- c(
  "CHI3L1", "CCL2", "PLCZ1", "ARID5A", "HGF", "RERGL", "ITGBL1", "LUZP2", "PXDN", "LRAT",
  "CD44", "CHI3L2", "THBS2", "NTNG1", "IFI16", "FABP5", "IRAK2", "TNFSF13B", "PRRX1", "F3",
  "GBP2", "ARHGEF3", "IGDCC4", "PLCXD3", "PMP2", "GRIK1", "KCNT2", "FAM20C", "SORCS3",
  "CCDC85A", "THSD7A", "MEIKIN", "MOXD1", "HELB", "FRMPD3", "STK17A", "RGS20", "ESR1",
  "ADRA1A", "SERPINA3", "SRPX", "CORIN", "HRH1", "PDPN", "CNDP1", "TNC", "UBASH3B", "APLNR",
  "IL33", "SPP1", "NIBAN1", "AHR", "SLC2A3", "CADPS", "DDB2", "QDPR", "ADCY8", "TPST1",
  "NTRK3", "LDLR", "OSMR", "PDE1A", "BRINP3", "PRR16", "ACTN2", "CLIC4", "SLC44A3",
  "SLC9A9", "PIRT", "NLN", "ZFYVE28", "LHFPL6", "DCLK1", "ERBIN", "GALNT15", "FAM126A",
  "NFIL3", "ARSB", "ARNTL2", "FAM189A2", "RBMS3", "RUNX1", "PRUNE2", "PAG1", "AQP4",
  "MYO5B", "PLSCR1", "AEBP2", "SNAP91", "TF", "PFKFB3", "PADI2", "SAMD4A", "MAN2A1",
  "SLC1A2", "PDE4B", "NTM", "PGM2L1", "CTTNBP2", "SGCD", "ADGRL3", "EDIL3"
)

# Ligands Upregulated in MS CSF and Correlated with Periventricular Damage
specified_ligands <- c(
  "TNFRSF1A", "TNF", "CCL21", "FGA", "FGB", "FGG", "CCL19", "IL2", "TNFSF13B", 
  "CXCL10", "IFNG", "CHI3L1", "PVALB", "NEFL", "CD163"
)

#====================================================
# 4. Creating Heatmap of Ligand-Target Interactions
#====================================================

# Create data table with tight junction genes and specified ligands
data.dt <- ligand_target_matrix[top_upregulated_genes, intersect(specified_ligands, colnames(ligand_target_matrix))] |> 
  as.data.table(keep.rownames = "target") |> 
  melt(id.vars = "target", variable.name = "ligand", value.name = "weight")

# Create heatmap with target genes ordered by weight
Heatmap_bar <- wrap_plots(
  ggplot(data.dt) +
    aes(
      x = fct_reorder(target, weight, .fun = sum),
      y = fct_reorder(ligand, weight, .fun = sum),
      fill = weight
    ) +
    geom_tile() +
    scale_fill_viridis_c(option = "turbo") +
    scale_x_discrete(position = "top") +
    theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
    labs(y = "Ligand") + 
    theme(legend.position = "left"),
  
  ggplot(data.dt) +
    aes(x = weight, y = fct_reorder(ligand, weight, .fun = sum)) + 
    geom_col(position = position_stack()) +
    theme_void(),
  
  widths = c(1, 0.15)
) & theme(plot.margin = margin(0, 0, 0, 0))

# Save the heatmap plot
ggsave("Heatmap_bar.png", Heatmap_bar, width = 14, height = 4)

#====================================================
# 5. Creating Heatmap of Ligands Most Predicted to Drive Target Genes
#====================================================

###################### Supp Fig 4B ######################

# ligands most predicted to drive the MS ependymal transcriptomic signature 
# found by summing the weighted interaction between all ligands in ligand_target_matrix and target weight


specified_ligands <- c(
  "TGFB1", "IL4", "IL1B", "TNF", "BMP2", "IL10", "IFNG", "NOG", "IL27", "TGFB3",
  "IL21", "IFNL1", "IL15", "PDGFD", "BDNF", "OSM", "TSLP", "EGF", "FGF2", "IL36A",
  "IL13", "IGF2", "IL26", "IL22", "LIF", "BMP6", "CSF2", "ELFN1", "IL24", "IL6",
  "TNFSF11", "IL25", "IGF1", "IFNA1", "EBI3", "PDGFB", "IL36G", "TGFB2", "VEGFA",
  "IL12A", "IL2", "IL11", "IFNB1", "IL23A", "IL36B", "IL20", "PDGFC", "CCL2",
  "HMGB1", "INHBA"
)

# Create data table with tight junction genes and specified ligands
data.dt <- ligand_target_matrix[top_upregulated_genes, intersect(specified_ligands, colnames(ligand_target_matrix))] |> 
  as.data.table(keep.rownames = "target") |> 
  melt(id.vars = "target", variable.name = "ligand", value.name = "weight")

# Create heatmap for predicted ligands
Heatmap_bar <- wrap_plots(
  ggplot(data.dt) +
    aes(
      x = fct_reorder(target, weight, .fun = sum),
      y = fct_reorder(ligand, weight, .fun = sum),
      fill = weight
    ) +
    geom_tile() +
    scale_fill_viridis_c(option = "turbo") +
    scale_x_discrete(position = "top") +
    theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
    labs(y = "Ligand") + 
    theme(legend.position = "left"),
  
  ggplot(data.dt) +
    aes(x = weight, y = fct_reorder(ligand, weight, .fun = sum)) + 
    geom_col(position = position_stack()) +
    theme_void(),
  
  widths = c(1, 0.15)
) & theme(plot.margin = margin(0, 0, 0, 0))

# Save the heatmap plot
ggsave("Heatmap_bar.png", Heatmap_bar, width = 14, height = 4)
