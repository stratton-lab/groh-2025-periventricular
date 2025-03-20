#====================================================
# Outline - Ependymal Cell Plotting
#====================================================
# 1. Loading Packages
# 2. Visualizing RePACT Results (3D PCA, Violin, Density)
# 3. Plotting Heatmaps & Treeplots
# 4. UMAP Plots for TF Activity-Based Subsets

#====================================================
# 1. Loading Packages
#====================================================
library(Seurat)
library(ggplot2)
library(dplyr)
library(enrichplot)
library(clusterProfiler)
library(plot3D)    
library(viridis)
library(pheatmap)

# Adjust as needed for your environment
setwd("your_wd")

#====================================================
# 2. Visualizing RePACT Results (3D PCA, Violin, Density)
#====================================================
# Load RePACT object from 'analysis.R'
T2D.scRNA.RePACT <- readRDS("E:/multiOmics/integrated object/snRNAseq/Subset_ependymal/out.rds")

#--------------------- 3D PCA Plot ---------------------
###################### Fig 3D ######################

Three_dim_PCA_plot <- scatter3D(
  T2D.scRNA.RePACT$BetaPCA[,"PC_1"],
  T2D.scRNA.RePACT$BetaPCA[,"PC_2"],
  T2D.scRNA.RePACT$BetaPCA[,"PC_3"],
  ticktype = "detailed",
  pch = 20,
  theta = -70,
  phi = 180,
  colvar = ifelse(
    T2D.scRNA.RePACT$BetaPCA[,"diseaseStat"] ==
      unique(T2D.scRNA.RePACT$BetaPCA[,"diseaseStat"])[1], 1, 0
  ),
  bty = "b2",
  cex = 0.6,
  col = alpha.col(col = c("blue", "brown"), 0.6)
)

ggsave(
  filename = "scatter3D_plot_custom_angle.png",
  plot     = Three_dim_PCA_plot,
  width    = 8,
  height   = 6,
  units    = "in",
  dpi      = 300
)

#--------------------- Violin Plot ---------------------
Violin_pseudo_index_per_pateint <- ggplot(
  T2D.scRNA.RePACT$BetaPCA %>% .[complete.cases(.),]
) +
  aes(Sample, pseudo.index.balanced, fill = diseaseStat) +
  geom_violin() +
  geom_boxplot(
    width = 0.2,
    outlier.shape = NA,
    notch = FALSE,
    coef = 0,
    fill = "grey25",
    color = "grey75"
  ) +
  coord_flip() +
  theme_classic() +
  scale_fill_manual(values = c("blue", "brown")) +
  theme_bw() +
  theme(legend.position = "none")

Violin_pseudo_index_per_pateint_bigger_text <- Violin_pseudo_index_per_pateint +
  theme(
    axis.text.x   = element_text(size = 15),
    axis.text.y   = element_text(size = 15),
    axis.title.x  = element_text(size = 15),
    axis.title.y  = element_text(size = 15),
    legend.text   = element_text(size = 12),
    legend.title  = element_text(size = 14)
  )

ggsave(
  filename = "Violin_pseudo_index_per_pateint_bigger_text.png",
  plot     = Violin_pseudo_index_per_pateint_bigger_text,
  width    = 8,
  height   = 6,
  units    = "in",
  dpi      = 300
)

#--------------------- Density Plot ---------------------
Density_plot <- ggplot(T2D.scRNA.RePACT$BetaPCA) +
  aes(pseudo.index.balanced, fill = diseaseStat) +
  geom_density(alpha = 0.75) +
  scale_fill_manual(values = c("blue", "brown")) +
  theme_classic() +
  theme(legend.position = "none")

ggsave(
  "density_plot_pseudo_index.png",
  plot   = Density_plot,
  width  = 8,
  height = 6,
  dpi    = 300
)

#====================================================
# 3. Plotting Heatmaps & Treeplots
#====================================================
#--------------------- RePACT Heatmap ---------------------
###################### Fig 3E ######################

# Example usage:
plot_scRNA_RePACT_heatmap_obj <- plot_scRNA_RePACT_heatmap(
  geneset    = c(T2D.scRNA.RePACT$RePACT_call$UP[1:16], T2D.scRNA.RePACT$RePACT_call$DN[1:15]),
  RePACT_OBJ = T2D.scRNA.RePACT,
  ifShowGene = TRUE
)

plot_scRNA_RePACT_heatmap_obj <- plot_scRNA_RePACT_heatmap_obj +
  theme(
    axis.text.x   = element_text(size = 15),
    axis.text.y   = element_text(size = 15),
    axis.title.x  = element_text(size = 15),
    axis.title.y  = element_text(size = 15),
    legend.text   = element_text(size = 12),
    legend.title  = element_text(size = 14)
  )

ggsave(
  filename = "plot_scRNA_RePACT_heatmap_15.png",
  plot     = plot_scRNA_RePACT_heatmap_obj,
  width    = 10,
  height   = 8,
  dpi      = 300
)

#--------------------- Treeplot for GO Enrichment ---------------------
###################### Fig 3F ######################


# Load the GO object saved from analysis.R
go_enrich <- readRDS("go_enrich_obj.RDS")

p2 <- treeplot(go_enrich)
ggsave(
  "go_enrichment_tree_1.png",
  plot   = p2,
  width  = 10,
  height = 8,
  dpi    = 300
)

#--------------------- Junction and Transporter  Heatmap ---------------------
###################### Fig 3G ######################

# Load average expression data & gene list from the analysis script
avg_df <- readRDS("average_expression_ependymal.RDS")
gene_list_morph <- readRDS("gene_list_morph.RDS")

ependymal_filtered_MS <- avg_df[rownames(avg_df) %in% gene_list_morph, ]
ependymal_filtered_MS <- ependymal_filtered_MS[gene_list_morph, ]

HM_morphology <- pheatmap(
  ependymal_filtered_MS,
  scale         = "row",
  cluster_cols  = FALSE,
  cluster_rows  = FALSE,
  main          = "Gene Expression Heatmap"
)

ggsave(
  filename = "Strucutre_ependymal.JPEG",
  plot     = HM_morphology,
  width    = 10,
  height   = 8,
  dpi      = 300
)

#====================================================
# 4. UMAP Plots for TF Activity-Based Subsets
#====================================================

# We'll load the final object saved with TF metadata
ependymal_cells <- readRDS("ependymal_cells_TF_above_015.RDS")

# Plot UMAP for STAT3_BACH2
UMAP_STAT3_BACH2 <- DimPlot(
  ependymal_cells,
  reduction = "umap",
  label = TRUE,
  group.by = "STAT3_BACH2",
  pt.size = 0.5
) + NoLegend()

ggsave(
  "UMAP_STAT3_BACH2.png",
  plot   = UMAP_STAT3_BACH2,
  width  = 8,
  height = 6,
  dpi    = 300
)

# Plot UMAP for ZBTB20_ZNF704
UMAP_ZBTB20_ZNF704 <- DimPlot(
  ependymal_cells,
  reduction = "umap",
  label = TRUE,
  group.by = "ZBTB20_ZNF704",
  pt.size = 0.5
) + NoLegend()

ggsave(
  "UMAP_ZBTB20_ZNF704.png",
  plot   = UMAP_ZBTB20_ZNF704,
  width  = 8,
  height = 6,
  dpi    = 300
)
