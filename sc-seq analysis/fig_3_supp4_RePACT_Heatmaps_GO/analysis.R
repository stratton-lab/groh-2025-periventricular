#====================================================
# Outline - Ependymal Cell Analysis
#====================================================
# 1. Loading Packages and Setting Seed
# 2. Loading Integrated Ependymal Cells & Preprocessing
# 3. RePACT Analysis Setup
# 4. Differential Expression Between MS and Control
# 5. Gene Ontology on Ependymal Upregulated Genes
# 6. Morphology/Tight Junction Heatmaps
# 7. TF Activity-Based Cell Subsetting & GO Analysis

#====================================================
# 1. Loading Packages and Setting Seed
#====================================================
library(Seurat)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(RePACT)

setwd("your_wd")
set.seed(1)

#====================================================
# 2. Loading Integrated Ependymal Cells & Preprocessing
#====================================================
ependymal_cells <- readRDS("E:/multiOmics/integrated object/snRNAseq/Subset_ependymal/ependymal_cells.RDS")

#--------------------- RePACT Setup ---------------------
raw_counts <- GetAssayData(object = ependymal_cells, layer = "counts")
ependymal_reclustered <- CreateSeuratObject(counts = raw_counts)

# Transfer metadata
ependymal_reclustered$cluster_ids <- ependymal_cells$seurat_clusters[colnames(ependymal_reclustered)]
ependymal_reclustered$Group       <- ependymal_cells$Group[colnames(ependymal_reclustered)]
ependymal_reclustered$Donor       <- ependymal_cells$Donor[colnames(ependymal_reclustered)]
ependymal_reclustered$percent.mt  <- ependymal_cells$percent.mt[colnames(ependymal_reclustered)]

# Rename metadata columns
meta_data <- ependymal_reclustered@meta.data
names(meta_data)[names(meta_data) == "Donor"] <- "Sample"
names(meta_data)[names(meta_data) == "Group"] <- "diseaseStat"
ependymal_reclustered@meta.data <- meta_data

# Basic scRNA-Seq Steps
ependymal_reclustered <- NormalizeData(ependymal_reclustered)
ependymal_reclustered <- FindVariableFeatures(ependymal_reclustered, selection.method = "vst", nfeatures = 2000)
ependymal_reclustered <- ScaleData(ependymal_reclustered, features = rownames(ependymal_reclustered))
ependymal_reclustered <- RunPCA(ependymal_reclustered, features = VariableFeatures(object = ependymal_reclustered))

#====================================================
# 3. RePACT Analysis Setup
#====================================================

###################### Fig 3E ######################
###################### Fig 3D ######################

# (Typically run on Ubuntu if there are multi-core issues on Windows. Use the scRNA.RePACT function see on the RePACT github )
# dont forget to convert you seurat object to V3 before you run scRNA.RePACT
# Example file: "out.rds" is your RePACT results

T2D.scRNA.RePACT <- readRDS("out.rds")

#====================================================
# 4. Differential Expression Between MS and Control
#====================================================
Idents(ependymal_cells) <- "Group"
cluster_Light.markers <- FindMarkers(
  ependymal_cells,
  ident.1 = "MS",
  ident.2 = "CTRL",
  avg_log2FC = 0.25
)

degs_MS_vs_CTRL <- cluster_Light.markers[cluster_Light.markers$avg_log2FC > 0, ]
degs_MS_vs_CTRL <- degs_MS_vs_CTRL[order(degs_MS_vs_CTRL$avg_log2FC, decreasing = TRUE), ]

write.csv(degs_MS_vs_CTRL, file = "MS_CTRL_degs_MS_vs_CTRL.csv")

#====================================================
# 5. Gene Ontology on Ependymal Upregulated Genes
#====================================================
###################### Fig 3F ######################

# Filter for p_val_adj < 0.05 & avg_log2FC > 1
degs_MS_vs_CTRL_SIG_1 <- subset(degs_MS_vs_CTRL, p_val_adj < 0.05 & avg_log2FC > 1)

gene_list <- rownames(degs_MS_vs_CTRL_SIG_1)
gene_list_entrez <- mapIds(
  org.Hs.eg.db,
  keys      = gene_list,
  column    = "ENTREZID",
  keytype   = "SYMBOL",
  multiVals = "first"
)
gene_list_entrez <- na.omit(gene_list_entrez)

go_enrich <- enrichGO(
  gene          = gene_list_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "ALL",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)
go_enrich <- pairwise_termsim(go_enrich)

saveRDS(go_enrich, file = "go_enrich_obj.RDS")  # Save for plotting

#====================================================
# 6. Morphology/Tight Junction Heatmaps
#====================================================
###################### Fig 3G ######################

# Save the needed expression data for later plotting
gene_list_morph <- c(
  "DST","CTNNA3","APC","CLIC4","DLG1","CNN3","TNC","CD44","ACTN2","LDLRAD3","CTNNA2",
  "CADM1","SAMD4A","SHISA9","DLG2","CTNND2","GRIK1","GRIN2A","PSD3","CADPS","GRID1",
  "CADM2","CAMK2D","PCDH9","THBS2","NRCAM","MACF1","SORBS1","CCDC85A","FRMD5","PDPN",
  "CTNNA1","NCAM1","ADCY8","ERBIN","HGF","ARHGEF3","AQP4","VCAN","GFAP","MAPT","MAP2",
  "TUBB3","CHI3L1","SLC1A2"
)

avg_exp_data <- AverageExpression(ependymal_cells, group.by = "Donor", return.seurat = FALSE)
# Convert to a data frame so it can be saved
avg_df <- as.data.frame(as.matrix(avg_exp_data[["RNA"]]))

saveRDS(avg_df, "average_expression_ependymal.RDS")
saveRDS(gene_list_morph, "gene_list_morph.RDS")

#====================================================
# 7. TF Activity-Based Cell Subsetting & GO Analysis
#====================================================

DefaultAssay(ependymal_cells) <- "RNA"

#--------------------- BACH2 & STAT3 ---------------------
# Took all cells wiith eRSS score above 0.15 

BACH2_cells <- read.csv("E:/multiOmics/cleaned_code/Elias_try/REPACT_HEATMAPS_GO/BACH2_direct__51g_barcodes_over_0.15.csv")
STAT3_cells <- read.csv("E:/multiOmics/cleaned_code/Elias_try/REPACT_HEATMAPS_GO/STAT3_direct__15g_barcodes_over_0.15.csv")

BACH2_cells <- BACH2_cells[-1, , drop = FALSE]
STAT3_cells <- STAT3_cells[-1, , drop = FALSE]

ependyml_barcodes         <- rownames(ependymal_cells@meta.data)
cleaned_ependyml_barcodes <- stringr::str_replace(ependyml_barcodes, "_.*$", "")

combined_barcodes         <- unique(c(STAT3_cells$Cell, BACH2_cells$Cell))
matching_barcodes         <- cleaned_ependyml_barcodes %in% combined_barcodes
subsetted_ependymal_cells <- ependymal_cells[, matching_barcodes]

ependymal_cells@meta.data$STAT3_BACH2 <- ifelse(
  cleaned_ependyml_barcodes %in% combined_barcodes,
  "above_0.15_auc",
  "below_0.15_auc"
)

Idents(ependymal_cells) <- "STAT3_BACH2"
STAT3_BACH2_markers <- FindMarkers(
  object   = ependymal_cells,
  ident.1  = "above_0.15_auc",
  ident.2  = "below_0.15_auc"
)
write.csv(STAT3_BACH2_markers, file = "STAT3_BACH2_markers.csv", row.names = TRUE)

saveRDS(STAT3_BACH2_markers, file = "STAT3_BACH2_markers_obj.RDS")

#--------------------- GO Enrichment for STAT3_BACH2 Markers ---------------------
###################### Supp Fig 4A ######################

# Filter significant DEGs
degs_MS_vs_CTRL_SIG_1 <- subset(STAT3_BACH2_markers, p_val_adj < 0.05 & avg_log2FC > 0.25)

# Convert gene symbols to ENTREZ IDs
gene_list <- rownames(degs_MS_vs_CTRL_SIG_1)
gene_list_entrez <- mapIds(
  org.Hs.eg.db,
  keys      = gene_list,
  column    = "ENTREZID",
  keytype   = "SYMBOL",
  multiVals = "first"
)
gene_list_entrez <- na.omit(gene_list_entrez)

# Enrich GO
go_enrich <- enrichGO(
  gene          = gene_list_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "ALL",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)
go_enrich <- pairwise_termsim(go_enrich)

p2 <- treeplot(go_enrich)
p2

ggsave(
  "BACH2_STAT3_GO.png",
  plot   = p2,
  width  = 8,
  height = 6,
  dpi    = 300
)

#--------------------- ZNF707 & ZBTB20 ---------------------
ZBTB20_cells <- read.csv("E:/multiOmics/cleaned_code/Elias_try/REPACT_HEATMAPS_GO/ZBTB20_direct__64g_barcodes_over_0.15.csv")
ZNF704_cells <- read.csv("E:/multiOmics/cleaned_code/Elias_try/REPACT_HEATMAPS_GO/ZNF704_direct__92g_barcodes_over_0.15.csv")

ZBTB20_cells <- ZBTB20_cells[-1, , drop = FALSE]
ZNF704_cells <- ZNF704_cells[-1, , drop = FALSE]

combined_barcodes <- unique(c(ZBTB20_cells$Cell, ZNF704_cells$Cell))
ependymal_cells@meta.data$ZBTB20_ZNF704 <- ifelse(
  cleaned_ependyml_barcodes %in% combined_barcodes,
  "above_0.15_auc",
  "below_0.15_auc"
)

Idents(ependymal_cells) <- "ZBTB20_ZNF704"
ZBTB20_ZNF704_markers <- FindMarkers(
  object   = ependymal_cells,
  ident.1  = "above_0.15_auc",
  ident.2  = "below_0.15_auc"
)

write.csv(ZBTB20_ZNF704_markers, file = "ZBTB20_ZNF704_markers.csv", row.names = TRUE)
saveRDS(ZBTB20_ZNF704_markers, file = "ZBTB20_ZNF704_markers_obj.RDS")

#--------------------- GO Enrichment for ZBTB20_ZNF704 Markers ---------------------
###################### Supp Fig 4A ######################

# Step 11: Filter significant DEGs
degs_MS_vs_CTRL_SIG_1 <- subset(ZBTB20_ZNF704_markers, p_val_adj < 0.05 & avg_log2FC > 0.25)

# Step 12: Convert gene symbols to ENTREZ IDs for enrichment analysis
gene_list <- rownames(degs_MS_vs_CTRL_SIG_1)
gene_list_entrez <- mapIds(
  org.Hs.eg.db,
  keys      = gene_list,
  column    = "ENTREZID",
  keytype   = "SYMBOL",
  multiVals = "first"
)
gene_list_entrez <- na.omit(gene_list_entrez)

# Step 13: Order ENTREZ IDs by log-fold change
degs_MS_vs_CTRL_SIG_1 <- degs_MS_vs_CTRL_SIG_1[order(degs_MS_vs_CTRL_SIG_1$avg_log2FC, decreasing = TRUE), ]
gene_list_entrez_ordered <- gene_list_entrez[match(rownames(degs_MS_vs_CTRL_SIG_1), names(gene_list_entrez))]

# Step 14: Run GO enrichment analysis
go_enrich <- enrichGO(
  gene          = gene_list_entrez_ordered,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "ALL",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)
go_enrich <- pairwise_termsim(go_enrich)

# Step 15: Create a tree plot of the enrichment results
p2 <- treeplot(go_enrich)
p2

ggsave(
  "ZBTB20_ZNF704_GO.png",
  plot   = p2,
  width  = 8,
  height = 6,
  dpi    = 300
)

saveRDS(ependymal_cells, "ependymal_cells_TF_above_015.RDS")