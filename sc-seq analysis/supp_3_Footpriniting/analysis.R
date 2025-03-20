#====================================================
# Outline of the Script
#====================================================
# 1. Loading Packages 
#    - Load necessary libraries for ATAC-seq processing, visualization, and motif analysis.
#
# 2. Making Ependymal Object
#    - Load integrated single-cell dataset.
#    - Extract ependymal cells and prepare ATAC count matrix.
#    - Create a Chromatin Assay and Seurat object.
#
# 3. Adding Gene Annotations
#    - Retrieve gene annotations from EnsDb.
#    - Convert chromosome names to UCSC format.
#    - Assign annotations to the ATAC object.
#
# 4. Computing QC Metrics
#    - Compute nucleosome signal and transcription start site (TSS) enrichment.
#
# 5. Dimensionality Reduction & Clustering
#    - Normalize data, perform PCA/UMAP, and identify cell clusters.
#
# 6. Adding Metadata for Disease
#    - Extract and integrate disease metadata into the ATAC-seq object.
#
# 7. Computing Motif Statistics
#    - Identify motifs, add motif information, and perform motif footprinting.
#
# 8. Making Accessibility Plots
#    - Compute motif activity using ChromVAR.
#    - Identify differentially accessible motifs in MS vs. control.
#    - Save results to a CSV file.

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
# Making Ependymal Object
#====================================================

# Load integrated multiomic object
multiomic_ventricle_integrated <- readRDS("path_to_your_integrated_ventricle_object")

# Extract only ependymal cells
ependymal_cells <- subset(multiomic_ventricle_integrated, subset = CellType == "ependymal")

# Get raw ATAC count matrix
mtx <- GetAssayData(object = ependymal_cells, assay = "ATAC", slot = "counts")

# Trim each barcode to the first 18 characters
colnames(mtx) <- substr(colnames(mtx), 1, 18)

# ------------------- Creating Chromatin Assay -------------------
test_assay <- CreateChromatinAssay(
  counts = mtx,  
  min.cells = 5,  
  min.features = 200,
  fragments = "E:/multiOmics/integrated object/multiome/filtered_fragments_sorted.tsv.gz",
  sep = c(":", "-"),  
  genome = "hg38"
)

# Create Seurat object
ATAC_ependymal <- CreateSeuratObject(
  counts = test_assay,
  assay = 'ATAC'
)

#====================================================
# Adding Gene Annotations
#====================================================

# Extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# Convert Ensembl chromosome names to UCSC format (e.g., "1" → "chr1")
seqlevelsStyle(annotations) <- "UCSC"

# Assign annotations to ATAC object
Annotation(ATAC_ependymal) <- annotations

#====================================================
# Computing QC Metrics
#====================================================

# Compute nucleosome signal score per cell
ATAC_ependymal <- NucleosomeSignal(ATAC_ependymal)

# Compute TSS enrichment score per cell
ATAC_ependymal <- TSSEnrichment(object = ATAC_ependymal, fast = FALSE)

View(ATAC_ependymal@meta.data)

#====================================================
# Dimensionality Reduction & Clustering
#====================================================

# Normalization and feature selection
ATAC_ependymal <- RunTFIDF(ATAC_ependymal)  # Normalization
ATAC_ependymal <- FindTopFeatures(ATAC_ependymal, min.cutoff = 'q0')  # Selecting top features
ATAC_ependymal <- RunSVD(ATAC_ependymal)  # Dimensionality reduction

# Non-linear dimensional reduction and clustering
ATAC_ependymal <- RunUMAP(object = ATAC_ependymal, reduction = 'lsi', dims = 1:30)
ATAC_ependymal <- FindNeighbors(object = ATAC_ependymal, reduction = 'lsi', dims = 2:30)
ATAC_ependymal <- FindClusters(object = ATAC_ependymal, algorithm = 3)

#====================================================
# Adding Metadata for Disease
#====================================================

# Extract the metadata
metadata <- ependymal_cells@meta.data
current_barcodes <- rownames(metadata)

# Trim each barcode to the first 18 characters
filtered_barcodes <- substr(current_barcodes, 1, 18)

# Extract relevant disease metadata
selected_metadata <- ependymal_cells@meta.data[, c("disease"), drop = FALSE]

# Update barcode names in metadata
barcodes <- rownames(selected_metadata)
stripped_barcodes <- substr(barcodes, 1, 18)
rownames(selected_metadata) <- stripped_barcodes

# Match barcodes and add metadata to ATAC object
atac_barcodes <- rownames(ATAC_ependymal@meta.data)
matched_metadata <- selected_metadata[match(atac_barcodes, stripped_barcodes), ]

# Integrate metadata
ATAC_ependymal@meta.data <- cbind(ATAC_ependymal@meta.data, matched_metadata)

# Rename the 'matched_metadata' column to 'disease'
colnames(ATAC_ependymal@meta.data)[colnames(ATAC_ependymal@meta.data) == "matched_metadata"] <- "disease"

#====================================================
# Computing Motif Statistics
#====================================================

DefaultAssay(ATAC_ependymal) <- "ATAC"
Idents(ATAC_ependymal) <- "disease"

# Extract position frequency matrices (PFM) for motifs
pwm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

# Add motif information to the object
ATAC_ependymal <- AddMotifs(
  object = ATAC_ependymal,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pwm
)

# Perform motif footprinting
ATAC_ependymal <- Footprint(
  object = ATAC_ependymal,
  motif.name = c("STAT3", "BACH2", "STAT1", "RFX2", "RFX3"),
  genome = BSgenome.Hsapiens.UCSC.hg38,
  fragments = "E:/multiOmics/Final_demultiplexed object/footprinting/ependymal_cells.tsv.gz"
)

#====================================================
# Making Accessibility Plots
#====================================================

# Compute motif activity
ATAC_ependymal <- RunChromVAR(
  object = ATAC_ependymal,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# Extract chromVAR data
chromvar_data <- ATAC_ependymal@assays[["chromvar"]]@data

# Extract motif names
motif_map <- ATAC_ependymal@assays[["ATAC"]]@motifs@motif.names

# Identify differential motif accessibility
deviations_markers <- FindMarkers(
  object = ATAC_ependymal,
  ident.1 = "MS",  
  ident.2 = "CTRL",  
  assay = "chromvar",
  test.use = "wilcox",      
  logfc.threshold = 0.25
)

# Map motif names
deviations_markers$MotifName <- sapply(rownames(deviations_markers), function(motif_id) {
  motif_map[[motif_id]] %||% motif_id
})

# Save results
write.csv(
  deviations_markers, 
  file = "TF_names_differential_motif_accessibility_MS_vs_CTRL.csv",
  row.names = TRUE
)


