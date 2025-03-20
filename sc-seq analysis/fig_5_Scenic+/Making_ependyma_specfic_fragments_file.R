#====================================================
# Outline - Subsetting Ependymal Cells and Creating a
# Fragment File for SCENICplus
#====================================================
# 1. Loading Packages
# 2. Loading Data & Exploring Structure
# 3. Subsetting Ependymal Cells
# 4. Extracting and Filtering Fragment Files
# 5. Writing Filtered Fragments to Disk
# 6. (Ubuntu Commands for Sorting & Indexing)

#====================================================
# 1. Loading Packages
#====================================================
library(Signac)
library(Seurat)
library(ggplot2)
library(dplyr)
library(data.table)

#====================================================
# 2. Loading Data & Exploring Structure
#====================================================
multiomic_ventricle_integrated <- readRDS(
  file = "E:/multiOmics/integrated object/multiome/sn_multiomic_samples.rds"
)

#====================================================
# 3. Subsetting Ependymal Cells
#====================================================
ependymal_cells <- readRDS(
  "E:/multiOmics/integrated object/snRNAseq/Subset_ependymal/ependymal_cells.RDS"
)

#====================================================
# 4. Extracting and Filtering Fragment Files
#====================================================
# Get fragment information from the ATAC assay
fragments_list <- Fragments(ependymal_cells[["ATAC"]])

# Verify fragment paths
lapply(fragments_list, function(fragment) fragment@path)

# Prepare a list of cell barcodes (trimmed to 18 characters)
subset_cells <- colnames(ependymal_cells)
filtered_subset_cells <- substr(subset_cells, 1, 18)

# Initialize a list to hold filtered fragments
filtered_fragments_list <- list()

# Loop through each fragment in the list and filter to desired cells
for (fragment in fragments_list) {
  # Read the fragment file into a data.table
  fragment_dt <- fread(fragment@path, header = FALSE)
  
  # Filter rows by matching the trimmed barcodes in column V4
  filtered_dt <- fragment_dt[V4 %in% filtered_subset_cells]
  
  # Append to the combined list
  filtered_fragments_list <- append(filtered_fragments_list, list(filtered_dt))
}

# Combine into a single data.table
filtered_fragments <- rbindlist(filtered_fragments_list, use.names = TRUE)

#====================================================
# 5. Writing Filtered Fragments to Disk
#====================================================
fwrite(filtered_fragments, file = "filtered_fragments.tsv", sep = "\t")

#====================================================
# 6. (Ubuntu Commands for Sorting & Indexing)
#====================================================
# The following commands are run externally (in a terminal),
# once you're in the directory containing 'filtered_fragments.tsv':
#
# sed '1d' filtered_fragments.tsv | sort -k1,1 -k2,2n -k3,3n | bgzip > filtered_fragments_sorted.tsv.gz
#
# tabix -p bed filtered_fragments_sorted.tsv.gz
#
# These commands remove any header, sort the file, then compress and index it.
