library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)

setwd("E:/005---ThirdProject/ThirdObject/0.RealData/")


# Load full matrix and annotation
mtx <- readMM(gzfile("GSE234129/GSE234129_count_matrix.mtx.gz"))
features <- read.delim(gzfile("GSE234129/GSE234129_features.tsv.gz"), header = FALSE)
barcodes <- read.delim(gzfile("GSE234129/GSE234129_barcodes.tsv.gz"), header = FALSE)
meta <- read.delim(gzfile("GSE234129/GSE234129_meta.tsv.gz"), header = TRUE)

# Assign row and column names
rownames(mtx) <- features$V1
colnames(mtx) <- barcodes$V1

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = mtx, project = "Gastric_all")

# Ensure metadata is ordered correctly
meta <- meta[match(colnames(seurat_obj), meta$cell_barcodes), ]
seurat_obj <- AddMetaData(seurat_obj, metadata = meta)

# Optional cleaning if needed:
# Remove any NA or mismatched rows
seurat_obj <- subset(seurat_obj, cells = which(!is.na(seurat_obj$celltype)))

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)


# 
# # Define new column for T cells
# seurat_obj$T_cell_status <- ifelse(seurat_obj$celltype %in% c("CD4T", "CD8T", "CD8_C7", "CD8_C2","CD8_C0", "CD4_C4" ), "T_cell", "Other")
# 
# p1 <- DimPlot(seurat_obj, group.by = "T_cell_status", cols = c("grey", "#4DBBD5FF"), pt.size = 0.4) +
#   ggtitle("UMAP: T cells highlighted") +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# p1
# 
# 
# p2 <- FeaturePlot(seurat_obj, features = "CD8A", cols = c("grey90", "red"), pt.size = 0.4) +
#   ggtitle("CD8A Expression") +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# p2
# 
# p2 <- FeaturePlot(seurat_obj, features = "CD4_C4", cols = c("grey90", "blue"), pt.size = 0.4) +
#   ggtitle("CD8A Expression") +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# p2
# 
# 








c("#A2C2E2", "#FEAF8A", "#F6B7C6", "#A2DADE", "#D8CBF0")


# Create new metadata column for T cell labeling
seurat_obj$T_cell_status <- ifelse(grepl("^CD4|^CD8", seurat_obj$celltype), "T_cell", "Other")

# Create new metadata column for CD8 T cell labeling
seurat_obj$CD8_status <- ifelse(grepl("^CD8", seurat_obj$celltype), "CD8_T", "Other")
seurat_obj$CD4_status <- ifelse(grepl("^CD4", seurat_obj$celltype), "CD4_T", "Other")

p_T <- DimPlot(seurat_obj, group.by = "T_cell_status", cols = c("grey80", "#A2DADE"), pt.size = 0.4) +
  ggtitle("T cells vs Other") +
  theme(plot.title = element_text(hjust = 0.5))

p_T

# Save as PDF
pdf("PLOT_UMAP/Plot_GSE234129_Tcells.pdf", width = 4, height = 4)
p_T
dev.off()

# Save as SVG
svg("PLOT_UMAP/Plot_GSE234129_Tcells.svg", width = 4, height = 4)
p_T
dev.off()



p_CD8 <- DimPlot(seurat_obj, group.by = "CD8_status", cols = c("#FEAF8A", "grey80"), pt.size = 0.4) +
  ggtitle("CD8 T cells vs Other") +
  theme(plot.title = element_text(hjust = 0.5))


p_CD8





p_CD4 <- DimPlot(seurat_obj, group.by = "CD4_status", cols = c("#A2C2E2", "grey80"), pt.size = 0.4) +
  ggtitle("CD4 T cells vs Other") +
  theme(plot.title = element_text(hjust = 0.5))


p_CD4



# Save as PDF
pdf("PLOT_UMAP/Plot_GSE234129_CD4.pdf", width = 4, height = 4)
p_CD4
dev.off()

# Save as SVG
svg("PLOT_UMAP/Plot_GSE234129_CD4.svg", width = 4, height = 4)
p_CD4
dev.off()














library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)

# Set your working directory containing only the relevant folders
base_dir <- "Sathe/gastric_scRNAseq_filtered"  # change to your path
target_folders <- c("5866_t1", "5931_t1", "5866_t2", "6342_t1", "5931_t2", 
                    "6592_t1", "5846_t1", "6207_t1", "6709_t1")

seurat_list <- list()

for (folder in target_folders) {
  path <- file.path(base_dir, folder)
  
  mtx <- readMM(gzfile(file.path(path, "matrix.mtx")))
  features <- read.delim(gzfile(file.path(path, "features.tsv")), header = FALSE)
  barcodes <- read.delim(gzfile(file.path(path, "barcodes.tsv")), header = FALSE)
  
  rownames(mtx) <- features$V1
  # Make column names unique with folder prefix
  clean_barcodes <- sub("-\\d+$", "", barcodes$V1)  # remove suffix
  colnames(mtx) <- paste0(folder, "_", clean_barcodes)
  
  obj <- CreateSeuratObject(counts = mtx, project = folder)
  obj$sample <- folder
  
  seurat_list[[folder]] <- obj
}

# Merge all samples
seurat_obj <- merge(seurat_list[[1]], y = seurat_list[-1])

# Standard Seurat workflow
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# Load metadata
cell_labels <- read.csv("Sathe/cell_labels.csv")

# Remove suffix from cell_labels barcodes
cell_labels$clean_barcode <- sub("-\\d+$", "", cell_labels$cell_barcode)

# Construct full barcodes matching Seurat object naming
cell_labels$barcode_full <- paste0(cell_labels$orig.ident, "_", cell_labels$clean_barcode)

# Subset only cells with matching barcodes
common_cells <- intersect(Cells(seurat_obj), cell_labels$barcode_full)
message("Number of common cells found: ", length(common_cells))

# Proceed if some cells match
if (length(common_cells) > 0) {
  rownames(cell_labels) <- cell_labels$barcode_full
  seurat_obj <- subset(seurat_obj, cells = common_cells)
  seurat_obj$celltype <- cell_labels[Cells(seurat_obj), "final_celltype"]
} else {
  stop("No matching cells found after barcode cleaning.")
}

# Label T cells and CD8 T cells
seurat_obj$T_cell_status <- ifelse(grepl("^CD4|^CD8", seurat_obj$celltype), "T_cell", "Other")
seurat_obj$CD8_status <- ifelse(grepl("^CD8", seurat_obj$celltype), "CD8_T", "Other")
seurat_obj$CD4_status <- ifelse(grepl("^CD4", seurat_obj$celltype), "CD4_T", "Other")

# Plot T cells
p_T <- DimPlot(seurat_obj, group.by = "T_cell_status", cols = c("grey80", "#A2DADE"), pt.size = 0.4) +
  ggtitle("T cells vs Other") +
  theme(plot.title = element_text(hjust = 0.5))


p_T


# Plot CD8 T cells
p_CD8 <- DimPlot(seurat_obj, group.by = "CD8_status", cols = c("#FEAF8A", "grey80"), pt.size = 0.4) +
  ggtitle("CD8 T cells vs Other") +
  theme(plot.title = element_text(hjust = 0.5))

p_CD8


# Plot CD8 T cells
p_CD4 <- DimPlot(seurat_obj, group.by = "CD4_status", cols = c("#A2C2E2", "grey80"), pt.size = 0.4) +
  ggtitle("CD4 T cells vs Other") +
  theme(plot.title = element_text(hjust = 0.5))

p_CD4





# Create output folder if not exists
#dir.create("PLOT_UMAP", showWarnings = FALSE)

# Save plots
pdf("PLOT_UMAP/Plot_Sathe_Tcells.pdf", width = 4, height = 4)
print(p_T)
dev.off()

svg("PLOT_UMAP/Plot_Sathe_Tcells.svg", width = 4, height = 4)
print(p_T)
dev.off()

pdf("PLOT_UMAP/Plot_Sathe_CD8.pdf", width = 4, height = 4)
print(p_CD8)
dev.off()

svg("PLOT_UMAP/Plot_Sathe_CD8.svg", width = 4, height = 4)
print(p_CD8)
dev.off()


pdf("PLOT_UMAP/Plot_Sathe_CD4.pdf", width = 4, height = 4)
print(p_CD4)
dev.off()

svg("PLOT_UMAP/Plot_Sathe_CD4.svg", width = 4, height = 4)
print(p_CD4)
dev.off()









library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)

setwd("E:/005---ThirdProject/ThirdObject/0.RealData/")

# Define base folder
base_folder <- "GSE134520_RAW"

# Find all dt.*.txt files inside subfolders
file_paths <- list.files(base_folder, pattern = "^dt\\..*\\.txt$", recursive = TRUE, full.names = TRUE)

# Function to load and create Seurat object from a single file
load_seurat_from_txt <- function(filepath) {
  message("Loading: ", filepath)
  raw_matrix <- read.table(filepath, header = TRUE, row.names = 1)
  sample_name <- tools::file_path_sans_ext(basename(filepath))  # e.g., dt.EGC
  seurat_obj <- CreateSeuratObject(counts = raw_matrix, min.features = 200)
  seurat_obj$sample <- sample_name
  return(seurat_obj)
}

# Load all files into Seurat objects
seurat_list <- lapply(file_paths, load_seurat_from_txt)

# Merge all into one object
combined <- merge(seurat_list[[1]], y = seurat_list[-1])

# ---- Seurat Preprocessing ----
# 1. Quality control
combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-")
combined <- subset(combined, subset = nFeature_RNA > 200 & nFeature_RNA < 6500 & percent.mt < 15)

# 2. Remove genes expressed in fewer than 3 cells
#counts <- GetAssayData(combined, slot = "counts", layer = "counts")
genes_use <- rowSums(counts > 0) >= 3
combined <- subset(combined, features = names(genes_use[genes_use]))

# 3. Normalize, scale, PCA, clustering
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
combined <- ScaleData(combined)
combined <- RunPCA(combined, features = VariableFeatures(combined))
combined <- FindNeighbors(combined, dims = 1:20)
combined <- FindClusters(combined, resolution = 0.5)
combined <- RunUMAP(combined, dims = 1:20)

# ---- Identify T/CD4/CD8 Cells by Marker Expression ----
# Add CD4/CD8 marker expressions manually
cd4_cells <- subset(combined, subset = CD4 > 1 & IL7R > 1)
cd8_cells <- subset(combined, subset = CD8A > 1 | CD8B > 1)

barcodes_cd4 <- colnames(cd4_cells)
barcodes_cd8 <- colnames(cd8_cells)

# Add labels to combined object
combined$CD4_label <- ifelse(Cells(combined) %in% barcodes_cd4, "CD4_T", "Other")
combined$CD8_label <- ifelse(Cells(combined) %in% barcodes_cd8, "CD8_T", "Other")
combined$T_label <- ifelse(combined$CD4_label == "CD4_T" | combined$CD8_label == "CD8_T", "T_cell", "Other")

# ---- Plotting ----
#dir.create("PLOT_UMAP", showWarnings = FALSE)

# # CD4
# combined$CD4_label <- factor(combined$CD4_label, levels = c("Other", "CD4_T"))
# 
# p_cd4 <- DimPlot(combined, group.by = "CD4_label", cols = c("grey80", "#A2C2E2"), pt.size = 0.4) +
#   ggtitle("CD4 T cells") + theme(plot.title = element_text(hjust = 0.5))
# 
# p_cd4

# CD8
combined$CD8_label <- factor(combined$CD8_label, levels = c("Other", "CD8_T"))
p_cd8 <- DimPlot(combined, group.by = "CD8_label", cols = c("grey80", "#FEAF8A"), pt.size = 0.4) +
  ggtitle("CD8 T cells") + theme(plot.title = element_text(hjust = 0.5))

p_cd8

# T cells
combined$T_label <- factor(combined$T_label, levels = c("Other", "T_cell"))
p_tcell <- DimPlot(combined, group.by = "T_label", cols = c("grey80", "#A2DADE"), pt.size = 0.4) +
  ggtitle("T cells (CD4 or CD8)") + theme(plot.title = element_text(hjust = 0.5))

p_tcell


# Save plots
pdf("PLOT_UMAP/Plot_GSE134520_Tcells.pdf", width = 4, height = 4)
print(p_tcell)
dev.off()

svg("PLOT_UMAP/Plot_GSE134520_Tcells.svg", width = 4, height = 4)
print(p_tcell)
dev.off()




# Manually extract CD8_T and Other cells
# cd8_cells <- WhichCells(combined, expression = CD8_label == "CD8_T")
# other_cd8_cells <- WhichCells(combined, expression = CD8_label == "Other")
# 
# # Plot CD8_T cells last (on top)
# p_cd8 <- DimPlot(combined,
#                  group.by = "CD8_label",
#                  cells = c(other_cd8_cells, cd8_cells),
#                  cols = c("Other" = "grey80", "CD8_T" = "#FEAF8A"),
#                  pt.size = 0.4) +
#   ggtitle("CD8 T cells") +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# p_cd8


pdf("PLOT_UMAP/Plot_GSE134520_CD8.pdf", width = 4, height = 4)
print(p_cd8)
dev.off()

svg("PLOT_UMAP/Plot_GSE134520_CD8.svg", width = 4, height = 4)
print(p_cd8)
dev.off()








# First, manually split cells based on CD4 label
cd4_cells <- WhichCells(combined, expression = CD4_label == "CD4_T")
other_cells <- WhichCells(combined, expression = CD4_label == "Other")

# Then, supply cells in desired order to DimPlot
p_cd4 <- DimPlot(combined, group.by = "CD4_label", 
                 cells = c(other_cells, cd4_cells),  # CD4 cells come last (on top)
                 cols = c("Other" = "grey80", "CD4_T" = "#A2C2E2"),
                 pt.size = 0.4) +
  ggtitle("CD4 T cells") +
  theme(plot.title = element_text(hjust = 0.5))

p_cd4

pdf("PLOT_UMAP/Plot_GSE134520_CD4.pdf", width = 4, height = 4)
print(p_cd4)
dev.off()

svg("PLOT_UMAP/Plot_GSE134520_CD4.svg", width = 4, height = 4)
print(p_cd4)
dev.off()




