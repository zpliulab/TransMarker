install.packages("Seurat")
install.packages("harmony")
install.packages("SingleCellExperiment")

library(Seurat)
library(dplyr)
library(harmony)
library(ggplot2)
library(SingleCellExperiment)

setwd("E:/005---ThirdProject/ThirdObject/0.RealData/")

# 1. Load raw matrix
raw_matrix <- read.table("ESCC/GSE199654/GSE199654_scTDN_UMI_matrix_epithelial_cells.txt/scTDN_UMI_matrix_epithelial_cells.txt", header = T)
raw_matrix <- read.table("ESCC/GSE160269/GSE160269_UMI_matrix_Tcell.txt/UMI_matrix_Tcell.txt", header = T)


seurat_obj <- CreateSeuratObject(counts = raw_matrix, min.features = 200)

# 2. Quality Control: remove cells with <200 genes or >6500 genes (possible doublets), and high mito%
# Quality Control Filtering
# - Filter out cells with <200 genes
# - Filter cells with >15% mitochondrial genes
# - Filter cells with >6500 genes (doublets)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & 
                       nFeature_RNA < 6500 & 
                       percent.mt < 15)

# 3. Remove genes expressed in fewer than 3 cells (alternative method)
counts <- GetAssayData(seurat_obj, slot = "counts")
genes_use <- rowSums(counts > 0) >= 3
seurat_obj <- subset(seurat_obj, features = names(genes_use[genes_use]))

# 4. Normalize the data
seurat_obj <- NormalizeData(seurat_obj)

# 5. Find variable features
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# 6. Scale data and run PCA
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
#seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj), npcs = 8)

# 7. (Optional) Harmony batch correction
#seurat_obj <- RunHarmony(seurat_obj, group.by.vars = "batch_column")


# 8. Clustering
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
#seurat_obj <- FindNeighbors(seurat_obj, dims = 1:8)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
#seurat_obj <- RunUMAP(seurat_obj, dims = 1:8)


# 9. Visualize known markers for CD4 and CD8 T cells
FeaturePlot(seurat_obj, features = c("CD4", "IL7R"), reduction = "umap")
FeaturePlot(seurat_obj, features = c("CD8A", "CD8B"), reduction = "umap")

# 10. Find markers to annotate clusters
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(markers)

# 11. Subset CD4+ and CD8+ cells
cd4_cells <- subset(seurat_obj, subset = CD4 > 1 & IL7R > 1)
cd8_cells <- subset(seurat_obj, subset = CD8A > 1 | CD8B > 1)

# myeloid_cells <- subset(seurat_obj, subset = LYZ > 1 & CD14 > 1)
# b_cells <- subset(seurat_obj, subset = MS4A1 > 1 | CD79A > 1)
# nk_cells <- subset(seurat_obj, subset = NKG7 > 1 & GNLY > 1)


# 12. Save barcodes of CD4 and CD8 cells
barcodes_cd4 <- colnames(cd4_cells)
barcodes_cd8 <- colnames(cd8_cells)

barcodes_cd4
barcodes_cd8

Mat_CD4 <- raw_matrix[ , colnames(raw_matrix) %in% barcodes_cd4, drop = FALSE]
Mat_CD8 <- raw_matrix[ , colnames(raw_matrix) %in% barcodes_cd8, drop = FALSE]



# 
# write.table(barcodes_cd4, file = "ESCC/GSE199654/ESCC_CD4_cell_barcodes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
# write.table(barcodes_cd8, file = "ESCC/GSE199654/ESCC_CD8_cell_barcodes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
# 
# 
# 
# write.table(Mat_CD4, file = "ESCC/GSE199654/ESCC_CD4.txt", quote = FALSE, row.names = T, col.names = T)
# write.table(Mat_CD4, file = "ESCC/GSE199654/ESCC_CD8.txt", quote = FALSE, row.names = T, col.names = T)
# 

barcodes_b_cells <- colnames(b_cells)
write.table(barcodes_b_cells, file = "ESCC/GSE199654/ESCC_barcodes_b_cells_barcodes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)



barcodes_myeloid_cells <- colnames(myeloid_cells)
write.table(barcodes_myeloid_cells, file = "ESCC/GSE199654/ESCC_barcodes_myeloid_cells_barcodes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)


barcodes_nk_cells <- colnames(nk_cells)
write.table(barcodes_nk_cells, file = "ESCC/GSE199654/ESCC_barcodes_nk_cells_barcodes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

