# install.packages("Seurat")
# install.packages("harmony")
# install.packages("SingleCellExperiment")

library(Seurat)
library(dplyr)
library(harmony)
library(ggplot2)
library(SingleCellExperiment)

setwd("E:/005---ThirdProject/ThirdObject/0.RealData/")

# 1. Load raw matrix
raw_matrix <- read.table("GSE134520_RAW/GSM3954958_processed_EGC.txt/dt.EGC.txt", header = TRUE, row.names = 1)
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

# 7. (Optional) Harmony batch correction
#seurat_obj <- RunHarmony(seurat_obj, group.by.vars = "batch_column")

# 8. Clustering
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

# 9. Visualize known markers for CD4 and CD8 T cells
FeaturePlot(seurat_obj, features = c("CD4", "IL7R"), reduction = "umap")
FeaturePlot(seurat_obj, features = c("CD8A", "CD8B"), reduction = "umap")

# 10. Find markers to annotate clusters
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#View(markers)

# 11. Subset CD4+ and CD8+ cells
cd4_cells <- subset(seurat_obj, subset = CD4 > 1 & IL7R > 1)
cd8_cells <- subset(seurat_obj, subset = CD8A > 1 | CD8B > 1)

# 12. Save barcodes of CD4 and CD8 cells
barcodes_cd4 <- colnames(cd4_cells)
barcodes_cd8 <- colnames(cd8_cells)

barcodes_cd4
barcodes_cd8

write.table(barcodes_cd4, file = "GSE134520_RAW/EGC_CD4_cell_barcodes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(barcodes_cd8, file = "GSE134520_RAW/EGC_CD8_cell_barcodes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

# 13. Inspect metadata (optional)
head(cd4_cells@meta.data)
head(cd8_cells@meta.data)




################### INTEGRATE RESULTS Per STAGES ####################

# NAG CD4
NAG1_CD4 <- read.table("GSE134520_RAW/NAG1_CD4_cell_barcodes.txt", header = F)
NAG1_matrix <- read.table("GSE134520_RAW/GSM3954946_processed_NAG1.txt/dt.NAG1.txt", header = TRUE, row.names = 1)
NAG1_matrix<- as.data.frame(NAG1_matrix)
NAG1_C_4 <- NAG1_matrix[ , colnames(NAG1_matrix) %in% NAG1_CD4$V1, drop = FALSE]


NAG2_CD4 <- read.table("GSE134520_RAW/NAG2_CD4_cell_barcodes.txt", header = F)
NAG2_matrix <- read.table("GSE134520_RAW/GSM3954947_processed_NAG2.txt/dt.NAG2.txt", header = TRUE, row.names = 1)
NAG2_matrix<- as.data.frame(NAG2_matrix)
NAG2_C_4 <- NAG2_matrix[ , colnames(NAG2_matrix) %in% NAG2_CD4$V1, drop = FALSE]


NAG3_CD4 <- read.table("GSE134520_RAW/NAG3_CD4_cell_barcodes.txt", header = F)
NAG3_matrix <- read.table("GSE134520_RAW/GSM3954948_processed_NAG3.txt/dt.NAG3.txt", header = TRUE, row.names = 1)
NAG3_matrix<- as.data.frame(NAG3_matrix)
NAG3_C_4 <- NAG3_matrix[ , colnames(NAG3_matrix) %in% NAG3_CD4$V1, drop = FALSE]


NAG_CD4<-cbind(NAG1_C_4,NAG2_C_4,NAG3_C_4)
dim(NAG_CD4)
colnames(NAG_CD4) <- paste0("Zhang_", colnames(NAG_CD4))
write.table(NAG_CD4, file = "GSE134520_RAW/NAG_CD4.txt", quote = FALSE, row.names = T, col.names = T)


#NAG CD8
NAG1_CD8 <- read.table("GSE134520_RAW/NAG1_CD8_cell_barcodes.txt", header = F)
NAG1_matrix <- read.table("GSE134520_RAW/GSM3954946_processed_NAG1.txt/dt.NAG1.txt", header = TRUE, row.names = 1)
NAG1_matrix<- as.data.frame(NAG1_matrix)
NAG1_C_8 <- NAG1_matrix[ , colnames(NAG1_matrix) %in% NAG1_CD8$V1, drop = FALSE]


NAG2_CD8 <- read.table("GSE134520_RAW/NAG2_CD8_cell_barcodes.txt", header = F)
NAG2_matrix <- read.table("GSE134520_RAW/GSM3954947_processed_NAG2.txt/dt.NAG2.txt", header = TRUE, row.names = 1)
NAG2_matrix<- as.data.frame(NAG2_matrix)
NAG2_C_8 <- NAG2_matrix[ , colnames(NAG2_matrix) %in% NAG2_CD8$V1, drop = FALSE]
dim(NAG2_C_8)

NAG3_CD8 <- read.table("GSE134520_RAW/NAG3_CD8_cell_barcodes.txt", header = F)
NAG3_matrix <- read.table("GSE134520_RAW/GSM3954948_processed_NAG3.txt/dt.NAG3.txt", header = TRUE, row.names = 1)
NAG3_matrix<- as.data.frame(NAG3_matrix)
NAG3_C_8 <- NAG3_matrix[ , colnames(NAG3_matrix) %in% NAG3_CD8$V1, drop = FALSE]
dim(NAG3_C_8)

NAG_CD8 <-cbind(NAG1_C_8,NAG2_C_8,NAG3_C_8)
dim(NAG_CD8)
colnames(NAG_CD8) <- paste0("Zhang_", colnames(NAG_CD8))
write.table(NAG_CD8, file = "GSE134520_RAW/NAG_CD8.txt", quote = FALSE, row.names = T, col.names = T)



# CAG CD4
CAG1_CD4 <- read.table("GSE134520_RAW/CAG1_CD4_cell_barcodes.txt", header = F)
CAG1_matrix <- read.table("GSE134520_RAW/GSM3954949_processed_CAG1.txt/dt.CAG1.txt", header = TRUE, row.names = 1)
CAG1_matrix<- as.data.frame(CAG1_matrix)
CAG1_C_4 <- CAG1_matrix[ , colnames(CAG1_matrix) %in% CAG1_CD4$V1, drop = FALSE]


CAG2_CD4 <- read.table("GSE134520_RAW/CAG2_CD4_cell_barcodes.txt", header = F)
CAG2_matrix <- read.table("GSE134520_RAW/GSM3954950_processed_CAG2.txt/dt.CAG2.txt", header = TRUE, row.names = 1)
CAG2_matrix<- as.data.frame(CAG2_matrix)
CAG2_C_4 <- CAG2_matrix[ , colnames(CAG2_matrix) %in% CAG2_CD4$V1, drop = FALSE]


CAG3_CD4 <- read.table("GSE134520_RAW/CAG3_CD4_cell_barcodes.txt", header = F)
CAG3_matrix <- read.table("GSE134520_RAW/GSM3954951_processed_CAG3.txt/dt.CAG3.txt", header = TRUE, row.names = 1)
CAG3_matrix<- as.data.frame(CAG3_matrix)
CAG3_C_4 <- CAG3_matrix[ , colnames(CAG3_matrix) %in% CAG3_CD4$V1, drop = FALSE]


CAG_CD4<-cbind(CAG1_C_4,CAG2_C_4,CAG3_C_4)
CAG_CD4<-cbind(CAG2_C_4,CAG3_C_4)
dim(CAG_CD4)
colnames(CAG_CD4) <- paste0("Zhang_", colnames(CAG_CD4))
write.table(CAG_CD4, file = "GSE134520_RAW/CAG_CD4.txt", quote = FALSE, row.names = T, col.names = T)


#CAG CD8
CAG1_CD8 <- read.table("GSE134520_RAW/CAG1_CD8_cell_barcodes.txt", header = F)
CAG1_matrix <- read.table("GSE134520_RAW/GSM3954949_processed_CAG1.txt/dt.CAG1.txt", header = TRUE, row.names = 1)
CAG1_matrix<- as.data.frame(CAG1_matrix)
CAG1_C_8 <- CAG1_matrix[ , colnames(CAG1_matrix) %in% CAG1_CD8$V1, drop = FALSE]


CAG2_CD8 <- read.table("GSE134520_RAW/CAG2_CD8_cell_barcodes.txt", header = F)
CAG2_matrix <- read.table("GSE134520_RAW/GSM3954950_processed_CAG2.txt/dt.CAG2.txt", header = TRUE, row.names = 1)
CAG2_matrix<- as.data.frame(CAG2_matrix)
CAG2_C_8 <- CAG2_matrix[ , colnames(CAG2_matrix) %in% CAG2_CD8$V1, drop = FALSE]
dim(CAG2_C_8)

CAG3_CD8 <- read.table("GSE134520_RAW/CAG3_CD8_cell_barcodes.txt", header = F)
CAG3_matrix <- read.table("GSE134520_RAW/GSM3954951_processed_CAG3.txt/dt.CAG3.txt", header = TRUE, row.names = 1)
CAG3_matrix<- as.data.frame(CAG3_matrix)
CAG3_C_8 <- CAG3_matrix[ , colnames(CAG3_matrix) %in% CAG3_CD8$V1, drop = FALSE]
dim(CAG3_C_8)

CAG_CD8 <-cbind(CAG1_C_8,CAG2_C_8,CAG3_C_8)
dim(CAG_CD8)
colnames(CAG_CD8) <- paste0("Zhang_", colnames(CAG_CD8))
write.table(CAG_CD8, file = "GSE134520_RAW/CAG_CD8.txt", quote = FALSE, row.names = T, col.names = T)



# IM CD4
IMW1_CD4 <- read.table("GSE134520_RAW/IMW1_CD4_cell_barcodes.txt", header = F)
IMW1_matrix <- read.table("GSE134520_RAW/GSM3954952_processed_IMW1.txt/dt.IMW1.txt", header = TRUE, row.names = 1)
IMW1_matrix<- as.data.frame(IMW1_matrix)
IMW1_C_4 <- IMW1_matrix[ , colnames(IMW1_matrix) %in% IMW1_CD4$V1, drop = FALSE]


IMW2_CD4 <- read.table("GSE134520_RAW/IMW2_CD4_cell_barcodes.txt", header = F)
IMW2_matrix <- read.table("GSE134520_RAW/GSM3954953_processed_IMW2.txt/dt.IMW2.txt", header = TRUE, row.names = 1)
IMW2_matrix<- as.data.frame(IMW2_matrix)
IMW2_C_4 <- IMW2_matrix[ , colnames(IMW2_matrix) %in% IMW2_CD4$V1, drop = FALSE]


IMS1_CD4 <- read.table("GSE134520_RAW/IMs1_CD4_cell_barcodes.txt", header = F)
IMS1_matrix <- read.table("GSE134520_RAW/GSM3954954_processed_IMS1.txt/dt.IMS1.txt", header = TRUE, row.names = 1)
IMS1_matrix<- as.data.frame(IMS1_matrix)
IMS1_C_4 <- IMS1_matrix[ , colnames(IMS1_matrix) %in% IMS1_CD4$V1, drop = FALSE]


IMS2_CD4 <- read.table("GSE134520_RAW/IMS2_CD4_cell_barcodes.txt", header = F)
IMS2_matrix <- read.table("GSE134520_RAW/GSM3954955_processed_IMS2.txt/dt.IMS2.txt", header = TRUE, row.names = 1)
IMS2_matrix<- as.data.frame(IMS2_matrix)
IMS2_C_4 <- IMS2_matrix[ , colnames(IMS2_matrix) %in% IMS2_CD4$V1, drop = FALSE]

IMS3_CD4 <- read.table("GSE134520_RAW/IMS3_CD4_cell_barcodes.txt", header = F)
IMS3_matrix <- read.table("GSE134520_RAW/GSM3954956_processed_IMS3.txt/dt.IMS3.txt", header = TRUE, row.names = 1)
IMS3_matrix<- as.data.frame(IMS3_matrix)
IMS3_C_4 <- IMS3_matrix[ , colnames(IMS3_matrix) %in% IMS3_CD4$V1, drop = FALSE]

IMS4_CD4 <- read.table("GSE134520_RAW/IMS4_CD4_cell_barcodes.txt", header = F)
IMS4_matrix <- read.table("GSE134520_RAW/GSM3954957_processed_IMS4.txt/dt.IMS4.txt", header = TRUE, row.names = 1)
IMS4_matrix<- as.data.frame(IMS4_matrix)
IMS4_C_4 <- IMS4_matrix[ , colnames(IMS4_matrix) %in% IMS4_CD4$V1, drop = FALSE]


IM_CD4<-cbind(IMW1_C_4,IMW2_C_4,IMS1_C_4,IMS2_C_4,IMS3_C_4,IMS4_C_4)
IM_CD4<-cbind(IMW1_C_4,IMS1_C_4,IMS2_C_4,IMS3_C_4,IMS4_C_4)
dim(IM_CD4)
colnames(IM_CD4) <- paste0("Zhang_", colnames(IM_CD4))
write.table(IM_CD4, file = "GSE134520_RAW/IM_CD4.txt", quote = FALSE, row.names = T, col.names = T)


# IM CD8
IMW1_CD8 <- read.table("GSE134520_RAW/IMW1_CD8_cell_barcodes.txt", header = F)
IMW1_matrix <- read.table("GSE134520_RAW/GSM3954952_processed_IMW1.txt/dt.IMW1.txt", header = TRUE, row.names = 1)
IMW1_matrix<- as.data.frame(IMW1_matrix)
IMW1_C_8 <- IMW1_matrix[ , colnames(IMW1_matrix) %in% IMW1_CD8$V1, drop = FALSE]
dim(IMW1_C_8)

IMW2_CD8 <- read.table("GSE134520_RAW/IMW2_CD8_cell_barcodes.txt", header = F)
IMW2_matrix <- read.table("GSE134520_RAW/GSM3954953_processed_IMW2.txt/dt.IMW2.txt", header = TRUE, row.names = 1)
IMW2_matrix<- as.data.frame(IMW2_matrix)
IMW2_C_8 <- IMW2_matrix[ , colnames(IMW2_matrix) %in% IMW2_CD8$V1, drop = FALSE]
dim(IMW2_C_8)

IMS1_CD8 <- read.table("GSE134520_RAW/IMs1_CD8_cell_barcodes.txt", header = F)
IMS1_matrix <- read.table("GSE134520_RAW/GSM3954954_processed_IMS1.txt/dt.IMS1.txt", header = TRUE, row.names = 1)
IMS1_matrix<- as.data.frame(IMS1_matrix)
IMS1_C_8 <- IMS1_matrix[ , colnames(IMS1_matrix) %in% IMS1_CD8$V1, drop = FALSE]
dim(IMS1_C_8)

IMS2_CD8 <- read.table("GSE134520_RAW/IMS2_CD8_cell_barcodes.txt", header = F)
IMS2_matrix <- read.table("GSE134520_RAW/GSM3954955_processed_IMS2.txt/dt.IMS2.txt", header = TRUE, row.names = 1)
IMS2_matrix<- as.data.frame(IMS2_matrix)
IMS2_C_8 <- IMS2_matrix[ , colnames(IMS2_matrix) %in% IMS2_CD8$V1, drop = FALSE]
dim(IMS2_C_8)

IMS3_CD8 <- read.table("GSE134520_RAW/IMS3_CD8_cell_barcodes.txt", header = F)
IMS3_matrix <- read.table("GSE134520_RAW/GSM3954956_processed_IMS3.txt/dt.IMS3.txt", header = TRUE, row.names = 1)
IMS3_matrix<- as.data.frame(IMS3_matrix)
IMS3_C_8 <- IMS3_matrix[ , colnames(IMS3_matrix) %in% IMS3_CD8$V1, drop = FALSE]
dim(IMS3_C_8)

IMS4_CD8 <- read.table("GSE134520_RAW/IMS4_CD8_cell_barcodes.txt", header = F)
IMS4_matrix <- read.table("GSE134520_RAW/GSM3954957_processed_IMS4.txt/dt.IMS4.txt", header = TRUE, row.names = 1)
IMS4_matrix<- as.data.frame(IMS4_matrix)
IMS4_C_8 <- IMS4_matrix[ , colnames(IMS4_matrix) %in% IMS4_CD8$V1, drop = FALSE]
dim(IMS4_C_8)

IM_CD8<-cbind(IMW1_C_8,IMW2_C_8,IMS1_C_8,IMS2_C_8,IMS3_C_8,IMS4_C_8)
dim(IM_CD8)
colnames(IM_CD8) <- paste0("Zhang_", colnames(IM_CD8))
write.table(IM_CD8, file = "GSE134520_RAW/IM_CD8.txt", quote = FALSE, row.names = T, col.names = T)



#EGC
EGC_CD4 <- read.table("GSE134520_RAW/EGC_CD4_cell_barcodes.txt", header = F)
EGC_matrix <- read.table("GSE134520_RAW/GSM3954958_processed_EGC.txt/dt.EGC.txt", header = TRUE, row.names = 1)
EGC_matrix<- as.data.frame(EGC_matrix)
EGC_C_4 <- EGC_matrix[ , colnames(EGC_matrix) %in% EGC_CD4$V1, drop = FALSE]
dim(EGC_C_4)

EGC_CD4<- EGC_C_4
dim(EGC_CD4)
colnames(EGC_CD4) <- paste0("Zhang_", colnames(EGC_CD4))
write.table(EGC_CD4, file = "GSE134520_RAW/EGC_CD4.txt", quote = FALSE, row.names = T, col.names = T)


EGC_CD8 <- read.table("GSE134520_RAW/EGC_CD8_cell_barcodes.txt", header = F)
EGC_matrix <- read.table("GSE134520_RAW/GSM3954958_processed_EGC.txt/dt.EGC.txt", header = TRUE, row.names = 1)
EGC_matrix<- as.data.frame(EGC_matrix)
EGC_C_8 <- EGC_matrix[ , colnames(EGC_matrix) %in% EGC_CD8$V1, drop = FALSE]
dim(EGC_C_8)

EGC_CD8<- EGC_C_8
dim(EGC_CD8)
colnames(EGC_CD8) <- paste0("Zhang_", colnames(EGC_CD8))
write.table(EGC_CD8, file = "GSE134520_RAW/EGC_CD8.txt", quote = FALSE, row.names = T, col.names = T)
