library(Matrix)
library(dplyr)


setwd("/home/Fatemeh/0--ThirdProject/")



##### Load Meta Data #####

load(file='DATA/Gastric/meta.Rdata')

ls()
head(meta)
dim(meta)

#write.csv(meta, file = "/home/fatemeh/Fatemeh/0--ThirdProject/GastricTME-GastricTME/input_data/cell_metadata_with_stage.csv", row.names = TRUE)

str(meta)
unique(meta$tissue)
unique(meta$celltype.big.L1)
length(unique(meta$sample))

tissue_counts <- as.data.frame(table(meta$tissue))
tissue_counts

meta_Sample_Tissue<- unique(meta[,c(1,2)])
meta_Sample_Tissue




##### Load Gastric-GSE234129 Data #####
mtx <- readMM(gzfile("DATA/Gastric/GSE234129/GSE234129_count_matrix.mtx.gz"))

features <- read.delim(gzfile("DATA/Gastric/GSE234129/GSE234129_features.tsv.gz"), header = FALSE)
barcodes <- read.delim(gzfile("DATA/Gastric/GSE234129/GSE234129_barcodes.tsv.gz"), header = FALSE)
meta_sub <- read.delim(gzfile("DATA/Gastric/GSE234129/GSE234129_meta.tsv.gz"), header = TRUE)


dim(mtx)
head(features)
head(barcodes)
head(meta_sub)

##### Prepare meta data for stages' analysis #####
meta_sub_C<-meta_sub
meta_sub_C$sample_clean <- gsub("^MDA_", "", meta_sub_C$sample)
head(meta_sub_C)
dim(meta_sub_C)

meta_sub_Ext<-meta_sub_C

meta_sub_Ext <- meta_sub_Ext %>%
  left_join(meta_Sample_Tissue[, c("sample", "tissue.level.8")], 
            by = c("sample_clean" = "sample"))

dim(meta_sub_Ext)
head(meta_sub_Ext)


### Select meta section data for this section
meta_pt <- meta[grepl("^Pt", meta$sample), ]


meta_sub_Ext_Ty<- meta_sub_Ext
meta_sub_Ext_Ty$celltype.big.L1<- meta_pt$celltype.big.L1

##### Prepare an ordinary matrix from sparse matrix ####
dense_mtx <- as.matrix(mtx)
dim(dense_mtx)


##### Prepare separated feature matrix for each stage #####

### Check the cell orders
# Check if they are exactly in the same order
Barcodes_order<- identical(meta_sub_Ext$cell_barcodes, barcodes$V1)

### Reorder meta_sub_Ext to match the order of barcodes$V1
#if(Barcodes_order==FALSE)
#  meta_sub_Ext <- meta_sub_Ext[match(barcodes$V1, meta_sub_Ext$cell_barcodes), ]




##### Select NAT cells #####

## NAT & CD4T
nat_CD4T_rows <- which(meta_sub_Ext_Ty$tissue.level.8 == "NAT"& meta_sub_Ext_Ty$celltype.big.L1 == "CD4T")
length(nat_CD4T_rows)

nat_CD4T_matrix <- dense_mtx[, nat_CD4T_rows]
dim(nat_CD4T_matrix)

colnames(nat_CD4T_matrix) <- paste0("MDA_", meta_sub_Ext_Ty[nat_CD4T_rows, ]$cell_barcodes)
rownames(nat_CD4T_matrix)<- features$V1
nat_CD4T_matrix[1:5,1:5]

write.table(nat_CD4T_matrix, file = "DATA/Gastric/CellTypeStage/MAD_nat_CD4T_matrix.txt", quote = FALSE, row.names = T, col.names = T)

## NAT & CD8T
nat_CD8T_rows <- which(meta_sub_Ext_Ty$tissue.level.8 == "NAT"& meta_sub_Ext_Ty$celltype.big.L1 == "CD8T")
length(nat_CD8T_rows)

nat_CD8T_matrix <- dense_mtx[, nat_CD8T_rows]
dim(nat_CD8T_matrix)

colnames(nat_CD8T_matrix) <- paste0("MDA_", meta_sub_Ext_Ty[nat_CD8T_rows, ]$cell_barcodes)
rownames(nat_CD8T_matrix)<- features$V1
nat_CD8T_matrix[1:5,1:5]

write.table(nat_CD8T_matrix, file = "DATA/Gastric/CellTypeStage/MAD_nat_CD8T_matrix.txt", quote = FALSE, row.names = T, col.names = T)


## NAT & Myeloid
# nat_Myeloid_rows <- which(meta_sub_Ext_Ty$tissue.level.8 == "NAT"& meta_sub_Ext_Ty$celltype.big.L1 == "Myeloid")
# length(nat_Myeloid_rows)
# 
# nat_Myeloid_matrix <- dense_mtx[, nat_Myeloid_rows]
# dim(nat_Myeloid_matrix)
# 




##### Select GAC_Primary cells #####

## GAC_Primary & CD4T
GAC_Primary_CD4T_rows <- which(meta_sub_Ext_Ty$tissue.level.8 == "GAC_Primary"& meta_sub_Ext_Ty$celltype.big.L1 == "CD4T")
length(GAC_Primary_CD4T_rows)

GAC_Primary_CD4T_matrix <- dense_mtx[, GAC_Primary_CD4T_rows]
dim(GAC_Primary_CD4T_matrix)

colnames(GAC_Primary_CD4T_matrix) <- paste0("MDA_", meta_sub_Ext_Ty[GAC_Primary_CD4T_rows, ]$cell_barcodes)
rownames(GAC_Primary_CD4T_matrix)<- features$V1
GAC_Primary_CD4T_matrix[1:5,1:5]

write.table(GAC_Primary_CD4T_matrix, file = "DATA/Gastric/CellTypeStage/MAD_GAC_Primary_CD4T_matrix.txt", quote = FALSE, row.names = T, col.names = T)

## GAC_Primary & CD8T
GAC_Primary_CD8T_rows <- which(meta_sub_Ext_Ty$tissue.level.8 == "GAC_Primary"& meta_sub_Ext_Ty$celltype.big.L1 == "CD8T")
length(GAC_Primary_CD8T_rows)

GAC_Primary_CD8T_matrix <- dense_mtx[, GAC_Primary_CD8T_rows]
dim(GAC_Primary_CD8T_matrix)

colnames(GAC_Primary_CD8T_matrix) <- paste0("MDA_", meta_sub_Ext_Ty[GAC_Primary_CD8T_rows, ]$cell_barcodes)
rownames(GAC_Primary_CD8T_matrix)<- features$V1
GAC_Primary_CD8T_matrix[1:5,1:5]

write.table(GAC_Primary_CD8T_matrix, file = "DATA/Gastric/CellTypeStage/MAD_GAC_Primary_CD8T_matrix.txt", quote = FALSE, row.names = T, col.names = T)

## GAC_Primary & Myeloid
# GAC_Primary_Myeloid_rows <- which(meta_sub_Ext_Ty$tissue.level.8 == "GAC_Primary"& meta_sub_Ext_Ty$celltype.big.L1 == "Myeloid")
# length(GAC_Primary_Myeloid_rows)
# 
# GAC_Primary_Myeloid_matrix <- dense_mtx[, GAC_Primary_Myeloid_rows]
# dim(GAC_Primary_Myeloid_matrix)





##### Select Metastasis cells #####

## Metastasis & CD4T
Metastasis_CD4T_rows <- which(meta_sub_Ext_Ty$tissue.level.8 == "Metastasis"& meta_sub_Ext_Ty$celltype.big.L1 == "CD4T")
length(Metastasis_CD4T_rows)

Metastasis_CD4T_matrix <- dense_mtx[, Metastasis_CD4T_rows]
dim(Metastasis_CD4T_matrix)

colnames(Metastasis_CD4T_matrix) <- paste0("MDA_", meta_sub_Ext_Ty[Metastasis_CD4T_rows, ]$cell_barcodes)
rownames(Metastasis_CD4T_matrix)<- features$V1
Metastasis_CD4T_matrix[1:5,1:5]

write.table(Metastasis_CD4T_matrix, file = "DATA/Gastric/CellTypeStage/MAD_Metastasis_CD4T_matrix.txt", quote = FALSE, row.names = T, col.names = T)

## Metastasis & CD8T
Metastasis_CD8T_rows <- which(meta_sub_Ext_Ty$tissue.level.8 == "Metastasis"& meta_sub_Ext_Ty$celltype.big.L1 == "CD8T")
length(Metastasis_CD8T_rows)

Metastasis_CD8T_matrix <- dense_mtx[, Metastasis_CD8T_rows]
dim(Metastasis_CD8T_matrix)

colnames(Metastasis_CD8T_matrix) <- paste0("MDA_", meta_sub_Ext_Ty[Metastasis_CD8T_rows, ]$cell_barcodes)
rownames(Metastasis_CD8T_matrix)<- features$V1
Metastasis_CD8T_matrix[1:5,1:5]

write.table(Metastasis_CD8T_matrix, file = "DATA/Gastric/CellTypeStage/MAD_Metastasis_CD8T_matrix.txt", quote = FALSE, row.names = T, col.names = T)

## Metastasis & Myeloid
# Metastasis_Myeloid_rows <- which(meta_sub_Ext_Ty$tissue.level.8 == "Metastasis"& meta_sub_Ext_Ty$celltype.big.L1 == "Myeloid")
# length(Metastasis_Myeloid_rows)
# 
# Metastasis_Myeloid_matrix <- dense_mtx[, Metastasis_Myeloid_rows]
# dim(Metastasis_Myeloid_matrix)



### Count number of patient in these 3 stages
n_unique_patients <- meta_sub_Ext %>%
  filter(tissue.level.8 != "GAC_PBMC") %>%
  pull(patient) %>%
  unique() %>%
  length()

n_unique_patients





##### Load Gastric-Sathe Data #####

##Each sample is named in the following format: patientID_condition_replicateID
##patientID = 4-digit identifier used throughout the manuscript
##condition = n: normal, t: tumor, pbmc: peripheral blood mononuclear cells. Note that 6649 metaplasia sample is labelled as '6649_t1'. 
##replicateID = 1 for samples with single replicate, 2 for 2nd replicate from the same biological sample.
##condition: normal/tumor/metaplasia/pbmc


### CELL_LABEL file ###
cell_labels <- read.csv("DATA/Gastric/Sathe/cell_labels.csv")

unique(cell_labels$orig.ident)
unique(cell_labels$final_celltype)

##### Select NAT cells #####

unique_n <- unique(grep("n[12]$", cell_labels$orig.ident, value = TRUE))
print(unique_n)


##### 5846_n1 #####
n_Label<- cell_labels[cell_labels$orig.ident=="5846_n1",]


mtx_Sathe <- readMM(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/5846_n1/matrix.mtx.gz"))
features_Sathe <- read.delim(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/5846_n1/features.tsv.gz"), header = FALSE)
barcodes_Sathe <- read.delim(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/5846_n1/barcodes.tsv.gz"), header = FALSE)

dim(mtx_Sathe)
dim(features_Sathe)
head(features_Sathe)
dim(barcodes_Sathe)
head(barcodes_Sathe)

## Prepare an ordinary matrix from sparse matrix 
dense_mtx_Sathe <- as.matrix(mtx_Sathe)
dim(mtx_Sathe)


# Remove barcode suffix 
clean_barcode_Sathe <- sub("-\\d+$", "", barcodes_Sathe$V1)

### CD4
S_CD4T<- n_Label[n_Label$final_celltype == "CD4", ]
# Remove barcode suffix 
clean_barcode_S_CD4T <- sub("-\\d+$", "", S_CD4T$cell_barcode)
# Compare without suffix
S_CD4T_rows <- which(clean_barcode_Sathe %in% clean_barcode_S_CD4T)
#S_CD4T_rows <- which(S_CD4T$cell_barcode %in% barcodes_Sathe$V1)
length(S_CD4T_rows)
S_CD4T_matrix <- as.matrix(mtx_Sathe[, S_CD4T_rows])
dim(S_CD4T_matrix)




### CD8
S_CD8T<- n_Label[n_Label$final_celltype == "CD8", ]
# Remove barcode suffix 
clean_barcode_S_CD8T <- sub("-\\d+$", "", S_CD8T$cell_barcode)
# Compare without suffix
S_CD8T_rows <- which(clean_barcode_Sathe %in% clean_barcode_S_CD8T)
#S_CD8T_rows <- which(S_CD8T$cell_barcode %in% barcodes_Sathe$V1)
length(S_CD8T_rows)
S_CD8T_matrix <- as.matrix(mtx_Sathe[, S_CD8T_rows])
dim(S_CD8T_matrix)


colnames(S_CD4T_matrix) <- paste0("Sathe_5846_n1_", clean_barcode_S_CD4T)
rownames(S_CD4T_matrix)<- features_Sathe$V2
S_CD4T_matrix[1:5,1:5]


colnames(S_CD8T_matrix) <- paste0("Sathe_5846_n1_", clean_barcode_S_CD8T)
rownames(S_CD8T_matrix)<- features_Sathe$V2
S_CD8T_matrix[1:5,1:5]

## Save
S_5846_n1_CD4 <- S_CD4T_matrix
S_5846_n1_CD8<- S_CD8T_matrix




##### 5866_n1 #####
n_Label<- cell_labels[cell_labels$orig.ident=="5866_n1",]


mtx_Sathe <- readMM(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/5866_n1/matrix.mtx.gz"))
features_Sathe <- read.delim(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/5866_n1/features.tsv.gz"), header = FALSE)
barcodes_Sathe <- read.delim(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/5866_n1/barcodes.tsv.gz"), header = FALSE)

dim(mtx_Sathe)
dim(features_Sathe)
head(features_Sathe)
dim(barcodes_Sathe)
head(barcodes_Sathe)

## Prepare an ordinary matrix from sparse matrix 
dense_mtx_Sathe <- as.matrix(mtx_Sathe)
dim(mtx_Sathe)



# Remove barcode suffix 
clean_barcode_Sathe <- sub("-\\d+$", "", barcodes_Sathe$V1)

### CD4
S_CD4T<- n_Label[n_Label$final_celltype == "CD4", ]
# Remove barcode suffix 
clean_barcode_S_CD4T <- sub("-\\d+$", "", S_CD4T$cell_barcode)
# Compare without suffix
S_CD4T_rows <- which(clean_barcode_Sathe %in% clean_barcode_S_CD4T)
#S_CD4T_rows <- which(S_CD4T$cell_barcode %in% barcodes_Sathe$V1)
length(S_CD4T_rows)
S_CD4T_matrix <- as.matrix(mtx_Sathe[, S_CD4T_rows])
dim(S_CD4T_matrix)

### CD8
S_CD8T<- n_Label[n_Label$final_celltype == "CD8", ]
# Remove barcode suffix 
clean_barcode_S_CD8T <- sub("-\\d+$", "", S_CD8T$cell_barcode)
# Compare without suffix
S_CD8T_rows <- which(clean_barcode_Sathe %in% clean_barcode_S_CD8T)
#S_CD8T_rows <- which(S_CD8T$cell_barcode %in% barcodes_Sathe$V1)
length(S_CD8T_rows)
S_CD8T_matrix <- as.matrix(mtx_Sathe[, S_CD8T_rows])
dim(S_CD8T_matrix)

colnames(S_CD4T_matrix) <- paste0("Sathe_5866_n1_", clean_barcode_S_CD4T)
rownames(S_CD4T_matrix)<- features_Sathe$V2
S_CD4T_matrix[1:5,1:5]

colnames(S_CD8T_matrix) <- paste0("Sathe_5866_n1_", clean_barcode_S_CD8T)
rownames(S_CD8T_matrix)<- features_Sathe$V2
S_CD8T_matrix[1:5,1:5]

## Save
S_5866_n1_CD4 <- S_CD4T_matrix
S_5866_n1_CD8<- S_CD8T_matrix


##### 5866_n2 #####
n_Label<- cell_labels[cell_labels$orig.ident=="5866_n2",]


mtx_Sathe <- readMM(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/5866_n2/matrix.mtx.gz"))
features_Sathe <- read.delim(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/5866_n2/features.tsv.gz"), header = FALSE)
barcodes_Sathe <- read.delim(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/5866_n2/barcodes.tsv.gz"), header = FALSE)

dim(mtx_Sathe)
dim(features_Sathe)
head(features_Sathe)
dim(barcodes_Sathe)
head(barcodes_Sathe)

## Prepare an ordinary matrix from sparse matrix 
dense_mtx_Sathe <- as.matrix(mtx_Sathe)
dim(mtx_Sathe)



# Remove barcode suffix 
clean_barcode_Sathe <- sub("-\\d+$", "", barcodes_Sathe$V1)

### CD4
S_CD4T<- n_Label[n_Label$final_celltype == "CD4", ]
# Remove barcode suffix 
clean_barcode_S_CD4T <- sub("-\\d+$", "", S_CD4T$cell_barcode)
# Compare without suffix
S_CD4T_rows <- which(clean_barcode_Sathe %in% clean_barcode_S_CD4T)
#S_CD4T_rows <- which(S_CD4T$cell_barcode %in% barcodes_Sathe$V1)
length(S_CD4T_rows)
S_CD4T_matrix <- as.matrix(mtx_Sathe[, S_CD4T_rows])
dim(S_CD4T_matrix)

### CD8
S_CD8T<- n_Label[n_Label$final_celltype == "CD8", ]
# Remove barcode suffix 
clean_barcode_S_CD8T <- sub("-\\d+$", "", S_CD8T$cell_barcode)
# Compare without suffix
S_CD8T_rows <- which(clean_barcode_Sathe %in% clean_barcode_S_CD8T)
#S_CD8T_rows <- which(S_CD8T$cell_barcode %in% barcodes_Sathe$V1)
length(S_CD8T_rows)
S_CD8T_matrix <- as.matrix(mtx_Sathe[, S_CD8T_rows])
dim(S_CD8T_matrix)

colnames(S_CD4T_matrix) <- paste0("Sathe_5866_n2_", clean_barcode_S_CD4T)
rownames(S_CD4T_matrix)<- features_Sathe$V2
S_CD4T_matrix[1:5,1:5]

colnames(S_CD8T_matrix) <- paste0("Sathe_5866_n2_", clean_barcode_S_CD8T)
rownames(S_CD8T_matrix)<- features_Sathe$V2
S_CD8T_matrix[1:5,1:5]
## Save
S_5866_n2_CD4 <- S_CD4T_matrix
S_5866_n2_CD8<- S_CD8T_matrix


##### 5931_n1 #####
n_Label<- cell_labels[cell_labels$orig.ident=="5931_n1",]


mtx_Sathe <- readMM("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/5931_n1/matrix.mtx")
features_Sathe <- read.delim("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/5931_n1/features.tsv", header = FALSE)
barcodes_Sathe <- read.delim("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/5931_n1/barcodes.tsv", header = FALSE)

dim(mtx_Sathe)
dim(features_Sathe)
head(features_Sathe)
dim(barcodes_Sathe)
head(barcodes_Sathe)

## Prepare an ordinary matrix from sparse matrix 
dense_mtx_Sathe <- as.matrix(mtx_Sathe)
dim(mtx_Sathe)



# Remove barcode suffix 
clean_barcode_Sathe <- sub("-\\d+$", "", barcodes_Sathe$V1)

### CD4
S_CD4T<- n_Label[n_Label$final_celltype == "CD4", ]
# Remove barcode suffix 
clean_barcode_S_CD4T <- sub("-\\d+$", "", S_CD4T$cell_barcode)
# Compare without suffix
S_CD4T_rows <- which(clean_barcode_Sathe %in% clean_barcode_S_CD4T)
#S_CD4T_rows <- which(S_CD4T$cell_barcode %in% barcodes_Sathe$V1)
length(S_CD4T_rows)
S_CD4T_matrix <- as.matrix(mtx_Sathe[, S_CD4T_rows])
dim(S_CD4T_matrix)

### CD8
S_CD8T<- n_Label[n_Label$final_celltype == "CD8", ]
# Remove barcode suffix 
clean_barcode_S_CD8T <- sub("-\\d+$", "", S_CD8T$cell_barcode)
# Compare without suffix
S_CD8T_rows <- which(clean_barcode_Sathe %in% clean_barcode_S_CD8T)
#S_CD8T_rows <- which(S_CD8T$cell_barcode %in% barcodes_Sathe$V1)
length(S_CD8T_rows)
S_CD8T_matrix <- as.matrix(mtx_Sathe[, S_CD8T_rows])
dim(S_CD8T_matrix)


colnames(S_CD4T_matrix) <- paste0("Sathe_5931_n1_", clean_barcode_S_CD4T)
rownames(S_CD4T_matrix)<- features_Sathe$V2
S_CD4T_matrix[1:5,1:5]

colnames(S_CD8T_matrix) <- paste0("Sathe_5931_n1_", clean_barcode_S_CD8T)
rownames(S_CD8T_matrix)<- features_Sathe$V2
S_CD8T_matrix[1:5,1:5]

## Save
S_5931_n1_CD4 <- S_CD4T_matrix
S_5931_n1_CD8<- S_CD8T_matrix


##### 5931_n2 #####
n_Label<- cell_labels[cell_labels$orig.ident=="5931_n2",]


mtx_Sathe <- readMM(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/5931_n2/matrix.mtx.gz"))
features_Sathe <- read.delim(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/5931_n2/features.tsv.gz"), header = FALSE)
barcodes_Sathe <- read.delim(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/5931_n2/barcodes.tsv.gz"), header = FALSE)

dim(mtx_Sathe)
dim(features_Sathe)
head(features_Sathe)
dim(barcodes_Sathe)
head(barcodes_Sathe)

## Prepare an ordinary matrix from sparse matrix 
dense_mtx_Sathe <- as.matrix(mtx_Sathe)
dim(mtx_Sathe)



# Remove barcode suffix 
clean_barcode_Sathe <- sub("-\\d+$", "", barcodes_Sathe$V1)

### CD4
S_CD4T<- n_Label[n_Label$final_celltype == "CD4", ]
# Remove barcode suffix 
clean_barcode_S_CD4T <- sub("-\\d+$", "", S_CD4T$cell_barcode)
# Compare without suffix
S_CD4T_rows <- which(clean_barcode_Sathe %in% clean_barcode_S_CD4T)
#S_CD4T_rows <- which(S_CD4T$cell_barcode %in% barcodes_Sathe$V1)
length(S_CD4T_rows)
S_CD4T_matrix <- as.matrix(mtx_Sathe[, S_CD4T_rows])
dim(S_CD4T_matrix)

### CD8
S_CD8T<- n_Label[n_Label$final_celltype == "CD8", ]
# Remove barcode suffix 
clean_barcode_S_CD8T <- sub("-\\d+$", "", S_CD8T$cell_barcode)
# Compare without suffix
S_CD8T_rows <- which(clean_barcode_Sathe %in% clean_barcode_S_CD8T)
#S_CD8T_rows <- which(S_CD8T$cell_barcode %in% barcodes_Sathe$V1)
length(S_CD8T_rows)
S_CD8T_matrix <- as.matrix(mtx_Sathe[, S_CD8T_rows])
dim(S_CD8T_matrix)

colnames(S_CD4T_matrix) <- paste0("Sathe_5931_n2_", clean_barcode_S_CD4T)
rownames(S_CD4T_matrix)<- features_Sathe$V2
S_CD4T_matrix[1:5,1:5]

colnames(S_CD8T_matrix) <- paste0("Sathe_5931_n2_", clean_barcode_S_CD8T)
rownames(S_CD8T_matrix)<- features_Sathe$V2
S_CD8T_matrix[1:5,1:5]

## Save
S_5931_n2_CD4 <- S_CD4T_matrix
S_5931_n2_CD8<- S_CD8T_matrix


##### 6709_n1 #####
n_Label<- cell_labels[cell_labels$orig.ident=="6709_n1",]


mtx_Sathe <- readMM(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/6709_n1/matrix.mtx.gz"))
features_Sathe <- read.delim(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/6709_n1/features.tsv.gz"), header = FALSE)
barcodes_Sathe <- read.delim(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/6709_n1/barcodes.tsv.gz"), header = FALSE)

dim(mtx_Sathe)
dim(features_Sathe)
head(features_Sathe)
dim(barcodes_Sathe)
head(barcodes_Sathe)

## Prepare an ordinary matrix from sparse matrix 
dense_mtx_Sathe <- as.matrix(mtx_Sathe)
dim(mtx_Sathe)



# Remove barcode suffix 
clean_barcode_Sathe <- sub("-\\d+$", "", barcodes_Sathe$V1)

### CD4
S_CD4T<- n_Label[n_Label$final_celltype == "CD4", ]
# Remove barcode suffix 
clean_barcode_S_CD4T <- sub("-\\d+$", "", S_CD4T$cell_barcode)
# Compare without suffix
S_CD4T_rows <- which(clean_barcode_Sathe %in% clean_barcode_S_CD4T)
#S_CD4T_rows <- which(S_CD4T$cell_barcode %in% barcodes_Sathe$V1)
length(S_CD4T_rows)
S_CD4T_matrix <- as.matrix(mtx_Sathe[, S_CD4T_rows])
dim(S_CD4T_matrix)

### CD8
S_CD8T<- n_Label[n_Label$final_celltype == "CD8", ]
# Remove barcode suffix 
clean_barcode_S_CD8T <- sub("-\\d+$", "", S_CD8T$cell_barcode)
# Compare without suffix
S_CD8T_rows <- which(clean_barcode_Sathe %in% clean_barcode_S_CD8T)
#S_CD8T_rows <- which(S_CD8T$cell_barcode %in% barcodes_Sathe$V1)
length(S_CD8T_rows)
S_CD8T_matrix <- as.matrix(mtx_Sathe[, S_CD8T_rows])
dim(S_CD8T_matrix)

colnames(S_CD4T_matrix) <- paste0("Sathe_6709_n1_", clean_barcode_S_CD4T)
rownames(S_CD4T_matrix)<- features_Sathe$V2
S_CD4T_matrix[1:5,1:5]

colnames(S_CD8T_matrix) <- paste0("Sathe_6709_n1_", clean_barcode_S_CD8T)
rownames(S_CD8T_matrix)<- features_Sathe$V2
S_CD8T_matrix[1:5,1:5]

## Save
S_6709_n1_CD4 <- S_CD4T_matrix
S_6709_n1_CD8<- S_CD8T_matrix


##### 6592_n1 #####
n_Label<- cell_labels[cell_labels$orig.ident=="6592_n1",]


mtx_Sathe <- readMM(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/6592_n1/matrix.mtx.gz"))
features_Sathe <- read.delim(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/6592_n1/features.tsv.gz"), header = FALSE)
barcodes_Sathe <- read.delim(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/6592_n1/barcodes.tsv.gz"), header = FALSE)

dim(mtx_Sathe)
dim(features_Sathe)
head(features_Sathe)
dim(barcodes_Sathe)
head(barcodes_Sathe)

## Prepare an ordinary matrix from sparse matrix 
dense_mtx_Sathe <- as.matrix(mtx_Sathe)
dim(mtx_Sathe)



# Remove barcode suffix 
clean_barcode_Sathe <- sub("-\\d+$", "", barcodes_Sathe$V1)

### CD4
S_CD4T<- n_Label[n_Label$final_celltype == "CD4", ]
# Remove barcode suffix 
clean_barcode_S_CD4T <- sub("-\\d+$", "", S_CD4T$cell_barcode)
# Compare without suffix
S_CD4T_rows <- which(clean_barcode_Sathe %in% clean_barcode_S_CD4T)
#S_CD4T_rows <- which(S_CD4T$cell_barcode %in% barcodes_Sathe$V1)
length(S_CD4T_rows)
S_CD4T_matrix <- as.matrix(mtx_Sathe[, S_CD4T_rows])
dim(S_CD4T_matrix)

### CD8
S_CD8T<- n_Label[n_Label$final_celltype == "CD8", ]
# Remove barcode suffix 
clean_barcode_S_CD8T <- sub("-\\d+$", "", S_CD8T$cell_barcode)
# Compare without suffix
S_CD8T_rows <- which(clean_barcode_Sathe %in% clean_barcode_S_CD8T)
#S_CD8T_rows <- which(S_CD8T$cell_barcode %in% barcodes_Sathe$V1)
length(S_CD8T_rows)
S_CD8T_matrix <- as.matrix(mtx_Sathe[, S_CD8T_rows])
dim(S_CD8T_matrix)

colnames(S_CD4T_matrix) <- paste0("Sathe_6592_n1_", clean_barcode_S_CD4T)
rownames(S_CD4T_matrix)<- features_Sathe$V2
S_CD4T_matrix[1:5,1:5]

colnames(S_CD8T_matrix) <- paste0("Sathe_6592_n1_", clean_barcode_S_CD8T)
rownames(S_CD8T_matrix)<- features_Sathe$V2
S_CD8T_matrix[1:5,1:5]

## Save
S_6592_n1_CD4 <- S_CD4T_matrix
S_6592_n1_CD8<- S_CD8T_matrix


##### 6342_n1 #####
n_Label<- cell_labels[cell_labels$orig.ident=="6342_n1",]


mtx_Sathe <- readMM(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/6342_n1/matrix.mtx.gz"))
features_Sathe <- read.delim(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/6342_n1/features.tsv.gz"), header = FALSE)
barcodes_Sathe <- read.delim(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/6342_n1/barcodes.tsv.gz"), header = FALSE)

dim(mtx_Sathe)
dim(features_Sathe)
head(features_Sathe)
dim(barcodes_Sathe)
head(barcodes_Sathe)

## Prepare an ordinary matrix from sparse matrix 
dense_mtx_Sathe <- as.matrix(mtx_Sathe)
dim(mtx_Sathe)



# Remove barcode suffix 
clean_barcode_Sathe <- sub("-\\d+$", "", barcodes_Sathe$V1)

### CD4
S_CD4T<- n_Label[n_Label$final_celltype == "CD4", ]
# Remove barcode suffix 
clean_barcode_S_CD4T <- sub("-\\d+$", "", S_CD4T$cell_barcode)
# Compare without suffix
S_CD4T_rows <- which(clean_barcode_Sathe %in% clean_barcode_S_CD4T)
#S_CD4T_rows <- which(S_CD4T$cell_barcode %in% barcodes_Sathe$V1)
length(S_CD4T_rows)
S_CD4T_matrix <- as.matrix(mtx_Sathe[, S_CD4T_rows])
dim(S_CD4T_matrix)

### CD8
S_CD8T<- n_Label[n_Label$final_celltype == "CD8", ]
# Remove barcode suffix 
clean_barcode_S_CD8T <- sub("-\\d+$", "", S_CD8T$cell_barcode)
# Compare without suffix
S_CD8T_rows <- which(clean_barcode_Sathe %in% clean_barcode_S_CD8T)
#S_CD8T_rows <- which(S_CD8T$cell_barcode %in% barcodes_Sathe$V1)
length(S_CD8T_rows)
S_CD8T_matrix <- as.matrix(mtx_Sathe[, S_CD8T_rows])
dim(S_CD8T_matrix)

colnames(S_CD4T_matrix) <- paste0("Sathe_6342_n1_", clean_barcode_S_CD4T)
rownames(S_CD4T_matrix)<- features_Sathe$V2
S_CD4T_matrix[1:5,1:5]

colnames(S_CD8T_matrix) <- paste0("Sathe_6342_n1_", clean_barcode_S_CD8T)
rownames(S_CD8T_matrix)<- features_Sathe$V2
S_CD8T_matrix[1:5,1:5]

## Save
S_6342_n1_CD4 <- S_CD4T_matrix
S_6342_n1_CD8<- S_CD8T_matrix

##### 6649_n1 #####
n_Label<- cell_labels[cell_labels$orig.ident=="6649_n1",]


mtx_Sathe <- readMM(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/6649_n1/matrix.mtx.gz"))
features_Sathe <- read.delim(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/6649_n1/features.tsv.gz"), header = FALSE)
barcodes_Sathe <- read.delim(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/6649_n1/barcodes.tsv.gz"), header = FALSE)

dim(mtx_Sathe)
dim(features_Sathe)
head(features_Sathe)
dim(barcodes_Sathe)
head(barcodes_Sathe)

## Prepare an ordinary matrix from sparse matrix 
dense_mtx_Sathe <- as.matrix(mtx_Sathe)
dim(mtx_Sathe)



# Remove barcode suffix 
clean_barcode_Sathe <- sub("-\\d+$", "", barcodes_Sathe$V1)

### CD4
S_CD4T<- n_Label[n_Label$final_celltype == "CD4", ]
# Remove barcode suffix 
clean_barcode_S_CD4T <- sub("-\\d+$", "", S_CD4T$cell_barcode)
# Compare without suffix
S_CD4T_rows <- which(clean_barcode_Sathe %in% clean_barcode_S_CD4T)
#S_CD4T_rows <- which(S_CD4T$cell_barcode %in% barcodes_Sathe$V1)
length(S_CD4T_rows)
S_CD4T_matrix <- as.matrix(mtx_Sathe[, S_CD4T_rows])
dim(S_CD4T_matrix)

### CD8
S_CD8T<- n_Label[n_Label$final_celltype == "CD8", ]
# Remove barcode suffix 
clean_barcode_S_CD8T <- sub("-\\d+$", "", S_CD8T$cell_barcode)
# Compare without suffix
S_CD8T_rows <- which(clean_barcode_Sathe %in% clean_barcode_S_CD8T)
#S_CD8T_rows <- which(S_CD8T$cell_barcode %in% barcodes_Sathe$V1)
length(S_CD8T_rows)
S_CD8T_matrix <- as.matrix(mtx_Sathe[, S_CD8T_rows])
dim(S_CD8T_matrix)

colnames(S_CD4T_matrix) <- paste0("Sathe_6649_n1_", clean_barcode_S_CD4T)
rownames(S_CD4T_matrix)<- features_Sathe$V2
S_CD4T_matrix[1:5,1:5]

colnames(S_CD8T_matrix) <- paste0("Sathe_6649_n1_", clean_barcode_S_CD8T)
rownames(S_CD8T_matrix)<- features_Sathe$V2
S_CD8T_matrix[1:5,1:5]

## Save
S_6649_n1_CD4 <- S_CD4T_matrix
S_6649_n1_CD8<- S_CD8T_matrix


##### 6207_n1 #####
n_Label<- cell_labels[cell_labels$orig.ident=="6207_n1",]


mtx_Sathe <- readMM(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/6207_n1/matrix.mtx.gz"))
features_Sathe <- read.delim(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/6207_n1/features.tsv.gz"), header = FALSE)
barcodes_Sathe <- read.delim(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/6207_n1/barcodes.tsv.gz"), header = FALSE)

dim(mtx_Sathe)
dim(features_Sathe)
head(features_Sathe)
dim(barcodes_Sathe)
head(barcodes_Sathe)

## Prepare an ordinary matrix from sparse matrix 
dense_mtx_Sathe <- as.matrix(mtx_Sathe)
dim(mtx_Sathe)



# Remove barcode suffix 
clean_barcode_Sathe <- sub("-\\d+$", "", barcodes_Sathe$V1)

### CD4
S_CD4T<- n_Label[n_Label$final_celltype == "CD4", ]
# Remove barcode suffix 
clean_barcode_S_CD4T <- sub("-\\d+$", "", S_CD4T$cell_barcode)
# Compare without suffix
S_CD4T_rows <- which(clean_barcode_Sathe %in% clean_barcode_S_CD4T)
#S_CD4T_rows <- which(S_CD4T$cell_barcode %in% barcodes_Sathe$V1)
length(S_CD4T_rows)
S_CD4T_matrix <- as.matrix(mtx_Sathe[, S_CD4T_rows])
dim(S_CD4T_matrix)

### CD8
S_CD8T<- n_Label[n_Label$final_celltype == "CD8", ]
# Remove barcode suffix 
clean_barcode_S_CD8T <- sub("-\\d+$", "", S_CD8T$cell_barcode)
# Compare without suffix
S_CD8T_rows <- which(clean_barcode_Sathe %in% clean_barcode_S_CD8T)
#S_CD8T_rows <- which(S_CD8T$cell_barcode %in% barcodes_Sathe$V1)
length(S_CD8T_rows)
S_CD8T_matrix <- as.matrix(mtx_Sathe[, S_CD8T_rows])
dim(S_CD8T_matrix)

colnames(S_CD4T_matrix) <- paste0("Sathe_6207_n1_", clean_barcode_S_CD4T)
rownames(S_CD4T_matrix)<- features_Sathe$V2
S_CD4T_matrix[1:5,1:5]

colnames(S_CD8T_matrix) <- paste0("Sathe_6207_n1_", clean_barcode_S_CD8T)
rownames(S_CD8T_matrix)<- features_Sathe$V2
S_CD8T_matrix[1:5,1:5]

## Save
S_6207_n1_CD4 <- S_CD4T_matrix
S_6207_n1_CD8<- S_CD8T_matrix




#

##### MAKE NORMAL MATRIX #####

S_NAT_CD4<- cbind(S_5846_n1_CD4, S_5866_n1_CD4, S_5866_n2_CD4, S_5931_n1_CD4, S_5931_n2_CD4, S_6709_n1_CD4, S_6592_n1_CD4, S_6342_n1_CD4, S_6649_n1_CD4, S_6207_n1_CD4)
rownames(S_NAT_CD4)<- features_Sathe$V2
# Assign gene symbols directly (with duplicates)
gene_names <- rownames(S_NAT_CD4)
# Collapse rows with the same gene name by rowMeans
S_NAT_CD4_dedup <- rowsum(S_NAT_CD4, group = gene_names) / as.vector(table(gene_names))


S_NAT_CD8<- cbind(S_5846_n1_CD8, S_5866_n1_CD8, S_5866_n2_CD8, S_5931_n1_CD8, S_5931_n2_CD8, S_6709_n1_CD8, S_6592_n1_CD8, S_6342_n1_CD8, S_6649_n1_CD8, S_6207_n1_CD8 )
rownames(S_NAT_CD8)<- features_Sathe$V2
# Assign gene symbols directly (with duplicates)
gene_names <- rownames(S_NAT_CD8)
# Collapse rows with the same gene name by rowMeans
S_NAT_CD8_dedup <- rowsum(S_NAT_CD8, group = gene_names) / as.vector(table(gene_names))



write.table(S_NAT_CD4_dedup, file = "DATA/Gastric/CellTypeStage/Sathe_NAT_CD4T_matrix.txt", quote = FALSE, row.names = T, col.names = T)

write.table(S_NAT_CD8_dedup, file = "DATA/Gastric/CellTypeStage/Sathe_NAT_CD8T_matrix.txt", quote = FALSE, row.names = T, col.names = T)





#


##### Select IM cells #####
##### 6649_t1 #####
n_Label<- cell_labels[cell_labels$orig.ident=="6649_t1",]


mtx_Sathe <- readMM("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/6649_t1/matrix.mtx")
features_Sathe <- read.delim("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/6649_t1/features.tsv", header = FALSE)
barcodes_Sathe <- read.delim("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/6649_t1/barcodes.tsv", header = FALSE)

dim(mtx_Sathe)
dim(features_Sathe)
head(features_Sathe)
dim(barcodes_Sathe)
head(barcodes_Sathe)

## Prepare an ordinary matrix from sparse matrix 
dense_mtx_Sathe <- as.matrix(mtx_Sathe)
dim(mtx_Sathe)



# Remove barcode suffix 
clean_barcode_Sathe <- sub("-\\d+$", "", barcodes_Sathe$V1)

### CD4
S_CD4T<- n_Label[n_Label$final_celltype == "CD4", ]
# Remove barcode suffix 
clean_barcode_S_CD4T <- sub("-\\d+$", "", S_CD4T$cell_barcode)
# Compare without suffix
S_CD4T_rows <- which(clean_barcode_Sathe %in% clean_barcode_S_CD4T)
#S_CD4T_rows <- which(S_CD4T$cell_barcode %in% barcodes_Sathe$V1)
length(S_CD4T_rows)
S_CD4T_matrix <- as.matrix(mtx_Sathe[, S_CD4T_rows])
dim(S_CD4T_matrix)

### CD8
S_CD8T<- n_Label[n_Label$final_celltype == "CD8", ]
# Remove barcode suffix 
clean_barcode_S_CD8T <- sub("-\\d+$", "", S_CD8T$cell_barcode)
# Compare without suffix
S_CD8T_rows <- which(clean_barcode_Sathe %in% clean_barcode_S_CD8T)
#S_CD8T_rows <- which(S_CD8T$cell_barcode %in% barcodes_Sathe$V1)
length(S_CD8T_rows)
S_CD8T_matrix <- as.matrix(mtx_Sathe[, S_CD8T_rows])
dim(S_CD8T_matrix)


colnames(S_CD4T_matrix) <- paste0("Sathe_6649_t1_", clean_barcode_S_CD4T)
rownames(S_CD4T_matrix)<- features_Sathe$V2
S_CD4T_matrix[1:2,1:2]

colnames(S_CD8T_matrix) <- paste0("Sathe_6649_t1_", clean_barcode_S_CD8T)
rownames(S_CD8T_matrix)<- features_Sathe$V2
S_CD8T_matrix[1:5,1:5]

## Save
S_6649_t1_CD4 <- S_CD4T_matrix
S_6649_t1_CD8<- S_CD8T_matrix

rownames(S_6649_t1_CD4)<- features_Sathe$V2
# Assign gene symbols directly (with duplicates)
gene_names <- rownames(S_6649_t1_CD4)
# Collapse rows with the same gene name by rowMeans
S_6649_t1_CD4_dedup <- rowsum(S_6649_t1_CD4, group = gene_names) / as.vector(table(gene_names))

rownames(S_6649_t1_CD8)<- features_Sathe$V2
# Assign gene symbols directly (with duplicates)
gene_names <- rownames(S_6649_t1_CD8)
# Collapse rows with the same gene name by rowMeans
S_6649_t1_CD8_dedup <- rowsum(S_6649_t1_CD8, group = gene_names) / as.vector(table(gene_names))


write.table(S_6649_t1_CD4_dedup, file = "DATA/Gastric/CellTypeStage/Sathe_IM_CD4T_matrix.txt", quote = FALSE, row.names = T, col.names = T)
write.table(S_6649_t1_CD8_dedup, file = "DATA/Gastric/CellTypeStage/Sathe_IM_CD8T_matrix.txt", quote = FALSE, row.names = T, col.names = T)





##### Select GAC_Primary cells #####

unique_t <- unique(grep("t[12]$", cell_labels$orig.ident, value = TRUE))
print(unique_t)

## "5866_t1" "5931_t1" "5866_t2" "6342_t1" "5931_t2" "6592_t1" "5846_t1" "6207_t1" "6709_t1"


##### 5866_t1 #####
n_Label<- cell_labels[cell_labels$orig.ident=="5866_t1",]


mtx_Sathe <- readMM(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/5866_t1/matrix.mtx.gz"))
features_Sathe <- read.delim(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/5866_t1/features.tsv.gz"), header = FALSE)
barcodes_Sathe <- read.delim(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/5866_t1/barcodes.tsv.gz"), header = FALSE)

dim(mtx_Sathe)
dim(features_Sathe)
head(features_Sathe)
dim(barcodes_Sathe)
head(barcodes_Sathe)

## Prepare an ordinary matrix from sparse matrix 
dense_mtx_Sathe <- as.matrix(mtx_Sathe)
dim(mtx_Sathe)



# Remove barcode suffix 
clean_barcode_Sathe <- sub("-\\d+$", "", barcodes_Sathe$V1)

### CD4
S_CD4T<- n_Label[n_Label$final_celltype == "CD4", ]
# Remove barcode suffix 
clean_barcode_S_CD4T <- sub("-\\d+$", "", S_CD4T$cell_barcode)
# Compare without suffix
S_CD4T_rows <- which(clean_barcode_Sathe %in% clean_barcode_S_CD4T)
#S_CD4T_rows <- which(S_CD4T$cell_barcode %in% barcodes_Sathe$V1)
length(S_CD4T_rows)
S_CD4T_matrix <- as.matrix(mtx_Sathe[, S_CD4T_rows])
dim(S_CD4T_matrix)

### CD8
S_CD8T<- n_Label[n_Label$final_celltype == "CD8", ]
# Remove barcode suffix 
clean_barcode_S_CD8T <- sub("-\\d+$", "", S_CD8T$cell_barcode)
# Compare without suffix
S_CD8T_rows <- which(clean_barcode_Sathe %in% clean_barcode_S_CD8T)
#S_CD8T_rows <- which(S_CD8T$cell_barcode %in% barcodes_Sathe$V1)
length(S_CD8T_rows)
S_CD8T_matrix <- as.matrix(mtx_Sathe[, S_CD8T_rows])
dim(S_CD8T_matrix)


colnames(S_CD4T_matrix) <- paste0("Sathe_5866_t1_", clean_barcode_S_CD4T)
rownames(S_CD4T_matrix)<- features_Sathe$V2
S_CD4T_matrix[1:5,1:5]

colnames(S_CD8T_matrix) <- paste0("Sathe_5866_t1_", clean_barcode_S_CD8T)
rownames(S_CD8T_matrix)<- features_Sathe$V2
S_CD8T_matrix[1:5,1:5]


## Save
S_5866_t1_CD4 <- S_CD4T_matrix
S_5866_t1_CD8<- S_CD8T_matrix




##### 5931_t1 #####
n_Label<- cell_labels[cell_labels$orig.ident=="5931_t1",]


mtx_Sathe <- readMM(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/5931_t1/matrix.mtx.gz"))
features_Sathe <- read.delim(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/5931_t1/features.tsv.gz"), header = FALSE)
barcodes_Sathe <- read.delim(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/5931_t1/barcodes.tsv.gz"), header = FALSE)

dim(mtx_Sathe)
dim(features_Sathe)
head(features_Sathe)
dim(barcodes_Sathe)
head(barcodes_Sathe)

## Prepare an ordinary matrix from sparse matrix 
dense_mtx_Sathe <- as.matrix(mtx_Sathe)
dim(mtx_Sathe)



# Remove barcode suffix 
clean_barcode_Sathe <- sub("-\\d+$", "", barcodes_Sathe$V1)

### CD4
S_CD4T<- n_Label[n_Label$final_celltype == "CD4", ]
# Remove barcode suffix 
clean_barcode_S_CD4T <- sub("-\\d+$", "", S_CD4T$cell_barcode)
# Compare without suffix
S_CD4T_rows <- which(clean_barcode_Sathe %in% clean_barcode_S_CD4T)
#S_CD4T_rows <- which(S_CD4T$cell_barcode %in% barcodes_Sathe$V1)
length(S_CD4T_rows)
S_CD4T_matrix <- as.matrix(mtx_Sathe[, S_CD4T_rows])
dim(S_CD4T_matrix)

### CD8
S_CD8T<- n_Label[n_Label$final_celltype == "CD8", ]
# Remove barcode suffix 
clean_barcode_S_CD8T <- sub("-\\d+$", "", S_CD8T$cell_barcode)
# Compare without suffix
S_CD8T_rows <- which(clean_barcode_Sathe %in% clean_barcode_S_CD8T)
#S_CD8T_rows <- which(S_CD8T$cell_barcode %in% barcodes_Sathe$V1)
length(S_CD8T_rows)
S_CD8T_matrix <- as.matrix(mtx_Sathe[, S_CD8T_rows])
dim(S_CD8T_matrix)


colnames(S_CD4T_matrix) <- paste0("Sathe_5931_t1_", clean_barcode_S_CD4T)
rownames(S_CD4T_matrix)<- features_Sathe$V2
S_CD4T_matrix[1:5,1:5]

colnames(S_CD8T_matrix) <- paste0("Sathe_5931_t1_", clean_barcode_S_CD8T)
rownames(S_CD8T_matrix)<- features_Sathe$V2
S_CD8T_matrix[1:5,1:5]

## Save
S_5931_t1_CD4 <- S_CD4T_matrix
S_5931_t1_CD8<- S_CD8T_matrix




##### 5866_t2 #####
n_Label<- cell_labels[cell_labels$orig.ident=="5866_t2",]


mtx_Sathe <- readMM(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/5866_t2/matrix.mtx.gz"))
features_Sathe <- read.delim(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/5866_t2/features.tsv.gz"), header = FALSE)
barcodes_Sathe <- read.delim(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/5866_t2/barcodes.tsv.gz"), header = FALSE)

dim(mtx_Sathe)
dim(features_Sathe)
head(features_Sathe)
dim(barcodes_Sathe)
head(barcodes_Sathe)

## Prepare an ordinary matrix from sparse matrix 
dense_mtx_Sathe <- as.matrix(mtx_Sathe)
dim(mtx_Sathe)



# Remove barcode suffix 
clean_barcode_Sathe <- sub("-\\d+$", "", barcodes_Sathe$V1)

### CD4
S_CD4T<- n_Label[n_Label$final_celltype == "CD4", ]
# Remove barcode suffix 
clean_barcode_S_CD4T <- sub("-\\d+$", "", S_CD4T$cell_barcode)
# Compare without suffix
S_CD4T_rows <- which(clean_barcode_Sathe %in% clean_barcode_S_CD4T)
#S_CD4T_rows <- which(S_CD4T$cell_barcode %in% barcodes_Sathe$V1)
length(S_CD4T_rows)
S_CD4T_matrix <- as.matrix(mtx_Sathe[, S_CD4T_rows])
dim(S_CD4T_matrix)

### CD8
S_CD8T<- n_Label[n_Label$final_celltype == "CD8", ]
# Remove barcode suffix 
clean_barcode_S_CD8T <- sub("-\\d+$", "", S_CD8T$cell_barcode)
# Compare without suffix
S_CD8T_rows <- which(clean_barcode_Sathe %in% clean_barcode_S_CD8T)
#S_CD8T_rows <- which(S_CD8T$cell_barcode %in% barcodes_Sathe$V1)
length(S_CD8T_rows)
S_CD8T_matrix <- as.matrix(mtx_Sathe[, S_CD8T_rows])
dim(S_CD8T_matrix)

colnames(S_CD4T_matrix) <- paste0("Sathe_5866_t2_", clean_barcode_S_CD4T)
rownames(S_CD4T_matrix)<- features_Sathe$V2
S_CD4T_matrix[1:5,1:5]

colnames(S_CD8T_matrix) <- paste0("Sathe_5866_t2_", clean_barcode_S_CD8T)
rownames(S_CD8T_matrix)<- features_Sathe$V2
S_CD8T_matrix[1:5,1:5]

## Save
S_5866_t2_CD4 <- S_CD4T_matrix
S_5866_t2_CD8<- S_CD8T_matrix




##### 6342_t1 #####
n_Label<- cell_labels[cell_labels$orig.ident=="6342_t1",]


mtx_Sathe <- readMM("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/6342_t1/matrix.mtx")
features_Sathe <- read.delim("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/6342_t1/features.tsv", header = FALSE)
barcodes_Sathe <- read.delim("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/6342_t1/barcodes.tsv", header = FALSE)

dim(mtx_Sathe)
dim(features_Sathe)
head(features_Sathe)
dim(barcodes_Sathe)
head(barcodes_Sathe)

## Prepare an ordinary matrix from sparse matrix 
dense_mtx_Sathe <- as.matrix(mtx_Sathe)
dim(mtx_Sathe)



# Remove barcode suffix 
clean_barcode_Sathe <- sub("-\\d+$", "", barcodes_Sathe$V1)

### CD4
S_CD4T<- n_Label[n_Label$final_celltype == "CD4", ]
# Remove barcode suffix 
clean_barcode_S_CD4T <- sub("-\\d+$", "", S_CD4T$cell_barcode)
# Compare without suffix
S_CD4T_rows <- which(clean_barcode_Sathe %in% clean_barcode_S_CD4T)
#S_CD4T_rows <- which(S_CD4T$cell_barcode %in% barcodes_Sathe$V1)
length(S_CD4T_rows)
S_CD4T_matrix <- as.matrix(mtx_Sathe[, S_CD4T_rows])
dim(S_CD4T_matrix)

### CD8
S_CD8T<- n_Label[n_Label$final_celltype == "CD8", ]
# Remove barcode suffix 
clean_barcode_S_CD8T <- sub("-\\d+$", "", S_CD8T$cell_barcode)
# Compare without suffix
S_CD8T_rows <- which(clean_barcode_Sathe %in% clean_barcode_S_CD8T)
#S_CD8T_rows <- which(S_CD8T$cell_barcode %in% barcodes_Sathe$V1)
length(S_CD8T_rows)
S_CD8T_matrix <- as.matrix(mtx_Sathe[, S_CD8T_rows])
dim(S_CD8T_matrix)

colnames(S_CD4T_matrix) <- paste0("Sathe_6342_t1_", clean_barcode_S_CD4T)
rownames(S_CD4T_matrix)<- features_Sathe$V2
S_CD4T_matrix[1:5,1:5]

colnames(S_CD8T_matrix) <- paste0("Sathe_6342_t1_", clean_barcode_S_CD8T)
rownames(S_CD8T_matrix)<- features_Sathe$V2
S_CD8T_matrix[1:5,1:5]

## Save
S_6342_t1_CD4 <- S_CD4T_matrix
S_6342_t1_CD8<- S_CD8T_matrix




##### 5931_t2 #####
n_Label<- cell_labels[cell_labels$orig.ident=="5931_t2",]


mtx_Sathe <- readMM(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/5931_t2/matrix.mtx.gz"))
features_Sathe <- read.delim(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/5931_t2/features.tsv.gz"), header = FALSE)
barcodes_Sathe <- read.delim(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/5931_t2/barcodes.tsv.gz"), header = FALSE)

dim(mtx_Sathe)
dim(features_Sathe)
head(features_Sathe)
dim(barcodes_Sathe)
head(barcodes_Sathe)

## Prepare an ordinary matrix from sparse matrix 
dense_mtx_Sathe <- as.matrix(mtx_Sathe)
dim(mtx_Sathe)



# Remove barcode suffix 
clean_barcode_Sathe <- sub("-\\d+$", "", barcodes_Sathe$V1)

### CD4
S_CD4T<- n_Label[n_Label$final_celltype == "CD4", ]
# Remove barcode suffix 
clean_barcode_S_CD4T <- sub("-\\d+$", "", S_CD4T$cell_barcode)
# Compare without suffix
S_CD4T_rows <- which(clean_barcode_Sathe %in% clean_barcode_S_CD4T)
#S_CD4T_rows <- which(S_CD4T$cell_barcode %in% barcodes_Sathe$V1)
length(S_CD4T_rows)
S_CD4T_matrix <- as.matrix(mtx_Sathe[, S_CD4T_rows])
dim(S_CD4T_matrix)

### CD8
S_CD8T<- n_Label[n_Label$final_celltype == "CD8", ]
# Remove barcode suffix 
clean_barcode_S_CD8T <- sub("-\\d+$", "", S_CD8T$cell_barcode)
# Compare without suffix
S_CD8T_rows <- which(clean_barcode_Sathe %in% clean_barcode_S_CD8T)
#S_CD8T_rows <- which(S_CD8T$cell_barcode %in% barcodes_Sathe$V1)
length(S_CD8T_rows)
S_CD8T_matrix <- as.matrix(mtx_Sathe[, S_CD8T_rows])
dim(S_CD8T_matrix)



colnames(S_CD4T_matrix) <- paste0("Sathe_5931_t2_", clean_barcode_S_CD4T)
rownames(S_CD4T_matrix)<- features_Sathe$V2
S_CD4T_matrix[1:5,1:5]

colnames(S_CD8T_matrix) <- paste0("Sathe_5931_t2_", clean_barcode_S_CD8T)
rownames(S_CD8T_matrix)<- features_Sathe$V2
S_CD8T_matrix[1:5,1:5]

## Save
S_5931_t2_CD4 <- S_CD4T_matrix
S_5931_t2_CD8<- S_CD8T_matrix




##### 6592_t1 #####
n_Label<- cell_labels[cell_labels$orig.ident=="6592_t1",]


mtx_Sathe <- readMM(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/6592_t1/matrix.mtx.gz"))
features_Sathe <- read.delim(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/6592_t1/features.tsv.gz"), header = FALSE)
barcodes_Sathe <- read.delim(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/6592_t1/barcodes.tsv.gz"), header = FALSE)

dim(mtx_Sathe)
dim(features_Sathe)
head(features_Sathe)
dim(barcodes_Sathe)
head(barcodes_Sathe)

## Prepare an ordinary matrix from sparse matrix 
dense_mtx_Sathe <- as.matrix(mtx_Sathe)
dim(mtx_Sathe)



# Remove barcode suffix 
clean_barcode_Sathe <- sub("-\\d+$", "", barcodes_Sathe$V1)

### CD4
S_CD4T<- n_Label[n_Label$final_celltype == "CD4", ]
# Remove barcode suffix 
clean_barcode_S_CD4T <- sub("-\\d+$", "", S_CD4T$cell_barcode)
# Compare without suffix
S_CD4T_rows <- which(clean_barcode_Sathe %in% clean_barcode_S_CD4T)
#S_CD4T_rows <- which(S_CD4T$cell_barcode %in% barcodes_Sathe$V1)
length(S_CD4T_rows)
S_CD4T_matrix <- as.matrix(mtx_Sathe[, S_CD4T_rows])
dim(S_CD4T_matrix)

### CD8
S_CD8T<- n_Label[n_Label$final_celltype == "CD8", ]
# Remove barcode suffix 
clean_barcode_S_CD8T <- sub("-\\d+$", "", S_CD8T$cell_barcode)
# Compare without suffix
S_CD8T_rows <- which(clean_barcode_Sathe %in% clean_barcode_S_CD8T)
#S_CD8T_rows <- which(S_CD8T$cell_barcode %in% barcodes_Sathe$V1)
length(S_CD8T_rows)
S_CD8T_matrix <- as.matrix(mtx_Sathe[, S_CD8T_rows])
dim(S_CD8T_matrix)



colnames(S_CD4T_matrix) <- paste0("Sathe_6592_t1_", clean_barcode_S_CD4T)
rownames(S_CD4T_matrix)<- features_Sathe$V2
S_CD4T_matrix[1:5,1:5]

colnames(S_CD8T_matrix) <- paste0("Sathe_6592_t1_", clean_barcode_S_CD8T)
rownames(S_CD8T_matrix)<- features_Sathe$V2
S_CD8T_matrix[1:5,1:5]

## Save
S_6592_t1_CD4 <- S_CD4T_matrix
S_6592_t1_CD8<- S_CD8T_matrix



##### 5846_t1 #####
n_Label<- cell_labels[cell_labels$orig.ident=="5846_t1",]


mtx_Sathe <- readMM(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/5846_t1/matrix.mtx.gz"))
features_Sathe <- read.delim(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/5846_t1/features.tsv.gz"), header = FALSE)
barcodes_Sathe <- read.delim(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/5846_t1/barcodes.tsv.gz"), header = FALSE)

dim(mtx_Sathe)
dim(features_Sathe)
head(features_Sathe)
dim(barcodes_Sathe)
head(barcodes_Sathe)

## Prepare an ordinary matrix from sparse matrix 
dense_mtx_Sathe <- as.matrix(mtx_Sathe)
dim(mtx_Sathe)



# Remove barcode suffix 
clean_barcode_Sathe <- sub("-\\d+$", "", barcodes_Sathe$V1)

### CD4
S_CD4T<- n_Label[n_Label$final_celltype == "CD4", ]
# Remove barcode suffix 
clean_barcode_S_CD4T <- sub("-\\d+$", "", S_CD4T$cell_barcode)
# Compare without suffix
S_CD4T_rows <- which(clean_barcode_Sathe %in% clean_barcode_S_CD4T)
#S_CD4T_rows <- which(S_CD4T$cell_barcode %in% barcodes_Sathe$V1)
length(S_CD4T_rows)
S_CD4T_matrix <- as.matrix(mtx_Sathe[, S_CD4T_rows])
dim(S_CD4T_matrix)

### CD8
S_CD8T<- n_Label[n_Label$final_celltype == "CD8", ]
# Remove barcode suffix 
clean_barcode_S_CD8T <- sub("-\\d+$", "", S_CD8T$cell_barcode)
# Compare without suffix
S_CD8T_rows <- which(clean_barcode_Sathe %in% clean_barcode_S_CD8T)
#S_CD8T_rows <- which(S_CD8T$cell_barcode %in% barcodes_Sathe$V1)
length(S_CD8T_rows)
S_CD8T_matrix <- as.matrix(mtx_Sathe[, S_CD8T_rows])
dim(S_CD8T_matrix)

colnames(S_CD4T_matrix) <- paste0("Sathe_5846_t1_", clean_barcode_S_CD4T)
rownames(S_CD4T_matrix)<- features_Sathe$V2
S_CD4T_matrix[1:5,1:5]

colnames(S_CD8T_matrix) <- paste0("Sathe_5846_t1_", clean_barcode_S_CD8T)
rownames(S_CD8T_matrix)<- features_Sathe$V2
S_CD8T_matrix[1:5,1:5]

## Save
S_5846_t1_CD4 <- S_CD4T_matrix
S_5846_t1_CD8<- S_CD8T_matrix



##### 6207_t1 #####
n_Label<- cell_labels[cell_labels$orig.ident=="6207_t1",]


mtx_Sathe <- readMM(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/6207_t1/matrix.mtx.gz"))
features_Sathe <- read.delim(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/6207_t1/features.tsv.gz"), header = FALSE)
barcodes_Sathe <- read.delim(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/6207_t1/barcodes.tsv.gz"), header = FALSE)

dim(mtx_Sathe)
dim(features_Sathe)
head(features_Sathe)
dim(barcodes_Sathe)
head(barcodes_Sathe)

## Prepare an ordinary matrix from sparse matrix 
dense_mtx_Sathe <- as.matrix(mtx_Sathe)
dim(mtx_Sathe)



# Remove barcode suffix 
clean_barcode_Sathe <- sub("-\\d+$", "", barcodes_Sathe$V1)

### CD4
S_CD4T<- n_Label[n_Label$final_celltype == "CD4", ]
# Remove barcode suffix 
clean_barcode_S_CD4T <- sub("-\\d+$", "", S_CD4T$cell_barcode)
# Compare without suffix
S_CD4T_rows <- which(clean_barcode_Sathe %in% clean_barcode_S_CD4T)
#S_CD4T_rows <- which(S_CD4T$cell_barcode %in% barcodes_Sathe$V1)
length(S_CD4T_rows)
S_CD4T_matrix <- as.matrix(mtx_Sathe[, S_CD4T_rows])
dim(S_CD4T_matrix)

### CD8
S_CD8T<- n_Label[n_Label$final_celltype == "CD8", ]
# Remove barcode suffix 
clean_barcode_S_CD8T <- sub("-\\d+$", "", S_CD8T$cell_barcode)
# Compare without suffix
S_CD8T_rows <- which(clean_barcode_Sathe %in% clean_barcode_S_CD8T)
#S_CD8T_rows <- which(S_CD8T$cell_barcode %in% barcodes_Sathe$V1)
length(S_CD8T_rows)
S_CD8T_matrix <- as.matrix(mtx_Sathe[, S_CD8T_rows])
dim(S_CD8T_matrix)



colnames(S_CD4T_matrix) <- paste0("Sathe_6207_t1_", clean_barcode_S_CD4T)
rownames(S_CD4T_matrix)<- features_Sathe$V2
S_CD4T_matrix[1:5,1:5]

colnames(S_CD8T_matrix) <- paste0("Sathe_6207_t1_", clean_barcode_S_CD8T)
rownames(S_CD8T_matrix)<- features_Sathe$V2
S_CD8T_matrix[1:5,1:5]

## Save
S_6207_t1_CD4 <- S_CD4T_matrix
S_6207_t1_CD8<- S_CD8T_matrix




##### 6709_t1 #####
n_Label<- cell_labels[cell_labels$orig.ident=="6709_t1",]


mtx_Sathe <- readMM(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/6709_t1/matrix.mtx.gz"))
features_Sathe <- read.delim(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/6709_t1/features.tsv.gz"), header = FALSE)
barcodes_Sathe <- read.delim(gzfile("DATA/Gastric/Sathe/gastric_scRNAseq_filtered/6709_t1/barcodes.tsv.gz"), header = FALSE)

dim(mtx_Sathe)
dim(features_Sathe)
head(features_Sathe)
dim(barcodes_Sathe)
head(barcodes_Sathe)

## Prepare an ordinary matrix from sparse matrix 
dense_mtx_Sathe <- as.matrix(mtx_Sathe)
dim(mtx_Sathe)



# Remove barcode suffix 
clean_barcode_Sathe <- sub("-\\d+$", "", barcodes_Sathe$V1)

### CD4
S_CD4T<- n_Label[n_Label$final_celltype == "CD4", ]
# Remove barcode suffix 
clean_barcode_S_CD4T <- sub("-\\d+$", "", S_CD4T$cell_barcode)
# Compare without suffix
S_CD4T_rows <- which(clean_barcode_Sathe %in% clean_barcode_S_CD4T)
#S_CD4T_rows <- which(S_CD4T$cell_barcode %in% barcodes_Sathe$V1)
length(S_CD4T_rows)
S_CD4T_matrix <- as.matrix(mtx_Sathe[, S_CD4T_rows])
dim(S_CD4T_matrix)

### CD8
S_CD8T<- n_Label[n_Label$final_celltype == "CD8", ]
# Remove barcode suffix 
clean_barcode_S_CD8T <- sub("-\\d+$", "", S_CD8T$cell_barcode)
# Compare without suffix
S_CD8T_rows <- which(clean_barcode_Sathe %in% clean_barcode_S_CD8T)
#S_CD8T_rows <- which(S_CD8T$cell_barcode %in% barcodes_Sathe$V1)
length(S_CD8T_rows)
S_CD8T_matrix <- as.matrix(mtx_Sathe[, S_CD8T_rows])
dim(S_CD8T_matrix)


colnames(S_CD4T_matrix) <- paste0("Sathe_6709_t1_", clean_barcode_S_CD4T)
rownames(S_CD4T_matrix)<- features_Sathe$V2
S_CD4T_matrix[1:5,1:5]

colnames(S_CD8T_matrix) <- paste0("Sathe_6709_t1_", clean_barcode_S_CD8T)
rownames(S_CD8T_matrix)<- features_Sathe$V2
S_CD8T_matrix[1:5,1:5]

## Save
S_6709_t1_CD4 <- S_CD4T_matrix
S_6709_t1_CD8<- S_CD8T_matrix



##### MAKE GAC-Primary MATRIX #####

S_GAC_CD4<- cbind(S_5866_t1_CD4, S_5931_t1_CD4, S_5866_t2_CD4, S_6342_t1_CD4, S_5931_t2_CD4, S_6592_t1_CD4, S_5846_t1_CD4, S_6207_t1_CD4, S_6709_t1_CD4)
S_GAC_CD8<- cbind(S_5866_t1_CD8, S_5931_t1_CD8, S_5866_t2_CD8, S_6342_t1_CD8, S_5931_t2_CD8, S_6592_t1_CD8, S_5846_t1_CD8, S_6207_t1_CD8, S_6709_t1_CD8)


rownames(S_GAC_CD4)<- features_Sathe$V2
# Assign gene symbols directly (with duplicates)
gene_names <- rownames(S_GAC_CD4)
# Collapse rows with the same gene name by rowMeans
S_GAC_CD4_dedup <- rowsum(S_GAC_CD4, group = gene_names) / as.vector(table(gene_names))


rownames(S_GAC_CD8)<- features_Sathe$V2
# Assign gene symbols directly (with duplicates)
gene_names <- rownames(S_GAC_CD8)
# Collapse rows with the same gene name by rowMeans
S_GAC_CD8_dedup <- rowsum(S_GAC_CD8, group = gene_names) / as.vector(table(gene_names))



write.table(S_GAC_CD4_dedup, file = "DATA/Gastric/CellTypeStage/Sathe_GAC_Primary_CD4T_matrix.txt", quote = FALSE, row.names = T, col.names = T)
write.table(S_GAC_CD8_dedup, file = "DATA/Gastric/CellTypeStage/Sathe_GAC_Primary_CD8T_matrix.txt", quote = FALSE, row.names = T, col.names = T)




 


 





 


 


 


 


 



  
  