### Integrate stages' data from Four data-sets where consider gene names in combination



setwd("/home/Fatemeh/0--ThirdProject")

##### Load Data #####


##### Integrate NAT Data #####

##### NAT CD4 #####
MAD_NAT_CD4T <- read.table("DATA/Gastric/CellTypeStage/MAD_nat_CD4T_matrix.txt", header = T)
dim(MAD_NAT_CD4T)


Sathe_NAT_CD4T <- read.table("DATA/Gastric/CellTypeStage/Sathe_NAT_CD4T_matrix.txt", header = T)
dim(Sathe_NAT_CD4T)


Zhang_NAG_CD4 <- read.table("DATA/Gastric/CellTypeStage/Zhang_NAG_CD4.txt", header = T)
dim(Zhang_NAG_CD4)


# Add gene names as a column
df1<- MAD_NAT_CD4T
df2<- Zhang_NAG_CD4
df3<- Sathe_NAT_CD4T

df1$gene <- rownames(MAD_NAT_CD4T)
df2$gene <- rownames(Zhang_NAG_CD4)
df3$gene <- rownames(Sathe_NAT_CD4T)

# Merge step by step
merged_df <- merge(df1, df2, by = "gene", all = TRUE)
merged_df <- merge(merged_df, df3, by = "gene", all = TRUE)

# Set rownames and remove gene column
rownames(merged_df) <- merged_df$gene
merged_df$gene <- NULL

merged_df[is.na(merged_df)] <- 0
dim(merged_df)


#WRITE
write.table(merged_df, file = "DATA/Gastric/CellTypeStage/Integrated/NAT_CD4T_matrix.txt", quote = FALSE, row.names = T, col.names = T)


##### NAT CD8 #####
MAD_NAT_CD8T <- read.table("DATA/Gastric/CellTypeStage/MAD_nat_CD8T_matrix.txt", header = T)
dim(MAD_NAT_CD8T)


Sathe_NAT_CD8T <- read.table("DATA/Gastric/CellTypeStage/Sathe_NAT_CD8T_matrix.txt", header = T)
dim(Sathe_NAT_CD8T)


Zhang_NAG_CD8 <- read.table("DATA/Gastric/CellTypeStage/Zhang_NAG_CD8.txt", header = T)
dim(Zhang_NAG_CD8)



# Add gene names as a column
df1<- MAD_NAT_CD8T
df2<- Zhang_NAG_CD8
df3<- Sathe_NAT_CD8T

df1$gene <- rownames(MAD_NAT_CD8T)
df2$gene <- rownames(Zhang_NAG_CD8)
df3$gene <- rownames(Sathe_NAT_CD8T)

# Merge step by step
merged_df <- merge(df1, df2, by = "gene", all = TRUE)
merged_df <- merge(merged_df, df3, by = "gene", all = TRUE)

# Set rownames and remove gene column
rownames(merged_df) <- merged_df$gene
merged_df$gene <- NULL

merged_df[is.na(merged_df)] <- 0
dim(merged_df)

#WRITE
write.table(merged_df, file = "DATA/Gastric/CellTypeStage/Integrated/NAT_CD8T_matrix.txt", quote = FALSE, row.names = T, col.names = T)





##### Integrate CAG Data #####
##### CAG CD4 #####
Zhang_CAG_CD4 <- read.table("DATA/Gastric/CellTypeStage/Zhang_CAG_CD4.txt", header = T)
dim(Zhang_CAG_CD4)


#WRITE
write.table(Zhang_CAG_CD4, file = "DATA/Gastric/CellTypeStage/Integrated/CAG_CD4T_matrix.txt", quote = FALSE, row.names = T, col.names = T)


##### CAG CD8 #####
Zhang_CAG_CD8 <- read.table("DATA/Gastric/CellTypeStage/Zhang_CAG_CD8.txt", header = T)
dim(Zhang_CAG_CD8)


#WRITE
write.table(Zhang_CAG_CD8, file = "DATA/Gastric/CellTypeStage/Integrated/CAG_CD8T_matrix.txt", quote = FALSE, row.names = T, col.names = T)






##### Integrate IM Data #####
##### IM CD4 #####
Sathe_IM_CD4T <- read.table("DATA/Gastric/CellTypeStage/Sathe_IM_CD4T_matrix.txt", header = T)
dim(Sathe_IM_CD4T)


Zhang_IM_CD4 <- read.table("DATA/Gastric/CellTypeStage/Zhang_IM_CD4.txt", header = T)
dim(Zhang_IM_CD4)


# Add gene names as a column
df1<- Sathe_IM_CD4T
df2<- Zhang_IM_CD4

df1$gene <- rownames(Sathe_IM_CD4T)
df2$gene <- rownames(Zhang_IM_CD4)


# Merge step by step
merged_df <- merge(df1, df2, by = "gene", all = TRUE)

# Set rownames and remove gene column
rownames(merged_df) <- merged_df$gene
merged_df$gene <- NULL

merged_df[is.na(merged_df)] <- 0
dim(merged_df)


#WRITE
write.table(merged_df, file = "DATA/Gastric/CellTypeStage/Integrated/IM_CD4T_matrix.txt", quote = FALSE, row.names = T, col.names = T)




##### IM CD8 #####
Sathe_IM_CD8T <- read.table("DATA/Gastric/CellTypeStage/Sathe_IM_CD8T_matrix.txt", header = T)
dim(Sathe_IM_CD8T)


Zhang_IM_CD8 <- read.table("DATA/Gastric/CellTypeStage/Zhang_IM_CD8.txt", header = T)
dim(Zhang_IM_CD8)


# Add gene names as a column
df1<- Sathe_IM_CD8T
df2<- Zhang_IM_CD8

df1$gene <- rownames(Sathe_IM_CD8T)
df2$gene <- rownames(Zhang_IM_CD8)


# Merge step by step
merged_df <- merge(df1, df2, by = "gene", all = TRUE)

# Set rownames and remove gene column
rownames(merged_df) <- merged_df$gene
merged_df$gene <- NULL

merged_df[is.na(merged_df)] <- 0
dim(merged_df)


#WRITE
write.table(merged_df, file = "DATA/Gastric/CellTypeStage/Integrated/IM_CD8T_matrix.txt", quote = FALSE, row.names = T, col.names = T)





##### Integrate PGAC Data #####
##### P_GAC CD4 #####
MAD_PGAC_CD4T <- read.table("DATA/Gastric/CellTypeStage/MAD_GAC_Primary_CD4T_matrix.txt", header = T)
dim(MAD_PGAC_CD4T)


Sathe_PGAC_CD4T <- read.table("DATA/Gastric/CellTypeStage/Sathe_GAC_Primary_CD4T_matrix.txt", header = T)
dim(Sathe_PGAC_CD4T)


Zhang_PGAC_CD4 <- read.table("DATA/Gastric/CellTypeStage/Zhang_EGC_CD4.txt", header = T)
dim(Zhang_PGAC_CD4)


# Add gene names as a column
df1<- MAD_PGAC_CD4T
df2<- Sathe_PGAC_CD4T
df3<- Zhang_PGAC_CD4

df1$gene <- rownames(MAD_PGAC_CD4T)
df2$gene <- rownames(Sathe_PGAC_CD4T)
df3$gene <- rownames(Zhang_PGAC_CD4)

# Merge step by step
merged_df <- merge(df1, df2, by = "gene", all = TRUE)
merged_df <- merge(merged_df, df3, by = "gene", all = TRUE)

# Set rownames and remove gene column
rownames(merged_df) <- merged_df$gene
merged_df$gene <- NULL

merged_df[is.na(merged_df)] <- 0
dim(merged_df)


#WRITE
write.table(merged_df, file = "DATA/Gastric/CellTypeStage/Integrated/PGAC_CD4T_matrix.txt", quote = FALSE, row.names = T, col.names = T)


##### P_GAC CD8 #####
MAD_PGAC_CD8T <- read.table("DATA/Gastric/CellTypeStage/MAD_GAC_Primary_CD8T_matrix.txt", header = T)
dim(MAD_PGAC_CD8T)


Sathe_PGAC_CD8T <- read.table("DATA/Gastric/CellTypeStage/Sathe_GAC_Primary_CD8T_matrix.txt", header = T)
dim(Sathe_PGAC_CD8T)


Zhang_PGAC_CD8 <- read.table("DATA/Gastric/CellTypeStage/Zhang_EGC_CD8.txt", header = T)
dim(Zhang_PGAC_CD8)


# Add gene names as a column
df1<- MAD_PGAC_CD8T
df2<- Sathe_PGAC_CD8T
df3<- Zhang_PGAC_CD8

df1$gene <- rownames(MAD_PGAC_CD8T)
df2$gene <- rownames(Sathe_PGAC_CD8T)
df3$gene <- rownames(Zhang_PGAC_CD8)

# Merge step by step
merged_df <- merge(df1, df2, by = "gene", all = TRUE)
merged_df <- merge(merged_df, df3, by = "gene", all = TRUE)

# Set rownames and remove gene column
rownames(merged_df) <- merged_df$gene
merged_df$gene <- NULL

merged_df[is.na(merged_df)] <- 0
dim(merged_df)


#WRITE
write.table(merged_df, file = "DATA/Gastric/CellTypeStage/Integrated/PGAC_CD8T_matrix.txt", quote = FALSE, row.names = T, col.names = T)






##### Integrate Metastasis Data #####
##### Metastasis CD4 #####
MAD_Met_CD4 <- read.table("DATA/Gastric/CellTypeStage/MAD_Metastasis_CD4T_matrix.txt", header = T)
dim(MAD_Met_CD4)


#WRITE
write.table(MAD_Met_CD4, file = "DATA/Gastric/CellTypeStage/Integrated/Metastasis_CD4T_matrix.txt", quote = FALSE, row.names = T, col.names = T)


MAD_Met_CD8 <- read.table("DATA/Gastric/CellTypeStage/MAD_Metastasis_CD8T_matrix.txt", header = T)
dim(MAD_Met_CD8)


#WRITE
write.table(MAD_Met_CD8, file = "DATA/Gastric/CellTypeStage/Integrated/Metastasis_CD8T_matrix.txt", quote = FALSE, row.names = T, col.names = T)



##### Make Full Matrices #####

### CD4 ###
# Read matrices
NAT <- read.table("DATA/Gastric/CellTypeStage/Integrated/NAT_CD4T_matrix.txt", header = TRUE, row.names = 1)
CAG <- read.table("DATA/Gastric/CellTypeStage/Integrated/CAG_CD4T_matrix.txt", header = TRUE, row.names = 1)
IM <- read.table("DATA/Gastric/CellTypeStage/Integrated/IM_CD4T_matrix.txt", header = TRUE, row.names = 1)
PGAC <- read.table("DATA/Gastric/CellTypeStage/Integrated/PGAC_CD4T_matrix.txt", header = TRUE, row.names = 1)
Metastasis <- read.table("DATA/Gastric/CellTypeStage/Integrated/Metastasis_CD4T_matrix.txt", header = TRUE, row.names = 1)

# Get union of all gene names
all_genes <- Reduce(union, list(rownames(NAT), rownames(CAG), rownames(IM), rownames(PGAC), rownames(Metastasis)))

# Function to add missing genes and set expression to 0
fill_missing_genes <- function(mat, all_genes) {
  missing_genes <- setdiff(all_genes, rownames(mat))
  if (length(missing_genes) > 0) {
    zeros <- matrix(0, nrow = length(missing_genes), ncol = ncol(mat))
    rownames(zeros) <- missing_genes
    colnames(zeros) <- colnames(mat)
    mat <- rbind(mat, zeros)
  }
  mat[all_genes, , drop = FALSE]  # Reorder rows
}

# Apply function to all stages
NAT_full <- fill_missing_genes(NAT, all_genes)
CAG_full <- fill_missing_genes(CAG, all_genes)
IM_full <- fill_missing_genes(IM, all_genes)
PGAC_full <- fill_missing_genes(PGAC, all_genes)
Metastasis_full <- fill_missing_genes(Metastasis, all_genes)

# Save the aligned matrices
write.table(NAT_full, "DATA/Gastric/CellTypeStage/Integrated/Full/NAT_CD4T_matrix_Full.txt", sep = "\t", quote = FALSE)
write.table(CAG_full, "DATA/Gastric/CellTypeStage/Integrated/Full/CAG_CD4T_matrix_Full.txt", sep = "\t", quote = FALSE)
write.table(IM_full, "DATA/Gastric/CellTypeStage/Integrated/Full/IM_CD4T_matrix_Full.txt", sep = "\t", quote = FALSE)
write.table(PGAC_full, "DATA/Gastric/CellTypeStage/Integrated/Full/PGAC_CD4T_matrix_Full.txt", sep = "\t", quote = FALSE)
write.table(Metastasis_full, "DATA/Gastric/CellTypeStage/Integrated/Full/Metastasis_CD4T_matrix_Full.txt", sep = "\t", quote = FALSE)



### CD8 ###
# Read matrices
NAT <- read.table("DATA/Gastric/CellTypeStage/Integrated/NAT_CD8T_matrix.txt", header = TRUE, row.names = 1)
CAG <- read.table("DATA/Gastric/CellTypeStage/Integrated/CAG_CD8T_matrix.txt", header = TRUE, row.names = 1)
IM <- read.table("DATA/Gastric/CellTypeStage/Integrated/IM_CD8T_matrix.txt", header = TRUE, row.names = 1)
PGAC <- read.table("DATA/Gastric/CellTypeStage/Integrated/PGAC_CD8T_matrix.txt", header = TRUE, row.names = 1)
Metastasis <- read.table("DATA/Gastric/CellTypeStage/Integrated/Metastasis_CD8T_matrix.txt", header = TRUE, row.names = 1)

# Get union of all gene names
all_genes <- Reduce(union, list(rownames(NAT), rownames(CAG), rownames(IM), rownames(PGAC), rownames(Metastasis)))

# Function to add missing genes and set expression to 0
fill_missing_genes <- function(mat, all_genes) {
  missing_genes <- setdiff(all_genes, rownames(mat))
  if (length(missing_genes) > 0) {
    zeros <- matrix(0, nrow = length(missing_genes), ncol = ncol(mat))
    rownames(zeros) <- missing_genes
    colnames(zeros) <- colnames(mat)
    mat <- rbind(mat, zeros)
  }
  mat[all_genes, , drop = FALSE]  # Reorder rows
}

# Apply function to all stages
NAT_full <- fill_missing_genes(NAT, all_genes)
CAG_full <- fill_missing_genes(CAG, all_genes)
IM_full <- fill_missing_genes(IM, all_genes)
PGAC_full <- fill_missing_genes(PGAC, all_genes)
Metastasis_full <- fill_missing_genes(Metastasis, all_genes)

# Save the aligned matrices
write.table(NAT_full, "DATA/Gastric/CellTypeStage/Integrated/Full/NAT_CD8T_matrix_Full.txt", sep = "\t", quote = FALSE)
write.table(CAG_full, "DATA/Gastric/CellTypeStage/Integrated/Full/CAG_CD8T_matrix_Full.txt", sep = "\t", quote = FALSE)
write.table(IM_full, "DATA/Gastric/CellTypeStage/Integrated/Full/IM_CD8T_matrix_Full.txt", sep = "\t", quote = FALSE)
write.table(PGAC_full, "DATA/Gastric/CellTypeStage/Integrated/Full/PGAC_CD8T_matrix_Full.txt", sep = "\t", quote = FALSE)
write.table(Metastasis_full, "DATA/Gastric/CellTypeStage/Integrated/Full/Metastasis_CD8T_matrix_Full.txt", sep = "\t", quote = FALSE)



