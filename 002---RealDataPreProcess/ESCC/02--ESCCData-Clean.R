if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")
library(edgeR)


setwd("E:/005---ThirdProject/ThirdObject/0.RealData/")


##### Read Data #####
ESCC_UMI<- read.table("ESCC/GSE199654/GSE199654_scTDN_UMI_matrix_epithelial_cells.txt/scTDN_UMI_matrix_epithelial_cells.txt")
dim(ESCC_UMI)

#ESCC_CD4<- read.table("ESCC/GSE199654/CD4_cells_barcodes.txt")

ESCC_CD8<- read.table("ESCC/GSE199654/CD8_cells_barcodes.txt")

#ESCC_b_cells<- read.table("ESCC/GSE199654/ESCC_barcodes_b_cells_barcodes.txt")

ESCC_myeloid_cells<- read.table("ESCC/GSE199654/ESCC_barcodes_myeloid_cells_barcodes.txt")
ESCC_nk_cells<- read.table("ESCC/GSE199654/ESCC_barcodes_nk_cells_barcodes.txt")




ESCC_Metadata <- read.table("ESCC/GSE199654/GSE199654_scTDN_METADATA_epithelial_cells.txt/scTDN_METADATA_epithelial_cells.txt", header = T)

ESCC_KEGG1<- read.table("ESCC/ESCC1_cancer_kegg_genes.csv", sep = ",", header = T)
ESCC_KEGG2<- read.table("ESCC/ESCC2_cancer_kegg_genes.csv", sep = ",", header = T)

ESCC_KEGG<- unique(rbind(data.frame(GeneSymbol=ESCC_KEGG1[,2]),data.frame(GeneSymbol=ESCC_KEGG2[,2])))
dim(ESCC_KEGG)

##### Prepare Data #####

unique(ESCC_Metadata$stage)
#These values represent a progression from healthy to cancerous:
#N (Normal) → LGIN → HGIN → INF → Stage_I → Stage_II → Stage_III → Stage_IV


# Count number of samples for each unique stage
stage_counts <- table(ESCC_Metadata$stage)
# Print the counts
print(stage_counts)


#dim(ESCC_Metadata[ESCC_Metadata$cellID %in% ESCC_CD8$V1,])


# Convert '.1' to '-1'
ESCC_CD8_fixed <- gsub("\\.1$", "-1", ESCC_CD8$V1)

# Now subset metadata
filtered_metadata <- ESCC_Metadata[ESCC_Metadata$cellID %in% ESCC_CD8_fixed, ]
dim(filtered_metadata)

sum(ESCC_Metadata$cellID %in% ESCC_CD8_fixed)  # Or the other way around

ESCC_Meta_CD8 <- ESCC_Metadata[ESCC_Metadata$cellID %in% ESCC_CD8_fixed,]
dim(ESCC_Meta_CD8)

print(table(ESCC_Meta_CD8$stage))



# Create output directory (optional)
dir.create("ESCC/cellIDs_by_stage", showWarnings = FALSE)

# Loop through each unique stage
for (stage in unique(ESCC_Metadata$stage)) {
  
  # Subset metadata for this stage
  subset_stage <- ESCC_Metadata[ESCC_Metadata$stage == stage, ]
  
  # Extract cellIDs
  cell_ids <- subset_stage$cellID
  
  # Create a valid file name (remove any special characters)
  file_name <- paste0("ESCC/cellIDs_by_stage/cellIDs_", gsub("[^a-zA-Z0-9]", "_", stage), ".txt")
  
  # Write to file
  writeLines(cell_ids, file_name)
}


N_cellID<- filtered_metadata[filtered_metadata$stage == "N",]$cellID
HGIN_cellID<- filtered_metadata[filtered_metadata$stage == "HGIN",]$cellID
INF_cellID<- filtered_metadata[filtered_metadata$stage == "INF",]$cellID
Stage_I_cellID<- filtered_metadata[filtered_metadata$stage == "Stage_I",]$cellID
Stage_II_cellID<- filtered_metadata[filtered_metadata$stage == "Stage_II",]$cellID
Stage_III_cellID<- filtered_metadata[filtered_metadata$stage == "Stage_III",]$cellID


N_cellID <- gsub("\\-1$", ".1", N_cellID)
HGIN_cellID<- gsub("\\-1$", ".1", HGIN_cellID)
INF_cellID<- gsub("\\-1$", ".1", INF_cellID)
Stage_I_cellID<- gsub("\\-1$", ".1", Stage_I_cellID)
Stage_II_cellID<- gsub("\\-1$", ".1", Stage_II_cellID)
Stage_III_cellID<- gsub("\\-1$", ".1", Stage_III_cellID)
  
fea_N <- ESCC_UMI[ , colnames(ESCC_UMI) %in% N_cellID, drop = FALSE]
fea_HGIN <- ESCC_UMI[ , colnames(ESCC_UMI) %in% HGIN_cellID, drop = FALSE]
fea_INF <- ESCC_UMI[ , colnames(ESCC_UMI) %in% INF_cellID, drop = FALSE]
fea_Stage_I <- ESCC_UMI[ , colnames(ESCC_UMI) %in% Stage_I_cellID, drop = FALSE]
fea_Stage_II <- ESCC_UMI[ , colnames(ESCC_UMI) %in% Stage_II_cellID, drop = FALSE]
fea_Stage_III <- ESCC_UMI[ , colnames(ESCC_UMI) %in% Stage_III_cellID, drop = FALSE]


write.table(fea_N, file = "ESCC/cellIDs_by_stage/ESCC_CD8_N.txt", quote = FALSE)
write.table(fea_HGIN, file = "ESCC/cellIDs_by_stage/ESCC_CD8_HGIN.txt", quote = FALSE)
write.table(fea_INF, file = "ESCC/cellIDs_by_stage/ESCC_CD8_INF.txt", quote = FALSE)
write.table(fea_Stage_I, file = "ESCC/cellIDs_by_stage/ESCC_CD8_Stage_I.txt", quote = FALSE)
write.table(fea_Stage_II, file = "ESCC/cellIDs_by_stage/ESCC_CD8_Stage_II.txt", quote = FALSE)
write.table(fea_Stage_III, file = "ESCC/cellIDs_by_stage/ESCC_CD8_Stage_III.txt", quote = FALSE)



clean_names <- function(x) {
  trimws(toupper(gsub("[[:space:]]+", "", x)))
}
clean_n_rownames <- clean_names(rownames(fea_N))


fea_N_cleaned <- rowsum(fea_N, group = clean_n_rownames)
fea_HGIN_cleaned <- rowsum(fea_HGIN, group = clean_n_rownames)
fea_INF_cleaned <- rowsum(fea_INF, group = clean_n_rownames)
fea_Stage_I_cleaned <- rowsum(fea_Stage_I, group = clean_n_rownames)
fea_Stage_II_cleaned <- rowsum(fea_Stage_II, group = clean_n_rownames)
fea_Stage_III_cleaned <- rowsum(fea_Stage_III, group = clean_n_rownames)



n_genes <- rownames(fea_N_cleaned)
n_genes <- toupper(n_genes)
rownames(fea_HGIN_cleaned)<- n_genes
rownames(fea_INF_cleaned)<- n_genes
rownames(fea_Stage_I_cleaned)<- n_genes
rownames(fea_Stage_II_cleaned)<- n_genes
rownames(fea_Stage_III_cleaned)<- n_genes



##### Step 2: Store stage matrices in a list #####
stage_list <- list(
  HGIN = fea_HGIN_cleaned,
  INF = fea_INF_cleaned,
  Stage_I = fea_Stage_I_cleaned,
  Stage_II = fea_Stage_II_cleaned,
  Stage_III = fea_Stage_III_cleaned
)


# Step 3: Define NAT matrix
nat <- fea_N_cleaned

# Step 4: DEG function

get_DEGs <- function(stage_matrix, nat_matrix, stage_name, logFC_cutoff = 1.0, pval_cutoff = 0.05) {
  # Step 1: Match genes
  common_genes <- intersect(rownames(stage_matrix), rownames(nat_matrix))
  stage_matrix <- stage_matrix[common_genes, ]
  nat_matrix <- nat_matrix[common_genes, ]
  
  # Step 2: Combine
  combined <- cbind(stage_matrix, nat_matrix)
  
  # Step 3: Convert all to numeric safely
  combined <- as.data.frame(lapply(combined, function(col) as.numeric(as.character(col))))
  rownames(combined) <- common_genes
  
  # Step 4: Remove columns with all NA or zero
  combined <- combined[, colSums(is.na(combined)) < nrow(combined)]        # remove all-NA columns
  combined <- combined[, colSums(combined, na.rm = TRUE) > 0]              # remove all-zero columns
  
  # Step 5: Remove rows with any NA
  combined <- combined[complete.cases(combined), ]
  
  # Step 6: Group label
  n_stage <- sum(colnames(combined) %in% colnames(stage_matrix))
  n_nat <- sum(colnames(combined) %in% colnames(nat_matrix))
  group <- factor(c(rep("Stage", n_stage), rep("NAT", n_nat)))
  
  # Step 7: Run edgeR
  dge <- DGEList(counts = combined)
  dge <- calcNormFactors(dge)
  design <- model.matrix(~group)
  dge <- estimateDisp(dge, design)
  fit <- glmFit(dge, design)
  lrt <- glmLRT(fit, coef = 2)
  
  # Step 8: Filter DEGs
  degs <- topTags(lrt, n = Inf)$table
  degs <- degs[abs(degs$logFC) > logFC_cutoff & degs$PValue < pval_cutoff, ]
  degs$gene <- rownames(degs)
  
  print(paste("Finished:", stage_name))
  return(degs)
}


# Step 5: Apply DEG analysis to each stage
deg_list <- lapply(names(stage_list), function(stage) {
  get_DEGs(stage_list[[stage]], nat, stage)
})
names(deg_list) <- names(stage_list)

#Step 6: Optionally, write to files
for(stage in names(deg_list)) {
  write.csv(deg_list[[stage]], paste0("DEGs_", stage, "_vs_NAT_CD8_20250520.csv"), row.names = FALSE)
}



#Step 6: Optionally, write to files
for(stage in names(deg_list)) {
  write.csv(deg_list[[stage]], paste0("ESCC/cellIDs_by_stage/DEGs_", stage, "_vs_Nor_CD8_20250527.csv"), row.names = FALSE)
}





##### Build Networks with RegNetwork as known knowledge #####
# Load Gastric KEGG genes
ESCC_KEGG

kegg_genes<- (unique(ESCC_KEGG))

# 1. Get union of all DEGs across all stages
all_deg_genes <- unique(unlist(lapply(deg_list, function(df) df$gene)))
cat("Total unique DEG genes across all stages:", length(all_deg_genes), "\n")

nat_genes <- rownames(nat)

# kegg_genes <- toupper(kegg_genes)
# all_deg_genes <- toupper(all_deg_genes)
# nat_genes <- toupper(nat_genes)

sum(kegg_genes %in% nat_genes)
sum(kegg_genes %in% all_deg_genes)

common_kegg_in_nat <- kegg_genes[kegg_genes$GeneSymbol %in% nat_genes,]

selected_genes <- unique(c(common_kegg_in_nat, all_deg_genes))
cat("Total unique Selected genes across all stages:", length(selected_genes), "\n")


##### Load RegNetwork Edges
RegNet<- read.table("RegNetwork-2022.human.source", head=F)
RegNet<- RegNet[,c(1,3)]


# Apply toupper to both columns
RegNet <- data.frame(lapply(RegNet, toupper), stringsAsFactors = FALSE)



# Filter RegNet where both source and target are in selected_genes
filtered_RegNet <- RegNet[RegNet[,1] %in% selected_genes & RegNet[,2] %in% selected_genes, ]

# Check result
cat("Number of edges in filtered RegNet:", nrow(filtered_RegNet), "\n")


# Get unique genes (nodes) from both columns of filtered_RegNet
remaining_genes_in_RegNet <- unique(c(filtered_RegNet[,1], filtered_RegNet[,2]))

# Check how many unique genes are in the filtered network
cat("Number of unique genes in filtered RegNet:", length(remaining_genes_in_RegNet), "\n")


# Find KEGG genes that are still present in the filtered RegNet
kegg_genes_in_RegNet <- intersect(common_kegg_in_nat, remaining_genes_in_RegNet)
# Count how many KEGG genes remain
cat("Number of KEGG genes in filtered RegNet:", length(kegg_genes_in_RegNet), "\n")



# Create a list of cleaned matrices
cleaned_matrices <- list(
  fea_N_cleaned = fea_N_cleaned,
  fea_HGIN_cleaned = fea_HGIN_cleaned,
  fea_INF_cleaned = fea_INF_cleaned,
  fea_Stage_I_cleaned = fea_Stage_I_cleaned,
  fea_Stage_II_cleaned = fea_Stage_II_cleaned,
  fea_Stage_III_cleaned = fea_Stage_III_cleaned
)



#  Create output folder
dir.create("ESCC/cellIDs_by_stage/Filtered_Cleaned_Inputs", showWarnings = FALSE)

# Filter and save each cleaned matrix with quoted colnames and rownames
for (name in names(cleaned_matrices)) {
  mat <- cleaned_matrices[[name]]
  filtered_mat <- mat[intersect(rownames(mat), remaining_genes_in_RegNet), ]
  
  # Add rownames as a separate column
  filtered_mat_with_rownames <- cbind(Gene = rownames(filtered_mat), filtered_mat)
  
  # Write CSV with quoted colnames and rownames
  write.table(
    filtered_mat_with_rownames,
    file = paste0("ESCC/cellIDs_by_stage/Filtered_Cleaned_Inputs/", name, "_Filtered.csv"),
    sep = ",",
    row.names = FALSE,
    col.names = TRUE,
    quote = TRUE
  )
}




# Step 1: Ensure directory exists
dir.create("ESCC/cellIDs_by_stage/Filtered_Cleaned_Inputs", showWarnings = FALSE)

# Step 2: Write header manually with quotes
writeLines('"Gene1","Gene2"', "ESCC/cellIDs_by_stage/Filtered_Cleaned_Inputs/filtered_RegNet_CD8.csv")

# Step 3: Append data without quotes, no row names, no column names
write.table(filtered_RegNet,
            file = "ESCC/cellIDs_by_stage/Filtered_Cleaned_Inputs/filtered_RegNet_CD8.csv",
            sep = ",",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE,
            append = TRUE)

cat("✅ Filtered cleaned matrices saved in 'Filtered_Cleaned_Inputs/' folder.\n")




dim(RegNet[RegNet[,1] %in% selected_genes & RegNet[,2] %in% selected_genes, ])
dim(RegNet[RegNet[,1] %in% deg_list$HGIN$gene & RegNet[,2] %in% deg_list$HGIN$gene, ])
dim(RegNet[RegNet[,1] %in% deg_list$INF$gene & RegNet[,2] %in% deg_list$INF$gene, ])
dim(RegNet[RegNet[,1] %in% deg_list$Stage_I$gene & RegNet[,2] %in% deg_list$Stage_I$gene, ])
dim(RegNet[RegNet[,1] %in% deg_list$Stage_II$gene & RegNet[,2] %in% deg_list$Stage_II$gene, ])
dim(RegNet[RegNet[,1] %in% deg_list$Stage_III$gene & RegNet[,2] %in% deg_list$Stage_III$gene, ])













