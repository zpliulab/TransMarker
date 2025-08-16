if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")
library(edgeR)

setwd("E:/005---ThirdProject/ThirdObject/0.RealData/")


##### DEG Identification with edgeR #####

# Step 1: Read matrices (genes as rownames, samples as columns)
NAT_CD8T_M <- read.table("Full/NAT_CD8T_matrix_Full.txt", header = TRUE, row.names = 1)
CAG_CD8T_M <- read.table("Full/CAG_CD8T_matrix_Full.txt", header = TRUE, row.names = 1)
IM_CD8T_M <- read.table("Full/IM_CD8T_matrix_Full.txt", header = TRUE, row.names = 1)
PGAC_CD8T_M <- read.table("Full/PGAC_CD8T_matrix_Full.txt", header = TRUE, row.names = 1)
Metastasis_CD8T_M <- read.table("Full/Metastasis_CD8T_matrix_Full.txt", header = TRUE, row.names = 1)


clean_names <- function(x) {
  trimws(toupper(gsub("[[:space:]]+", "", x)))
}


clean_nat_rownames <- clean_names(rownames(NAT_CD8T_M))
NAT_CD8T_M_cleaned <- rowsum(NAT_CD8T_M, group = clean_nat_rownames)
CAG_CD8T_M_cleaned <- rowsum(CAG_CD8T_M, group = clean_nat_rownames)
IM_CD8T_M_cleaned <- rowsum(IM_CD8T_M, group = clean_nat_rownames)
PGAC_CD8T_M_cleaned <- rowsum(PGAC_CD8T_M, group = clean_nat_rownames)
Metastasis_CD8T_M_cleaned <- rowsum(Metastasis_CD8T_M, group = clean_nat_rownames)

n_genes <- rownames(NAT_CD8T_M_cleaned)
n_genes <- toupper(n_genes)
rownames(NAT_CD8T_M_cleaned)<- n_genes
rownames(CAG_CD8T_M_cleaned)<- n_genes
rownames(IM_CD8T_M_cleaned)<- n_genes
rownames(PGAC_CD8T_M_cleaned)<- n_genes
rownames(Metastasis_CD8T_M_cleaned)<- n_genes




write.table(data.frame(rownames(NAT_CD8T_M_cleaned)), "gene_Names_CD8.txt", sep = "\t", col.names= FALSE, row.names = FALSE, quote = FALSE)




##### Step 2: Store stage matrices in a list #####
stage_list <- list(
  CAG = CAG_CD8T_M_cleaned,
  IM = IM_CD8T_M_cleaned,
  PGAC = PGAC_CD8T_M_cleaned,
  Metastasis = Metastasis_CD8T_M_cleaned
)

# Step 3: Define NAT matrix
nat <- NAT_CD8T_M_cleaned

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
  write.csv(deg_list[[stage]], paste0("DEGs_", stage, "_vs_NAT_CD8.csv"), row.names = FALSE)
}




##### Build Networks with RegNetwork as known knowledge #####
# Load Gastric KEGG genes
kegg_genes <- read.csv("gastric_cancer_kegg_genes.csv", stringsAsFactors = FALSE)
kegg_genes <- kegg_genes[,2]
kegg_genes<- (unique(kegg_genes))

# 1. Get union of all DEGs across all stages
all_deg_genes <- unique(unlist(lapply(deg_list, function(df) df$gene)))
cat("Total unique DEG genes across all stages:", length(all_deg_genes), "\n")

nat_genes <- rownames(nat)

kegg_genes <- toupper(kegg_genes)
all_deg_genes <- toupper(all_deg_genes)
nat_genes <- toupper(nat_genes)

sum(kegg_genes %in% nat_genes)
sum(kegg_genes %in% all_deg_genes)

common_kegg_in_nat <- kegg_genes[kegg_genes %in% nat_genes]

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


# # Find KEGG genes that are still present in the filtered RegNet
# kegg_genes_in_RegNet <- intersect(kegg_genes, remaining_genes_in_RegNet)
# # Count how many KEGG genes remain
# cat("Number of KEGG genes in filtered RegNet:", length(kegg_genes_in_RegNet), "\n")



# Create a list of cleaned matrices
cleaned_matrices <- list(
  NAT_CD8T_M_cleaned = NAT_CD8T_M_cleaned,
  CAG_CD8T_M_cleaned = CAG_CD8T_M_cleaned,
  IM_CD8T_M_cleaned = IM_CD8T_M_cleaned,
  PGAC_CD8T_M_cleaned = PGAC_CD8T_M_cleaned,
  Metastasis_CD8T_M_cleaned = Metastasis_CD8T_M_cleaned
)

#  Create output folder

dir.create("Filtered_Cleaned_Inputs", showWarnings = FALSE)

# Filter and save each cleaned matrix with quoted colnames and rownames
for (name in names(cleaned_matrices)) {
  mat <- cleaned_matrices[[name]]
  filtered_mat <- mat[intersect(rownames(mat), remaining_genes_in_RegNet), ]
  
  # Add rownames as a separate column
  filtered_mat_with_rownames <- cbind(Gene = rownames(filtered_mat), filtered_mat)
  
  # Write CSV with quoted colnames and rownames
  write.table(
    filtered_mat_with_rownames,
    file = paste0("Filtered_Cleaned_Inputs/", name, "_Filtered.csv"),
    sep = ",",
    row.names = FALSE,
    col.names = TRUE,
    quote = TRUE
  )
}

# Step 1: Ensure directory exists
dir.create("Filtered_Cleaned_Inputs", showWarnings = FALSE)

# Step 2: Write header manually with quotes
writeLines('"Gene1","Gene2"', "Filtered_Cleaned_Inputs/filtered_RegNet_CD8.csv")

# Step 3: Append data without quotes, no row names, no column names
write.table(filtered_RegNet,
            file = "Filtered_Cleaned_Inputs/filtered_RegNet_CD8.csv",
            sep = ",",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE,
            append = TRUE)

cat("✅ Filtered cleaned matrices saved in 'Filtered_Cleaned_Inputs/' folder.\n")













##### Use ANOVA to find DEGs #####
# NAT<- NAT_CD4T_M_cleaned
# CAG <- CAG_CD4T_M_cleaned
# IM <- IM_CD4T_M_cleaned
# PGAC <- PGAC_CD4T_M_cleaned
# Metastasis <- Metastasis_CD4T_M_cleaned
# 
# 
# 
# # Vector indicating stage of each sample
# group <- factor(c(
#   rep("NAT", ncol(NAT)),
#   rep("CAG", ncol(CAG)),
#   rep("IM", ncol(IM)),
#   rep("PGAC", ncol(PGAC)),
#   rep("Metastasis", ncol(Metastasis))
# ))
# 
# # Combine all matrices into one expression matrix
# combined_matrix <- cbind(NAT, CAG, IM, PGAC, Metastasis)
# 
# # Perform ANOVA for each gene (row)
# p_values <- apply(combined_matrix, 1, function(expr) {
#   fit <- aov(expr ~ group)
#   summary(fit)[[1]][["Pr(>F)"]][1]
# })
# 
# # Adjust p-values
# adj_p_values <- p.adjust(p_values, method = "BH")
# 
# # Identify DEGs
# deg_anova <- names(adj_p_values[adj_p_values < 0.0001])
# cat("Number of DEGs identified by ANOVA:", length(deg_anova), "\n")
# 
# # Optional: write result
# write.csv(data.frame(Gene = names(adj_p_values), Adjusted_P = adj_p_values),
#           "anova_deg_results.csv", row.names = FALSE)
# 


# 
# ##### Build Co-expression Networks with PCC and Resampling #####
# 
# library(matrixStats)  # For row/column stats
# #install.packages("Hmisc")
# library(Hmisc)
# 
# # Load Gastric KEGG genes
# kegg_genes <- read.csv("gastric_cancer_kegg_genes.csv", stringsAsFactors = FALSE)
# #kegg_genes <- kegg_genes[,2]
# #kegg_genes<- (unique(kegg_genes))
# 
# ## Common KEGG and All Genes
# # All_G<- data.frame(rownames(NAT_CD4T_M))
# # common_All_KEGG<- All_G[All_G$rownames.NAT_CD4T_M. %in% kegg_genes,]
# # 
# 
# 
# 
# # 1. Get union of all DEGs across all stages
# all_deg_genes <- unique(unlist(lapply(deg_list, function(df) df$gene)))
# cat("Total unique DEG genes across all stages:", length(all_deg_genes), "\n")
# 
# 
# 
# kegg_genes <- unique(trimws(kegg_genes[, 2]))  # remove leading/trailing spaces
# all_deg_genes <- unique(unlist(lapply(deg_list, function(df) trimws(df$gene))))
# nat_genes <- trimws(rownames(NAT_CD4T_M))
# 
# kegg_genes <- toupper(kegg_genes)
# all_deg_genes <- toupper(all_deg_genes)
# nat_genes <- toupper(nat_genes)
# 
# # Find common genes pairwise
# common_kegg_deg <- intersect(kegg_genes, all_deg_genes)
# common_kegg_nat <- intersect(kegg_genes, nat_genes)
# common_deg_nat <- intersect(all_deg_genes, nat_genes)
# 
# # Print results
# cat("KEGG ∩ DEG:", length(common_kegg_deg), "genes\n")
# cat("KEGG ∩ NAT:", length(common_kegg_nat), "genes\n")
# cat("DEG ∩ NAT:", length(common_deg_nat), "genes\n")
# 
# 
# 
# # Ensure both vectors are cleaned and standardized
# kegg_genes_clean <- unique(toupper(trimws(kegg_genes)))           # KEGG genes
# all_deg_genes_clean <- unique(toupper(trimws(all_deg_genes)))     # DEG genes
# 
# # Create the union
# deg_kegg_union <- unique(union(all_deg_genes_clean, common_kegg_nat))
# 
# 
# # Optional: check the number of genes
# cat("Total unique genes in DEG ∪ KEGG:", length(deg_kegg_union), "\n")
# 
# dim(nat[rownames(nat) %in% deg_kegg_union,])
# 
# 
# # 2. Filter all matrices to keep deg_kegg_union genes
# filtered_stage_list <- lapply(stage_list, function(mat) {
#   mat[deg_kegg_union,]
# })
# 
# # Include NAT layer too
# filtered_nat <- nat[deg_kegg_union, ]
# 
# 
# 
# 
# 
# clean_names <- function(x) {
#   trimws(toupper(gsub("[[:space:]]+", "", x)))
# }
# 
# clean_deg_kegg_union <- clean_names(deg_kegg_union)
# clean_nat_rownames <- clean_names(rownames(nat))
# 
# # Now try again
# sum(clean_deg_kegg_union %in% clean_nat_rownames)  # should now be closer to 2494
# setdiff(clean_deg_kegg_union, clean_nat_rownames)  # should be much smaller
# 
# nat_cleaned <- rowsum(nat, group = clean_nat_rownames)
# 
# nat[clean_deg_kegg_union,]
# 
# CAG_CD4T_M [clean_deg_kegg_union,]
# 
# CAG_CD4T_M1<- CAG_CD4T_M
# rownames(CAG_CD4T_M1)<- clean_nat_rownames
# 
# 
# 
# # Combine all for network construction
# all_layers <- c(filtered_stage_list, list(NAT = filtered_nat))
# layer_names <- names(all_layers)
# 
# 
# filtered_gene_list<- rownames(filtered_nat)
# 
# 
# 
# 
# ## Load RegNetwork Edges
# RegNet<- read.table("RegNetwork-2022.human.source", head=F)
# RegNet<- RegNet[,c(1,3)]
# 
# 
# 
# 
# 
# 
# 
# # 3. Define function to compute PCC and resampling threshold
# build_coexp_network <- function(expr_matrix, n_resample = 1000, pcc_thresh = 0.7) {
#   # Convert to numeric matrix
#   gene_expr <- as.matrix(expr_matrix)
#   gene_expr <- apply(gene_expr, 2, as.numeric)
#   rownames(gene_expr) <- rownames(expr_matrix)  # Restore rownames
#   
#   # Remove rows with NA after coercion
#   gene_expr <- gene_expr[complete.cases(gene_expr), ]
#   
#   # Remove constant rows
#   gene_expr <- gene_expr[rowSds(gene_expr) > 0, ]
#   
#   # Compute PCC
#   cor_mat <- cor(t(gene_expr), method = "pearson", use = "pairwise.complete.obs")
#   
#   # Resampling-based null distribution
#   null_pccs <- replicate(n_resample, {
#     i <- sample(1:nrow(gene_expr), 1)
#     j <- sample(1:nrow(gene_expr), 1)
#     while (i == j) j <- sample(1:nrow(gene_expr), 1)
#     cor(gene_expr[i, ], gene_expr[j, ])
#   })
#   
#   cutoff <- quantile(null_pccs, 0.99, na.rm = TRUE)
#   cat("Resampling PCC threshold:", round(cutoff, 3), "\n")
#   final_thresh <- max(cutoff, pcc_thresh)
#   
#   adj <- (abs(cor_mat) >= final_thresh) * 1
#   diag(adj) <- 0
#   return(adj)
# }
# 
# 
# 
# # 4. Build networks for each layer
# coexp_networks <- list()
# for (layer in layer_names) {
#   cat("Building network for:", layer, "\n")
#   coexp_networks[[layer]] <- build_coexp_network(all_layers[[layer]])
# }
# 
# 
# 
# # Convert adjacency matrices to igraph graphs and visualize
# library(igraph)
# #install.packages("ggraph")
# library(ggraph)
# library(ggplot2)
# 
# 
# 
# # Store graphs and counts
# graph_list <- list()
# node_edge_stats <- data.frame(Layer=character(), Nodes=integer(), Edges=integer(), stringsAsFactors=FALSE)
# 
# for (layer in names(coexp_networks)) {
#   adj <- coexp_networks[[layer]]
#   
#   # Create igraph object
#   g <- graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
#   graph_list[[layer]] <- g
#   
#   # Count nodes and edges
#   node_count <- vcount(g)
#   edge_count <- ecount(g)
#   node_edge_stats <- rbind(node_edge_stats, data.frame(Layer=layer, Nodes=node_count, Edges=edge_count))
#   
#   # Plot
#   cat(paste0("Layer: ", layer, " | Nodes: ", node_count, " | Edges: ", edge_count, "\n"))
#   
#   print(
#     ggraph(g, layout = 'fr') +  # force-directed layout
#       geom_edge_link(alpha = 0.5) +
#       geom_node_point(size = 2, color = "steelblue") +
#       ggtitle(paste("Co-expression Network -", layer)) +
#       theme_void()
#   )
# }
# 
# 
# 
# 
# 
# 
# 
# # Get full union of all nodes from all layers
# all_genes_union <- unique(unlist(lapply(coexp_networks, function(mat) rownames(mat))))
# 
# # Re-align all adjacency matrices to same size (square: union x union)
# aligned_networks <- list()
# for (layer in names(coexp_networks)) {
#   adj <- coexp_networks[[layer]]
#   
#   # Add missing genes as zero rows/columns
#   missing_genes <- setdiff(all_genes_union, rownames(adj))
#   if (length(missing_genes) > 0) {
#     zero_block <- matrix(0, nrow = length(missing_genes), ncol = ncol(adj),
#                          dimnames = list(missing_genes, colnames(adj)))
#     adj <- rbind(adj, zero_block)
#     zero_block <- matrix(0, nrow = nrow(adj), ncol = length(missing_genes),
#                          dimnames = list(rownames(adj), missing_genes))
#     adj <- cbind(adj, zero_block)
#   }
#   
#   # Final ordering
#   adj <- adj[all_genes_union, all_genes_union]
#   aligned_networks[[layer]] <- adj
# }
# 
# 
# 
# 
# graph_list <- list()
# node_edge_stats <- data.frame(Layer = character(), Nodes = integer(), Edges = integer(), stringsAsFactors = FALSE)
# 
# for (layer in names(aligned_networks)) {
#   adj <- aligned_networks[[layer]]
#   
#   g <- graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
#   graph_list[[layer]] <- g
#   
#   node_count <- vcount(g)
#   edge_count <- ecount(g)
#   node_edge_stats <- rbind(node_edge_stats, data.frame(Layer = layer, Nodes = node_count, Edges = edge_count))
#   
#   cat(paste0("Layer: ", layer, " | Nodes: ", node_count, " | Edges: ", edge_count, "\n"))
#   
#   print(
#     ggraph(g, layout = 'fr') +
#       geom_edge_link(alpha = 0.3, color = "gray") +
#       geom_node_point(size = 2, color = "steelblue") +
#       ggtitle(paste("Co-expression Network -", layer)) +
#       theme_void()
#   )
# }
# 
# 
# 
# 
# 
# # Directory to save outputs
# dir.create("Layer_Features", showWarnings = FALSE)
# dir.create("Adjacency_Matrices", showWarnings = FALSE)
# dir.create("Edge_Lists", showWarnings = FALSE)
# 
# # 1. Save feature matrix per layer (rows: genes, columns: samples)
# for (layer in names(all_layers)) {
#   expr_mat <- all_layers[[layer]]
#   
#   # Align rows to all_genes_union
#   missing_genes <- setdiff(all_genes_union, rownames(expr_mat))
#   if (length(missing_genes) > 0) {
#     zero_block <- matrix(0, nrow = length(missing_genes), ncol = ncol(expr_mat),
#                          dimnames = list(missing_genes, colnames(expr_mat)))
#     expr_mat <- rbind(expr_mat, zero_block)
#   }
#   
#   # Final ordering
#   expr_mat <- expr_mat[all_genes_union, ]
#   
#   # Save to file
#   write.table(expr_mat,
#               file = paste0("Layer_Features/FeatureMatrix_", layer, ".txt"),
#               sep = "\t", quote = FALSE, col.names = NA)
# }
# 
# 
# # 2. Save aligned adjacency matrices
# for (layer in names(aligned_networks)) {
#   adj_mat <- aligned_networks[[layer]]
#   write.table(adj_mat,
#               file = paste0("Adjacency_Matrices/AdjMatrix_", layer, ".txt"),
#               sep = "\t", quote = FALSE, col.names = NA)
# }
# 
# # 3. Save edge lists from coexp_networks
# for (layer in names(coexp_networks)) {
#   adj_mat <- coexp_networks[[layer]]
#   
#   edge_list <- which(adj_mat == 1, arr.ind = TRUE)
#   edge_df <- data.frame(
#     Source = rownames(adj_mat)[edge_list[, 1]],
#     Target = colnames(adj_mat)[edge_list[, 2]]
#   )
#   
#   # Remove self-loops and duplicate edges (undirected)
#   edge_df <- edge_df[edge_df$Source != edge_df$Target, ]
#   edge_df <- edge_df[!duplicated(t(apply(edge_df, 1, sort))), ]
#   
#   # Save
#   write.table(edge_df,
#               file = paste0("Edge_Lists/EdgeList_", layer, ".txt"),
#               sep = "\t", row.names = FALSE, quote = FALSE)
# }
# 
# cat("All feature matrices, adjacency matrices, and edge lists saved.\n")
# 
# write.table(data.frame(all_genes_union), "Node_List.txt", sep = "\t", row.names = FALSE,col.names = TRUE ,quote = FALSE)
# 
# #write.table(data.frame(all_genes_union,1), "Node_List_1.txt", sep = "\t", row.names = FALSE,col.names = TRUE ,quote = FALSE)
# 
