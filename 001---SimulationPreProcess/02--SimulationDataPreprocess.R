
# Load necessary library
library(igraph)

# Define the full set of nodes
all_nodes <- paste0("G", 1:100)

# Define the file names
input_files <- paste0("DATA/SERGIO-Data/Interaction1_100_", 1:4, ".txt")
output_files <- paste0("RESULTS/SimulationDataPreprocess/adj_matrix_G100_", 1:4, ".txt")

# Process each file
for (i in 1:4) {
  # Read the weighted adjacency matrix
  W <- read.table(input_files[i], header = FALSE)
  colnames(W) <- c("TF", "Gene", "weight")
  
  # Create a directed graph
  g <- graph_from_data_frame(W, directed = TRUE, vertices = data.frame(name = all_nodes))
  
  # Generate adjacency matrix with weights
  adj_matrix <- as.matrix(as_adjacency_matrix(g, attr = "weight", sparse = FALSE))
  
  # Convert to binary adjacency matrix (1 if edge exists, 0 otherwise)
  adj_matrix_pure <- (adj_matrix != 0) * 1
  
  print(adj_matrix_pure[1:10, 1:10])
  
  # Write the binary adjacency matrix to a new file
  write.table(adj_matrix_pure, file = output_files[i], row.names = FALSE, col.names = FALSE)
  
  # Print confirmation
  cat("Processed:", input_files[i], "->", output_files[i], "\n")
}

cat("All files have been processed and saved.\n")



##########################################################


expr_Gene50_StageI <- read.csv("DATA/SERGIO-Data/Gene100_Stage_IS.csv", header = TRUE, row.names = 1)
expr_Gene50_StageII <- read.csv("DATA/SERGIO-Data/Gene100_Stage_IIS.csv", header = TRUE, row.names = 1)
expr_Gene50_StageIII <- read.csv("DATA/SERGIO-Data/Gene100_Stage_IIIS.csv", header = TRUE, row.names = 1)
expr_Gene50_StageIV <- read.csv("DATA/SERGIO-Data/Gene100_Stage_IVS.csv", header = TRUE, row.names = 1)




Log_expr_Gene50_StageI <- log2(expr_Gene50_StageI + 1)
Log_expr_Gene50_StageII <- log2(expr_Gene50_StageII + 1)
Log_expr_Gene50_StageIII <- log2(expr_Gene50_StageIII + 1)
Log_expr_Gene50_StageIV <- log2(expr_Gene50_StageIV + 1)


write.table(Log_expr_Gene50_StageI, "RESULTS/SimulationDataPreprocess/Fea_G100_Stage_1.txt", row.names = FALSE, col.names = FALSE)
write.table(Log_expr_Gene50_StageII, "RESULTS/SimulationDataPreprocess/Fea_G100_Stage_2.txt", row.names = FALSE, col.names = FALSE)
write.table(Log_expr_Gene50_StageIII, "RESULTS/SimulationDataPreprocess/Fea_G100_Stage_3.txt", row.names = FALSE, col.names = FALSE)
write.table(Log_expr_Gene50_StageIV, "RESULTS/SimulationDataPreprocess/Fea_G100_Stage_4.txt", row.names = FALSE, col.names = FALSE)


write.table(data.frame(all_nodes),"RESULTS/SimulationDataPreprocess/all_nodes_G100.txt", row.names = FALSE, col.names = FALSE, quote = F)
