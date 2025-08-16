import numpy as np
import networkx as nx
import os
import pandas as pd

# Step 1: Load the adjacency matrix from a file
cell = 'CD8'  # Graph name (example)
stage = 'NAT'    # State (example) NAT, CAG,  IM, PGAC, Metastasis

base_dir = "/home/fatemeh/ThirdObject/0.RealData/Filtered_Cleaned_Inputs"

net_stage_path = os.path.join(base_dir, f"CMI_Stage_{stage}_{cell}_Adj_UnionGenes.csv")
#adj_matrix = np.loadtxt(f'./Input/Adjacency_Matrices/{GN}_AdjMatrix_{ST}.txt').astype(int)  # Ensure integer values
df = pd.read_csv(net_stage_path, index_col=0)
print(df.shape)
print(df.head())
adj_matrix = df.values.astype(int)  # Convert to NumPy array of integers

# Ensure the adjacency matrix is square
assert adj_matrix.shape[0] == adj_matrix.shape[1], "Adjacency matrix must be square!"

# Step 2: Create a networkx graph from the adjacency matrix
G = nx.from_numpy_matrix(adj_matrix)

# Step 3: Compute Local Structure (Shortest Path Distance)
shortest_paths = dict(nx.all_pairs_shortest_path_length(G))

local_weights = np.zeros((len(G.nodes), len(G.nodes)))
for i in G.nodes:
    for j in G.nodes:
        # For nodes that are disconnected, assign a large number (infinity)
        local_weights[i, j] = shortest_paths[i].get(j, 1000)  # Use 1000 for disconnected nodes

# Normalize the local weights to make them comparable
local_weights = np.exp(-local_weights)  # Use an exponential decay to convert distance to similarity

# Avoid NaN by setting the maximum value of the weights (e.g., 1000) to a reasonable number
local_weights = np.nan_to_num(local_weights, nan=0.0)  # Replace NaNs with 0.0

# Normalize rows
local_weights = local_weights / np.sum(local_weights, axis=1, keepdims=True)

# Check if NaN values appear after local weights calculation
if np.isnan(local_weights).any():
    print("NaN detected in local weights!")

# Step 4: Compute Global Structure using PageRank
pagerank_scores = nx.pagerank(G)  # Compute PageRank

# Convert PageRank scores into a global weight matrix
global_weights = np.array([pagerank_scores[node] for node in G.nodes])
global_weights = np.outer(global_weights, global_weights)  # Create a matrix of global relations

# Normalize the global weights
global_weights = global_weights / np.sum(global_weights, axis=1, keepdims=True)

# Check if NaN values appear after global weights calculation
if np.isnan(global_weights).any():
    print("NaN detected in global weights!")

# Step 5: Combine Local and Global Structures using a weighted sum
alpha = 0.3  # Weighting parameter
combined_weights = alpha * local_weights + (1 - alpha) * global_weights

# Check if NaN values appear in the combined weights
if np.isnan(combined_weights).any():
    print("NaN detected in combined weights!")

# Print out the shapes
print("Local Weights Shape:", local_weights.shape)
print("Global Weights Shape:", global_weights.shape)
print("Combined Weights Shape:", combined_weights.shape)

# Now you can use the combined_weights as the final adjacency matrix


# Now you can use the combined_weights as the final adjacency matrix
print("Local Weights:", local_weights)
print("Global Weights:", global_weights)
print("Combined Weights:", combined_weights)

# Save embeddings
output_dir = "./StructuralEncoding"
file_path_export = os.path.join(output_dir, stage + '_' + 'StructuralEncoding_Stage' + cell + '_Reg.txt')

# Create the directory if it does not exist
os.makedirs(output_dir, exist_ok=True)

# Save embeddings
np.savetxt(file_path_export, combined_weights, fmt='%.6f')

print("Saved StructuralEncoding")

