import pandas as pd
import os
from PCA_CMI import pca_cmi
import numpy as np
import networkx as nx

cell = 'CD8'
stage = 'Metastasis'  # 'NAT', 'CAG', 'IM', 'PGAC', 'Metastasis'

# Define the base directory
base_dir = "/home/fatemeh/ThirdObject/0.RealData/Filtered_Cleaned_Inputs"

# Construct file paths dynamically
net_stage1_path = os.path.join(base_dir, f"filtered_RegNet_{cell}.csv")
exp_stage1_path = os.path.join(base_dir, f"{stage}_{cell}T_M_cleaned_Filtered.csv")

# Read the CSV files
Net_Stage1 = pd.read_csv(net_stage1_path)
Exp_Stage1 = pd.read_csv(exp_stage1_path, index_col=0)

# Display structure and preview
print(Net_Stage1.head())
print(Net_Stage1.info())
print(Exp_Stage1.head())
print(Exp_Stage1.info())

# Run PCA-CMI
max_order = 1
show = False
node_feature = Exp_Stage1
net_bit = Net_Stage1

binary_adjMatrix_pmi, predicted_graph_pmi = pca_cmi(node_feature, net_bit, theta=0.1,
                                                    max_order=max_order, L=-1, show=show)

print('PCA-CMI finished.')
binary_adjMatrix_pmi = binary_adjMatrix_pmi.toarray()

# ✅ Remove self-loops
predicted_graph_pmi.remove_edges_from(nx.selfloop_edges(predicted_graph_pmi))
print("Self-loops removed.")

# Convert to adjacency matrix
adj_matrix = nx.to_numpy_array(predicted_graph_pmi)

# Save adjacency matrix
output_file = os.path.join(base_dir, f"CMI_Stage_{stage}_{cell}"+"NoL"+".txt")
np.savetxt(output_file, adj_matrix, fmt="%.2f", delimiter=",")
print(f"Adjacency matrix saved successfully at: {output_file}")

# Extract and save edge list
edges_with_weights = []
for gene1, gene2, data in predicted_graph_pmi.edges(data=True):
    weight = data.get('weight', 1)  # Default weight = 1
    edges_with_weights.append((gene1, gene2, weight))

edges_df = pd.DataFrame(edges_with_weights, columns=["Gene1", "Gene2", "Weight"])
edges_file = os.path.join(base_dir, f"CMI_Stage_{stage}_{cell}_edges_NoL.csv")
edges_df.to_csv(edges_file, index=False, quotechar='"')
print(f"Edge list with weights saved successfully at: {edges_file}")
