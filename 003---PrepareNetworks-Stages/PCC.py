
import os
import pandas as pd
import numpy as np
import networkx as nx
from scipy.stats import pearsonr

# Parameters
cell = 'CD8'
stage = 'Metastasis'  # 'NAT' 'CAG' 'IM' 'PGAC' 'Metastasis'
base_dir = "/home/fatemeh/ThirdObject/0.RealData/Filtered_Cleaned_Inputs"

# File paths
net_stage1_path = os.path.join(base_dir, f"CMI_Stage_{stage}_{cell}_edges.csv")
exp_stage1_path = os.path.join(base_dir, f"{stage}_{cell}T_M_cleaned_Filtered.csv")

# Read data
print("Reading input files...")
net_df = pd.read_csv(net_stage1_path)
exp_df = pd.read_csv(exp_stage1_path, index_col=0)

print(f"Network edges: {net_df.shape[0]}")
print(f"Expression matrix: {exp_df.shape}")

# Create graph
print("Building graph...")
G = nx.from_pandas_edgelist(net_df, source='Gene1', target='Gene2')

# Compute Pearson correlation and assign as weights
print("Calculating Pearson correlation for each edge...")
for u, v in G.edges():
    if u in exp_df.index and v in exp_df.index:
        expr_u = exp_df.loc[u]
        expr_v = exp_df.loc[v]
        try:
            # Compute PCC
            corr, _ = pearsonr(expr_u, expr_v)
            G[u][v]['weight'] = round(corr, 4)
        except ValueError:
            # In case of PCC calculation error, assign weight 0
            G[u][v]['weight'] = 0
            print(f"Edge ({u}, {v}) failed to compute PCC, assigning weight 0")
    else:
        # If gene is missing in the expression data, assign weight 0
        G[u][v]['weight'] = 0
        print(f"Edge ({u}, {v}) skipped due to missing gene data, assigning weight 0")

# Print final edge count (should be the same as original)
print(f"Total edges in the graph after adding weights: {G.number_of_edges()}")

# Save weighted edge list
edge_list_df = nx.to_pandas_edgelist(G)
edge_list_path = os.path.join(base_dir, f"Weighted_EdgeList_Stage_{stage}_{cell}.csv")
edge_list_df.to_csv(edge_list_path, index=False)
print(f"Weighted edge list saved to: {edge_list_path}")

# Save weighted adjacency matrix
adj_matrix = nx.to_numpy_array(G, nodelist=sorted(G.nodes()), weight='weight')
adj_df = pd.DataFrame(adj_matrix, index=sorted(G.nodes()), columns=sorted(G.nodes()))
adj_matrix_path = os.path.join(base_dir, f"Weighted_AdjMatrix_Stage_{stage}_{cell}.csv")
adj_df.to_csv(adj_matrix_path)
print(f"Weighted adjacency matrix saved to: {adj_matrix_path}")



