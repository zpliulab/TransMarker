import os
import pandas as pd
import numpy as np

# Define base directory and parameters
base_dir = "/home/fatemeh/ThirdObject/0.RealData/Filtered_Cleaned_Inputs"
cell = 'CD8'
stages = ['NAT', 'CAG', 'IM', 'PGAC', 'Metastasis']

# Step 1: Collect all genes (union) from edge lists
print("Collecting genes from edge list files...")
all_genes = set()

for stage in stages:
    edge_path = os.path.join(base_dir, f"CMI_Stage_{stage}_{cell}_edges_NoL.csv")
    edges_df = pd.read_csv(edge_path)
    genes_stage = set(edges_df['Gene1']) | set(edges_df['Gene2'])
    all_genes.update(genes_stage)
    print(f"Stage {stage} - Nodes: {len(genes_stage)}, Edges: {len(edges_df)}")

print(f"Total unique genes from all stages: {len(all_genes)}")
sorted_genes = sorted(all_genes)  # For consistent order
print(f"unique genes: {len(sorted_genes)}")

# ✅ Save union genes
union_genes_path = os.path.join(base_dir, f"Union_Genes_{cell}.csv")
pd.Series(sorted_genes).to_csv(union_genes_path, index=False, header=['Gene'])
print(f"Union genes saved to: {union_genes_path}")


# Step 2: Filter expression matrices and save filtered versions
for stage in stages:
    exp_path = os.path.join(base_dir, f"{stage}_{cell}T_M_cleaned_Filtered.csv")
    exp_df = pd.read_csv(exp_path, index_col=0)

    # Filter expression matrix
    filtered_exp_df = exp_df.loc[exp_df.index.intersection(all_genes)]

    # Save filtered matrix
    output_path = os.path.join(base_dir, f"{stage}_{cell}T_M_cleaned_Filtered_UnionGenes.csv")
    filtered_exp_df.to_csv(output_path)

    print(f"Filtered expression matrix saved for stage {stage}: {filtered_exp_df.shape[0]} genes")

# Step 3: Build and save square adjacency matrices for all stages
print("\nCreating square adjacency matrices...")

gene_index = {gene: idx for idx, gene in enumerate(sorted_genes)}

for stage in stages:
    edge_path = os.path.join(base_dir, f"CMI_Stage_{stage}_{cell}_edges_NoL.csv")
    edges_df = pd.read_csv(edge_path)

    # Initialize zero matrix
    size = len(sorted_genes)
    adj_matrix = np.zeros((size, size))

    for _, row in edges_df.iterrows():
        g1, g2, weight = row['Gene1'], row['Gene2'], row.get('Weight', 1)
        if g1 in gene_index and g2 in gene_index:
            i, j = gene_index[g1], gene_index[g2]
            adj_matrix[i, j] = weight
            adj_matrix[j, i] = weight  # Symmetric (undirected)

    # Save matrix
    adj_output_path = os.path.join(base_dir, f"CMI_Stage_{stage}_{cell}_Adj_UnionGenes.csv")
    adj_df = pd.DataFrame(adj_matrix, index=sorted_genes, columns=sorted_genes)
    adj_df.to_csv(adj_output_path)

    print(f"Adjacency matrix saved for stage {stage}: {adj_output_path}")

print("\n✅ All adjacency matrices generated and saved using union genes.")


# Step 4: Read RegNetwork and build adjacency matrix using union genes
print("\nProcessing RegNetwork edge list...")

# Read RegNetwork edge list (columns 1 and 3 are 0-based columns 0 and 2 in Python)
regnet_path = "/home/fatemeh/ThirdObject/0.RealData/RegNetwork-2022.human.source"
regnet_df = pd.read_csv(regnet_path, sep="\t", header=None, usecols=[0, 2], names=["Regulator", "Target"])

# Initialize zero matrix
size = len(sorted_genes)
regnet_adj_matrix = np.zeros((size, size))

# Build index mapping
gene_index = {gene: idx for idx, gene in enumerate(sorted_genes)}

# Fill in adjacency matrix
edge_count = 0
for _, row in regnet_df.iterrows():
    g1, g2 = row["Regulator"], row["Target"]
    if g1 in gene_index and g2 in gene_index:
        i, j = gene_index[g1], gene_index[g2]
        regnet_adj_matrix[i, j] = 1  # Directed edge (Regulator → Target)
        edge_count += 1

print(f"Filtered RegNetwork edges retained: {edge_count}")

# Save the RegNetwork adjacency matrix
regnet_adj_df = pd.DataFrame(regnet_adj_matrix, index=sorted_genes, columns=sorted_genes)
regnet_adj_path = os.path.join(base_dir, f"RegNetwork_Adj_UnionGenes_"+cell+".csv")
regnet_adj_df.to_csv(regnet_adj_path)
print(f"✅ RegNetwork adjacency matrix saved: {regnet_adj_path}")
'''

import os
import pandas as pd

# Define base directory and parameters
base_dir = "/home/fatemeh/ThirdObject/0.RealData/Filtered_Cleaned_Inputs"
cell = 'CD8'
stages = ['NAT', 'CAG', 'IM', 'PGAC', 'Metastasis']

# Step 1: Collect all genes (union) from edge lists
print("Collecting genes from edge list files...")
all_genes = set()

for stage in stages:
    edge_path = os.path.join(base_dir, f"CMI_Stage_{stage}_{cell}_edges_NoL.csv")
    edges_df = pd.read_csv(edge_path)
    genes_stage = set(edges_df['Gene1']) | set(edges_df['Gene2'])
    all_genes.update(genes_stage)
    # Print number of nodes and edges for the current stage
    print(f"Stage {stage} - Nodes: {len(genes_stage)}, Edges: {len(edges_df)}")

print(f"Total unique genes from all stages: {len(all_genes)}")

# ✅ Save union of all genes
union_genes_path = os.path.join(base_dir, f"Union_Genes_{cell}.csv")
pd.Series(sorted(all_genes)).to_csv(union_genes_path, index=False, header=['Gene'])
print(f"Union genes saved to: {union_genes_path}")

# Step 2: Filter expression matrices and save filtered versions
for stage in stages:
    exp_path = os.path.join(base_dir, f"{stage}_{cell}T_M_cleaned_Filtered.csv")
    exp_df = pd.read_csv(exp_path, index_col=0)

    # Filter expression matrix
    filtered_exp_df = exp_df.loc[exp_df.index.intersection(all_genes)]

    # Save filtered matrix
    output_path = os.path.join(base_dir, f"{stage}_{cell}T_M_cleaned_Filtered_UnionGenes.csv")
    filtered_exp_df.to_csv(output_path)

    print(f"Filtered expression matrix saved for stage {stage}: {filtered_exp_df.shape[0]} genes")

print("All expression matrices processed and saved.")
'''