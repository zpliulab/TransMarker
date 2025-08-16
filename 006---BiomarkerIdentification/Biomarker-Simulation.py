import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from kneed import KneeLocator
import math



GN = 'G50'


# Paths
path = '/home/fatemeh/ThirdObject'

# Load gene names and alignment scores
gene_list = pd.read_csv('./Input/all_nodes_'+GN+'.txt', header=None, names=["Gene"])  # Assuming gene names
GW_scores = np.loadtxt(path + '/3.GWAlignment/gw_cumulative_alignment_'+GN+'.csv', delimiter=",")
alignment_scores = np.array(GW_scores)
gene_names = gene_list["Gene"].values

print(alignment_scores)
print(gene_names)

# Ensure same length
assert len(alignment_scores) == len(gene_names), "Mismatch between gene names and alignment scores"

# ----------------------------
# Filtering Option 1: Top-k%
# ----------------------------
def filter_top_k_percent(alignment_scores, gene_names, percent=1):
    k = math.ceil(len(alignment_scores) * percent / 100)
    sorted_indices = np.argsort(alignment_scores)[::-1]
    top_k_indices = sorted_indices[:k]
    return gene_names[top_k_indices], alignment_scores[top_k_indices]


# ----------------------------
# Filtering Option 2: Mean + N*Std
# ----------------------------
def filter_stat_threshold(alignment_scores, gene_names, std_multiplier=2):
    threshold = np.mean(alignment_scores) + std_multiplier * np.std(alignment_scores)
    indices = np.where(alignment_scores >= threshold)[0]
    return gene_names[indices], alignment_scores[indices]


# ----------------------------
# Filtering Option 3: Elbow Method
# ----------------------------
def filter_elbow_method(alignment_scores, gene_names, plot=True):
    sorted_scores = np.sort(alignment_scores)
    indices_sorted = np.argsort(alignment_scores)

    x = np.arange(len(sorted_scores))
    kneedle = KneeLocator(x, sorted_scores, curve="convex", direction="increasing")
    elbow_idx = kneedle.knee

    if plot:
        plt.figure(figsize=(8, 5))
        plt.plot(x, sorted_scores, label='Sorted Alignment Scores')
        if elbow_idx is not None:
            plt.axvline(elbow_idx, color='red', linestyle='--', label=f'Elbow Point: {elbow_idx}')
        plt.title("Elbow Method on Alignment Scores")
        plt.xlabel("Ranked Genes")
        plt.ylabel("Alignment Score")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.savefig('elbow_plot_'+GN+'.pdf')  # Save as PDF
        plt.savefig('elbow_plot_'+GN+'.svg')  # Save as SVG
        plt.show()


    if elbow_idx is None:
        print("Elbow point not detected. Returning top 1% as fallback.")
        return filter_top_k_percent(alignment_scores, gene_names, percent=1)

    selected_indices = indices_sorted[elbow_idx:]
    return gene_names[selected_indices], alignment_scores[selected_indices]


# Option 1: Top 1%
genes_top_k, scores_top_k = filter_top_k_percent(alignment_scores, gene_names, percent=30)
print(genes_top_k)
print(scores_top_k)



# Option 2: Mean + 2*STD
genes_stat, scores_stat = filter_stat_threshold(alignment_scores, gene_names, std_multiplier=2)
print(genes_stat)
print(scores_stat)

# Option 3: Elbow method (with plot)
genes_elbow, scores_elbow = filter_elbow_method(alignment_scores, gene_names, plot=True)
print(genes_elbow)
print(scores_elbow)




# ----------------------------
# Construct Union Network of Candidate Genes
# ----------------------------

import networkx as nx

# Load adjacency matrix instead of edge list
A1 = np.loadtxt(path+'/2.GraphEmbedding/Input/adj_matrix_'+GN+'_1.txt').astype(int)
A2 = np.loadtxt(path+'/2.GraphEmbedding/Input/adj_matrix_'+GN+'_2.txt').astype(int)
A3 = np.loadtxt(path+'/2.GraphEmbedding/Input/adj_matrix_'+GN+'_3.txt').astype(int)
A4 = np.loadtxt(path+'/2.GraphEmbedding/Input/adj_matrix_'+GN+'_4.txt').astype(int)




def build_union_network(candidate_genes, gene_names, adjacency_matrices):
    gene_to_idx = {gene: idx for idx, gene in enumerate(gene_names)}
    candidate_indices = [gene_to_idx[gene] for gene in candidate_genes if gene in gene_to_idx]

    union_graph = nx.Graph()
    union_graph.add_nodes_from(candidate_genes)  # Add gene nodes by name

    for A in adjacency_matrices:
        for i in candidate_indices:
            for j in candidate_indices:
                if i != j and A[i, j] == 1:
                    gene_i = gene_names[i]
                    gene_j = gene_names[j]
                    union_graph.add_edge(gene_i, gene_j)

    return union_graph

adj_matrices = [A1, A2, A3, A4]

# Use any of the selected gene sets:
candidate_genes = genes_elbow  # or genes_stat or genes_top_k
#candidate_genes = genes_top_k
#candidate_genes = genes_stat


'''
# Save the node list of the best component to a text file
with open('best_component_nodes_' + GN + '.txt', 'w') as f:
    for node in candidate_genes:
        f.write(f"{node}\n")

print(candidate_genes)

'''

union_net = build_union_network(candidate_genes, gene_names, adj_matrices)

# Print basic stats
print(f"Union Network has {union_net.number_of_nodes()} nodes and {union_net.number_of_edges()} edges")



# ----------------------------
# Visualize Union Network
# ----------------------------
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

def visualize_union_network(graph, output_prefix="union_network"):
    plt.figure(figsize=(5, 5))

    # Generate spring layout positions
    pos = nx.spring_layout(graph, iterations=40)

    # Draw nodes
    degrees = dict(graph.degree())
    node_sizes = [700 + 5 * degrees[n] for n in graph.nodes()]
    node_colors = [degrees[n] for n in graph.nodes()]

    nx.draw_networkx_nodes(graph, pos, node_size=node_sizes, node_color=node_colors,
                           cmap=plt.cm.viridis, alpha=0.6)

    # Draw edges using LineCollection manually
    edge_segments = [(pos[u], pos[v]) for u, v in graph.edges()]
    edge_collection = LineCollection(edge_segments, linewidths=0.6, alpha=0.6, color='gray')
    ax = plt.gca()
    ax.add_collection(edge_collection)

    # Draw labels
    nx.draw_networkx_labels(graph, pos, font_size=10, font_color='black')

    plt.title("Union Network of Candidate Genes", fontsize=14)
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(f"{output_prefix+'_'+GN}.pdf")
    plt.savefig(f"{output_prefix+'_'+GN}.svg")
    plt.show()



#visualize_union_network(union_net, output_prefix="union_network_elbow")





# ----------------------------
# Extract Connected Components
# ----------------------------

components = list(nx.connected_components(union_net))
filtered_components = [c for c in components if len(c) > 1]  # Filter out components with just one node
component_subgraphs = [union_net.subgraph(c).copy() for c in filtered_components]

print(f"Found {len(components)} connected components.")
for i, comp in enumerate(component_subgraphs):
    print(f"Component {i+1}: {comp.number_of_nodes()} nodes, {comp.number_of_edges()} edges")

#visualize_union_network(component_subgraphs[0], output_prefix="union_network_elbow")



# ----------------------------
# Evaluate Components with DN Index
# ----------------------------

# DN index
def compute_module_cohesion_density(subgraph):
    n = subgraph.number_of_nodes()
    m = subgraph.number_of_edges()
    if n <= 1:
        return 0
    return (2 * m) / (n * (n - 1))  # edge density

###????
def compute_alignment_std(component_nodes, alignment_dict):
    scores = [alignment_dict[node] for node in component_nodes if node in alignment_dict]
    if len(scores) == 0:
        return 1  # maximum uncertainty if no scores
    return np.std(scores)

def compute_alignment_sum(component_nodes, alignment_dict):
    scores = [alignment_dict[node] for node in component_nodes if node in alignment_dict]

    if len(scores) == 0:
        return 0  # No alignment scores found, return 0

    return np.sum(scores)  # Return the summation of the alignment scores

def compute_dni(subgraph, alignment_scores):
    alignment_dict = {gene: score for gene, score in zip(gene_names, alignment_scores)}
    nodes = list(subgraph.nodes())

    std_score = compute_alignment_std(nodes, alignment_dict)
    #std_score = compute_alignment_sum(nodes, alignment_dict)
    cohesion = compute_module_cohesion_density(subgraph)
    print(std_score)

    gcs = (1 - std_score) * cohesion
    dni = np.exp(-gcs)

    return dni, gcs, std_score, cohesion


component_scores = []
for i, sg in enumerate(component_subgraphs):
    dni, gcs, std_val, coh = compute_dni(sg, alignment_scores)
    component_scores.append((i, dni, gcs, std_val, coh, sg))

# Sort by DNI (lower is better)
sorted_components = sorted(component_scores, key=lambda x: x[1])

# Print top results
for i, (idx, dni, gcs, std_val, coh, sg) in enumerate(sorted_components[:3]):
    print(f"Component #{idx+1} — DNI: {dni:.4f}, GCS: {gcs:.4f}, STD: {std_val:.4f}, Cohesion: {coh:.4f}")

# Extract node names of best component (based on DNI)
best_component_idx = sorted_components[0][0]  # Get the index of the best component
best_component = component_subgraphs[best_component_idx]
best_component_nodes = list(best_component.nodes())

print(f"Best component nodes: {best_component_nodes}")

visualize_union_network(best_component, output_prefix="best_component_network")


nx.write_edgelist(best_component, 'best_component_edges_'+GN+'.txt')

# Save the node list of the best component to a text file
with open('best_component_nodes_' + GN + '.txt', 'w') as f:
    for node in best_component_nodes:
        f.write(f"{node}\n")

print(f"Node list for the best component saved to 'best_component_nodes_{GN}.txt'")

