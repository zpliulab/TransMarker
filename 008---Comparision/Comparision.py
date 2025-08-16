import pandas as pd
import numpy as np
import networkx as nx
import os
from collections import Counter
from collections import defaultdict

#First Method
# ---- Paper-based Method 1: Entropy-Based Centrality for Multilayer Networks ---- #
def multilayer_entropy_centrality(graphs, gene_list):
    from collections import defaultdict
    import numpy as np

    L = len(graphs)
    entropy_layers = defaultdict(float)

    # Step 1: Compute betweenness centrality and degree for each layer
    centralities = []
    degrees = []
    for G in graphs:
        bc = nx.betweenness_centrality(G)
        deg = dict(G.degree())
        centralities.append(bc)
        degrees.append(deg)

    for gene in gene_list:
        total_entropy = 0
        for l in range(L):
            bc_l = centralities[l].get(gene, 0)
            deg_l = degrees[l].get(gene, 0)
            score = bc_l + deg_l

            # Prepare probability distribution over neighbors
            G = graphs[l]
            neighbors = list(G.neighbors(gene)) if G.has_node(gene) else []
            p_vals = []

            for nbr in neighbors:
                nbr_bc = centralities[l].get(nbr, 0)
                nbr_deg = degrees[l].get(nbr, 0)
                val = nbr_bc + nbr_deg
                p_vals.append(val)

            p_sum = sum(p_vals)
            if p_sum == 0:
                continue
            probs = [p / p_sum for p in p_vals]
            entropy = -sum(p * np.log(p + 1e-10) for p in probs)
            total_entropy += entropy

        entropy_layers[gene] = total_entropy

    return entropy_layers


#Secound method
# Add ground node to ensure strong connectivity
def add_ground_node(G):
    G = G.copy().to_directed()
    ground_node = "__ground__"
    G.add_node(ground_node)

    in_deg = dict(G.in_degree(weight='weight'))
    for n in G.nodes():
        if n != ground_node:
            G.add_edge(ground_node, n, weight=1)
            G.add_edge(n, ground_node, weight=in_deg.get(n, 1))
    return G


# Weighted PageRank on a single layer with ground node
def weighted_pagerank_layer(G, alpha=0.85, max_iter=100, tol=1e-6):
    G = add_ground_node(G)
    pr = nx.pagerank(G, alpha=alpha, max_iter=max_iter, tol=tol, weight='weight')
    pr.pop("__ground__", None)
    return pr


# Compute AHP weights (assuming ordinal layer weights)
def compute_ahp_weights(num_layers):
    # Simple geometric scale (ordinal weights): 1, 2, ..., num_layers
    scores = np.arange(1, num_layers + 1, dtype=float)
    weights = scores / scores.sum()
    return weights


# ElementRank: apply weighted PageRank + AHP weighting
def elementrank(graphs):
    num_layers = len(graphs)
    pagerank_per_layer = []

    for G in graphs:
        # Convert edge weight if not present
        if not nx.get_edge_attributes(G, "weight"):
            for u, v in G.edges():
                G[u][v]['weight'] = 1.0
        pr = weighted_pagerank_layer(G)
        pagerank_per_layer.append(pr)

    # Get consistent set of nodes across all layers
    all_nodes = set(pagerank_per_layer[0].keys())
    for pr in pagerank_per_layer[1:]:
        all_nodes &= set(pr.keys())

    ahp_weights = compute_ahp_weights(num_layers)

    # Aggregate using AHP weights
    final_scores = defaultdict(float)
    for node in all_nodes:
        for l in range(num_layers):
            final_scores[node] += ahp_weights[l] * pagerank_per_layer[l].get(node, 0.0)

    return final_scores

#Third Method
def versatility_centrality(graphs, gene_list):
    """
    Compute node versatility based on eigenvector centrality across multiple layers.
    Versatility is defined as the L2 norm of centrality scores across layers.

    Args:
        graphs (List[nx.Graph]): List of NetworkX graphs, one per layer.
        gene_list (List[str]): List of consistent node names.

    Returns:
        dict: Node -> versatility score
    """
    from collections import defaultdict
    import numpy as np

    # Step 1: Compute eigenvector centrality per layer
    centrality_per_layer = []
    for G in graphs:
        try:
            centrality = nx.eigenvector_centrality_numpy(G)
        except nx.NetworkXException:
            # Fallback for disconnected graphs
            centrality = {node: 0 for node in G.nodes()}
        centrality_per_layer.append(centrality)

    # Step 2: Collect centrality vectors per node
    node_vectors = defaultdict(list)
    for centrality in centrality_per_layer:
        for node in gene_list:
            node_vectors[node].append(centrality.get(node, 0))

    # Step 3: Compute versatility (L2 norm)
    versatility = {node: np.linalg.norm(scores) for node, scores in node_vectors.items()}
    return versatility


#Forth method
def versatility_degree_centrality(graphs, gene_list):
    """
    Compute versatility score for each node based on degree centrality across multiple layers.

    Args:
        graphs (List[nx.Graph]): List of NetworkX graphs (layers).
        gene_list (List[str]): List of all node names.

    Returns:
        dict: Node -> versatility score
    """
    import numpy as np
    from collections import defaultdict

    total_degree = defaultdict(float)
    layer_count = defaultdict(int)

    for G in graphs:
        degrees = dict(G.degree())
        for gene in gene_list:
            deg = degrees.get(gene, 0)
            if deg > 0:
                total_degree[gene] += deg
                layer_count[gene] += 1

    versatility_score = {}
    for gene in gene_list:
        if layer_count[gene] > 0:
            versatility_score[gene] = total_degree[gene] / layer_count[gene]
        else:
            versatility_score[gene] = 0.0

    return versatility_score

#Fifth method
import numpy as np
import networkx as nx


def eigenvector_multicentrality(graphs, max_iter=100, tol=1e-6):
    """
    Compute eigenvector multicentrality scores for multilayer network.

    Parameters:
        graphs: list of nx.Graph objects (layers), all with same nodes
        max_iter: number of iterations for power method
        tol: tolerance for convergence

    Returns:
        dict: node -> multicentrality score
    """
    num_layers = len(graphs)
    nodes = sorted(graphs[0].nodes())
    num_nodes = len(nodes)
    node_idx = {node: i for i, node in enumerate(nodes)}

    # Build adjacency tensor A[i, l, j, m] where i/j are node idx, l/m are layers
    A = np.zeros((num_nodes, num_layers, num_nodes, num_layers))

    for l, G in enumerate(graphs):
        for u, v in G.edges():
            i, j = node_idx[u], node_idx[v]
            A[i, l, j, l] = 1  # intra-layer edges

    # Add interlayer links between the same node across layers (weight=1)
    for i in range(num_nodes):
        for l in range(num_layers):
            for m in range(num_layers):
                if l != m:
                    A[i, l, i, m] = 1

    # Flatten tensor to 2D matrix for power iteration
    N = num_nodes * num_layers
    A_matrix = A.reshape((N, N))

    # Power iteration
    x = np.ones(N)
    x /= np.linalg.norm(x)

    for _ in range(max_iter):
        x_new = A_matrix @ x
        x_new /= np.linalg.norm(x_new)
        if np.linalg.norm(x_new - x) < tol:
            break
        x = x_new

    # Aggregate scores across layers for each node
    x_reshaped = x.reshape((num_nodes, num_layers))
    node_scores = x_reshaped.sum(axis=1)

    return {node: score for node, score in zip(nodes, node_scores)}


'''
centrality_methods = {
    "multilayer_entropy": lambda Gs: multilayer_entropy_centrality(Gs, gene_list),
    "elementrank": lambda Gs: elementrank(Gs)
}
'''

# Load a single network from an adjacency matrix CSV file
def load_network(adj_path, gene_list):
    adj_df = pd.read_csv(adj_path, index_col=0)
    adj_df = adj_df.loc[gene_list, gene_list]
    G = nx.from_pandas_adjacency(adj_df)
    return G
'''
# Rank nodes in each graph using the selected method
def get_ranked_nodes(graphs, method_func, top_k=20):
    rankings = []
    for G in graphs:
        scores = method_func(G)
        ranked = sorted(scores.items(), key=lambda x: -x[1])
        rankings.append([n for n, _ in ranked[:top_k]])
    return rankings


def get_ranked_nodes(graphs, method_func, top_k=30):
    rankings = []
    for G in graphs:
        scores = method_func(G)
        sorted_scores = sorted(scores.items(), key=lambda x: -x[1])
        if len(sorted_scores) <= top_k:
            rankings.append([n for n, _ in sorted_scores])
        else:
            threshold_score = sorted_scores[top_k - 1][1]
            top_nodes = [n for n, s in sorted_scores if s >= threshold_score]
            rankings.append(top_nodes)
    return rankings
'''
'''
def get_ranked_nodes(graphs, method_func, top_k=30):
    rankings = []
    for G in graphs:
        scores = method_func(G)
        # Sort nodes by descending score
        ranked = sorted(scores.items(), key=lambda x: -x[1])
        # Keep only the first top_k nodes
        top_nodes = [n for n, _ in ranked[:top_k]]
        rankings.append(top_nodes)
    return rankings
'''


def get_ranked_nodes(graphs, method_func, top_k=30):
    # If method is multilayer-based, use all graphs at once
    try:
        scores = method_func(graphs)
        ranked = sorted(scores.items(), key=lambda x: -x[1])
        return [[n for n, _ in ranked[:top_k]]]
    except:
        rankings = []
        for G in graphs:
            scores = method_func(G)
            ranked = sorted(scores.items(), key=lambda x: -x[1])
            rankings.append([n for n, _ in ranked[:top_k]])
        return rankings


'''
# Aggregate top-k results across layers
def aggregate_rankings(rank_lists, top_k=20):
    count = Counter()
    for ranks in rank_lists:
        for i, node in enumerate(ranks[:top_k]):
            count[node] += (top_k - i)  # weighted rank
    return [n for n, _ in count.most_common()]
'''

def aggregate_rankings(rank_lists, top_k=30):
    count = Counter()
    for ranks in rank_lists:
        for i, node in enumerate(ranks[:top_k]):
            count[node] += (top_k - i)  # weighted scoring
    top_combined = count.most_common(top_k)
    return [n for n, _ in top_combined]


# ---- Main execution ---- #
def run_centrality_pipeline(cell, stages, GN, top_k=30):
    base_dir = "/home/fatemeh/ThirdObject/0.RealData/Filtered_Cleaned_Inputs"
    gene_path = f"/home/fatemeh/ThirdObject//0.RealData/Filtered_Cleaned_Inputs/Union_Genes_"+GN+".csv"
    gene_list = pd.read_csv(gene_path, header=0, names=["Gene"])["Gene"].values

    graphs = []
    for stage in stages:
        adj_path = os.path.join(base_dir, f"CMI_Stage_{stage}_{cell}_Adj_UnionGenes.csv")
        print(f"Loading: {adj_path}")
        G = load_network(adj_path, gene_list)
        graphs.append(G)

    # Define centrality methods (including the one from paper) here AFTER gene_list is available
    centrality_methods = {
        "multilayer_entropy": lambda Gs: multilayer_entropy_centrality(Gs, gene_list),
        "elementrank": lambda Gs: elementrank(Gs),
         "versatility": lambda Gs: versatility_centrality(Gs, gene_list),
         "versatility_deg": lambda Gs: versatility_degree_centrality(Gs, gene_list),
         "eigenvector_multicentrality":  lambda Gs: eigenvector_multicentrality(Gs)

    }

    results = {}
    for method, func in centrality_methods.items():
        print(f"Running {method} centrality...")
        ranks = get_ranked_nodes(graphs, func, top_k)
        results[method] = aggregate_rankings(ranks, top_k)

    return results


# Example call
if __name__ == "__main__":
    cell = "CD8"
    GN = "CD8"
    stages = ["NAT", "CAG", "IM", "PGAC", "Metastasis"]
    result_dict = run_centrality_pipeline(cell, stages, GN, top_k=30)

    # Save to CSV
    out_path = f"/home/fatemeh/ThirdObject/7.RankingTest/Comparison_methods_rankings_{cell}.csv"
    df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in result_dict.items()]))
    df.to_csv(out_path, index=False)
    print(f"Combined rankings saved to {out_path}")

for method, nodes in result_dict.items():
    print(f"{method}: {len(nodes)} nodes")


