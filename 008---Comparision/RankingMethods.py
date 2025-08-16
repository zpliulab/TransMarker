import pandas as pd
import numpy as np
import networkx as nx
import os
from collections import Counter

# Define all centrality methods

def betweenness_centrality(G): return nx.betweenness_centrality(G)
def bottleneck_centrality(G): return nx.betweenness_centrality(G)  # Approximate
def degree_centrality(G): return nx.degree_centrality(G)
def diffusion_degree(G): return {n: sum(1 / G.degree(nb) for nb in G.neighbors(n)) for n in G}
def latora_closeness(G):
    paths = dict(nx.all_pairs_shortest_path_length(G))
    return {n: sum(1 / d for t, d in paths[n].items() if d > 0) / (len(G) - 1) for n in G}
'''
def lin_centrality(G):
    ecc = nx.eccentricity(G)
    return {n: G.degree(n) / (1 + ecc[n]) for n in G}
'''

def lin_centrality(G):
    centrality = {}
    if not nx.is_connected(G):
        # compute on largest connected component
        largest_cc = max(nx.connected_components(G), key=len)
        G_cc = G.subgraph(largest_cc).copy()
        ecc = nx.eccentricity(G_cc)
        for n in G.nodes():
            centrality[n] = 1 / ecc[n] if n in ecc else 0
    else:
        ecc = nx.eccentricity(G)
        centrality = {n: 1 / e if e > 0 else 0 for n, e in ecc.items()}
    return centrality
'''
def laplacian_centrality(G):
    base = sum(nx.laplacian_matrix(G).todense().diagonal())
    return {n: base - sum(nx.laplacian_matrix(nx.restricted_view(G, {n}, {})).todense().diagonal()) for n in G}
'''

def laplacian_centrality(G):
    import numpy as np
    from networkx.linalg.laplacianmatrix import laplacian_matrix

    # Base Laplacian energy of full graph
    L = laplacian_matrix(G).todense()
    base_energy = np.sum(np.square(np.linalg.eigvals(L)).real)

    centrality = {}
    for n in G.nodes():
        G_minus_n = G.copy()
        G_minus_n.remove_node(n)
        if G_minus_n.number_of_nodes() == 0:
            centrality[n] = 0
            continue
        L_minus = laplacian_matrix(G_minus_n).todense()
        energy_minus = np.sum(np.square(np.linalg.eigvals(L_minus)).real)
        centrality[n] = (base_energy - energy_minus) / base_energy if base_energy != 0 else 0

    return centrality


def local_centrality(G):
    return {n: G.degree(n) + sum(G.degree(nb) for nb in G.neighbors(n)) for n in G}
def leaderrank(G):
    G = G.to_directed()
    G.add_node('ground')
    for n in G.nodes():
        if n != 'ground':
            G.add_edge('ground', n)
            G.add_edge(n, 'ground')
    score = {n: 1.0 for n in G.nodes()}
    for _ in range(50):
        new_score = {n: 0 for n in G.nodes()}
        for n in G:
            if G.out_degree(n) == 0: continue
            share = score[n] / G.out_degree(n)
            for nb in G.successors(n):
                new_score[nb] += share
        score = new_score
    score.pop('ground')
    return score
def leverage_centrality(G):
    return {n: (G.degree(n) - np.mean([G.degree(nb) for nb in G.neighbors(n)])) /
               (G.degree(n) + np.mean([G.degree(nb) for nb in G.neighbors(n)]) + 1e-6)
            for n in G if G.degree(n) > 0}
def residual_closeness(G):
    close = nx.closeness_centrality(G)
    return {n: close[n] / G.degree(n) if G.degree(n) > 0 else 0 for n in G}
def radiality_centrality(G):
    sp = dict(nx.all_pairs_shortest_path_length(G))
    max_dist = max(max(d.values()) for d in sp.values())
    return {n: sum(max_dist + 1 - d for t, d in sp[n].items() if t != n) / (len(G) - 1) for n in G}
def pagerank_centrality(G): return nx.pagerank(G)

centrality_methods = {
    "betweenness": betweenness_centrality,
    "bottleneck": bottleneck_centrality,
    "degree": degree_centrality,
    "diffusion": diffusion_degree,
    "latora": latora_closeness,
    "lin": lin_centrality,
    "laplacian": laplacian_centrality,
    "local": local_centrality,
    "leaderrank": leaderrank,
    "leverage": leverage_centrality,
    "residual": residual_closeness,
    "radiality": radiality_centrality,
    "pagerank": pagerank_centrality
}

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
    #gene_list = pd.read_csv(gene_path, header=None, names=["Gene"])["Gene"].tolist()
    gene_list = pd.read_csv(gene_path, header=0,
                            names=["Gene"])
    gene_list = gene_list["Gene"].values

    graphs = []
    for stage in stages:
        adj_path = os.path.join(base_dir, f"CMI_Stage_{stage}_{cell}_Adj_UnionGenes.csv")
        print(f"Loading: {adj_path}")
        G = load_network(adj_path, gene_list)
        graphs.append(G)

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
    out_path = f"/home/fatemeh/ThirdObject/centrality_combined_rankings_{cell}.csv"
    df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in result_dict.items()]))
    df.to_csv(out_path, index=False)
    print(f"Combined rankings saved to {out_path}")

for method, nodes in result_dict.items():
    print(f"{method}: {len(nodes)} nodes")


