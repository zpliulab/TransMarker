import pandas as pd
import mygene

# Step 1: Load gene symbols
with open("gene_Names.txt") as f:
    g_symbols = [line.strip() for line in f if line.strip()]

# Step 2: Load KEGG Entrez IDs from the first column and gene symbols from the second
kegg_df = pd.read_csv("gastric_cancer_kegg_genes.csv", dtype=str)
kegg_entrez_ids = set(kegg_df.iloc[:, 0].dropna().astype(str))
kegg_gene_symbols = set(kegg_df.iloc[:, 1].dropna().astype(str))

# Step 3: Map DEG symbols → Entrez IDs
mg = mygene.MyGeneInfo()
entrez_results = mg.querymany(g_symbols, scopes="symbol", fields="entrezgene", species="human")
symbol_to_entrez = {
    r['query']: str(r['entrezgene'])
    for r in entrez_results if not r.get('notfound') and 'entrezgene' in r
}
g_entrez_ids = set(symbol_to_entrez.values())

# Step 4: Map Entrez IDs → Gene Symbols
symbol_results = mg.querymany(list(g_entrez_ids), scopes="entrezgene", fields="symbol", species="human")
entrez_to_symbol = {
    str(r['query']): r['symbol']
    for r in symbol_results if not r.get('notfound') and 'symbol' in r
}
g_mapped_symbols = set(entrez_to_symbol.values())

# Step 5: Compare to KEGG symbols
common_symbols = g_mapped_symbols.intersection(kegg_gene_symbols)
common_entrez = g_entrez_ids.intersection(kegg_entrez_ids)
# Report
print(f"Mapped DEG Entrez IDs → Gene Symbols: {len(g_mapped_symbols)}")
print(f"KEGG Gene Symbols: {len(kegg_gene_symbols)}")
print(f"Common Gene Symbols: {len(common_symbols)}")
print("Common Gene Symbols:")
print(sorted(common_symbols))
print(f"Common Gene IDs: {len(common_entrez)}")


