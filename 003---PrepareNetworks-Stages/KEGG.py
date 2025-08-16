from bs4 import BeautifulSoup
import pandas as pd
import mygene

# Step 1: Load KGML (KEGG XML) file
with open("hsa05226.xml", "r") as f:
    soup = BeautifulSoup(f, "xml")

# Step 2: Extract Entrez IDs from gene entries
gene_entries = soup.find_all("entry", {"type": "gene"})

entrez_ids = set()
for entry in gene_entries:
    for item in entry["name"].split():
        if item.startswith("hsa:"):
            entrez_ids.add(item.replace("hsa:", ""))

print(f"✅ Extracted {len(entrez_ids)} Entrez IDs.")

# Step 3: Use MyGene to map Entrez ID → Gene Symbol
mg = mygene.MyGeneInfo()
results = mg.querymany(list(entrez_ids), scopes="entrezgene", fields="symbol", species="human")

# Step 4: Collect mappings
mapping = {}
for r in results:
    if "symbol" in r:
        mapping[r["query"]] = r["symbol"]
    else:
        mapping[r["query"]] = "UNKNOWN"

# Step 5: Save to CSV
df = pd.DataFrame(list(mapping.items()), columns=["EntrezID", "GeneSymbol"])
df.to_csv("gastric_cancer_kegg_genes.csv", index=False)

print("✅ Final gene symbol list saved to 'gastric_cancer_kegg_genes.csv'")


