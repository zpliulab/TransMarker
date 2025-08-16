# [TransMarker: Unveiling dynamic network biomarkers in cancer progression through cross-state graph alignment and optimal transport](https://github.com/zpliulab/TransMarker)

**TransMarker** is a computational framework for detecting dynamic biomarkers in cancer progression using single-cell data. It models disease states as multilayer networks, aligns regulatory shifts across states, and identifies genes with evolving regulatory roles.

If you have any questions about NetWalkRank, please directly contact the corresponding author Prof. Zhi-Ping Liu with the E-mail: zpliu@sdu.edu.cn


## Highlights:
- **Dynamic biomarker discovery:** Identifies state-specific genes with shifting regulatory roles during cancer progression.

- **Multilayer network modeling:** Encodes each disease state as a separate layer, integrating gene expression with prior biological interactions.

- **Cross-state alignment:** Uses Gromov–Wasserstein optimal transport to quantify structural shifts of genes across disease states.

- **Dynamic Network Index (DNI):** Introduces a new metric to capture regulatory variability and prioritize dynamic network biomarkers.

- **Deep learning integration:** Combines Graph Attention Networks (GATs) and deep neural networks (DNNs) for contextual embedding and disease state classification.

- **Robust validation:** Demonstrated superior performance on simulated and real single-cell gastric adenocarcinoma datasets, outperforming existing multilayer ranking approaches.

- **Biological relevance:** Detects genes that reflect meaningful rewiring of gene regulation, providing insight into molecular drivers of cancer progression.


## Data:

In the Data directory, we provide only the essential files required to run the code. The datasets, which were obtained from publicly available sources, are introduced and referenced in the paper, along with the respective website links where they can be accessed.


## R and Python Code to TransMarker:

The serial number (1) (2) ... (9) represents the order in which the program runs in our work.

- (1) Code for prepare simulation data and preprocess on simulation data and simulation data files.
- (2) Preprocess code for real data and results of this step.
- (3) R and Python code to prepare states of disease as network's layers for both simulation and real data.
- (4) Embedding process code for embedding layers.
- (5) Alignment Process code for both simulation and real data.
- (6) Biomarker Identification code for both simulation and real data.
- (7) Multi-Class classification Python code.
- (8) Codes for comparision with other similar and baseline methods.
- (9) R code for visualize all results.



