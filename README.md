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

## R and Python Code for TransMarker

The numbered folders **(1)–(9)** correspond to the sequential workflow of TransMarker.  

- **(1)** Preparation of simulation data, including preprocessing scripts and simulation data files.  
- **(2)** Preprocessing of real datasets and storage of intermediate results.  
- **(3)** R and Python scripts to construct disease states as network layers for both simulation and real data.  
- **(4)** Embedding scripts for generating representations of each network layer.  
- **(5)** Alignment process for comparing layers across states (simulation and real data).  
- **(6)** Biomarker identification for detecting transition genes in both simulation and real data.  
- **(7)** Python code for multi-class disease state classification.  
- **(8)** Benchmarking scripts for comparison with baseline and state-of-the-art methods.  
- **(9)** R scripts for visualization of all results.  




