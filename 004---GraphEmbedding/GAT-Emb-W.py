import torch
from torch_geometric.data import Data
import torch.nn.functional as F
from torch_geometric.nn import GATConv
import torch.optim as optim
import pandas as pd
import os
import numpy as np

GN = 'G30'
ST = '3'


def adjacency_to_edge_list(adj_matrix):
    """
    Convert an adjacency matrix to an edge list format.
    """
    adj_matrix = np.array(adj_matrix)  # Ensure it's a NumPy array
    edge_index = np.array(np.nonzero(adj_matrix))  # Get non-zero (connected) indices
    return edge_index  # Shape should be [2, num_edges]


def prepare_data(expression_matrix, adjacency_matrix):
    # Convert gene expression matrix to a tensor
    node_features = torch.tensor(expression_matrix, dtype=torch.float)

    # Convert adjacency matrix to edge list
    edge_index = adjacency_to_edge_list(adjacency_matrix)

    # Convert edge list to PyTorch tensor (ensure long type for PyG)
    edge_index = torch.tensor(edge_index, dtype=torch.long)

    # Create PyG Data object
    data = Data(x=node_features, edge_index=edge_index)
    return data


# Load expression matrix
expr_stage1 = np.loadtxt('./Input/Fea_' + GN + '_Stage_' + ST + '.txt').astype(float)

# Load adjacency matrix instead of edge list
adj_matrix = np.loadtxt('./StructuralEncoding/StructuralEncoding_Stage' + ST + '_' + GN + '.txt').astype(float)

# Convert adjacency matrix to edge list and prepare data
data_stage1 = prepare_data(expr_stage1, adj_matrix)

# Debug print
print(f"Expression matrix shape: {expr_stage1.shape}")  # Expected: (30, 100)/(50, 100)/(100, 100)
print(f"Adjacency matrix shape: {adj_matrix.shape}")  # Expected: (30, 30)/(50, 50)/(100, 100)
print(f"Edge index shape: {data_stage1.edge_index.shape}")  # Expected: (2, num_edges)


# Define GAT Model with Dropout
class GAT(torch.nn.Module):
    def __init__(self, in_channels, hidden_channels, out_channels, heads=2, dropout=0.6):
        super(GAT, self).__init__()
        self.conv1 = GATConv(in_channels, hidden_channels, heads=heads, dropout=dropout)
        self.conv2 = GATConv(hidden_channels * heads, out_channels, heads=1, dropout=dropout)

    def forward(self, data):
        x, edge_index = data.x, data.edge_index
        x = self.conv1(x, edge_index)
        x = F.relu(x)
        x = self.conv2(x, edge_index)
        return x  # Output embeddings


# Initialize model
in_dim = expr_stage1.shape[1]  # Number of features per node
gat_model = GAT(in_dim, hidden_channels=64, out_channels=32, dropout=0.6)

# Move to GPU if available
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
gat_model = gat_model.to(device)
data_stage1 = data_stage1.to(device)

optimizer = optim.Adam(gat_model.parameters(), lr=0.001, weight_decay=5e-4)


# Early Stopping
class EarlyStopping:
    def __init__(self, patience=10, delta=0):
        self.patience = patience
        self.delta = delta
        self.counter = 0
        self.best_loss = None
        self.early_stop = False

    def __call__(self, val_loss):
        if self.best_loss is None:
            self.best_loss = val_loss
        elif val_loss < self.best_loss - self.delta:
            self.best_loss = val_loss
            self.counter = 0
        else:
            self.counter += 1
            if self.counter >= self.patience:
                self.early_stop = True


# Training function with Early Stopping
def train(model, data, epochs=2000, patience=50): #epochs=200, patience=10
    early_stopping = EarlyStopping(patience=patience)

    model.train()
    for epoch in range(epochs):
        optimizer.zero_grad()
        embeddings = model(data)

        # Contrastive loss for connected nodes
        loss = torch.norm(embeddings[data.edge_index[0]] - embeddings[data.edge_index[1]], p=2).mean()
        loss.backward()
        optimizer.step()

        # Check early stopping
        early_stopping(loss.item())

        if early_stopping.early_stop:
            print(f"Early stopping at epoch {epoch}")
            break

        if epoch % 2 == 0:
            print(f"Epoch {epoch}, Loss: {loss.item()}")


train(gat_model, data_stage1)

gat_model.eval()
embeddings_stage1 = gat_model(data_stage1).detach().cpu().numpy()

# Save embeddings
output_dir = "./StructuralEncoding"
file_path_export = os.path.join(output_dir, 'embeddings_W_stage' + ST + '_' + GN + '.txt')

# Create the directory if it does not exist
os.makedirs(output_dir, exist_ok=True)

# Save embeddings
np.savetxt(file_path_export, embeddings_stage1, fmt='%.6f')

print("Saved embeddings")
