import ot
import numpy as np
from sklearn.metrics.pairwise import cosine_distances


GN = 'CD4'

path='/home/fatemeh/ThirdObject'

# Load embeddings for four stages
E1 = np.loadtxt(path+'/2.GraphEmbedding/StructuralEncoding/embeddings_W_stage1_' + GN + '.txt', delimiter=" ")
E2 = np.loadtxt(path+'/2.GraphEmbedding/StructuralEncoding/embeddings_W_stage2_' + GN + '.txt', delimiter=" ")
E3 = np.loadtxt(path+'/2.GraphEmbedding/StructuralEncoding/embeddings_W_stage3_' + GN + '.txt', delimiter=" ")
E4 = np.loadtxt(path+'/2.GraphEmbedding/StructuralEncoding/embeddings_W_stage4_' + GN + '.txt', delimiter=" ")

print("Min and Max of C1:", np.min(E1), np.max(E1))
print("Min and Max of C2:", np.min(E2), np.max(E2))
print("Min and Max of C3:", np.min(E3), np.max(E3))
print("Min and Max of C4:", np.min(E4), np.max(E4))

# Compute distance matrices
#C1 = np.linalg.norm(E1[:, None] - E1, axis=2)
#C2 = np.linalg.norm(E2[:, None] - E2, axis=2)
#C3 = np.linalg.norm(E3[:, None] - E3, axis=2)
#C4 = np.linalg.norm(E4[:, None] - E4, axis=2)

# Compute cosine distances (1 - cosine similarity)
C1 = cosine_distances(E1)
C2 = cosine_distances(E2)
C3 = cosine_distances(E3)
C4 = cosine_distances(E4)

print("Min and Max of C1:", np.min(C1), np.max(C1))
print("Min and Max of C2:", np.min(C2), np.max(C2))
print("Min and Max of C3:", np.min(C3), np.max(C3))
print("Min and Max of C4:", np.min(C4), np.max(C4))

print("Normalize")
# Normalize distances
#C1, C2, C3, C4 = [C / (np.std(C) + 1e-6) for C in [C1, C2, C3, C4]]
##C1, C2, C3, C4 = [C / (np.max(C) + 1e-6) for C in [C1, C2, C3, C4]]

# Add a small constant to avoid division by zero during normalization
#epsilon = 1e-6
#C1, C2, C3, C4 = [C + epsilon for C in [C1, C2, C3, C4]]

# Normalize distances
#C1, C2, C3, C4 = [C / (np.std(C) + 1e-6) for C in [C1, C2, C3, C4]]



print("Min and Max of C1:", np.min(C1), np.max(C1))
print("Min and Max of C2:", np.min(C2), np.max(C2))
print("Min and Max of C3:", np.min(C3), np.max(C3))
print("Min and Max of C4:", np.min(C4), np.max(C4))


print("C1 sample:\n", C1[:5, :5])
print("C2 sample:\n", C2[:5, :5])
print("Mean of C1:", np.mean(C1), "Std of C1:", np.std(C1))
print("Mean of C2:", np.mean(C2), "Std of C2:", np.std(C2))

# Compute GW alignments


def compute_gromov_wasserstein(C_source, C_target, epsilon=5e-4, max_iter=500):
    """
    Computes the Gromov-Wasserstein transport plan between two cost matrices.

    Parameters:
    - C_source: np.ndarray, cost matrix of the source space
    - C_target: np.ndarray, cost matrix of the target space
    - epsilon: float, regularization parameter
    - max_iter: int, maximum number of iterations

    Returns:
    - np.ndarray, the computed transport plan
    """
    p = np.sum(C_source, axis=1)
    p /= np.sum(p)  # Normalize

    q = np.sum(C_target, axis=1)
    q /= np.sum(q)  # Normalize

    transport_plan = ot.gromov.gromov_wasserstein(
        C_source, C_target, p, q, 'square_loss', verbose=True, epsilon=epsilon, max_iter=max_iter
    )

    return transport_plan


# Compute transport plans for consecutive cost matrices
gw_matrices = []
C_matrices = [C1, C2, C3, C4]

for i in range(len(C_matrices) - 1):
    T = compute_gromov_wasserstein(C_matrices[i], C_matrices[i + 1])
    gw_matrices.append(T)
    print(f"T{i + 1}{i + 2} sample:\n", T[:5, :5])

# Compute cumulative change for each gene
#If you care about variability in alignment patterns (diverse interactions across stages) → Use Standard Deviation.
'''
alignment_scores = alignment_scores = (
    np.std(T12, axis=1) +
    np.std(T23, axis=1) +
    np.std(T34, axis=1)
)
'''

alignment_scores = sum(np.std(T, axis=1) for T in gw_matrices)

print(alignment_scores)

#If you care about large shifts in gene alignment (genes significantly changing their network position) → Use Absolute Deviation.
'''
alignment_scores = (
    np.abs(T12 - np.mean(T12, axis=1, keepdims=True)).sum(axis=1) +
    np.abs(T23 - np.mean(T23, axis=1, keepdims=True)).sum(axis=1) +
    np.abs(T34 - np.mean(T34, axis=1, keepdims=True)).sum(axis=1)
)
print(alignment_scores)
'''


# Save alignment scores
np.savetxt('./gw_cumulative_alignment_SE_'+GN+'.csv', alignment_scores, delimiter=",")
print("Cumulative GW alignment computed!")



'''

p = np.array(np.sum(C1, axis=1))  # Sum of distances per node
p = p / np.sum(p)  # Normalize

q = np.array(np.sum(C2, axis=1))
q = q / np.sum(q)
T12 = ot.gromov.gromov_wasserstein(
    C1, C2, p, q, 'square_loss', verbose=True, epsilon=5e-4, max_iter=500
)
print("T12 sample:\n", T12[:5, :5])

p = np.array(np.sum(C2, axis=1))  # Sum of distances per node
p = p / np.sum(p)  # Normalize

q = np.array(np.sum(C3, axis=1))
q = q / np.sum(q)
T23 = ot.gromov.gromov_wasserstein(
    C2, C3, p, q, 'square_loss', verbose=True, epsilon=5e-4, max_iter=500
)
print("T23 sample:\n", T23[:5, :5])

p = np.array(np.sum(C3, axis=1))  # Sum of distances per node
p = p / np.sum(p)  # Normalize

q = np.array(np.sum(C4, axis=1))
q = q / np.sum(q)
T34 = ot.gromov.gromov_wasserstein(
    C3, C4, p, q, 'square_loss', verbose=True, epsilon=5e-4, max_iter=500
)
print("T34 sample:\n", T34[:5, :5])
'''
