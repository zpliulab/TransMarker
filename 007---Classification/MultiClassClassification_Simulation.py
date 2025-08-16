import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, roc_auc_score, precision_score, recall_score, f1_score, average_precision_score
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.utils import to_categorical

GN = 'G50'

# Paths
path = '/home/fatemeh/ThirdObject'

# Load gene names and biomarker nodes
gene_list = pd.read_csv(path + '/4.BiomarkerIdentification/Input/all_nodes_' + GN + '.txt', header=None, names=["Gene"])  # All gene names
Biomarker_nodes = pd.read_csv(path + '/4.BiomarkerIdentification/best_component_nodes_SE_' + GN + '.txt', header=None, names=["Gene"])

# Create empty lists to store expression data for all stages
expression_matrices = []
labels = []

# Loop through each stage (1 to 4)
for ST in ['1', '2', '3', '4']:
    # Load expression matrix for the current stage
    expr_stage = np.loadtxt(path + '/2.GraphEmbedding/Input/Fea_' + GN + '_Stage_' + ST + '.txt').astype(float)

    # Find the indices of the biomarker genes in the gene list
    biomarker_indices = gene_list[gene_list["Gene"].isin(Biomarker_nodes["Gene"])].index.tolist()

    # Filter the expression matrix to include only rows corresponding to biomarker genes
    expr_filtered_stage = expr_stage[biomarker_indices, :].T  # Transpose to get samples as rows

    # Append the filtered expression matrix
    expression_matrices.append(expr_filtered_stage)

    # Create label vector for the current stage
    # Set the label for the current stage to 1 (one-hot encoding for the 4 classes)
    label_stage = np.zeros((expr_filtered_stage.shape[0], 4))  # 4 columns for 4 stages
    label_stage[:, int(ST) - 1] = 1  # Set the label for the current stage to 1

    # Append labels
    labels.append(label_stage)

# Concatenate all stages' expression matrices and labels
X = np.concatenate(expression_matrices, axis=0)  # Combine all stages' features (samples as rows, biomarker genes as columns)
y = np.concatenate(labels, axis=0)  # Combine all stages' labels

print(f"Feature matrix shape: {X.shape}")  # Samples as rows, biomarker genes as columns
print(f"Label matrix shape: {y.shape}")  # Samples as rows, 4 columns for stages

# Step 1: Split data into train, validation, and test sets (70%, 10%, 20%)
X_train, X_temp, y_train, y_temp = train_test_split(X, y, test_size=0.3, random_state=42)
X_val, X_test, y_val, y_test = train_test_split(X_temp, y_temp, test_size=0.66, random_state=42)  # 0.66 of 0.3 => 0.2 for test set

# Step 2: Standardize the features (important for neural networks)
scaler = StandardScaler()
X_train = scaler.fit_transform(X_train)
X_val = scaler.transform(X_val)
X_test = scaler.transform(X_test)

# Step 3: Convert labels to categorical format (for multi-class classification)
y_train_cat = to_categorical(np.argmax(y_train, axis=1), num_classes=4)
y_val_cat = to_categorical(np.argmax(y_val, axis=1), num_classes=4)
y_test_cat = to_categorical(np.argmax(y_test, axis=1), num_classes=4)

# Verify the shape of the labels (y_train_cat, y_val_cat, y_test_cat)
print(f"y_train_cat shape: {y_train_cat.shape}")
print(f"y_val_cat shape: {y_val_cat.shape}")
print(f"y_test_cat shape: {y_test_cat.shape}")

#  Ensure that y_train_cat, y_val_cat, and y_test_cat have the same number of samples as X_train, X_val, and X_test respectively.
print(f"X_train shape: {X_train.shape}")
print(f"X_val shape: {X_val.shape}")
print(f"X_test shape: {X_test.shape}")

# Step 4: Build the MLP model
model = Sequential()
model.add(Dense(128, input_dim=X_train.shape[1], activation='relu'))
model.add(Dense(64, activation='relu'))
model.add(Dense(32, activation='relu'))
model.add(Dense(4, activation='softmax'))  # 4 classes for stages


# Step 5: Compile the model
model.compile(optimizer=Adam(), loss='categorical_crossentropy', metrics=['accuracy'])

#from tensorflow.keras.metrics import AUC
#model.compile(optimizer=Adam(), loss='categorical_crossentropy', metrics=[AUC()])

# Step 6: Train the model
history = model.fit(X_train, y_train_cat, epochs=100, batch_size=32, validation_data=(X_val, y_val_cat))

# Step 7: Make predictions
y_pred_probs = model.predict(X_test)  # Probabilities for AUPRC, AUROC
y_pred_classes = np.argmax(y_pred_probs, axis=1)  # Predicted class labels

# Step 8: Evaluate the model
# Accuracy
accuracy = accuracy_score(np.argmax(y_test_cat, axis=1), y_pred_classes)

# AUROC (for multi-class, use macro-average)
auroc = roc_auc_score(y_test_cat, y_pred_probs, multi_class='ovr', average='macro')

# AUPRC (Precision-Recall curve for multi-class)
auprc = average_precision_score(y_test_cat, y_pred_probs, average='macro')

# F1, Precision, Recall (for multi-class, use average)
f1 = f1_score(np.argmax(y_test_cat, axis=1), y_pred_classes, average='macro')
precision = precision_score(np.argmax(y_test_cat, axis=1), y_pred_classes, average='macro')
recall = recall_score(np.argmax(y_test_cat, axis=1), y_pred_classes, average='macro')


from sklearn.metrics import confusion_matrix

# True and predicted labels
y_true = np.argmax(y_test_cat, axis=1)
y_pred = y_pred_classes

# Confusion matrix
cm = confusion_matrix(y_true, y_pred)
n_classes = cm.shape[0]

# Sensitivity (Recall) per class
sensitivity_per_class = []
specificity_per_class = []

for i in range(n_classes):
    TP = cm[i, i]
    FN = np.sum(cm[i, :]) - TP
    FP = np.sum(cm[:, i]) - TP
    TN = np.sum(cm) - (TP + FN + FP)

    sensitivity = TP / (TP + FN) if (TP + FN) > 0 else 0
    specificity = TN / (TN + FP) if (TN + FP) > 0 else 0

    sensitivity_per_class.append(sensitivity)
    specificity_per_class.append(specificity)

# Macro averages

macro_sensitivity = np.mean(sensitivity_per_class)
macro_specificity = np.mean(specificity_per_class)




# Step 9: Print the evaluation results
print(f"Accuracy: {accuracy:.4f}")
print(f"AUROC: {auroc:.4f}")
print(f"AUPRC: {auprc:.4f}")
print(f"F1 Score: {f1:.4f}")
print(f"Precision: {precision:.4f}")
print(f"Recall: {recall:.4f}")
print(f"Macro Sensitivity (Recall): {macro_sensitivity:.4f}")
print(f"Macro Specificity: {macro_specificity:.4f}")
print("Specificity per class:", [f"{s:.4f}" for s in sensitivity_per_class])
print("Specificity per class:", [f"{s:.4f}" for s in specificity_per_class])

formatted_sensitivity = [f"{s:.4f}" for s in sensitivity_per_class]
formatted_sensitivity = [f"{s:.4f}" for s in specificity_per_class]

# Save evaluation results to a text file
evaluation_results = f"""
Accuracy: {accuracy:.4f}
AUROC: {auroc:.4f}
AUPRC: {auprc:.4f}
F1 Score: {f1:.4f}
Precision: {precision:.4f}
Recall: {recall:.4f}
Sensitivity (Macro): {macro_sensitivity:.4f}
Specificity (Macro): {macro_specificity:.4f}
Sensitivity per class: {formatted_sensitivity}
Specificity per class: {formatted_sensitivity}
"""

with open('./evaluation_results_SE_'+GN+'.txt', "w") as f:
    f.write(evaluation_results)

print("Evaluation results saved to 'evaluation_results.txt'")

# Save true labels and predicted probabilities for ROC curve in a CSV file
roc_data = pd.DataFrame({
    'True_Label': np.argmax(y_test_cat, axis=1),
    'Pred_Prob_Class_1': y_pred_probs[:, 0],
    'Pred_Prob_Class_2': y_pred_probs[:, 1],
    'Pred_Prob_Class_3': y_pred_probs[:, 2],
    'Pred_Prob_Class_4': y_pred_probs[:, 3]
})

roc_data.to_csv('./roc_curve_data_SE_'+GN+'.csv', index=False)

print("ROC curve data saved to 'roc_curve_data.csv'")




import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc

# Compute ROC curve for each class
fpr = {}
tpr = {}
roc_auc = {}

for i in range(4):  # Since there are 4 classes (stages)
    fpr[i], tpr[i], _ = roc_curve(np.argmax(y_test_cat, axis=1), y_pred_probs[:, i], pos_label=i)
    roc_auc[i] = auc(fpr[i], tpr[i])

# Plot all ROC curves
plt.figure(figsize=(10, 8))
for i in range(4):
    plt.plot(fpr[i], tpr[i], label=f'Class {i+1} (AUC = {roc_auc[i]:.2f})')

# Plot the diagonal (random classifier) line
plt.plot([0, 1], [0, 1], color='gray', linestyle='--')

# Add labels and title
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curve for Multi-Class Classification')
plt.legend(loc='lower right')
plt.savefig('ROC_plot_'+GN+'.pdf')  # Save as PDF
# Show the plot
plt.show()



