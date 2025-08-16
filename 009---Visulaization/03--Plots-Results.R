# Load required libraries
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)




setwd("E:/005---ThirdProject/ThirdObject/6.Results/")


# Define the tags and corresponding file names
tags <- paste0("T", seq(10, 50, 10), "p")
file_names <- paste0("evaluation_metrics_50_runs_", tags, ".csv")

# Custom colors (modify to your color codes)
custom_colors <- c(
  "T10p" = "#A2C2E2",
  "T20p" = "#FEAF8A",
  "T30p" = "#F6B7C6",
  "T40p" = "#A2DADE",
  "T50p" = "#D8CBF0"
)

# Initialize empty list to collect AUROC data
auroc_data_list <- list()

# Read each file and extract AUROC values
for (tag in tags) {
  file_path <- paste0("evaluation_metrics_50_runs_", tag, ".csv")
  if (file.exists(file_path)) {
    df <- readr::read_csv(file_path, show_col_types = FALSE)
    df$Tag <- tag
    auroc_data_list[[tag]] <- df %>% dplyr::select(Tag, auroc)
  } else {
    warning(paste("File not found:", file_path))
  }
}

# Combine all AUROC values into one data frame
auroc_all <- dplyr::bind_rows(auroc_data_list)

# Define label mapping for x-axis
label_map <- setNames(paste0(seq(10, 50, 10), "%"), tags)

# Plot boxplot with updated x-axis labels
Auc_P<- ggplot(auroc_all, aes(x = Tag, y = auroc, fill = Tag)) +
  geom_boxplot(width = 0.6, outlier.shape = 16, outlier.size = 1.5) +
  scale_fill_manual(values = custom_colors) +
  scale_x_discrete(labels = label_map) +
  labs(
    title = "",
    x = "",
    y = "AUROC"
  ) +
  theme_bw(base_size = 14) + 
  theme(text = element_text(family = "Times"))+# First apply the bw theme
  theme(
    legend.position = "none",               # Then override to hide legend
    axis.text.x = element_text(angle = 0, hjust = 1)  # Optional: rotate x labels
  )

Auc_P



pdf("0---GAC/20250620/Plot_Auc_P1.pdf",width = 4,height =3)
Auc_P
dev.off()

svg("0---GAC/20250620/Plot_Auc_P1.svg", width = 4, height = 3)
# Print the plot
Auc_P
# Close the SVG device to save the file
dev.off()




# Load necessary libraries
library(ggplot2)
library(readr)
library(dplyr)

# Define tags and file names for the 'a' series
tags_a <- paste0("a", seq(10, 90, 10))
file_names_a <- paste0("evaluation_metrics_50_runs_", tags_a, ".csv")

# Define your base colors and interpolate
base_colors <- c("#A2C2E2", "#FEAF8A", "#F6B7C6", "#A2DADE", "#D8CBF0")
color_palette_a <- colorRampPalette(base_colors)(length(tags_a))
names(color_palette_a) <- tags_a

# Create mapping for display labels: a10 → α = 0.1
alpha_labels <- paste0("α = ", seq(0.1, 0.9, 0.1))
names(alpha_labels) <- tags_a

# Read and collect AUROC data
auroc_data_list_a <- list()
for (tag in tags_a) {
  file_path <- paste0("evaluation_metrics_50_runs_", tag, ".csv")
  if (file.exists(file_path)) {
    df <- read_csv(file_path, show_col_types = FALSE)
    df$Tag <- tag
    auroc_data_list_a[[tag]] <- df %>% select(Tag, auroc)
  } else {
    warning(paste("File not found:", file_path))
  }
}

# Combine into one data frame
auroc_all_a <- bind_rows(auroc_data_list_a)
auroc_all_a$Tag <- factor(auroc_all_a$Tag, levels = tags_a)

# Plot the boxplot
AUC_a<- ggplot(auroc_all_a, aes(x = Tag, y = auroc, fill = Tag)) +
  geom_boxplot(width = 0.6, outlier.shape = 16, outlier.size = 1.5) +
  scale_fill_manual(values = color_palette_a) +
  scale_x_discrete(labels = alpha_labels) +
  labs(
    title = "",
    x = "",
    y = "AUROC"
  ) +
  theme_bw(base_size = 14) +
  theme(text = element_text(family = "Times"))+
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

AUC_a

pdf("0---GAC/20250620/Plot_Auc_a2.pdf",width = 6,height =3)
AUC_a
dev.off()

svg("0---GAC/20250620/Plot_Auc_a2.svg", width = 6, height = 3)
# Print the plot
AUC_a
# Close the SVG device to save the file
dev.off()









library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# Define tags and percentage labels
tags_T <- paste0("T", seq(10, 50, 10), "p")
percent_labels <- paste0(seq(10, 50, 10), "%")
names(percent_labels) <- tags_T

# Extended color palette
color_palette_extended <- colorRampPalette(c("#A2C2E2", "#FEAF8A", "#F6B7C6", "#A2DADE", "#D8CBF0"))(7)

# Read files and compute mean metrics
mean_metrics_T <- data.frame()
for (tag in tags_T) {
  file_path <- paste0("evaluation_metrics_50_runs_", tag, ".csv")
  if (file.exists(file_path)) {
    df <- read_csv(file_path, show_col_types = FALSE)
    df_mean <- df %>%
      summarise(across(where(is.numeric), mean)) %>%
      mutate(Tag = tag)
    mean_metrics_T <- bind_rows(mean_metrics_T, df_mean)
  } else {
    warning(paste("File not found:", file_path))
  }
}

# Reshape to long format
long_T <- pivot_longer(mean_metrics_T, 
                       cols = c(accuracy, auroc, auprc, f1, precision, recall, macro_specificity),
                       names_to = "Metric", values_to = "Mean")
long_T$Tag <- factor(long_T$Tag, levels = tags_T)

# Plot
P_line <- ggplot(long_T, aes(x = Tag, y = Mean, group = Metric, color = Metric)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 3) +
  scale_color_manual(values = color_palette_extended) +
  scale_x_discrete(labels = percent_labels) +  # Replace tag with percentage
  labs(
    title = "",
    x = "",
    y = "Mean Value"
  ) +
  theme_bw(base_size = 14) +
  theme(text = element_text(family = "Times"))+
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

# Show plot
P_line

# Save as PDF
pdf("0---GAC/20250620/Plot_line1.pdf", width = 6, height = 3)
P_line
dev.off()

# Save as SVG
svg("0---GAC/20250620/Plot_line1.svg", width = 6, height = 3)
P_line
dev.off()





library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# Define tags and pretty alpha labels
tags_a <- paste0("a", seq(10, 90, 10))
alpha_labels <- paste0("α = ", seq(0.1, 0.9, 0.1))
names(alpha_labels) <- tags_a

# Generate interpolated color palette
color_palette_extended <- colorRampPalette(c("#A2C2E2", "#FEAF8A", "#F6B7C6", "#A2DADE", "#D8CBF0"))(length(tags_a))

# Read all files and calculate mean of each metric
mean_metrics_a <- data.frame()

for (tag in tags_a) {
  file_path <- paste0("evaluation_metrics_50_runs_", tag, ".csv")
  if (file.exists(file_path)) {
    df <- read_csv(file_path, show_col_types = FALSE)
    df_mean <- df %>%
      summarise(across(where(is.numeric), mean)) %>%
      mutate(Tag = tag)
    mean_metrics_a <- bind_rows(mean_metrics_a, df_mean)
  } else {
    warning(paste("File not found:", file_path))
  }
}

# Reshape for plotting
long_a <- pivot_longer(mean_metrics_a, 
                       cols = c(accuracy, auroc, auprc, f1, precision, recall, macro_specificity),
                       names_to = "Metric", values_to = "Mean")
long_a$Tag <- factor(long_a$Tag, levels = tags_a)

# Plot
P_alpha <- ggplot(long_a, aes(x = Tag, y = Mean, group = Metric, color = Metric)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 3) +
  scale_color_manual(values = color_palette_extended) +
  scale_x_discrete(labels = alpha_labels) +  # Replace tag with alpha labels
  labs(
    title = "",
    x = "",
    y = "Mean Value"
  ) +
  theme_bw(base_size = 14) +
  theme(text = element_text(family = "Times"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Display plot
P_alpha

# Save as PDF
pdf("0---GAC/20250620/Plot_alpha_line1.pdf", width = 8, height = 3)
P_alpha
dev.off()

# Save as SVG
svg("0---GAC/20250620/Plot_alpha_line1.svg", width = 8, height = 3)
P_alpha
dev.off()








library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# File tags and display labels
file_tags <- c("a_n", "a0", "a100", "a30")
#display_labels <- c("w/o local & global", "w/o local", "w/o global", "method")
display_labels <- c("None", "w/o local", "w/o global", "method")
file_paths <- paste0("evaluation_metrics_50_runs_", file_tags, ".csv")

# Read and combine data
all_data <- list()

for (i in seq_along(file_paths)) {
  if (file.exists(file_paths[i])) {
    df <- read_csv(file_paths[i], show_col_types = FALSE) %>%
      select(accuracy, auroc, auprc) %>%
      mutate(Method = display_labels[i])
    all_data[[i]] <- df
  } else {
    warning(paste("File not found:", file_paths[i]))
  }
}

# Combine all
df_all <- bind_rows(all_data) %>%
  pivot_longer(cols = c(accuracy, auroc, auprc), names_to = "Metric", values_to = "Value")

# Order metrics
df_all$Metric <- factor(df_all$Metric, levels = c("accuracy", "auroc", "auprc"))

# Summary: mean and SE
summary_df <- df_all %>%
  group_by(Method, Metric) %>%
  summarise(Mean = mean(Value), SE = sd(Value)/sqrt(n()), .groups = "drop")

summary_df$Metric <- factor(summary_df$Metric, levels = c("accuracy", "auroc", "auprc"))

# Order methods
method_levels <- display_labels
df_all$Method <- factor(df_all$Method, levels = method_levels)
summary_df$Method <- factor(summary_df$Method, levels = method_levels)

# Define colors
metric_colors <- c("accuracy" = "#A2C2E2", "auroc" = "#D8CBF0", "auprc" = "#A2DADE")

#c("#A2C2E2", "#FEAF8A", "#F6B7C6", "#A2DADE", "#D8CBF0"))

# Plot
p <- ggplot(summary_df, aes(x = Method, y = Mean, fill = Metric)) +
  geom_col(position = position_dodge(0.8), width = 0.7, color = "black") +
  geom_jitter(data = df_all, aes(x = Method, y = Value, color = Metric),
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              size = 1.5, alpha = 0.7, inherit.aes = FALSE) +
  geom_errorbar(aes(ymin = Mean - (2*SE), ymax = Mean + (2*SE)),
                position = position_dodge(0.8), width = 0.3) +
  geom_text(aes(label = sprintf("%.3f ± %.3f", Mean, (2*SE)), y = 0.21 + SE + 0.01), # Put near bottom
            position = position_dodge(0.8),
            angle = 90,
            vjust = 0.5,
            hjust = 0.15,
            size = 3.5,
            fontface = "bold",
            color = "black") +
  scale_fill_manual(values = metric_colors) +
  scale_color_manual(values = metric_colors) +
  labs(
    x = NULL,
    y = "Metric Value",
    fill = "Metric",
    color = "Metric"
  ) +
  coord_cartesian(ylim = c(0.2, NA)) +
  theme_bw(base_size = 14) +
  theme(text = element_text(family = "Times"))+
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )


# Print plot
p

# 
# Save as PDF
pdf("0---GAC/20250620/Plot_aaa_barplot.pdf", width = 7, height = 4)
p
dev.off()

# Save as SVG
svg("0---GAC/20250620/Plot_aaa_barplot.svg", width = 7, height = 4)
p
dev.off()








library(readr)
library(ggplot2)
library(dplyr)

# 1. Read single-column alignment score files
gw <- read_csv("AlignmentScores/gw_cumulative_alignment_Reg_CD8.csv", col_names = "Alignment", show_col_types = FALSE)
best <- read_csv("AlignmentScores/Best_alignment_scoresCD8.csv", col_names = "Alignment", show_col_types = FALSE)
dim(gw)
dim(best)
# 2. Add source labels
gw$Source <- "GW"
best$Source <- "Best"

# 3. Filter out zero values
gw <- gw %>% filter(Alignment != 0)
best <- best %>% filter(Alignment != 0)

# 4. Combine the two datasets
combined <- bind_rows(gw, best)
#combined<-gw


#install.packages("patchwork")
library(ggplot2)
library(ggridges)
library(patchwork)

# Add Source column for clarity if missing
gw$Source <- "GW"
best$Source <- "Best"

p1 <- ggplot(gw, aes(x = Alignment, y = Source, fill = Source)) +
  geom_density_ridges(scale = 1.5, alpha = 0.9) +
  scale_fill_manual(values = c("GW" = "#A2DADE")) +
  scale_y_discrete(expand = expansion(add = c(0.9, 0.9))) +
  labs(x = "Alignment Score", y = NULL) +
  theme_bw(base_size = 14) +
  theme(legend.position = "right")

p2 <- ggplot(best, aes(x = Alignment, y = Source, fill = Source)) +
  geom_density_ridges(scale = 1.5, alpha = 0.9) +
  scale_fill_manual(values = c("Best" = "#D8CBF0")) +
  scale_y_discrete(expand = expansion(add = c(0.9, 0.9))) +
  labs(x = "Alignment Score", y = NULL) +
  theme_bw(base_size = 14) +
  theme(legend.position = "right")

# Combine plots side by side
#p1 + p2



#p_align<- p1 + p2

#p_align

pdf("Plot_align-G-1.pdf", width = 6, height = 4)
p1
dev.off()

# Save as SVG
svg("Plot_align-G-1.svg", width = 6, height = 4)
p1
dev.off()


pdf("Plot_align-G-Best-1.pdf", width = 4, height = 6)
p2
dev.off()

# Save as SVG
svg("Plot_align-G-Best-1.svg", width = 4, height = 6)
p2
dev.off()





library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# Files and display labels
file_tags <- c("a30", "a30-ESCC")
disease_labels <- c("Gastric", "ESCC")
file_paths <- paste0("evaluation_metrics_50_runs_", file_tags, ".csv")

# Read and combine data
all_data <- list()

for (i in seq_along(file_paths)) {
  if (file.exists(file_paths[i])) {
    df <- read_csv(file_paths[i], show_col_types = FALSE) %>%
      select(accuracy, auroc, auprc, f1, precision, recall, macro_specificity) %>%
      mutate(Disease = disease_labels[i])
    all_data[[i]] <- df
  } else {
    warning(paste("File not found:", file_paths[i]))
  }
}

# Combine all
df_all <- bind_rows(all_data) %>%
  pivot_longer(cols = c(accuracy, auroc, auprc, f1, precision, recall, macro_specificity), 
               names_to = "Metric", 
               values_to = "Value")

# Order metrics (optional, for nice order on x axis)
df_all$Metric <- factor(df_all$Metric, 
                        levels = c("accuracy", "auroc", "auprc", "f1", "precision", "recall", "macro_specificity"))

# Summary: mean and SE
summary_df <- df_all %>%
  group_by(Disease, Metric) %>%
  summarise(Mean = mean(Value), SE = sd(Value)/sqrt(n()), .groups = "drop")

c("#A2C2E2", "#FEAF8A", "#F6B7C6", "#A2DADE", "#D8CBF0")
# Colors for Diseases
disease_colors <- c("Gastric" = "#A2DADE", "ESCC" = "#D8CBF0")

p2 <- ggplot(summary_df, aes(x = Metric, y = Mean, fill = Disease)) +
  geom_col(position = position_dodge(0.8), width = 0.7, color = "black") +
  geom_jitter(data = df_all, aes(x = Metric, y = Value, color = Disease),
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              size = 1.5, alpha = 0.7, inherit.aes = FALSE) +
  geom_errorbar(aes(ymin = Mean - (2*SE), ymax = Mean + (2*SE)),
                position = position_dodge(0.8), width = 0.3) +
  geom_text(aes(label = sprintf("%.3f ± %.3f", Mean, (2*SE)), y = 0.21 + SE + 0.01), # adjust if needed
            position = position_dodge(0.8),
            angle = 90,
            vjust = 0.5,
            hjust = 0.15,
            size = 3.5,
            fontface = "bold",
            color = "black") +
  scale_fill_manual(values = disease_colors) +
  scale_color_manual(values = disease_colors) +
  labs(
    x = NULL,
    y = "Metric Value",
    fill = "Disease",
    color = "Disease"
  ) +
  coord_cartesian(ylim = c(0.2, NA)) +
  theme_bw(base_size = 14) +
  theme(text = element_text(family = "Times"))+
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 0.9)
  )

# Print plot
p2

# Save as PDF
pdf("0---GAC/20250620/Plot_Gastric_ESCC_barplot.pdf", width = 6, height = 4)
p2
dev.off()

# Save as SVG
svg("0---GAC/20250620/Plot_Gastric_ESCC_barplot.svg", width = 6, height = 4)
p2
dev.off()












library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# Define your files and disease labels
file_paths <- c("evaluation_metrics_50_runs_a30.csv", "evaluation_metrics_50_runs_a30-ESCC.csv")
disease_labels <- c("Gastric", "ESCC")

# Define your metric levels and colors
metric_levels <- c("accuracy", "auroc", "auprc", "f1", "precision", "recall", "macro_specificity")
metric_colors <- c("#FBCC6B", "#FEAF8A", "#F6B7C6", "#A2C2E2",   "#A2DADE", "#D8CBF0", "#B2D2E2")
metric_colors <- colorRampPalette(c("#A2C2E2", "#FEAF8A", "#F6B7C6", "#A2DADE", "#D8CBF0"))(7)
metric_colors <- colorRampPalette(c("#A2C2E2", "#F6B7C6", "#A2DADE", "#D8CBF0"))(7)
# Read data and combine
all_data <- list()

for (i in seq_along(file_paths)) {
  if (file.exists(file_paths[i])) {
    df <- read_csv(file_paths[i], show_col_types = FALSE) %>%
      select(all_of(metric_levels)) %>%
      mutate(Disease = disease_labels[i])
    all_data[[i]] <- df
  } else {
    warning(paste("File not found:", file_paths[i]))
  }
}

# Combine all
df_all <- bind_rows(all_data) %>%
  pivot_longer(cols = all_of(metric_levels), names_to = "Metric", values_to = "Value")

# Order metrics and disease factor levels
df_all$Metric <- factor(df_all$Metric, levels = c("macro_specificity", "recall", "precision", "f1", "auprc", "auroc", "accuracy"))
df_all$Disease <- factor(df_all$Disease, levels = c("Gastric", "ESCC"))  # ensure correct order

# Summary statistics (mean + SE)
summary_df <- df_all %>%
  group_by(Disease, Metric) %>%
  summarise(Mean = mean(Value), SE = sd(Value)/sqrt(n()), .groups = "drop")


p <- ggplot(summary_df, aes(x = Disease, y = Mean, fill = Metric)) +
  geom_col(position = position_dodge(0.8), width = 0.7, color = "black") +
  geom_jitter(data = df_all, 
              aes(x = Disease, y = Value, color = Metric, fill = Metric),
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              size = 1.5, alpha = 0.7, inherit.aes = FALSE, shape = 21) +
  geom_errorbar(aes(ymin = Mean - (2*SE), ymax = Mean + (2*SE)),
                position = position_dodge(0.8), width = 0.3) +
  geom_text(aes(label = sprintf("%.3f ± %.3f", Mean, (2*SE)), y = 0.21 + SE + 0.01),# adjust if needed
            position = position_dodge(0.8),
            angle = 90,
            vjust = 0.5,
            hjust = 0.15,
            size = 3.5,
            fontface = "bold",
            color = "black") +
  scale_fill_manual(values = metric_colors) +
  scale_color_manual(values = metric_colors) +
  labs(
    x = NULL,
    y = "Metric Value",
    fill = "Metric"
  ) +
  coord_cartesian(ylim = c(0.2, NA)) +
  theme_bw(base_size = 14) +
  theme(text = element_text(family = "Times"))+
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

# Print plot
p

 
# Save as PDF
pdf("0---GAC/20250620/Plot_Gastric_ESCC_barplot_Group.pdf", width = 6, height = 4)
p
dev.off()

# Save as SVG
svg("0---GAC/20250620/Plot_Gastric_ESCC_barplot_Group.svg", width = 6, height = 4)
p
dev.off()





library(readr)
library(ggplot2)
library(dplyr)

# 1. Read single-column alignment score files
gw <- read_csv("AlignmentScores/ESCC_gw_cumulative_alignment_Reg_CD8.csv", col_names = "Alignment", show_col_types = FALSE)
best <- read_csv("AlignmentScores/Best-alignment_score_ESCC.csv", col_names = "Alignment", show_col_types = FALSE)
dim(gw)
dim(best)
# 2. Add source labels
gw$Source <- "GW"
best$Source <- "Best"

# 3. Filter out zero values
gw <- gw %>% filter(Alignment != 0)
best <- best %>% filter(Alignment != 0)

# 4. Combine the two datasets
combined <- bind_rows(gw, best)
#combined<-gw


#install.packages("patchwork")
library(ggplot2)
library(ggridges)
library(patchwork)

# Add Source column for clarity if missing
gw$Source <- "GW"
best$Source <- "Best"

p1 <- ggplot(gw, aes(x = Alignment, y = Source, fill = Source)) +
  geom_density_ridges(scale = 1.5, alpha = 0.9) +
  scale_fill_manual(values = c("GW" = "#A2DADE")) +
  scale_y_discrete(expand = expansion(add = c(0.9, 0.9))) +
  labs(x = "Alignment Score", y = NULL) +
  theme_bw(base_size = 14) +
  theme(legend.position = "right")

p2 <- ggplot(best, aes(x = Alignment, y = Source, fill = Source)) +
  geom_density_ridges(scale = 1.5, alpha = 0.9) +
  scale_fill_manual(values = c("Best" = "#D8CBF0")) +
  scale_y_discrete(expand = expansion(add = c(0.9, 0.9))) +
  labs(x = "Alignment Score", y = NULL) +
  theme_bw(base_size = 14) +
  theme(legend.position = "right")

# Combine plots side by side
p1 + p2



p_align<- p1 + p2

p_align

pdf("Plot_align_ESCC.pdf", width = 8, height = 6)
p1
dev.off()

# Save as SVG
svg("Plot_align_ESCC.svg", width = 8, height = 6)
p1
dev.off()


pdf("Plot_align_ESCC_Best.pdf", width = 4, height = 6)
p2
dev.off()

# Save as SVG
svg("Plot_align_ESCC_Best.svg", width = 4, height = 6)
p2
dev.off()








library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# File tags and display labels
file_tags <- c("NoFilter", "Correlation", "a30", "10p", "20p", "30p")

display_labels <- c("No Filter", "Correlation", "method",  "+10%",  "+20%",  "+30%")
file_paths <- paste0("Ablation2/evaluation_metrics_50_runs_", file_tags, ".csv")

# Read and combine data
all_data <- list()

for (i in seq_along(file_paths)) {
  if (file.exists(file_paths[i])) {
    df <- read_csv(file_paths[i], show_col_types = FALSE) %>%
      select(accuracy, auroc, auprc) %>%
      mutate(Method = display_labels[i])
    all_data[[i]] <- df
  } else {
    warning(paste("File not found:", file_paths[i]))
  }
}

# Combine all
df_all <- bind_rows(all_data) %>%
  pivot_longer(cols = c(accuracy, auroc, auprc), names_to = "Metric", values_to = "Value")

# Order metrics
df_all$Metric <- factor(df_all$Metric, levels = c("accuracy", "auroc", "auprc"))

# Summary: mean and SE
summary_df <- df_all %>%
  group_by(Method, Metric) %>%
  summarise(Mean = mean(Value), SE = sd(Value)/sqrt(n()), .groups = "drop")

summary_df$Metric <- factor(summary_df$Metric, levels = c("accuracy", "auroc", "auprc"))

# Order methods
method_levels <- display_labels
df_all$Method <- factor(df_all$Method, levels = method_levels)
summary_df$Method <- factor(summary_df$Method, levels = method_levels)

# Define colors
metric_colors <- c("accuracy" = "#A2C2E2", "auroc" = "#D8CBF0", "auprc" = "#A2DADE")

#c("#A2C2E2", "#FEAF8A", "#F6B7C6", "#A2DADE", "#D8CBF0"))

# Plot
p <- ggplot(summary_df, aes(x = Method, y = Mean, fill = Metric)) +
  geom_col(position = position_dodge(0.8), width = 0.7, color = "black") +
  geom_jitter(data = df_all, aes(x = Method, y = Value, color = Metric),
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              size = 1.5, alpha = 0.7, inherit.aes = FALSE) +
  geom_errorbar(aes(ymin = Mean - (2*SE), ymax = Mean + (2*SE)),
                position = position_dodge(0.8), width = 0.3) +
  geom_text(aes(label = sprintf("%.3f ± %.3f", Mean, (2*SE)), y = 0.21 + SE + 0.01), # Put near bottom
            position = position_dodge(0.8),
            angle = 90,
            vjust = 0.5,
            hjust = 0.15,
            size = 3.5,
            fontface = "bold",
            color = "black") +
  scale_fill_manual(values = metric_colors) +
  scale_color_manual(values = metric_colors) +
  labs(
    x = NULL,
    y = "Metric Value",
    fill = "Metric",
    color = "Metric"
  ) +
  coord_cartesian(ylim = c(0.2, NA)) +
  theme_bw(base_size = 14) +
  theme(text = element_text(family = "Times"))+
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )


# Print plot
p

 
# Save as PDF
pdf("0---GAC/20250620/Plot_DifferentFiltering.pdf", width = 7, height = 4)
p
dev.off()

# Save as SVG
svg("0---GAC/20250620/Plot_DifferentFiltering.svg", width = 7, height = 4)
p
dev.off()



# Load required libraries
library(ggplot2)
library(tidyr)

# Read the CSV file
df <- read.csv("Ablation2/TopSwitchers/evaluation_metrics_TOP.csv", header = TRUE)

# Reshape the data to long format
df_long <- pivot_longer(df, cols = c(ACC, AUROC, AUPRC), 
                        names_to = "Metric", values_to = "Value")

# Define custom colors
metric_colors <- c("ACC" = "#A2C2E2", 
                   "AUROC" = "#D8CBF0", 
                   "AUPRC" = "#A2DADE")



plot_top<- ggplot(df_long, aes(x = Top, y = Value, color = Metric)) +
  geom_line(size = 1) +
  scale_color_manual(values = metric_colors) +
  scale_x_continuous(
    breaks = seq(1, 30, by = 2),
    limits = c(0, 30),
    expand = c(0, 0)  # remove extra space at ends
  ) +
  theme_bw() +
  labs(
    title = "",
    x = "Top percents",
    y = "Metric Value",
    color = "Metric"
  ) +
  theme(
    text = element_text(size = 12, family = "Times"),
    legend.position = "right",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

plot_top

# Save as PDF
pdf("0---GAC/20250620/Plot_Top.pdf", width = 5, height = 3)
plot_top
dev.off()

# Save as SVG
svg("0---GAC/20250620/Plot_Top.svg", width = 5, height = 3)
plot_top
dev.off()








library(ggplot2)
library(tidyr)

# Your original data
df <- data.frame(
  Method = c("Betweenness", "Bottleneck", "Degree", "Diffusion", "Latora",
             "Lin", "Laplacian", "Local", "Leaderrank", "Leverage",
             "Residual", "Radiality", "pagerank", "Method"),
  ACC = c(0.5858, 0.54, 0.6049, 0.6048, 0.3987,
          0.5707, 0.5893, 0.5627, 0.5987, 0.5667,
          0.6707, 0.4218,0.6338, 0.8755),
  AUROC = c(0.7256, 0.722, 0.801, 0.7607, 0.6786,
            0.7302, 0.7286, 0.7241, 0.7371, 0.7386,
            0.8126, 0.6967, 0.8045, 0.923),
  AUPRC = c(0.4226, 0.4174, 0.4622, 0.4898, 0.3564,
            0.4762, 0.4401, 0.4622, 0.4543, 0.4642,
            0.6131, 0.4041, 0.5486, 0.8871)
)



# Long format
df_long <- pivot_longer(df, cols = c("ACC", "AUROC", "AUPRC"),
                        names_to = "Metric", values_to = "Score")
df_long$Metric <- factor(df_long$Metric, levels = c("AUPRC", "AUROC", "ACC"))
df_long$Method <- factor(df_long$Method, levels = c("Betweenness", "Bottleneck", "Degree", "Diffusion", "Latora",
                                                    "Lin", "Laplacian", "Local", "Leaderrank", "Leverage",
                                                    "Residual", "Radiality","pagerank", "Method"))
# Vertical heatmap
plot_Heat_Rank<-ggplot(df_long, aes(x = Method, y = Metric, fill = Score)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Score, 3)), size = 3, angle = 45) +
  scale_fill_gradient2(low = "#A2DADE", mid = "white", high = "#D8CBF0", midpoint = 0.6)+
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, family = "Times"),
    axis.text.y = element_text(family = "Times"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    title = "",
    x = "",
    y = "",
    fill = "Score"
  )
plot_Heat_Rank

# Save as PDF
pdf("0---GAC/20250620/plot_Heat_Rank2.pdf", width = 6, height = 3)
plot_Heat_Rank
dev.off()

# Save as SVG
svg("0---GAC/20250620/plot_Heat_Rank2.svg", width = 6, height = 3)
plot_Heat_Rank
dev.off()







# install.packages("UpSetR")
# install.packages("tidyverse")

library(UpSetR)
library(tidyverse)



# Create the same presence matrix
df <- read.csv("0---GAC/centrality_rankings_Pmethod.csv", stringsAsFactors = FALSE)


# Gather into long format: gene-method combinations
long_df <- df %>%
  pivot_longer(cols = everything(), names_to = "method", values_to = "gene") %>%
  filter(!is.na(gene) & gene != "") %>%
  distinct()

# Build binary matrix: rows = genes, columns = methods
binary_matrix <- long_df %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = method, values_from = value, values_fill = 0) %>%
  column_to_rownames("gene")


upset(binary_matrix,
      nsets = 13,              # Show all 13 methods
      nintersects = 30,        # Show top 30 intersection combinations
      order.by = "freq",       # Order by size of intersection
      sets.bar.color = "#A2C2E2",
      main.bar.color = "#D8CBF0",
      matrix.color = "#A2DADE",
      text.scale = 1.4)



pdf("0---GAC/20250620/gene_upset_plot_Ranking.pdf", width = 7, height = 8, family = "Times")
upset(binary_matrix,
      nsets = 13,              # Show all 13 methods
      nintersects = 30,        # Show top 30 intersection combinations
      order.by = "freq",       # Order by size of intersection
      sets.bar.color = "#A2C2E2",
      main.bar.color = "#D8CBF0",
      matrix.color = "#A2DADE",
      text.scale = 1.4)
dev.off()







library(tidyverse)
library(forcats)

# Read your data
df <- read.csv("0---GAC/ComparedMethodsResults.csv", check.names = FALSE)

# Reshape to long format
df_long <- df %>%
  pivot_longer(
    cols = -Method,
    names_to = c("Metric", "Type"),
    names_pattern = "([^_]+)_?(.*)"
  ) %>%
  mutate(Type = ifelse(Type == "", "Value", Type)) %>%
  pivot_wider(names_from = Type, values_from = value)

# Rename error column and replace NA with 0
df_long <- df_long %>% rename(Error = error)
df_long$Error[is.na(df_long$Error)] <- 0

df_long$Metric <- trimws(df_long$Metric)
# Factor levels for Metric (order of metrics)
metric_levels <- c("Specificity", "Recall", "Precision","F1", "AUPRC","AUROC","ACC")
df_long$Metric <- factor(df_long$Metric, levels = metric_levels)

# Factor levels for Method (order of methods) - change as needed
method_levels <- unique(df_long$Method)
df_long$Method <- factor(df_long$Method, levels = method_levels)

# Your color palette for bars (fill)
metric_colors <- colorRampPalette(c("#A2C2E2", "#F6B7C6", "#A2DADE", "#D8CBF0"))(length(metric_levels))

# Plot
p <- ggplot(df_long, aes(x = Method, y = Value, fill = fct_relevel(Metric, metric_levels))) +
  geom_col(position = position_dodge(0.8), width = 0.7, color = "black") +
  geom_errorbar(aes(ymin = Value - (Error), ymax = Value + ( Error)),
                position = position_dodge(0.8), width = 0.3, color = "black") +
  scale_fill_manual(values = metric_colors) +
  labs(
    x = NULL,
    y = "Metric Value",
    fill = "Metric"
  ) +
  coord_cartesian(ylim = c(0.2, NA)) +
  theme_bw(base_size = 14) +
  theme(text = element_text(family = "Times")) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 60, hjust = 1)
  )

print(p)


# Save as PDF
pdf("0---GAC/20250620/Comparision.pdf", width = 9, height = 4)
p
dev.off()

# Save as SVG
svg("0---GAC/20250620/Comparision.svg", width = 9, height = 4)
p
dev.off()





#install.packages("venn")
#install.packages("ggpolypath")
library(venn)
library(ggpolypath)
library(ggplot2)
# 1. Prepare the list (already done)
# deg_list <- list(...)  # your named list of 6 gene sets

# 2. Plot
venn_colors <- colorRampPalette(c("#A2C2E2", "#FEAF8A", "#F6B7C6", "#A2DADE", "#D8CBF0"))(6)

CompareVenn<- venn(
  deg_list,
  zcolor = venn_colors,
  opacity = 0.4,
  ilcs = 1.5,  # intersection label font size
  sncs = 1.5,  # set name font size
  box = FALSE,
  ggplot = FALSE  # <- base R plot, includes numbers
)

venn(
  deg_list,
  zcolor = venn_colors,
  opacity = 0.4,
  ilcs = 1.5,  # intersection label font size
  sncs = 1.5,  # set name font size
  box = FALSE,
  ggplot = FALSE  # <- base R plot, includes numbers
)



pdf("0---GAC/20250620/CompareVenn.pdf", width = 7, height = 8, family = "Times")
venn(
  deg_list,
  zcolor = venn_colors,
  opacity = 0.4,
  ilcs = 1.5,  # intersection label font size
  sncs = 1.5,  # set name font size
  box = FALSE,
  ggplot = FALSE  # <- base R plot, includes numbers
)
dev.off()

CompareVenn
write.csv(CompareVenn,"0---GAC/CompareVenn.csv",quote = F)









library(ggplot2)
library(tidyr)
library(dplyr)


# Sample dropout data (replace with actual results)
data_dropout <- data.frame(
  DropoutRate = c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8),
  ACC = c(0.7909, 0.8093, 0.8278, 0.8701, 0.8755, 0.861, 0.8455),
  AUROC = c(0.8066, 0.8354, 0.8593, 0.9077, 0.9230, 0.8966, 0.8766),
  AUPRC = c(0.8112, 0.8289, 0.8466, 0.8766, 0.8871, 0.864, 0.8593)
)


# Reshape to long format
data_long_dropout <- pivot_longer(
  data_dropout,
  cols = -DropoutRate,
  names_to = "Metric",
  values_to = "Value"
)

# Define custom colors
metric_colors <- c("ACC" = "#A2C2E2", "AUROC" ="#D8CBF0" , "AUPRC" = "#A2DADE")



pDropout <- ggplot(data_long_dropout, aes(x = as.numeric(as.character(DropoutRate)), y = Value, color = Metric)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = metric_colors) +
  scale_x_continuous(
    breaks = seq(0.2, 0.8, by = 0.1),
    limits = c(0.15, 0.85),
    expand = c(0, 0)  # remove extra space at ends
  ) +
  scale_y_continuous(limits = c(0.74, 0.97), expand = c(0, 0)) +
  theme_bw() +
  labs(
    title = "",
    x = "Dropout Rate",
    y = "Metric Value",
    color = "Metric"
  ) +
  theme(
    text = element_text(size = 12, family = "Times"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "right"
  )




pDropout

pdf("0---GAC/20250620/plot_Dropout.pdf",width = 5,height =2)
pDropout
dev.off()


svg("0---GAC/20250620/plot_Dropout.svg", width = 5, height = 2)
# Print the plot
print(pDropout)
# Close the SVG device to save the file
dev.off()







library(tidyverse)



# Define metric colors
metric_colors <- c("ACC" = "#A2C2E2", "AUROC" = "#D8CBF0", "AUPRC" = "#A2DADE")

# Example data (replace with actual values)
data <- data.frame(
  Embedding_size = c(32, 64, 128, 256, 512),
  ACC = c(0.841, 0.8755, 0.860, 0.857, 0.850),
  AUROC = c(0.875, 0.9230, 0.893, 0.890, 0.858),
  AUPRC = c(0.831, 0.8871, 0.856, 0.844, 0.850)
)

# Reshape data to long format
data_long <- pivot_longer(
  data,
  cols = -Embedding_size,
  names_to = "Metric",
  values_to = "Values"
)

# Ensure Embedding_size is treated as a factor for consistent labeling
data_long$Embedding_size <- factor(data_long$Embedding_size, levels = c(32, 64, 128, 256, 512))

# Plot
pEMBg <- ggplot(data_long, aes(x = as.numeric(as.factor(Embedding_size)), y = Values, color = Metric)) +
  geom_line(size = 1) +
  geom_point(size = 2, shape = 19) +  # All points are solid circles
  scale_x_continuous(
    breaks = 1:5,
    labels = c("32", "64", "128", "256", "512"),
    limits = c(0.5, 5.5),             # Add space before 32 and after 512
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0.79, 0.95),
    expand = c(0, 0)
  ) +
  scale_color_manual(values = metric_colors) +
  labs(
    x = "Embedding Size",
    y = "Metric Value",
    color = "Metric"
  ) +
  theme_bw() +
  theme(
    text = element_text(family = "Times", size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "right"
  )

pEMBg



pdf("0---GAC/20250620/Embedding.pdf",width = 5,height =2)
pEMBg
dev.off()


svg("0---GAC/20250620/Embedding.svg", width = 5, height = 2)
# Print the plot
print(pEMBg)
# Close the SVG device to save the file
dev.off()







# Load necessary library
library(readr)

# Step 1: Read the CSV file
# Replace 'your_file.csv' with your actual file name
df <- read_csv("0---GAC/Comparision_methods_30Top.csv")

# Step 2: Extract genes from each column
method_genes <- df$Method
other_genes <- unique(unlist(df[ , !names(df) %in% "Method"]))

# Step 3: Find genes in 'Method' that are not in any other column
unique_to_method <- setdiff(method_genes, other_genes)

#  Show results
cat("Genes unique to 'Method' column:\n")
print(unique_to_method)


write.csv(unique_to_method, "0---GAC/20250620/unique_genes_in_method.csv", row.names = FALSE)


# Step 4: Read the DEG/KEGG annotation txt file
deg_kegg <- read_delim("0---GAC/best_component_nodes_Reg_CD8.txt", delim = "\t", col_names = c("Gene", "Source"))
deg_kegg<- data.frame(deg_kegg)
# Step 5: Filter only DEGs
deg_only <- deg_kegg[deg_kegg$Source=="DEG",]

# Step 6: Find overlap between unique_to_method and DEG genes
unique_deg_genes <- intersect(unique_to_method, deg_only$Gene)

# Step 7: Show the result
cat("Unique genes from 'Method' column that are DEGs:\n")
print(unique_deg_genes)




# Load required libraries
library(tidyverse)

# Set working directory
setwd("E:/005---ThirdProject/ThirdObject/6.Results/")

# Set the gene name tag used in filenames (e.g., "GAC", "CD8")
GN <- "CD8"

# Define paths
base_path <- "EXPVioPlot"
expr_path <- file.path(base_path)
biomarker_file <- file.path(base_path, paste0("best_component_nodes_Reg_", GN, ".txt"))
gene_list_file <- file.path(expr_path, paste0("Union_Genes_", GN, ".csv"))

# Read biomarker gene list (replace with a manual subset if needed)
biomarkers <- read.table(biomarker_file, header = FALSE, col.names = "Gene")$Gene
#biomarkers <- c("TFF3", "KLF4", "GSPT1")  # Optional override

# Read full list of gene names
all_genes <- read.csv(gene_list_file)$Gene

# Disease stages
stages <- c("NAT", "CAG", "IM", "PGAC", "Metastasis")

# Initialize list to collect long-format expression data
expr_long_all <- list()

# Loop through stages and collect expression data
for (stage in stages) {
  file_name <- paste0(stage, "_", GN, "T_M_cleaned_Filtered_UnionGenes.csv")
  full_path <- file.path(expr_path, file_name)
  
  expr_matrix <- read.csv(full_path, row.names = 1, check.names = FALSE)
  rownames(expr_matrix) <- all_genes  # Ensure consistent rownames
  
  expr_filtered <- expr_matrix[rownames(expr_matrix) %in% biomarkers, ]
  expr_long <- as.data.frame(t(expr_filtered))  # Transpose
  
  expr_long$Sample <- rownames(expr_long)
  expr_long$Stage <- stage
  
  expr_long_all[[stage]] <- expr_long
}

# Merge all stages
expr_combined <- bind_rows(expr_long_all)

# Convert to long format
expr_long_format <- expr_combined %>%
  pivot_longer(cols = -c(Sample, Stage), names_to = "Gene", values_to = "Expression") %>%
  mutate(Stage = factor(Stage, levels = stages))








library(tidyverse)
library(ggrepel)

biomarkers <- c("GSPT1", "TFF3", "KLF4", "ITK", "IL8", "ETS2", "GNB2L1", "STAT1", "FLI1", "GNG5", "RHOA")

# biomarkers2 <- read.table(biomarker_file, header = FALSE, col.names = "Gene")$Gene
# 
# # Find genes in biomarkers2 but not in biomarkers
# biomarkers3 <- setdiff(biomarkers2, biomarkers)
# 
# # Print result
# print(biomarkers3)
# biomarkers<-biomarkers3

stages <- c("NAT", "CAG", "IM", "PGAC", "Metastasis")

# Filter data for biomarkers only
filtered_data <- expr_long_format %>% filter(Gene %in% biomarkers)

# Compute mean with condition:
# If all values zero → mean = 0
# Else mean of non-zero values only
expression_summary <- filtered_data %>%
  group_by(Gene, Stage) %>%
  summarise(
    MeanExpr = if (all(Expression == 0)) {
      0
    } else {
      mean(Expression[Expression > 0], na.rm = TRUE)
    },
    .groups = "drop"
  )

# Create full grid of gene-stage pairs to ensure completeness
full_grid <- expand.grid(Gene = biomarkers, Stage = stages)

# Left join and fill missing with 0
expression_complete <- full_grid %>%
  left_join(expression_summary, by = c("Gene", "Stage")) %>%
  mutate(MeanExpr = if_else(is.na(MeanExpr), 0, MeanExpr)) %>%
  mutate(Stage = factor(Stage, levels = stages))

# Prepare colors
metric_levels <- unique(expression_complete$Gene)
#metric_colors <- colorRampPalette(c("#5986AB", "#CC6784", "#48989D", "#9B81C3"))(length(metric_levels))
metric_colors <- colorRampPalette(c("#FBCC6B", "#A2C2E2", "#FEAF8A", "#F6B7C6", "#A2DADE", "#D8CBF0"))(length(metric_levels))
# Plot
p <- ggplot(expression_complete, aes(x = Stage, y = MeanExpr, group = Gene, color = Gene)) +
  geom_line(size = 1.25) +
  geom_point(size = 3) +
  geom_text_repel(
    data = expression_complete %>% group_by(Gene) %>% filter(Stage == last(Stage)),
    aes(label = Gene),
    hjust = -0.1, size = 3, family = "Times", segment.color = "gray70"
  ) +
  scale_y_continuous(limits = c(0, 4)) +
  scale_color_manual(values = metric_colors) +
  theme_bw() +
  theme(
    text = element_text(family = "Times", size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "none"
  ) +
  labs(
    title = "",
    x = "Stage",
    y = "Mean Expression"
  )

print(p)



# Save as PDF
pdf("EXPVioPlot/Changes1.pdf", width = 9, height = 4)
p
dev.off()

# Save as SVG
svg("EXPVioPlot/Changes1.svg", width = 9, height = 4)
p
dev.off()

