# Load required packages
# install.packages("scales")
# install.packages("ggplot2")
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)


##### Violin Plot #####
# Define stages
stages <- c("NAT", "CAG", "IM", "PGAC", "Metastasis")

# Read and reshape data
combined_data <- lapply(stages, function(stage) {
  df <- read.csv(paste0("DATA/Gastric/FeatureMat/", stage, "_CD8T_M_cleaned_Filtered_UnionGenes.csv"))
  if (!"Gene" %in% colnames(df)) {
    colnames(df)[1] <- "Gene"
  }
  df_long <- df %>%
    pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Expression")
  df_long$Stage <- stage
  return(df_long)
}) %>% bind_rows()

# Make sure it's numeric
combined_data$Expression <- as.numeric(as.character(combined_data$Expression))

# Filter and log2-transform
combined_data <- combined_data %>% filter(Expression > 0)
combined_data$Log2Expr <- log2(combined_data$Expression)


# Set factor levels
combined_data$Stage <- factor(combined_data$Stage, levels = stages)

# Define colors
stage_colors <- c("#A2C2E2", "#FEAF8A", "#F6B7C6", "#A2DADE", "#D8CBF0")

# Plot
ggplot(combined_data, aes(x = Stage, y = Log2Expr, fill = Stage)) +
  geom_violin(trim = TRUE, scale = "width", color = "black") +
  #geom_boxplot(width = 0.1, fill = "white", alpha = 0.3, outlier.size = 0.5)+
  scale_fill_manual(values = stage_colors) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(size = 12)
  ) +
  ylab("log2(Expression)")




##### DEGs Upset and venn PLOTS ######
# --- Load libraries ---
library(ggplot2)
library(dplyr)
library(UpSetR)

# --- Define stages and color palette ---
stages <- c("CAG", "IM", "PGAC", "Metastasis")
my_colors <- c("CAG" = "#FEAF8A", "IM" = "#F6B7C6", "PGAC" = "#A2DADE", "Metastasis" = "#D8CBF0")

# --- Read DEG files and extract gene sets ---
deg_list <- list()
deg_counts <- data.frame(Stage = character(), Count = integer(), stringsAsFactors = FALSE)

for (stage in stages) {
  file_path <- paste0("DATA/Gastric/DEGs/CD8_DEGs_", stage, "_vs_NAT.csv")
  df <- read.csv(file_path)
  deg_genes <- df$gene
  deg_list[[stage]] <- unique(deg_genes)
  deg_counts <- rbind(deg_counts, data.frame(Stage = stage, Count = length(deg_genes)))
}

# --- Set Stage as factor to preserve order ---
deg_counts$Stage <- factor(deg_counts$Stage, levels = stages)

# --- Prepare data for proportional horizontal bar plot ---
deg_counts <- deg_counts %>%
  mutate(
    fraction = Count / sum(Count),
    xmin = c(0, head(cumsum(fraction), -1)),
    xmax = cumsum(fraction),
    label_pos = (xmin + xmax) / 2
  )

# --- Draw proportional horizontal bar plot ---
ggplot(deg_counts) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1, fill = Stage), color = "white") +
  geom_text(aes(x = label_pos, y = 0.5, label = Count), color = "black", size = 5) +
  scale_fill_manual(values = my_colors) +
  theme_void() +
  theme(legend.position = "bottom") +
  labs(title = "Number of DEGs per Stage", fill = "Stage")

# --- Create binary matrix for UpSet plot ---
all_genes <- unique(unlist(deg_list))
binary_matrix <- sapply(stages, function(stage) {
  as.integer(all_genes %in% deg_list[[stage]])
})
colnames(binary_matrix) <- stages
rownames(binary_matrix) <- all_genes
binary_df <- as.data.frame(binary_matrix)

# --- Draw UpSet plot showing DEG overlap ---
CD8_upset <- upset(binary_df, 
      sets = stages,
      sets.bar.color = my_colors[stages],
      order.by = "freq",
      main.bar.color = "#4878A8",
      matrix.color = "#A2C2E2",
      text.scale = 1.3)

# Save to PDF
pdf("PLOTS/DEG_CD8_upset.pdf", width = 4, height = 4)  # You can adjust size
CD8_upset
dev.off()  # Close the device


# Save to SVG
svg("PLOTS/DEG_CD8_upset.svg", width = 4, height = 4)
CD8_upset
dev.off()



# --- Install  ---
#install.packages("VennDiagram")

library(VennDiagram)

# --- Create Venn diagram and save to file ---
venn.plot <- venn.diagram(
  x = deg_list,
  category.names = names(deg_list),
  filename = NULL,
  fill = c("#FEAF8A", "#F6B7C6", "#A2DADE", "#D8CBF0"),
  cex = 1.5,
  cat.cex = 1.5,
  main = "DEG Overlap Between Stages"
)

# --- Plot to screen ---
grid::grid.newpage()
grid::grid.draw(venn.plot)


# Save to PDF
pdf("PLOTS/DEG_CD8_VennDiagram.pdf", width = 4, height = 4)  # You can adjust size
grid::grid.newpage()
grid::grid.draw(venn.plot)
dev.off()  # Close the device


# Save to SVG
svg("PLOTS/DEG_CD8_VennDiagram.svg", width = 4, height = 4)
grid::grid.newpage()
grid::grid.draw(venn.plot)
dev.off()





##### Donut Plot #####
# Load required packages
# install.packages("scales")
# install.packages("ggplot2")
library(ggplot2)
library(dplyr)

# Sample data
data <- data.frame(
  Stage = c("NAT", "CAG", "IM", "PGAC", "Metastasis"),
  count = c(654, 907, 884, 386, 1611)
)

# Fix the order of groups
data$Stage <- factor(data$Stage, levels = c("NAT", "CAG", "IM", "PGAC", "Metastasis"))

# Prepare data
data <- data %>%
  arrange(Stage) %>%  # now arranges by the specified factor levels
  mutate(
    fraction = count / sum(count),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n = -1)),
    label_pos = (ymax + ymin) / 2,
    label = paste0(count)
  )

# Draw donut plot
P_Donut_CD8 <- ggplot(data, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = Stage)) +
  geom_rect(color = "black") +
  coord_polar(theta = "y") +
  xlim(c(2, 4)) +
  theme_void() +
  geom_text(aes(x = 3.5, y = label_pos, label = label), color = "black", size = 4) +
  annotate("text", x = 0, y = 0, label = paste0("Total = ", sum(data$count)), size = 5) +
  scale_fill_manual(values = c("#A2C2E2", "#FEAF8A", "#F6B7C6", "#A2DADE", "#D8CBF0"))


P_Donut_CD8



pdf("PLOTS/Plot_Donut_CD8.pdf",width = 4,height =4)
P_Donut_CD8
dev.off()

svg("PLOTS/Plot_Donut_CD8.svg", width = 4, height = 4)
# Print the plot
print(P_Donut_CD8)
# Close the SVG device to save the file
dev.off()

