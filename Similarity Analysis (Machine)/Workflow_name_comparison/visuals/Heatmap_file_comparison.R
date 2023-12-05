# Load necessary libraries
library(ggplot2)
library(reshape2)

# Read the CSV file
data <- read.csv("file_comparison_matrix.csv", row.names = 1)

# Convert similarity values to factors
data_factor <- apply(data, c(1,2), factor, levels = c("D", "SS", "S", "VS", "I"))

# Convert the matrix to a format suitable for ggplot
data_melted <- melt(as.matrix(data_factor))

# Define colors for the similarity levels
colors <- c("D" = "cornflowerblue", "SS" = "lightblue", "S" = "cyan", "VS" = "aquamarine", "I" = "green")

# Create the heatmap with labeled color scale
ggplot(data_melted, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_manual(values = colors, name = "Similarity Level", breaks = c("I", "VS", "S", "SS", "D")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Column") +
  ylab("Row") +
  ggtitle("Heatmap of Similarity Matrix of the QC and Trimming Modules")
