library(ggplot2)
library(reshape2)

# Function to read CSV and convert to factor matrix
read_and_factorize <- function(file) {
  data <- read.csv(file, row.names = 1)
  return(apply(data, c(1, 2), factor, levels = c("D", "SS", "S", "VS", "I")))
}

# Read and process each CSV file
data1 <- read_and_factorize("QC.csv")
data2 <- read_and_factorize("trimming.csv")
data3 <- read_and_factorize("Quantification.csv")
data4 <- read_and_factorize("Post_Align_QC.csv")
data5 <- read_and_factorize("mapping.csv")
data6 <- read_and_factorize("Normalisation.csv")
data7 <- read_and_factorize("Visualisation.csv")

# Combine the row and column names from all matrices
all_rows <- unique(c(rownames(data1), rownames(data2), rownames(data3), rownames(data4), rownames(data5), rownames(data6), rownames(data7)))
all_cols <- unique(c(colnames(data1), colnames(data2), colnames(data3), colnames(data4), colnames(data5), colnames(data6), colnames(data7)))

# Create an empty matrix with combined dimensions
combined_matrix <- matrix("D", nrow = length(all_rows), ncol = length(all_cols),
                          dimnames = list(all_rows, all_cols))

# Function to fill the combined matrix
fill_matrix <- function(data, combined_matrix) {
  for (row in rownames(data)) {
    for (col in colnames(data)) {
      combined_matrix[row, col] <- as.character(data[row, col])
    }
  }
  return(combined_matrix)
}

# Fill the combined matrix with each data matrix
combined_matrix <- fill_matrix(data1, combined_matrix)
combined_matrix <- fill_matrix(data2, combined_matrix)
combined_matrix <- fill_matrix(data3, combined_matrix)
combined_matrix <- fill_matrix(data4, combined_matrix)
combined_matrix <- fill_matrix(data5, combined_matrix)
combined_matrix <- fill_matrix(data6, combined_matrix)
combined_matrix <- fill_matrix(data7, combined_matrix)

# Convert the combined matrix to a factor matrix
combined_matrix_factor <- apply(combined_matrix, c(1, 2), factor, levels = c("D", "SS", "S", "VS", "I"))

# Melt the matrix for ggplot
data_melted <- melt(as.matrix(combined_matrix_factor))

# Define colors for the similarity levels
colors <- c("D" = "cornflowerblue", "SS" = "lightblue", "S" = "cyan", "VS" = "aquamarine", "I" = "green")

# Create the heatmap
ggplot(data_melted, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_manual(values = colors, name = "Similarity Level", breaks = c("I", "VS", "S", "SS", "D")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Column") +
  ylab("Row") +
  ggtitle("Heatmap of Combined Similarity Matrix")
