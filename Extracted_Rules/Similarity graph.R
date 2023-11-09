# Install and load required packages
if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")
library(igraph)

# Read the CSV file
df <- read.csv("comparison_results.csv")

# Filter to include only words with 'similar' similarity levels
edge_list <- subset(df, Similarity_Level == "similar", 
                    select = c("Word1", "Word2", "Similarity_Level"))

# Create an undirected graph from the edge list
graph <- graph_from_data_frame(edge_list, directed = FALSE)

# Set edge attributes based on similarity level
# Since we only have 'similar' edges, we can assign a uniform color and weight
E(graph)$weight <- 1
E(graph)$color <- "blue"

# Plot the graph with igraph
plot(graph, layout  = layout_with_fr, edge.color = E(graph)$color, vertex.size = 5, area = vcount(graph)^100)