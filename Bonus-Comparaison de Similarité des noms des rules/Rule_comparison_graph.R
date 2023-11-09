library(dplyr)
library(ggplot2)
library(readr)



# Load the data from the CSV file
comparison_data <- read_csv("comparison_results.csv")

# Count the number of occurrences in each similarity level for each file pair
comparison_counts <- comparison_data %>%
  group_by(File1, File2, Similarity_Level) %>%
  summarise(Count = n(), .groups = 'drop') # Added .groups = 'drop' to avoid a warning

# Create a bar plot
ggplot(comparison_counts, aes(x = interaction(File1, File2), y = Count, fill = Similarity_Level)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "File Pairs", y = "Count of Word Pairs", fill = "Similarity Level") +
  ggtitle("Comparison of Word Pairs by Similarity Level")

# Save the plot as an image
ggsave("word_similarity_comparison.png", width = 12, height = 8)
