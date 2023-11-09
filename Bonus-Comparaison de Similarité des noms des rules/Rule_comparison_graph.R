# Load necessary packages
library(ggplot2)
library(dplyr)

# Read the CSV file
data <- read.csv("comparison_results.csv")

# Identify words that only have 'dissimilar' as their similarity level
words_with_only_dissimilar <- data %>%
  group_by(Word1) %>%
  filter(all(Similarity_Level == "dissimilar")) %>%
  select(Word1) %>%
  unique() %>%
  ungroup() %>%
  bind_rows(
    data %>%
      group_by(Word2) %>%
      filter(all(Similarity_Level == "dissimilar")) %>%
      select(Word1 = Word2) %>%
      unique() %>%
      ungroup()
  ) %>%
  pull(Word1)

# Filter out 'dissimilar' rows and words that only have 'dissimilar' associations
filtered_data <- data %>%
  filter(Similarity_Level != "dissimilar", 
         !Word1 %in% words_with_only_dissimilar, 
         !Word2 %in% words_with_only_dissimilar)

# Get an ordered list of unique words
unique_words <- unique(c(filtered_data$Word1, filtered_data$Word2))
unique_words <- sort(unique_words)

# Convert Similarity_Level to an ordered factor
filtered_data$Similarity_Level <- factor(filtered_data$Similarity_Level, 
                                         levels = c("vaguely similar", "similar", "identical"),
                                         ordered = TRUE)

# Convert Word1 and Word2 to factors with the same levels
filtered_data$Word1 <- factor(filtered_data$Word1, levels = unique_words)
filtered_data$Word2 <- factor(filtered_data$Word2, levels = unique_words)

# Plot heatmap with discrete colors and sorted axes
ggplot(filtered_data, aes(x = Word1, y = Word2, fill = Similarity_Level)) +
  geom_tile() +
  scale_fill_manual(values = c("vaguely similar" = "#ADD8E6", "similar" = "#FDB45C", "identical" = "#2E8B57")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        axis.text.y = element_text(hjust = 1)) +
  labs(x = "Word1", y = "Word2", fill = "Similarity Level") +
  scale_x_discrete(limits = unique_words) +
  scale_y_discrete(limits = unique_words)
