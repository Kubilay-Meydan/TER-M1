import os
import difflib
import pandas as pd

# Adjust the directory path to where your .txt files are located
directory_path = 'Extracted_Rules'

# Load words from all files into a set
def load_all_words(directory_path):
    all_words = set()
    for filename in os.listdir(directory_path):
        if filename.endswith('.txt'):
            filepath = os.path.join(directory_path, filename)
            with open(filepath, 'r') as file:
                all_words.update(file.read().splitlines())
    return all_words

# Compare every word with every other word
def compare_all_words(words):
    comparison_results = []
    for word1 in words:
        for word2 in words:
            ratio = difflib.SequenceMatcher(None, word1, word2).ratio()
            if ratio == 1:
                similarity_level = 'identical'
            elif ratio >= 0.8:
                similarity_level = 'similar'
            elif ratio >= 0.6:
                similarity_level = 'vaguely similar'
            else:
                similarity_level = 'dissimilar'
            comparison_results.append({
                'Word1': word1,
                'Word2': word2,
                'Similarity_Level': similarity_level
            })
    return comparison_results

# Load all words from the directory
all_words = load_all_words(directory_path)

# Perform the comparison
comparison_results = compare_all_words(all_words)

# Convert the results into a DataFrame
df = pd.DataFrame(comparison_results)

# Sort the DataFrame for better readability
df.sort_values(by=['Similarity_Level', 'Word1', 'Word2'], inplace=True)

# Save the DataFrame to a CSV file
output_csv_path = os.path.join(directory_path, 'comparison_results.csv')
df.to_csv(output_csv_path, index=False)

# Print the path to the CSV file
print(f"The comparison results have been saved to: {output_csv_path}")
