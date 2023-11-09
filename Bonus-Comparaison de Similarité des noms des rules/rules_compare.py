import os
import difflib
import pandas as pd
from collections import defaultdict

# Adjust the directory path to where your .txt files are located
directory_path = 'Extracted_Rules'

# Load words from a file into a set
def load_words_from_file(file_path):
    with open(file_path, 'r') as file:
        return set(file.read().splitlines())

# Compare two sets of words and group them by similarity
def compare_words(set1, set2):
    groups = defaultdict(list)
    for word1 in set1:
        for word2 in set2:
            ratio = difflib.SequenceMatcher(None, word1, word2).ratio()
            if ratio == 1:
                groups['identical'].append((word1, word2))
            elif ratio >= 0.8:
                groups['similar'].append((word1, word2))
            elif ratio >= 0.6:
                groups['vaguely similar'].append((word1, word2))
            else:
                groups['dissimilar'].append((word1, word2))
    return groups

# List all txt files in the directory
txt_files = [f for f in os.listdir(directory_path) if f.endswith('.txt')]

# Load words from all files
words_per_file = {file: load_words_from_file(os.path.join(directory_path, file)) for file in txt_files}

# Compare words between all files and organize into a DataFrame
dataframes = []

for i, file1 in enumerate(txt_files):
    for file2 in txt_files[i+1:]:
        comparison = compare_words(words_per_file[file1], words_per_file[file2])
        for similarity_level, word_pairs in comparison.items():
            for word1, word2 in word_pairs:
                dataframes.append({
                    'File1': file1,
                    'File2': file2,
                    'Word1': word1,
                    'Word2': word2,
                    'Similarity_Level': similarity_level
                })

# Convert the results into a DataFrame
df = pd.DataFrame(dataframes)

# Sort the DataFrame for better readability
df.sort_values(by=['Similarity_Level', 'File1', 'File2'], inplace=True)

# Save the DataFrame to a CSV file
output_csv_path = os.path.join(directory_path, 'comparison_results.csv')
df.to_csv(output_csv_path, index=False)

# Print the path to the CSV file
print(f"The comparison results have been saved to: {output_csv_path}")
