import os
import pandas as pd
import Levenshtein

def read_text_files(folder):
    files = [f for f in os.listdir(folder) if f.endswith('.txt')]
    texts = []
    for file in files:
        with open(os.path.join(folder, file), 'r', encoding='utf-8') as file:
            texts.append(file.readline().strip())
    return texts, files

def calculate_levenshtein_similarity(text1, text2):
    distance = Levenshtein.distance(text1, text2)
    max_length = max(len(text1), len(text2))
    if max_length == 0:  # To avoid division by zero
        return 1
    return 1 - distance / max_length

def build_and_save_levenshtein_similarity_matrix(folder, csv_file):
    texts, files = read_text_files(folder)
    n = len(texts)
    similarity_matrix = [[0] * n for _ in range(n)]

    for i in range(n):
        for j in range(n):
            similarity_matrix[i][j] = calculate_levenshtein_similarity(texts[i], texts[j])
    
    df = pd.DataFrame(similarity_matrix, index=files, columns=files)
    df.to_csv(csv_file)

# Using the function
folder = 'Attributes/rule_name'  # Replace with the path to your folder
csv_file = 'levenshtein_similarity_matrix.csv'  # Output CSV file name
build_and_save_levenshtein_similarity_matrix(folder, csv_file)
