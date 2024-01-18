import os
import pandas as pd
import Levenshtein

def read_text_files(folder):
    files = [f for f in os.listdir(folder)]
    texts = []
    for file in files:
        with open(os.path.join(folder, file), 'r', encoding='utf-8') as file:
            texts.append(file.read().strip())  # Changed to read() to get the full text
    return texts, files

def calculate_levenshtein_similarity(text1, text2):
    distance = Levenshtein.distance(text1, text2)
    max_length = max(len(text1), len(text2))
    if max_length == 0:  # To avoid division by zero
        return 1
    return 1 - distance / max_length

def build_and_save_levenshtein_similarity_matrix(folder, csv_file):
    texts, files = read_text_files(folder)
    pairs = []
    scores = []

    for i in range(len(files)):
        for j in range(i+1, len(files)):  # Only comparing each pair once
            similarity_score = calculate_levenshtein_similarity(texts[i], texts[j])
            pairs.append(f"{files[i]},{files[j]}")
            scores.append(similarity_score)

    df = pd.DataFrame({'Pair': pairs, 'Similarity': scores})
    df.to_csv(csv_file, index=False)
