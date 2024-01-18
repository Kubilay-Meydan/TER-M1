import os
import pandas as pd

def read_text_files(folder):
    files = [f for f in os.listdir(folder) if f.endswith('.txt')]
    texts = []
    for file in files:
        with open(os.path.join(folder, file), 'r', encoding='utf-8') as file:
            texts.append(set(file.read().strip().split('\n')))
    return texts, files

def calculate_jaccard_similarity(set1, set2):
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    if union == 0:
        return 1
    return intersection / union

def build_and_save_jaccard_similarity_matrix(folder):
    sets, files = read_text_files(folder)
    pairs = []
    scores = []

    for i in range(len(files)):
        for j in range(i+1, len(files)):
            similarity_score = calculate_jaccard_similarity(sets[i], sets[j])
            pairs.append(f"{files[i]},{files[j]}")
            scores.append(similarity_score)

    # Déterminer le nom du fichier de sortie basé sur le nom du dossier
    folder_name = os.path.basename(folder)
    csv_file = f"{folder_name}_jaccard.csv"
    df = pd.DataFrame({'Pair': pairs, 'Similarity': scores})
    df.to_csv(csv_file, index=False)
    print(f"Les résultats ont été sauvegardés dans {csv_file}")

# Demander à l'utilisateur de saisir le chemin du dossier
folder = input("Veuillez saisir le chemin du dossier : ")
build_and_save_jaccard_similarity_matrix(folder)
