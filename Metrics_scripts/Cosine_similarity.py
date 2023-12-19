import os
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import cosine_similarity
import pandas as pd

def lire_fichiers_texte(dossier):
    fichiers = [f for f in os.listdir(dossier) if f.endswith('.txt')]
    textes = []
    for fichier in fichiers:
        with open(os.path.join(dossier, fichier), 'r', encoding='utf-8') as file:
            textes.append(file.readline())
    return textes, fichiers

def construire_matrice_similarite_et_sauvegarder(dossier, fichier_csv):
    textes, fichiers = lire_fichiers_texte(dossier)
    vecteur = TfidfVectorizer()
    X = vecteur.fit_transform(textes)
    matrice_similarite = cosine_similarity(X)
    
    # Création d'un DataFrame et sauvegarde en CSV
    df = pd.DataFrame(matrice_similarite, index=fichiers, columns=fichiers)
    df.to_csv(fichier_csv)

# Utiliser la fonction
dossier = 'Attributes/output_type'  # Remplacer par le chemin vers votre dossier
fichier_csv = 'output_similarity_cosine.csv'  # Nom du fichier CSV de sortie
construire_matrice_similarite_et_sauvegarder(dossier, fichier_csv)
