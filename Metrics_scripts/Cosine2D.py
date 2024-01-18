import os
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import cosine_similarity
import pandas as pd

def lire_fichiers_texte(dossier):
    fichiers = [f for f in os.listdir(dossier)]
    textes = []
    for fichier in fichiers:
        with open(os.path.join(dossier, fichier), 'r', encoding='utf-8') as file:
            textes.append(file.read())
    return textes, fichiers

def construire_matrice_similarite_et_sauvegarder(dossier, fichier_csv):
    textes, fichiers = lire_fichiers_texte(dossier)
    vecteur = TfidfVectorizer()
    X = vecteur.fit_transform(textes)
    matrice_similarite = cosine_similarity(X)

    # Préparation des données pour le DataFrame
    paires = []
    scores = []
    for i in range(len(fichiers)):
        for j in range(i+1, len(fichiers)):
            paires.append(f"{fichiers[i]},{fichiers[j]}")
            scores.append(matrice_similarite[i][j])
    
    # Création d'un DataFrame et sauvegarde en CSV
    df = pd.DataFrame({'Pair': paires, 'Similarity': scores})
    df.to_csv(fichier_csv, index=False)

