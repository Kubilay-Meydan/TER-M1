import os
import pandas as pd

def lire_fichiers_texte(dossier):
    fichiers = [f for f in os.listdir(dossier)]
    nombre_lignes = []
    for fichier in fichiers:
        with open(os.path.join(dossier, fichier), 'r', encoding='utf-8') as file:
            lignes = file.readlines()
            nombre_lignes.append(len(lignes))
    return nombre_lignes, fichiers

def construire_matrice_similarite_et_sauvegarder(dossier, fichier_csv):
    nombre_lignes, fichiers = lire_fichiers_texte(dossier)
    
    # Préparation des données pour le DataFrame
    paires = []
    similarites = []
    for i in range(len(fichiers)):
        for j in range(i+1, len(fichiers)):
            similarite = 1 if nombre_lignes[i] == nombre_lignes[j] else 0
            paires.append(f"{fichiers[i]},{fichiers[j]}")
            similarites.append(similarite)
    
    # Création d'un DataFrame et sauvegarde en CSV
    df = pd.DataFrame({'Pair': paires, 'Similarity': similarites})
    df.to_csv(fichier_csv, index=False)

# Utiliser la fonction
dossier = 'Attributes/input_type'  # Remplacer par le chemin vers votre dossier
fichier_csv = 'input_number_similarity.csv'  # Nom du fichier CSV de sortie
construire_matrice_similarite_et_sauvegarder(dossier, fichier_csv)
