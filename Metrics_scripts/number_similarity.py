import os
import pandas as pd

def lire_fichiers_texte(dossier):
    fichiers = [f for f in os.listdir(dossier) if f.endswith('.txt')]
    nombre_lignes = []
    for fichier in fichiers:
        with open(os.path.join(dossier, fichier), 'r', encoding='utf-8') as file:
            lignes = file.readlines()
            nombre_lignes.append(len(lignes))
    return nombre_lignes, fichiers

def construire_matrice_similarite_et_sauvegarder(dossier, fichier_csv):
    nombre_lignes, fichiers = lire_fichiers_texte(dossier)
    
    # Création de la matrice de similarité
    taille = len(fichiers)
    matrice_similarite = [[1 if nombre_lignes[i] == nombre_lignes[j] else 0 for j in range(taille)] for i in range(taille)]
    
    # Création d'un DataFrame et sauvegarde en CSV
    df = pd.DataFrame(matrice_similarite, index=fichiers, columns=fichiers)
    df.to_csv(fichier_csv)

# Utiliser la fonction
dossier = '/Users/kubilaymeydan/Desktop/M1 Bibs/S1/TER/TER-M1/Attributes/output_type'  # Remplacer par le chemin vers votre dossier
fichier_csv = 'output_number_similarity.csv'  # Nom du fichier CSV de sortie
construire_matrice_similarite_et_sauvegarder(dossier, fichier_csv)
