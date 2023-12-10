
import pandas as pd
import numpy as np

#Charger les matrices à partir des fichiers CSV
matrice_1 = pd.read_csv('/Users/natha/Desktop/Cours/TER-M1/TER-M1/Annotation_Nathan/8_Vizualisation/file_comparison_matrix.csv', index_col=0)
matrice_2 = pd.read_csv('/Users/natha/Desktop/Cours/TER-M1/TER-M1/Annotation_Kubilay/Visualisation/Visualisation.csv', index_col=0)

#Trouver les règles communes (intersection des indices et des colonnes)
indices_communs = matrice_1.index.intersection(matrice_2.index)
colonnes_communes = matrice_1.columns.intersection(matrice_2.columns)

#Extraire les sous-matrices pour les éléments communs
sous_matrice_1 = matrice_1.loc[indices_communs, colonnes_communes]
sous_matrice_2 = matrice_2.loc[indices_communs, colonnes_communes]

#Identifier où les sous-matrices diffèrent
differences = sous_matrice_1.values != sous_matrice_2.values
lignes_diff, colonnes_diff = np.where(differences)

#Compter le nombre d'éléments différents dans les sous-matrices
difference = np.sum(differences)

#Compter tous les éléments qui ne sont pas communs
elements_non_communs = (matrice_1.size + matrice_2.size) - (sous_matrice_1.size + sous_matrice_2.size)

#Calculer le pourcentage de différence
pourcentage_difference = ((difference + elements_non_communs) / (matrice_1.size + matrice_2.size)) * 100

#Afficher le pourcentage de différence
print(f"Pourcentage de différence : {pourcentage_difference:.2f}%")

# Afficher les différences dans les sous-matrices communes
for ligne, colonne in zip(lignes_diff, colonnes_diff):
    print(f"Différence trouvée à la ligne '{indices_communs[ligne]}' et à la colonne '{colonnes_communes[colonne]}'")

#Identifier et afficher les éléments non communs
indices_non_communs = matrice_1.index.symmetric_difference(matrice_2.index)
colonnes_non_communes = matrice_1.columns.symmetric_difference(matrice_2.columns)

print("Lignes non communes :", indices_non_communs.tolist())
print("Colonnes non communes :", colonnes_non_communes.tolist())
