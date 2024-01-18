import os

def ajouter_extension_txt(dossier):
    # Parcourir tous les fichiers dans le dossier spécifié
    for nom_fichier in os.listdir(dossier):
        chemin_ancien_fichier = os.path.join(dossier, nom_fichier)
        # Ignorer les dossiers
        if os.path.isfile(chemin_ancien_fichier):
            chemin_nouveau_fichier = chemin_ancien_fichier + ".txt"
            # Renommer le fichier pour ajouter l'extension .txt
            os.rename(chemin_ancien_fichier, chemin_nouveau_fichier)
            print(f"Le fichier {nom_fichier} a été renommé en {nom_fichier}.txt")

# Demander à l'utilisateur de saisir le chemin du dossier
dossier = input("Veuillez saisir le chemin du dossier : ")
ajouter_extension_txt(dossier)
