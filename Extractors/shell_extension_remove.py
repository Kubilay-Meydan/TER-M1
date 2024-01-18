import os

def renommer_fichiers(dossier):
    suffixe_a_supprimer = "_tools_in_shell"

    # Parcourir tous les fichiers dans le dossier spécifié
    for nom_fichier in os.listdir(dossier):
        if nom_fichier.endswith('.txt') and suffixe_a_supprimer in nom_fichier:
            chemin_ancien_fichier = os.path.join(dossier, nom_fichier)
            # Construire le nouveau nom de fichier
            nouveau_nom = nom_fichier.replace(suffixe_a_supprimer, '') 
            chemin_nouveau_fichier = os.path.join(dossier, nouveau_nom)

            # Renommer le fichier
            os.rename(chemin_ancien_fichier, chemin_nouveau_fichier)
            print(f"Le fichier {nom_fichier} a été renommé en {nouveau_nom}")

# Demander à l'utilisateur de saisir le chemin du dossier
dossier = input("Veuillez saisir le chemin du dossier : ")
renommer_fichiers(dossier)
