import os
import shutil

def deplacer_fichiers_txt(dossier_principal):
    # Parcourir tous les sous-dossiers et fichiers dans le dossier principal
    for dossier, sous_dossiers, fichiers in os.walk(dossier_principal):
        for nom_fichier in fichiers:
            # Vérifier si le fichier est un fichier .txt
            if nom_fichier.endswith('.txt'):
                chemin_complet = os.path.join(dossier, nom_fichier)
                nouvelle_destination = os.path.join(dossier_principal, nom_fichier)

                # Vérifier si un fichier avec le même nom existe déjà dans le dossier principal
                if os.path.exists(nouvelle_destination):
                    print(f"Le fichier {nom_fichier} existe déjà dans le dossier principal. Renommage nécessaire.")
                    base, extension = os.path.splitext(nom_fichier)
                    i = 1
                    nouveau_nom = f"{base}_{i}{extension}"
                    nouvelle_destination = os.path.join(dossier_principal, nouveau_nom)

                    # Trouver un nouveau nom de fichier qui n'existe pas encore
                    while os.path.exists(nouvelle_destination):
                        i += 1
                        nouveau_nom = f"{base}_{i}{extension}"
                        nouvelle_destination = os.path.join(dossier_principal, nouveau_nom)

                # Déplacer le fichier
                shutil.move(chemin_complet, nouvelle_destination)
                print(f"Le fichier {nom_fichier} a été déplacé vers {nouvelle_destination}")

# Demander à l'utilisateur de saisir le chemin du dossier principal
dossier_principal = input("Veuillez saisir le chemin du dossier principal : ")
deplacer_fichiers_txt(dossier_principal)
