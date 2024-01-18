import os

def supprimer_espaces_et_sauts_de_ligne(dossier):
    # Parcourir tous les fichiers dans le dossier spécifié
    for nom_fichier in os.listdir(dossier):
        if nom_fichier.endswith('.txt'):
            chemin_fichier = os.path.join(dossier, nom_fichier)

            # Lire le contenu du fichier
            with open(chemin_fichier, 'r') as fichier:
                contenu = fichier.read()

            # Supprimer les espaces et les sauts de ligne
            contenu_modifie = contenu.replace(" ", "").replace("\n", "")

            # Réécrire le fichier avec le contenu modifié
            with open(chemin_fichier, 'w') as fichier:
                fichier.write(contenu_modifie)
            print(f"Le fichier {nom_fichier} a été modifié.")

# Demander à l'utilisateur de saisir le chemin du dossier
dossier = input("Veuillez saisir le chemin du dossier contenant les fichiers .txt : ")
supprimer_espaces_et_sauts_de_ligne(dossier)
