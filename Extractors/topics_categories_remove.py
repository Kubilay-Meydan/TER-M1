import os

def nettoyer_contenu_fichiers(dossier):
    elements_a_supprimer = ['"', ',', 'operation_', 'topic_']

    # Parcourir tous les fichiers dans le dossier spécifié
    for nom_fichier in os.listdir(dossier):
        if nom_fichier.endswith('.txt'):
            chemin_fichier = os.path.join(dossier, nom_fichier)

            # Lire le contenu du fichier
            with open(chemin_fichier, 'r') as fichier:
                contenu = fichier.read()

            # Supprimer les éléments spécifiés
            for element in elements_a_supprimer:
                contenu = contenu.replace(element, "")

            # Réécrire le fichier avec le contenu modifié
            with open(chemin_fichier, 'w') as fichier:
                fichier.write(contenu)
            print(f"Le fichier {nom_fichier} a été nettoyé.")

# Demander à l'utilisateur de saisir le chemin du dossier
dossier = input("Veuillez saisir le chemin du dossier contenant les fichiers .txt : ")
nettoyer_contenu_fichiers(dossier)
