import os

def copy_structure_and_convert(src_folder, dest_folder):
    for root, dirs, files in os.walk(src_folder):
        # Construire le chemin de destination
        dest_path = root.replace(src_folder, dest_folder, 1)

        # Créer le dossier de destination s'il n'existe pas
        if not os.path.exists(dest_path):
            os.makedirs(dest_path)

        # Pour chaque fichier .py, créer un fichier .txt correspondant
        for file in files:
            if file.endswith('.py'):
                # Construire le chemin complet du nouveau fichier .txt
                new_file_path = os.path.join(dest_path, file[:-3] + '.txt')
                # Créer un fichier .txt vide
                open(new_file_path, 'w').close()

copy_structure_and_convert('Classification Nathan', 'Attributes/shell_tools_categories')
