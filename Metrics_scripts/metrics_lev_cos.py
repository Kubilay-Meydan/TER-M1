import os
import Cosine2D  # Replace with your actual script's name
import Levenshtein2D  # Replace with your actual script's name

def process_folders(folders):
    script_directory = os.path.dirname(os.path.realpath(__file__))  # Directory of the script

    for folder in folders:
        # Create results folder in the script directory
        results_folder_name = f"{os.path.basename(folder)}_results"
        results_folder = os.path.join(script_directory, results_folder_name)
        os.makedirs(results_folder, exist_ok=True)

        # Run cosine similarity script
        cosine_csv = os.path.join(results_folder, f"{os.path.basename(folder)}_cosine.csv")
        Cosine2D.construire_matrice_similarite_et_sauvegarder(folder, cosine_csv)

        # Run Levenshtein similarity script
        levenshtein_csv = os.path.join(results_folder, f"{os.path.basename(folder)}_levenshtein.csv")
        Levenshtein2D.build_and_save_levenshtein_similarity_matrix(folder, levenshtein_csv)

# List of folders to process
folders_to_process = ['Attributes/shell_tools_description_biotools_no_white_space', 'Attributes/shell_tools_description_biotools']  # Replace with your actual folder paths

process_folders(folders_to_process)
