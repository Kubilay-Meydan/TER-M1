import os

def rename_files(folder_path, common_ending):
    for filename in os.listdir(folder_path):
        if filename.endswith(common_ending):
            new_name = filename.replace(common_ending, '')
            os.rename(os.path.join(folder_path, filename), os.path.join(folder_path, new_name))
    print(f"'{common_ending}' deleted from all file names")

# Replace with your folder path and common ending
path_to_folder = 'Attributes/Whole_Text_Shell'
common_ending = '_shell_commands.txt'  # Example: '_old.txt'

rename_files(path_to_folder, common_ending)
