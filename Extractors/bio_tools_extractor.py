import os
import json
import requests

# Fonctions pour interagir avec l'API bio.tools
def get_terms_function(dico):
    tab = []
    for d in dico:
        for o in d['operation']:
            tab.append(o['term'].lower())
    return list(set(tab))

def get_terms_topic(dico):
    tab = []
    for d in dico:
        tab.append(d['term'].lower())
    return tab

def get_info_biotools(tool, archive):
    tool = tool.lower()  # Convertir le nom de l'outil en minuscules
    if tool in archive:
        return archive[tool]

    try:
        response = requests.get(f"https://bio.tools/api/{tool}/?format=json")
        if response.status_code == 200:
            dico = response.json()
            res = {'description': dico['description'].lower(), 'function': get_terms_function(dico['function']), 'topic': get_terms_topic(dico['topic'])}
            archive[tool] = res
            return res
        else:
            print(f"Failed to fetch tool info for {tool}. Status code: {response.status_code}")
            return None
    except Exception as e:
        print(f"Error fetching tool info: {e}")
        return None

# Script principal
def read_tools_from_file(file_path):
    with open(file_path, 'r') as file:
        tools = file.read().splitlines()
    return tools

def main(rules_folder_path):
    # Chemin vers le fichier bioweb_archive.json
    bioweb_archive_path = os.path.join(rules_folder_path, "bioweb_archive.json")

    # Chargement ou initialisation de l'archive
    try:
        with open(bioweb_archive_path) as json_file:
            archive = json.load(json_file)
    except FileNotFoundError:
        archive = {}

    tools_info = {}
    for file_name in os.listdir(rules_folder_path):
        if file_name.endswith('_tools_in_shell.txt'):
            file_path = os.path.join(rules_folder_path, file_name)
            tools = read_tools_from_file(file_path)
            for tool in tools:
                tool_info = get_info_biotools(tool, archive)
                if tool_info:
                    tools_info[tool] = tool_info

    # Sauvegarde des informations dans bioweb_archive.json
    with open(bioweb_archive_path, "w") as outfile:
        json.dump(archive, outfile, indent=4)

    # Sauvegarde des résultats dans resultat.json
    json_path = os.path.join(rules_folder_path, 'a_resultat.json')
    with open(json_path, 'w') as file:
        json.dump(tools_info, file, indent=4)
    print(f'Tools information has been saved in "{json_path}".')

if __name__ == '__main__':
    rules_folder_path = input('Enter the path to the Rules folder: ')
    main(rules_folder_path)
