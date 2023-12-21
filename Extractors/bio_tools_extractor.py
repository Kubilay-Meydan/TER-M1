import os
import json
import requests

# Fonctions pour interagir avec l'API bio.tools
def get_terms_function(dico):
    return [o['term'].lower() for d in dico for o in d['operation']]

def get_terms_topic(dico):
    return [d['term'].lower() for d in dico]

def get_info_biotools(tool):
    tool = tool.lower()  # Convertir le nom de l'outil en minuscules
    try:
        response = requests.get(f"https://bio.tools/api/{tool}/?format=json")
        if response.status_code == 200:
            dico = response.json()
            return {
                'description': dico['description'].lower(),
                'function': get_terms_function(dico['function']),
                'topic': get_terms_topic(dico['topic'])
            }
        else:
            print(f"Failed to fetch tool info for {tool}. Status code: {response.status_code}")
            return None
    except Exception as e:
        print(f"Error fetching tool info: {e}")
        return None

# Script principal
def read_tools_from_file(file_path):
    with open(file_path, 'r') as file:
        return file.read().splitlines()

def main(rules_folder_path):
    unique_tools = set()
    for file_name in os.listdir(rules_folder_path):
        if file_name.endswith('_tools_in_shell.txt'):
            file_path = os.path.join(rules_folder_path, file_name)
            unique_tools.update(read_tools_from_file(file_path))

    tools_info = {tool: get_info_biotools(tool) for tool in unique_tools if get_info_biotools(tool)}

    # Sauvegarde des résultats dans resultat.json
    json_path = os.path.join(rules_folder_path, 'a_resultat_biotools.json')
    with open(json_path, 'w') as file:
        json.dump(tools_info, file, indent=4)

    print(f'Tools information has been saved in "{json_path}".')

if __name__ == '__main__':
    rules_folder_path = input('Enter the path to the Rules folder: ')
    main(rules_folder_path)
