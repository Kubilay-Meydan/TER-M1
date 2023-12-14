import os
import re

def extract_input_file_types(file_path):
    file_types = []
    with open(file_path, 'r', encoding='utf-8', errors='ignore') as file:
        inside_input_block = False
        for line in file:
            line = line.strip()
            if line.startswith('output:'):
                inside_input_block = True
            elif inside_input_block and (line.startswith('input:') or line.startswith('params:') or line.startswith('log:') or line.startswith('benchmark:') or line.startswith('conda:') or line.startswith('shell:')):
                inside_input_block = False
            elif inside_input_block:
                matches = re.findall(r'\.\w+(?=["\'])', line)
                # Remove the leading '.' and add each occurrence to the list
                cleaned_types = [match[1:] for match in matches]
                file_types.extend(cleaned_types)
    return file_types

def ensure_directory_exists(directory_path):
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)

def save_file_types(extracted_types_dir, file_name, file_types):
    ensure_directory_exists(extracted_types_dir)
    with open(os.path.join(extracted_types_dir, f'{file_name}_file_types.txt'), 'w') as f:
        for file_type in file_types:
            f.write(f'{file_type}\n')

def main(repo_path):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    extracted_types_dir = os.path.join(script_dir, 'Extracted_File_Types')
    
    for root, dirs, files in os.walk(repo_path):
        for file in files:
            if file.endswith('.py'):
                file_path = os.path.join(root, file)
                file_types = extract_input_file_types(file_path)
                if file_types:
                    save_file_types(extracted_types_dir, os.path.splitext(file)[0], file_types)

    print(f'File types have been saved into the "{extracted_types_dir}" directory.')

if __name__ == '__main__':
    repo_path = input('Enter the path to the directory: ')
    main(repo_path)
