import os

def extract_conda_commands(file_path):
    conda_commands = []
    with open(file_path, 'r', encoding='utf-8', errors='ignore') as file:
        inside_conda_block = False
        current_command = ""
        for line in file:
            stripped_line = line.strip()
            if stripped_line.startswith('conda:'):
                inside_conda_block = True
                current_command = stripped_line[6:].strip()  # Remove 'conda:' part
            elif inside_conda_block:
                if stripped_line.startswith('rule') or stripped_line.startswith('input:') or stripped_line.startswith('output:') or stripped_line.startswith('params:') or stripped_line.startswith('log:') or stripped_line.startswith('benchmark:') or stripped_line.startswith('shell:'):
                    inside_conda_block = False
                    conda_commands.append(current_command.strip())
                    current_command = ""
                else:
                    current_command += line.strip()
        # Check if the last command was not added
        if current_command:
            conda_commands.append(current_command.strip())
    return conda_commands

def ensure_directory_exists(directory_path):
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)

def save_conda_commands(extracted_commands_dir, file_name, conda_commands):
    ensure_directory_exists(extracted_commands_dir)
    with open(os.path.join(extracted_commands_dir, f'{file_name}_conda_commands.txt'), 'w') as f:
        for command in conda_commands:
            f.write(f'{command}\n')

def main(repo_path):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    extracted_commands_dir = os.path.join(script_dir, 'Extracted_Conda_Commands')

    for root, dirs, files in os.walk(repo_path):
        for file in files:
            if file.endswith('.py'):
                file_path = os.path.join(root, file)
                conda_commands = extract_conda_commands(file_path)
                if conda_commands:
                    save_conda_commands(extracted_commands_dir, os.path.splitext(file)[0], conda_commands)

    print(f'Conda commands have been saved into the "{extracted_commands_dir}" directory.')

if __name__ == '__main__':
    repo_path = input('Enter the path to the directory: ')
    main(repo_path)
