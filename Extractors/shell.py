import os

def extract_shell_commands(file_path):
    shell_commands = []
    with open(file_path, 'r', encoding='utf-8', errors='ignore') as file:
        inside_shell_block = False
        current_command = ""
        for line in file:
            stripped_line = line.strip()
            if stripped_line.startswith('shell:'):
                inside_shell_block = True
                current_command = stripped_line[6:].strip()  # Remove 'shell:' part
                if current_command:  # If command starts on the same line
                    current_command += ' '
            elif inside_shell_block:
                if stripped_line.startswith('rule') or stripped_line.startswith('input:') or stripped_line.startswith('output:') or stripped_line.startswith('params:') or stripped_line.startswith('log:') or stripped_line.startswith('benchmark:') or stripped_line.startswith('conda:'):
                    inside_shell_block = False
                    shell_commands.append(current_command.strip())
                    current_command = ""
                else:
                    current_command += line.strip() + ' '
        # Check if the last command was not added
        if current_command:
            shell_commands.append(current_command.strip())
    return shell_commands

def remove_quotes(string):
    return string.replace('"', '').replace("'", "")

def ensure_directory_exists(directory_path):
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)

def save_shell_commands(extracted_commands_dir, file_name, shell_commands):
    ensure_directory_exists(extracted_commands_dir)
    with open(os.path.join(extracted_commands_dir, f'{file_name}_shell_commands.txt'), 'w') as f:
        for command in shell_commands:
            f.write(f'{remove_quotes(command)}\n')

def main(repo_path):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    extracted_commands_dir = os.path.join(script_dir, 'Extracted_Shell_Commands')
    
    for root, dirs, files in os.walk(repo_path):
        for file in files:
            if file.endswith('.py'):
                file_path = os.path.join(root, file)
                shell_commands = extract_shell_commands(file_path)
                if shell_commands:
                    save_shell_commands(extracted_commands_dir, os.path.splitext(file)[0], shell_commands)

    print(f'Shell commands have been saved into the "{extracted_commands_dir}" directory.')

if __name__ == '__main__':
    repo_path = input('Enter the path to the directory: ')
    main(repo_path)
