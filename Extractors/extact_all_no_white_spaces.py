import os

def save_file(content, output_dir, file_name):
    """Saves the content to a text file without white spaces."""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    with open(os.path.join(output_dir, f'{file_name}.txt'), 'w', encoding='utf-8') as f:
        for line in content:
            # Strip white spaces from each line before writing
            f.write(line.strip() + '\n')

def main(repo_path):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.join(script_dir, 'Extracted_Raw_Rules_No_White_Spaces')
    
    for root, dirs, files in os.walk(repo_path):
        for file in files:
            if file.endswith('.py'):
                file_path = os.path.join(root, file)
                with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                    # Read lines and strip white spaces
                    content = [line.strip() for line in f if line.strip()]
                save_file(content, output_dir, os.path.splitext(file)[0])

    print(f'Files have been saved into the "{output_dir}" directory.')

if __name__ == '__main__':
    repo_path = input('Enter the path to the directory: ')
    main(repo_path)
