import os

def extract_words_after_rule(file_path):
    words_after_rule = []
    with open(file_path, 'r', encoding='utf-8', errors='ignore') as file:
        lines = file.readlines()
        for line in lines:
            words = line.split()
            if words and words[0] == 'rule':
                if len(words) > 2 and (words[2] == ':' or words[2] == ' :'):
                    words_after_rule.append(words[1])
                elif len(words) > 1 and (words[1][-1] == ':'):
                    words_after_rule.append(words[1][:-1])
    return words_after_rule

def ensure_directory_exists(directory_path):
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)

def save_rules_by_workflow(extracted_rules_dir, workflow_name, rules):
    ensure_directory_exists(extracted_rules_dir)
    with open(os.path.join(extracted_rules_dir, f'{workflow_name}_rules.txt'), 'w') as f:
        for rule in rules:
            f.write(f'{rule}\n')

def main(repo_path):
    extracted_rules_dir = os.path.join(repo_path, 'Extracted_Rules')
    for workflow_folder in os.listdir(repo_path):
        workflow_path = os.path.join(repo_path, workflow_folder)
        if os.path.isdir(workflow_path):  # Check if it's a directory
            all_words = set()
            for root, dirs, files in os.walk(workflow_path):
                for file in files:
                    if file.endswith(('.py', '.java', '.c', '.cpp', '.js', '.ts', '.html', '.css', '.smk')) or file.lower() == "snakefile":
                        file_path = os.path.join(root, file)
                        words_after_rule = extract_words_after_rule(file_path)
                        all_words.update(words_after_rule)
            save_rules_by_workflow(extracted_rules_dir, workflow_folder, all_words)

    print(f'Rules have been saved into the "{extracted_rules_dir}" directory.')

if __name__ == '__main__':
    repo_path = input('Enter the path to the Workflows folder: ')
    main(repo_path)
