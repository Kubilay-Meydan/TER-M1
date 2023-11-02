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

def main(repo_path):
    words_by_file = {}
    all_words = set()
    for root, dirs, files in os.walk(repo_path):
        for file in files:
            if file.endswith(('.py', '.java', '.c', '.cpp', '.js', '.ts', '.html', '.css', '.smk')):
                file_path = os.path.join(root, file)
                words_after_rule = extract_words_after_rule(file_path)
                if words_after_rule:
                    words_by_file[file] = words_after_rule
                    all_words.update(words_after_rule)
    
    print('Words after "rule" sorted by file:')
    for file, words_after_rule in words_by_file.items():
        print(f'{file}: {words_after_rule}')
    
    print('\nAll words found in all the files:', sorted(all_words), "\nNombre de rules trouvés", len(all_words))

if __name__ == '__main__':
    repo_path = input('Enter the path to the local repository: ')
    main(repo_path)
