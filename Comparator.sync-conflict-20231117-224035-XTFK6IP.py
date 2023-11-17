import tkinter as tk
import itertools
import os
import csv
import glob

# Specify the directory containing the files
directory = '/Users/kubilaymeydan/Desktop/M1 Bibs/S1/TER/TER-M1/Rules_classified/QC qnd trimming'

# Function to read all .py files from the directory
def get_files_from_directory(directory):
    return [os.path.basename(file) for file in glob.glob(os.path.join(directory, '*.py'))]

# Function to load or create the comparison matrix
def load_or_create_matrix(files, directory):
    matrix_path = os.path.join(directory, 'comparison_matrix.csv')
    if os.path.exists(matrix_path):
        with open(matrix_path, mode='r', newline='') as file:
            reader = csv.reader(file)
            return {frozenset({row[0], row[1]}): row[2] for row in reader}
    else:
        return {frozenset({file1, file2}): 'SS' for file1, file2 in itertools.combinations(files, 2)}

# Function to save the matrix to a CSV file
def save_matrix_to_csv(matrix, directory):
    matrix_path = os.path.join(directory, 'comparison_matrix.csv')
    with open(matrix_path, mode='w', newline='') as file:
        writer = csv.writer(file)
        for pair, value in matrix.items():
            writer.writerow(list(pair) + [value])

# Function to update the matrix and display next pair
def update_matrix_and_display_next(value):
    global current_pair
    comparison_matrix[frozenset(current_pair)] = value
    save_matrix_to_csv(comparison_matrix, directory)
    display_next_pair()

# Function to display the next pair of files
def display_next_pair():
    global current_pair
    try:
        current_pair = next(pairs_iterator)
        file1_label.config(text=current_pair[0])
        file2_label.config(text=current_pair[1])
    except StopIteration:
        file1_label.config(text="Done")
        file2_label.config(text="Done")

# Main
files = get_files_from_directory(directory)
comparison_matrix = load_or_create_matrix(files, directory)

# Set up the main window
root = tk.Tk()
root.title("File Comparison")

# Labels to display file names
file1_label = tk.Label(root, text="", font=("Helvetica", 16))
file1_label.pack()
file2_label = tk.Label(root, text="", font=("Helvetica", 16))
file2_label.pack()

# Buttons for comparison choices
buttons_frame = tk.Frame(root)
buttons_frame.pack(pady=10)

for value in ["VS", "S", "SS", "D"]:
    button = tk.Button(buttons_frame, text=value, command=lambda v=value: update_matrix_and_display_next(v))
    button.pack(side=tk.LEFT, padx=5)

# Initialize the iterator and display the first pair
pairs_iterator = iter(itertools.combinations(files, 2))
current_pair = None
display_next_pair()

# Start the GUI event loop
root.mainloop()
