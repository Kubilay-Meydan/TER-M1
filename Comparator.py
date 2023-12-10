import tkinter as tk
from tkinter import scrolledtext, ttk
import os
import itertools
import pandas as pd
from pygments import lex
from pygments.lexers import PythonLexer
from pygments.styles import get_style_by_name

#test
def create_gui(directory):
    # Initialize the matrix and list of file pairs
    files = [f for f in os.listdir(directory) if f.endswith('.py')]
    matrix = pd.DataFrame('SS', index=files, columns=files)
    file_pairs = list(itertools.combinations(files, 2))
    for f in files:
        matrix.at[f, f] = 'I'
    if not file_pairs:
        print("Need at least two Python files in the directory for comparison.")
        return

    # Create the main window
    root = tk.Tk()
    root.title("Python File Comparer")

    # Set the window to fullscreen
    root.state('zoomed')  # For Windows, use 'zoomed'
    # root.attributes('-fullscreen', True)  # Uncomment this line for Linux/MacOS

    # Create scrolled text areas for the two files
    text_area_left = scrolledtext.ScrolledText(root, wrap=tk.WORD, width=40, height=10)
    text_area_right = scrolledtext.ScrolledText(root, wrap=tk.WORD, width=40, height=10)

    # Place the text areas in the window
    text_area_left.grid(column=0, row=0, padx=10, pady=10, sticky="nsew")
    text_area_right.grid(column=1, row=0, padx=10, pady=10, sticky="nsew")

    # Configure grid weights to allow resizing
    root.grid_rowconfigure(0, weight=1)
    root.grid_columnconfigure(0, weight=1)
    root.grid_columnconfigure(1, weight=1)
    
    progress_frame = tk.Frame(root)
    progress_frame.grid(column=0, row=2, columnspan=2, pady=10)
    progress_var = tk.IntVar()
    progress_bar = ttk.Progressbar(progress_frame, orient=tk.HORIZONTAL, length=100, mode='determinate', variable=progress_var)
    progress_bar.pack(fill=tk.X, expand=True)
    total_comparisons = len(file_pairs)
    progress_var.set(0)

    def update_progress_bar():
        completed_comparisons = total_comparisons - len(file_pairs)
        progress_var.set((completed_comparisons / total_comparisons) * 100)
    
    # Function to apply syntax highlighting
    def add_syntax_highlighting(text_area, file_content):
        style = get_style_by_name('default')
        lexer = PythonLexer()
        tokens = lex(file_content, lexer)
        start = '1.0'
        for token_type, value in tokens:
            end = text_area.index(f"{start}+{len(value)}c")
            color = style.style_for_token(token_type)['color']
            if color:
                # Ensure the color is in the correct format for Tkinter
                color = '#' + color
            text_area.tag_add(str(token_type), start, end)
            text_area.tag_config(str(token_type), foreground=color)
            start = end

    # Function to load file content into a text area and apply syntax highlighting
    def load_file_content(file_path, text_area):
        text_area.delete('1.0', tk.END)  # Clear existing content
        with open(file_path, 'r') as file:
            content = file.read()
            text_area.insert(tk.INSERT, content)
            add_syntax_highlighting(text_area, content)


    # Function to save the matrix to a CSV file
    def save_matrix_to_csv():
        csv_path = os.path.join(directory, 'file_comparison_matrix.csv')
        matrix.to_csv(csv_path)
        print(f"Matrix saved to {csv_path}")

    # Function to update the matrix and display the next file pair
    def update_and_advance(value):
        left_file, right_file = file_pairs[0]
        matrix.at[left_file, right_file] = value
        matrix.at[right_file, left_file] = value
        file_pairs.pop(0)  # Remove the current pair
        if file_pairs:
            # Load the next pair of files
            next_pair = file_pairs[0]
            load_file_content(os.path.join(directory, next_pair[0]), text_area_left)
            load_file_content(os.path.join(directory, next_pair[1]), text_area_right)
            update_progress_bar()
        else:
            print("All file comparisons completed.")
            save_matrix_to_csv()
            root.destroy()  # Close the GUI

    # Add buttons for updating the matrix
    buttons_frame = tk.Frame(root)
    buttons_frame.grid(column=0, row=1, columnspan=2, pady=10)
    for value in ["VS", "S", "SS", "D"]:
        button = tk.Button(buttons_frame, text=value, command=lambda v=value: update_and_advance(v))
        button.pack(side=tk.LEFT, padx=5)

    # Load the first pair of files
    first_pair = file_pairs[0]
    load_file_content(os.path.join(directory, first_pair[0]), text_area_left)
    load_file_content(os.path.join(directory, first_pair[1]), text_area_right)

    # Start the GUI event loop
    root.mainloop()

# Replace this path with the path to your directory containing Python files
directory_path = 'Classification Kubilay/1_Quality_Control'
create_gui(directory_path)

