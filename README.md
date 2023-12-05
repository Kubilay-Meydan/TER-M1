

     _______ ______ _____         __  __ __ 
    |__   __|  ____|  __ \       |  \/  /_ |
       | |  | |__  | |__) |      | \  / || |
       | |  |  __| |  _  /       | |\/| || |
       | |  | |____| | \ \       | |  | || |
       |_|  |______|_|  \_\      |_|  |_||_|
                                         
                                         
Kubilay MEYDAN - Nathan CARRE
              


<div id="badges">
  <a href="Benchling">
    <img src="https://img.shields.io/badge/Benchling-blue?style=flat"/>
  </a>

</div>


# Snakemake rules Extractor

This Python script extracts module names from Snakemake workflow files. It searches through all the files in a specified local repository and extracts words that follow the word "rule" at the beginning of a line, provided that the word is followed by a colon `:` or ` :`. This is the typical snakemake sythax for defining rules or modules

## How to Use

1. Download the script to your local machine.
2. Open your terminal or command prompt and navigate to the location where you downloaded the script.
3. Run the script
4. When prompted, enter the path to your local repository containing the Snakemake workflow files.

The script will then process the files in the specified repository, workflow by workflow and save a list of txt files for each one containing a list of module names. It will save those txt files in a directory called `Extracted_Rules` it creates.

Searches through '.py', '.java', '.c', '.cpp', '.js', '.ts', '.html', '.css', '.smk', '.snakefile' file types `OR` files that have no type but are named "Snakefile". (The capitalisation doesn't matter)


# Workflow Module Comparer GUI

This Python script creates a Graphical User Interface (GUI) for comparing workflow modules stored in .py format within a specified directory. It is an effective tool for anyone needing to analyze and compare different versions or variations of workflow modules.

## Requirements

Python 3.x

Tkinter (typically included with Python)

Pandas: pip install pandas

Pygments: pip install pygments

## Usage

Prepare Directory: Ensure the script is placed in a directory or the directory_path variable in the script correctly references the directory containing the workflow modules.

Execute the Script: Run the script to open the GUI window.

Module Comparison: Two text areas in the GUI will each display a workflow module for visual comparison.

Document Findings: Utilize the buttons ('VS', 'S', 'SS', 'D') beneath the text areas to categorize the similarity between the currently displayed module pair.

Progress Through Modules: After making a selection, the script will load the next pair of modules. Continue this process until you have reviewed all pairs.

Save Results: Once all comparisons are complete, the matrix will be saved as file_comparison_matrix.csv in your designated directory.




# License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.
